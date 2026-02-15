### Transcriptome Analysis Pipeline with R
### Project: Differential Expression Analysis of RNA-seq Data (E-MTAB-5338)

## Step 1. Setup Environment & Load Data
# 패키지 설치 및 로드 (Check and Install required packages)
if(!requireNamespace("BiocManager", quietly = T))
  install.packages("BiocManager")

packages <- c("ArrayExpress", "TCC", "reshape2", "ggplot2", "ggrepel", "ggdendro", "pheatmap", "clusterProfiler", "org.Mm.eg.db")

for(p in packages){
  if(!requireNamespace(p, quietly = T)) BiocManager::install(p)
}

library(ArrayExpress)

# 분석 데이터 로드
dat <- ArrayExpress::getAE("E-MTAB-5338", type="processed")
dataDF <- read.table(dat$processedFiles, sep='\t', row.names=1, header=T, check.names=F)
head(dataDF)

# 메타 데이터 로드
metadata <- read.table('metadata.txt', sep='\t', quote='', header=T)

# 메타 데이터 내 샘플 ID 순서로 정렬
dataDF <- dataDF[metadata$Sample.Name]


## Step 2. Preprocessing & Filtering
# 추출된 데이터프레임을 매트릭스로 전환
dataDF <- data.matrix(dataDF)

# Read count가 5 이상인 샘플이 전체의 75% 이상인 유전자만 유지
dataDF <- dataDF[rowSums(dataDF >= 5) >= ncol(dataDF)*0.75,]


## Step 3. Normalization (TMM)
library(TCC)

# TCC 객체 생성 및 TMM 정규화
tcc <- new('TCC', dataDF, group=as.numeric(factor(metadata$Group)))
tcc <- calcNormFactors(tcc, norm.method='tmm')

# 정규화 데이터 추출
normDF <- getNormalizedData(tcc)


## Step 4. Visualization
library(reshape2)
library(ggplot2)

# 비표준화 데이터 재구조화
dataDF_melt <- melt(dataDF)
colnames(dataDF_melt) <- c('gene', 'sample', 'count')
dataDF_melt <- merge(dataDF_melt, metadata[c('Sample.Name', 'Group')], by.x='sample', by.y='Sample.Name')

# 표준화 데이터 재구조화
normDF_melt <- melt(normDF) 
colnames(normDF_melt) <- c('gene', 'sample', 'count') 
normDF_melt <- merge(normDF_melt, metadata[c('Sample.Name', 'Group')], by.x='sample', by.y='Sample.Name') 

# Unnormalized value
p1 <- ggplot(dataDF_melt)+
  geom_boxplot(aes(x=sample, y=log2(count+1), fill=Group))+
  labs(title='counts before normalization')+
  theme_bw()+
  theme(axis.text=element_text(angle=90))
ggsave('Fig1.Boxplot-1.png', p1, width=5, height=4)

# Normalized value
p2 <- ggplot(normDF_melt)+
  geom_boxplot(aes(x=sample, y=log2(count+1), fill=Group))+
  labs(title='counts after normalization')+
  theme_bw()+
  theme(axis.text=element_text(angle=90))
ggsave('Fig2.Boxplot-2.png', p2, width=5, height=4)

# PCA plot
library(ggrepel)
log2DF <- log2(normDF+1)
PCA <- prcomp(t(log2DF))
pcaDF <- merge(PCA$x, metadata[c('Sample.Name', 'Group')], by.x=0, by.y='Sample.Name')

p3 <- ggplot(pcaDF)+
  geom_point(aes(x=PC1, y=PC2, color=Group), size=3)+
  geom_hline(yintercept=0, color='black', linetype='dashed')+
  geom_vline(xintercept=0, color='black', linetype='dashed')+
  geom_text_repel(aes(x=PC1, y=PC2, label=Row.names), size=3)+
  labs(x=paste('PC1(', round(PCA$sdev[1]/sum(PCA$sdev)*100, 2), '%)', sep=''),
       y=paste('PC2(', round(PCA$sdev[2]/sum(PCA$sdev)*100, 2), '%)', sep=''))+
  theme_bw()
ggsave('Fig3.PCAplot.png', p3, width=6, height=4)

# dendrogram
library (ggdendro)
hc <- hclust(dist(t(log2DF)), method='average')
ggsave('Fig4.dendrogram.png', ggdendrogram(hc), width=4, height=4)

## Step 5. DEG Analysis
# 분석 데이터를 TCC class로 변환
cellTCC <- new('TCC', dataDF, group=as.numeric(factor(metadata$Cell.Type)))
cellTCC <- calcNormFactors(cellTCC, norm.method='tmm')

# edger 패키지를 이용한 DEG 추정
cellTCC <- estimateDE(cellTCC, test.method = 'edger', FDR=0.05)
cellTCC <- getResult(cellTCC, sort=F)

# DEG 선별 및 Z-score
cellDEG <- cellTCC[which(abs(cellTCC$m.value) >= 2 & cellTCC$q.value <= 0.05),]
cellDEG <- log2DF[rownames(log2DF)%in% cellDEG$gene_id,]
cellDEG <- t(scale(t(cellDEG)))

# DEG Heatmap
if (!requireNamespace("pheatmap", quietly = TRUE))
  install.packages("pheatmap")
library(pheatmap)
pheatmap(cellDEG, 
         main='DEGs by cell types',
         color=colorRampPalette(c('green', 'black', 'red'))(100),
         show_rownames=F, cluster_cols=F,
         clustering_method='ward.D',
         treeheight_row=10, cutree_rows=4,
         gaps_col=c(3, 6, 9),
         filename='Fig5.DEG-celltype.png',
         width=5, height=5)


## Step 6. Functional Enrichment Analysis
library(clusterProfiler)

# CD4 세포에서 Down-regulation 된 유전자
CD4_up <- as.numeric(cellTCC$gene_id[cellTCC$m.value <= -2 & cellTCC$q.value <= 0.05])

# GO Enrichment
CD4_GSEA <- enrichGO(CD4_up, 'org.Mm.eg.db', ont='BP', keyType='ENTREZID', minGSSize = 5, maxGSSize = 200)
CD4_GSEA <- CD4_GSEA@result
CD4_GSEA <- CD4_GSEA[CD4_GSEA$qvalue <= 0.05,]
CD4_GSEA <- CD4_GSEA[order(CD4_GSEA$GeneRatio, decreasing=T),]

# Bubble Plot
CD4_GSEA$ID <- factor(CD4_GSEA$ID, levels = CD4_GSEA$ID)

p6 <- ggplot(CD4_GSEA[1:10,])+
  geom_point(aes(x=ID, y=GeneRatio, size=Count, color=qvalue))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90))
ggsave('Fig6.CD4-GOenrich.png', p6, width=5, height = 4)
write.table(CD4_GSEA, 'Table1.CD4-GOenrich.txt', sep='\t', quote=F, col.names=T)
