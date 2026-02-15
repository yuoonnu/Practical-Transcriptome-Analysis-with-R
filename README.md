# 🧬 Transcriptome Analysis Pipeline
R 언어와 Bioconductor 패키지를 활용하여 RNA-seq 데이터(Count Matrix)를 정규화하고, 차등 발현 유전자(DEG)를 선별하여 생물학적 기능을 해석하는 전사체 분석 파이프라인을 구축했습니다.

## 👤 Identity & Goals
- **Interest**: Computational Biology & Bioinformatics
- **Objective**: 대용량의 유전체 데이터를 IT 기술로 다루며, 통계적 기법을 통해 질병 메커니즘과 연관된 유의미한 바이오 마커를 발굴하는 역량을 기릅니다.

## 🛠 Tech Stack
- **Language**: `R`
- **Libraries**: `BiocManager`, `TCC` (Normalization), `edgeR` (Statistical Test), `ggplot2` & `pheatmap` (Visualization), `clusterProfiler` (Enrichment)
- **Core Skills**: RNA-seq Workflow, TMM Normalization, Dimensionality Reduction (PCA), DEG Analysis, GO Enrichment

## 📂 Project Structure & Workflow

### 1. Data Acquisition & Preprocessing
* **ArrayExpress 데이터 로드**: `E-MTAB-5338` 데이터셋의 Processed Count Matrix를 활용했습니다.
* **Filtering**: 저발현 유전자(Low count genes)를 제거하여 분석의 신뢰도를 높였습니다.
    * *Criterion*: Read count ≥ 5인 샘플이 전체의 75% 이상인 유전자 선별.
* **Normalization (TMM)**: `TCC` 패키지를 사용하여 샘플 간 라이브러리 크기 차이에 따른 편향을 보정(Trimmed Mean of M-values)했습니다.

### 2. Exploratory Data Analysis (EDA)
* **Box Plot**: 정규화 전/후의 데이터 분포를 비교 시각화하여 전처리의 적절성을 검증했습니다.
* **PCA (Principal Component Analysis)**: 고차원 데이터를 축소하여 PC1, PC2 상에서 샘플 그룹(Group) 간의 군집 경향을 확인했습니다.
* **Hierarchical Clustering**: Dendrogram을 통해 샘플 간의 거리를 계층적으로 시각화했습니다.

### 3. DEG Analysis (Differential Expression Gene)
* **Statistical Testing**: `edgeR` 알고리즘을 적용하여 세포 타입(Cell Type) 간의 발현 차이를 검정했습니다.
* **Gene Selection**:
    * **FDR (q-value) < 0.05**: 다중 검정 오류를 보정한 유의수준 적용.
    * **|log2FC| ≥ 1**: 2배 이상 발현 차이가 나는 유전자 선별.
* **Visualization**: `pheatmap`을 활용하여 선별된 DEG들의 발현 패턴을 Z-score 기반의 Heatmap으로 표현했습니다.

### 4. Functional Enrichment Analysis
* **GO (Gene Ontology) Analysis**: `clusterProfiler`를 사용하여 과발현된 유전자군(Up-regulated genes)이 관여하는 생물학적 프로세스(Biological Process)를 규명했습니다.
* **Interpretation**: Bubble Plot을 통해 유의한 GO Term과 Gene Ratio를 시각화하여 통계적 결과의 생물학적 의미를 도출했습니다.

<br>

## 💡 Key Experience
- **Bioinformatics Pipeline 구축**: Raw Count 데이터부터 기능 분석까지 이어지는 RNA-seq 분석의 전체 흐름을 코드로 구현했습니다.
- **Data-Driven Insight**: 단순히 수치를 산출하는 것을 넘어, PCA와 Heatmap 등 다양한 시각화 기법을 통해 데이터의 잠재적 구조를 파악하는 능력을 길렀습니다.

---
**Blog**: [밈뮴 블로그 링크]()
