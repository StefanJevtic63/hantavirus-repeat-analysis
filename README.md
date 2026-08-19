# Hantavirus Genomic Repeat Analysis

This repository contains a modular Python and Jupyter Notebook pipeline for extracting, preprocessing, and clustering repetitive sequences from Hantavirus genomic data.

The study focuses on four Hantavirus species:

* *Orthohantavirus hantanense* (Hantaan virus)
* *Orthohantavirus dobravaense* (Dobrava-Belgrade virus)
* *Orthohantavirus puumalaense* (Puumala virus)
* *Orthohantavirus sinnombreense* (Sin Nombre virus)

The main objective is to investigate whether genomic repetitive-sequence patterns can be used to identify structural similarities and naturally occurring groups among Hantavirus nucleotide sequences.

---

## Overview

The analysis is organized as a reproducible and modular pipeline consisting of:

- FASTA sequence preprocessing
-  Repetitive-sequence extraction using **StatRepeats**
- Construction of binary repeat-presence matrices
- Repeat-length filtering
- **TF-IDF** transformation
- Dimensionality reduction using **PCA**
- Hyperparameter optimization using **Grid Search**
- Unsupervised clustering
- Evaluation using clustering quality metrics
- Analysis of the most characteristic repetitive sequences

The complete workflow is designed to operate on highly dimensional and sparse genomic repeat representations.

---

## Features & Workflow

### Sequence preprocessing

Raw genomic sequences are provided in FASTA format.

The preprocessing stage:

* validates the nucleotide alphabet;
* retains only `A`, `C`, `G`, `T`, and `U`;
* removes invalid symbols and structural inconsistencies;
* generates consistent `Sequence_ID` to virus-species mappings.

The resulting cleaned FASTA files are subsequently processed by StatRepeats.

### Repetitive-sequence extraction

The main feature-extraction component is **StatRepeats**, a specialized tool for identifying statistically significant maximal repetitive sequences.

Four structural repeat types are analyzed:

* **Direct Non-Complementary (`dn`)**
* **Direct Complementary (`dc`)**
* **Inverse Non-Complementary (`in`)**
* **Inverse Complementary (`ic`)**

StatRepeats outputs are subsequently transformed into a unified representation of repetitive sequences.

The StatRepeats tool is available at:

http://bioinfo.matf.bg.ac.rs/home/downloads.wafl?cat=Software&project=statrepeats

### Binary repeat matrix construction

For every genomic sequence, the pipeline determines which repetitive sequences are present.

A binary matrix is constructed where:

* each row represents one genome;
* each column represents one unique repetitive sequence;
* `1` indicates that the repeat is present;
* `0` indicates that the repeat is absent.

The resulting matrices are extremely high-dimensional and sparse.

### Repeat-length filtering

Repeat-length distributions are analyzed before clustering.

For `all_sequences`, the selected repeat lengths are:

```text
12, 14, 16, 18, 20
```

For `complete_sequences`, the selected repeat lengths are:

```text
10, 12, 14, 16, 18, 20
```

In both datasets, repeats of length 16 are the most frequent.

### TF-IDF transformation

Because the binary repeat matrices are highly sparse, **TF-IDF** is applied before clustering.

TF-IDF reduces the influence of repetitive sequences that occur in most genomes while increasing the relative importance of less frequent and potentially more discriminative repeats.

The transformed representation remains sparse throughout the preprocessing pipeline.

### PCA dimensionality reduction

Principal Component Analysis (**PCA**) is used to further reduce the dimensionality of the transformed feature space.

The following cumulative explained-variance levels are evaluated:

```text
5%, 10%, 15%, 20%, 40%, 60%
```

The number of retained components is selected according to the minimum number required to reach each target variance level.

---

## Datasets

Two final datasets are used in the experiments.

| Dataset                  | Description                                         | Samples | Repeat Features |
| :----------------------- | :-------------------------------------------------- | ------: | --------------: |
| **`all_sequences`**      | Partial, fragmented, and complete genomic sequences |   7,087 |          47,135 |
| **`complete_sequences`** | Fully sequenced and verified complete genomes       |     435 |          62,364 |

The `all_sequences` dataset provides a larger and more diverse sample set, while `complete_sequences` contains fewer samples but a larger number of unique repeat features. Both datasets are processed as sparse representations.

### Dataset files

| File                        | Role                                         |
| :-------------------------- | :------------------------------------------- |
| `all_sequences.fasta`       | Raw genomic sequences                        |
| `all_sequences_clean.fasta` | Validated FASTA sequences                    |
| `all_dn.txt`                | StatRepeats direct non-complementary output  |
| `all_dc.txt`                | StatRepeats direct complementary output      |
| `all_in.txt`                | StatRepeats inverse non-complementary output |
| `all_ic.txt`                | StatRepeats inverse complementary output     |
| `all_sequences.csv`         | Binary repeat matrix for all sequences       |
| `complete_sequences.csv`    | Binary repeat matrix for complete sequences  |

---

## Clustering Models

Five unsupervised clustering algorithms are evaluated:

### K-Means

Partitions the data into a predefined number of clusters by minimizing within-cluster squared distances to centroids.

Parameters investigated:

```text
n_clusters = {2, 3, 4, 6, 8}
init       = k-means++
n_init     = 10
max_iter   = 300
```

### Agglomerative Clustering

A hierarchical bottom-up method that initially treats every observation as a separate cluster and progressively merges the closest clusters.

Parameters investigated:

```text
n_clusters = {2, 3, 4, 6, 8}
linkage    = {ward, average, complete}
metric     = {euclidean, cosine}
```

### BIRCH

A hierarchical clustering method based on a tree of summarized statistics, designed to efficiently process large datasets.

Parameters investigated:

```text
number_clusters = {2, 3, 4, 6, 8}
diameter        = {0.1, 0.3, 0.5}
```

### DBSCAN

A density-based clustering algorithm capable of identifying clusters of arbitrary shape and treating low-density observations as noise.

The `eps` parameter is estimated from the distribution of distances to the k-nearest neighbors.

```text
neighbors = {3, 5, 10}
```

### CLIQUE

A grid-based clustering method designed for high-dimensional datasets. It identifies dense regions in subspaces and combines neighboring dense regions.

Parameters investigated:

```text
amount_intervals    = {3, 4}
density_threshold   = {1, 2}
```

The hyperparameter spaces are systematically evaluated using Grid Search.

---

## Evaluation Metrics

The primary optimization criterion is the **Silhouette Coefficient**.

The following metrics are calculated for each clustering configuration:

### Silhouette Coefficient

Measures both intra-cluster compactness and inter-cluster separation.

Higher values indicate better-defined clusters.

### Davies-Bouldin Index

Measures the relationship between within-cluster dispersion and between-cluster separation.

Lower values indicate better clustering.

### Calinski-Harabasz Index

Measures the ratio between between-cluster dispersion and within-cluster dispersion.

Higher values indicate better-separated clusters.

Silhouette is used to select the optimal hyperparameter configuration, while Davies-Bouldin and Calinski-Harabasz provide additional evaluation of the selected solutions.

---

## Experimental Results

The best configurations obtained during the experiments are summarized below.

| Model         | Dataset              | Repeat Filtering | PCA Variance | Components | Silhouette |
| :------------ | :------------------- | :--------------: | -----------: | ---------: | ---------: |
| K-Means       | `all_sequences`      |        Yes       |           5% |          3 | **0.7829** |
| Agglomerative | `all_sequences`      |        Yes       |           5% |          3 | **0.7950** |
| **BIRCH**     | **`all_sequences`**  |      **Yes**     |       **5%** |      **3** | **0.7974** |
| DBSCAN        | `complete_sequences` |        No        |          40% |         19 | **0.6077** |
| CLIQUE        | `all_sequences`      |        Yes       |          10% |          5 | **0.7004** |

### Best K-Means configuration

```text
Dataset:       all_sequences
Repeat filter: enabled
PCA variance:  5%
Components:    3
n_clusters:    2
init:          k-means++
n_init:        10
max_iter:      300
Silhouette:    0.7829
```

### Best Agglomerative configuration

```text
Dataset:       all_sequences
Repeat filter: enabled
PCA variance:  5%
Components:    3
n_clusters:    3
linkage:       average
metric:        euclidean
Silhouette:    0.7950
```

### Best BIRCH configuration

```text
Dataset:       all_sequences
Repeat filter: enabled
PCA variance:  5%
Components:    3
number_clusters: 3
diameter:        0.3
Silhouette:      0.7974
```

### Best DBSCAN configuration

```text
Dataset:       complete_sequences
Repeat filter: disabled
PCA variance:  40%
Components:    19
eps:           0.0050278427
neighbors:     3
Silhouette:    0.6077
```

### Best CLIQUE configuration

```text
Dataset:       all_sequences
Repeat filter: enabled
PCA variance:  10%
Components:    5
amount_intervals:  4
density_threshold: 1
Silhouette:       0.7004
```

---

## Analysis of Characteristic Repeats

The most significant repetitive sequences are identified by comparing their local representation within a cluster with their global representation across the dataset.

Several repetitive sequences consistently appear among the most characteristic features identified by K-Means, Agglomerative, and BIRCH, including:

```text
AAAGAAAAGAAA
AAGAAAAAAGAA
AGAAGAAGAAGA
AAGGAAAAGGAA
TATTATTATTAT
AAAAAAAAAAAA
GAAAAGGAAAAG
AGAAAAAAAAGA
```

The agreement between K-Means, Agglomerative, and BIRCH indicates that the identified repeat patterns are relatively stable across different clustering approaches.

Agglomerative and BIRCH generally produce larger differences between local and global repeat representation than K-Means, suggesting stronger cluster-specific repeat patterns.

---

## Conclusions

The experiments demonstrate that repetitive genomic sequences can provide useful structural information for unsupervised analysis of Hantavirus genomes.

The main conclusions are:

1. **BIRCH achieved the best overall Silhouette coefficient**, reaching **0.7974** on the `all_sequences` dataset.
2. **Agglomerative Clustering** achieved a very similar result of **0.7950**, followed by **K-Means with 0.7829**.
3. The best-performing configurations predominantly use the `all_sequences` dataset with repeat-length filtering and strong dimensionality reduction.
4. **DBSCAN** performed substantially better on `complete_sequences` than on the `all_sequences` configuration, but its best Silhouette score of **0.6077** remained below the best hierarchical and centroid-based approaches.
5. **CLIQUE** achieved a Silhouette coefficient of **0.7004**.
6. The similarity of the most significant repeats identified by K-Means, Agglomerative, and BIRCH supports the stability of the detected genomic patterns.
7. The results demonstrate the importance of preprocessing and dimensionality reduction for clustering highly dimensional sparse repeat representations.

Overall, **BIRCH applied to the filtered `all_sequences` dataset after TF-IDF transformation and PCA reduction to three components represents the best-performing configuration in the presented experiments.**

---

## Prerequisites & Installation

Python **3.10+** is recommended.

Clone the repository:

```bash
git clone https://github.com/StefanJevtic63/hantavirus-repeat-analysis.git
cd hantavirus-repeat-analysis
```

Install the required dependencies:

```bash
pip install -r requirements.txt
```

The main experimental preprocessing is implemented in the Jupyter notebook:

```text
clustering_sparse.ipynb
```

