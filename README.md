# STEVE
A quantitative framework developed to evaluate the accuracy and reliability of cell-type annotations.

Version: 1.0.0

Updated date: 2026.12.06

Citation: Elijah Torbenson, Xiao Ma, Daniel J. Garry, Stephen C. Jameson, Zhengdong Zhang, Laura J. Niedernhofer, Lei Zhang, Meiyi Li, and Xiao Dong. STEVE: Single-cell Transcriptomics Expression Visualization and Evaluation. in review at Briefings in Bioinformatics.
#####
## Author and License

Author: Elijah Torbenson? 

Email: elijahtorbenson@creighton.edu

Licensed under the GNU Affero General Public License version 3 or later.

#####
## Dependencies
• R 4.4.1

• dplyr 1.1.4

• tidyr 1.3.1

• ggplot2 3.5.1

• lattice 0.22.6

• caret 6.0.94

• e1071 1.7.14


## Release Notes

• v1.0.0, 2024.12.06, 1st version uploades.



# STEVE — Conceptual Overview

STEVE is a dataset-specific evaluation framework for **cell-type annotation** in scRNA-seq.

Instead of declaring a single “best” annotation tool, STEVE quantifies **accuracy, robustness, and uncertainty** of annotation *within the context of a full analysis pipeline*:

> normalization → HVGs → scaling → PCA → integration → UMAP → probabilistic mapping

---

# Framework Overview

STEVE includes four modules:

1. **Subsampling Evaluation** — how stable annotation is as reference size (and partitions) change  
2. **Novel Cell Evaluation** — whether a cell type absent from the reference is correctly flagged as **unknown**  
3. **Annotation Benchmarking** — compares **candidate annotations** (e.g., SingleR / ScType / manual) to **ground truth**  
4. **Reference Transfer Annotation** — transfers ground-truth labels from an external reference dataset to a new dataset  

---

# Core Mapping Model

All modules use the same probabilistic mapping framework:

1. Integrate reference + query into a shared 2D embedding (UMAP by default)
2. Estimate per-cell-type KDE density in the reference
3. Compute posterior probabilities per query cell
4. Assign a label only if **posterior odds** exceed a threshold  
   Otherwise label as **unknown**

---

# Installation

STEVE runs from a single script:

```
steve_master.R
```

## Required Packages

- `Seurat`
- `SeuratDisk`
- `sctransform`
- `MASS` (for `kde2d`)
- `harmony` (if using HarmonyIntegration)
- `ComplexHeatmap`
- `circlize`
- `ggplot2`
- `tidyr`
- `magrittr`

Optional:
- `scCustomize` (improved plotting palettes)

### Installation Example

```r
install.packages(c(
  "optparse",
  "Seurat",
  "SeuratDisk",
  "sctransform",
  "MASS",
  "harmony",
  "ggplot2",
  "tidyr",
  "magrittr"
))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("ComplexHeatmap", "circlize"))

install.packages("scCustomize") # optional
```

---

# Input Format

STEVE expects:

```
.h5seurat
```

### Key Principle: Labels Live in `meta.data`

Each Seurat object contains:

```
obj@meta.data
```

- Each row = one cell  
- Each column = metadata annotation  

STEVE reads required label columns from `meta.data`.

---

# Required Metadata Columns

Depending on module:

---

## A) Ground-Truth Label Column  
(Required for evaluation modules)

This is the **best available cell-type identity**, ideally:

- Experimentally defined (e.g., FACS-sorted)
- Expert-curated

Default column name:
```
groundtruth
```

Override using:
```
--reference_annotation
```

---

## B) Candidate Annotation Column  
(Required for Annotation Benchmarking only)

This is the annotation you want to evaluate, e.g.:

- Labels from SingleR  
- Labels from ScType  
- Manual labels  
- Any other annotation pipeline output  

Default column name:
```
user_annotation
```

Override using:
```
--user_annotation
```

---

## C) Optional Truth Column for User Dataset  
(Transfer evaluation only)

If the user/query dataset also has known truth labels, STEVE can compute post-transfer metrics.

Provide using:
```
--user_truth_annotation
```

---

# Quickstart (CLI)

All modules run through the same script:

```bash
Rscript steve_master.R --module <subsampling|novel|benchmark|transfer> \
  --userdata <path/to/dataset.h5seurat> \
  --outputdir <output_folder>
```

## Common Options

```
--integration HarmonyIntegration
--dims 1:20
--cluster_resolution 0.8
--grid_n 200
--odds_cutoff 2
--min_cells 10
```

---

# How STEVE Decides "Unknown"

For each query cell, STEVE computes posterior probabilities:

```
P(type | x, y)
```

Derived from:

- KDE likelihoods
- Reference priors

It then computes posterior odds:

```
odds = (highest posterior) / (second-highest posterior)
```

Decision rule:

- If `odds ≥ odds_cutoff` (default = 2) → assign top type  
- If `odds < odds_cutoff` → label as **unknown**

This follows the logic described in the manuscript Methods.

---

# Module 1 — Subsampling Evaluation

## Goal
Quantify robustness when reference size varies.

### Concept

1. Split dataset into unequal partitions (10/90, 20/80, …, 90/10)
2. Use smaller partition as reference
3. Annotate larger partition (query)
4. Compute sensitivity/specificity vs original labels

### Run

```bash
Rscript steve_master.R --module subsampling \
  --userdata data/pbmc.h5seurat \
  --reference_annotation groundtruth \
  --ratios 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 \
  --n_repeats 10 \
  --outputdir output/subsampling_pbmc
```

### Outputs

- `summary_macro_metrics.csv`
- `macro_metrics_plot.pdf`

Per-run folders:
- `integrated_umap.pdf`
- `predictions_query.csv`
- `per_celltype_metrics.csv`
- `confusion_heatmap.pdf`
- `confusion_heatmap.csv`

### Interpretation

Performance drops at small reference sizes suggest instability in:

- Embedding
- Density estimation
- Partition sensitivity

---

# Module 2 — Novel Cell Evaluation

## Goal
Test whether truly absent cell types are correctly labeled **unknown**.

### Concept

1. Split dataset 50/50
2. Remove one cell type from reference
3. Annotate query
4. Ask:

> Are omitted-type cells labeled unknown  
> or incorrectly mapped to another type?

### Run

```bash
Rscript steve_master.R --module novel \
  --userdata data/stewart_bcells.h5seurat \
  --reference_annotation groundtruth \
  --n_repeats 5 \
  --outputdir output/novel_stewart
```

### Outputs

- `novel_summary.csv`

Per omitted-type folders:
- `predictions_query.csv`
- `novel_misclassification_counts.csv`

### Key Metric

```
Novel sensitivity = fraction of omitted-type cells predicted as unknown
```

### Interpretation

- High novel sensitivity → embedding separates that type well  
- Low novel sensitivity → type is absorbed by nearby clusters  

Misclassification counts reveal closest competing types.

---

# Module 3 — Annotation Benchmarking

## Goal
Compare a candidate annotation method to ground truth.

### Requirements

Your `.h5seurat` must contain:

- Ground truth column (default `groundtruth`)
- Candidate annotation column (default `user_annotation`)

Example benchmarking ScType:

```bash
Rscript steve_master.R --module benchmark \
  --userdata data/pbmc_with_annotations.h5seurat \
  --reference_annotation groundtruth \
  --user_annotation ScType_labels \
  --n_repeats 5 \
  --outputdir output/benchmark_sctype
```

### Outputs

- `benchmark_summary.csv`
- `rep_XX_tool_vs_truth_per_celltype.csv`
- `rep_XX_tool_confusion_heatmap.pdf`
- `rep_XX_tool_confusion_heatmap.csv`

Also computes:

```
rep_XX_steve_vs_truth_per_celltype.csv
```

This baseline helps distinguish:

- Embedding/mapping limitations
- Candidate annotation method limitations

---

# Module 4 — Reference Transfer Annotation

## Goal
Annotate a new dataset using a separate labeled reference dataset.

### Concept

1. Integrate labeled reference + unlabeled user dataset
2. Estimate KDE densities from reference
3. Transfer labels via posterior probabilities

### Run

```bash
Rscript steve_master.R --module transfer \
  --userdata data/user_dataset.h5seurat \
  --groundtruthdata data/reference_dataset.h5seurat \
  --reference_annotation groundtruth \
  --outputdir output/transfer_run1
```

### Outputs

- `transfer_predictions_user.csv`
- `integrated_umap.pdf`

---

## Optional Evaluation (If User Dataset Has Truth Labels)

```bash
Rscript steve_master.R --module transfer \
  --userdata data/user_dataset.h5seurat \
  --groundtruthdata data/reference_dataset.h5seurat \
  --reference_annotation groundtruth \
  --user_truth_annotation user_truth \
  --outputdir output/transfer_eval
```

Additional outputs:

- `transfer_metrics_user_vs_truth.csv`
- `transfer_confusion_heatmap.pdf`
- `transfer_confusion_heatmap.csv`

---

# Best Practices 

- Reference label quality matters (prefer FACS or strong expert curation)
- Subsampling strong + Novel weak → pipeline absorbs unseen types
- Batch effects and dataset complexity reduce performance 
- Adjust `--odds_cutoff` to control conservativeness:

  - Higher cutoff → more unknown calls (conservative)
  - Lower cutoff → more assignments (aggressive)

---

# Conceptual Positioning

STEVE does not declare a universal “best” tool.

It measures:

- Robustness  
- Calibration  
- Uncertainty  
- Context dependence  

Within a complete single-cell analysis workflow.
