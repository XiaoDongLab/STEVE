# STEVE Tutorial: Download Public Datasets + Prepare Inputs (Modules 1–4)

This tutorial walks through downloading **two publicly available datasets** and preparing them for use with **STEVE**. By the end, you’ll have:

* `combined_5subsets.seurat.rds`
  *(E-MTAB-9544: 5 B-cell subsets merged; used for Modules 1–2)*
* `Blood_TSP1_30.seurat.rds`
  *(Tabula Sapiens blood reference; used for SingleR + transfer demos)*
* `pbmc3k.SingleR_annotated.seurat.rds`
  *(PBMC3k annotated with SingleR using Tabula; used for Module 3)*

---

## Recommended folder layout (optional but helpful)

You can organize files however you want, but a simple layout is:

```
tutorial/
  data/
    emtab9544/
      dn_filtered_feature_bc_matrix.h5
      classical_filtered_feature_bc_matrix.h5
      naive_filtered_feature_bc_matrix.h5
      igmmem_filtered_feature_bc_matrix.h5
      trans_filtered_feature_bc_matrix.h5
    tabula/
      Blood_TSP1_30_version2d_10X_smartseq_scvi_Nov122024.h5ad
      Blood_TSP1_30_obs.csv
      Blood_TSP1_30_10x_raw_counts/     # created by export script
    pbmc3k/
      pbmc3k_filtered_gene_bc_matrices/ # downloaded from 10x
      pbmc3k.seurat.rds                 # if you already have it
  scripts/
    export_h5ad_to_10x.py
    master_STEVE_script.R
```

Throughout this tutorial, replace paths like `data/...` or `scripts/...` with wherever you actually put the files.

---

## 1) Download the datasets

### Dataset A: E-MTAB-9544 (ArrayExpress / BioStudies)

Go to:
[https://www.ebi.ac.uk/biostudies/ArrayExpress/studies/E-MTAB-9544?query=E-MTAB-9544](https://www.ebi.ac.uk/biostudies/ArrayExpress/studies/E-MTAB-9544?query=E-MTAB-9544)

Download these **five** files:

* `dn_filtered_feature_bc_matrix.h5`
* `classical_filtered_feature_bc_matrix.h5`
* `naive_filtered_feature_bc_matrix.h5`
* `igmmem_filtered_feature_bc_matrix.h5`
* `trans_filtered_feature_bc_matrix.h5`

---

### Dataset B: Tabula Sapiens v2 (Figshare)

Go to:
[https://figshare.com/articles/dataset/Tabula_Sapiens_v2/27921984](https://figshare.com/articles/dataset/Tabula_Sapiens_v2/27921984)

Download:

* `Blood_TSP1_30_version2d_10X_smartseq_scvi_Nov122024.h5ad`

You will also need the metadata CSV used in the code below:

* `Blood_TSP1_30_obs.csv`

*(If you already have this from your workflow/export, place it next to the `.h5ad`.)*

---

## 2) Convert Tabula Sapiens `.h5ad` → 10X-style counts (for Seurat)

Tabula is provided as an **AnnData/Scanpy** object (`.h5ad`). The script below exports a **10X-style matrix** so Seurat can load it cleanly.

### Install/upgrade required Python packages

```bash
python3 -m pip install --upgrade anndata scipy pandas numpy
```

### Run the export script

```bash
python3 scripts/export_h5ad_to_10x.py \
  --in "data/tabula/Blood_TSP1_30_version2d_10X_smartseq_scvi_Nov122024.h5ad" \
  --outdir "data/tabula/Blood_TSP1_30_10x_raw_counts" \
  --layer raw_counts \
  --meta_csv "data/tabula/Blood_TSP1_30_obs.csv"
```

**What this does:**

* Reads the `.h5ad`
* Exports **raw counts** into a 10X-like folder (matrix + genes/features + barcodes)
* Writes/uses a metadata CSV that will become Seurat `meta.data`

---

## 3) Load datasets into Seurat + save as `.seurat.rds`

Run the following **R code** to:

1. Load each E-MTAB-9544 `.h5` subset into Seurat
2. Merge them into one combined object
3. Load the exported Tabula 10X matrix + metadata
4. Add `groundtruth` labels
5. Save both as `*.seurat.rds`

> **Important:** This assumes you’re using Seurat v5 objects (layers).

```r
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)  # SaveSeuratRds / LoadSeuratRds
})

# Be explicit: build v5 assays
options(Seurat.object.assay.version = "v5")

# -----------------------
# 1) File paths
# -----------------------
h5_paths <- c(
  dn        = "data/emtab9544/dn_filtered_feature_bc_matrix.h5",
  classical = "data/emtab9544/classical_filtered_feature_bc_matrix.h5",
  naive     = "data/emtab9544/naive_filtered_feature_bc_matrix.h5",
  igmmem    = "data/emtab9544/igmmem_filtered_feature_bc_matrix.h5",
  trans     = "data/emtab9544/trans_filtered_feature_bc_matrix.h5"
)

tabula_10x_dir  <- "data/tabula/Blood_TSP1_30_10x_raw_counts"
tabula_meta_csv <- "data/tabula/Blood_TSP1_30_obs.csv"

combined_out <- "data/emtab9544/combined_5subsets.seurat.rds"
tabula_out   <- "data/tabula/Blood_TSP1_30.seurat.rds"

# -----------------------
# 2) Read each 10X .h5 -> Seurat object
# -----------------------
seurat_list_h5 <- lapply(names(h5_paths), function(nm) {
  mtx <- Read10X_h5(h5_paths[[nm]])
  obj <- CreateSeuratObject(counts = mtx, project = nm)  # defaults only
  obj$sample_id <- nm
  obj
})
names(seurat_list_h5) <- names(h5_paths)

# -----------------------
# 3) Combine all .h5 objects into one Seurat object
# -----------------------
combined_h5 <- merge(
  x = seurat_list_h5[[1]],
  y = seurat_list_h5[-1],
  add.cell.ids = names(seurat_list_h5),
  project = "E-MTAB-9544_combined"
)

combined_h5 <- JoinLayers(combined_h5, assay = "RNA", layers = "counts")

cat("combined_h5 RNA layers after JoinLayers:\n")
print(Layers(combined_h5[["RNA"]]))  # should print just "counts"

# -----------------------
# 4) Tabula Sapiens from exported 10X MTX + metadata CSV
# -----------------------
mtx <- Read10X(tabula_10x_dir)

meta <- read.csv(tabula_meta_csv, row.names = 1)
meta <- meta[colnames(mtx), , drop = FALSE]  # align metadata to barcodes

tabula <- CreateSeuratObject(counts = mtx, meta.data = meta, project = "Blood_TSP1_30")

cat("tabula RNA layers:\n")
print(Layers(tabula[["RNA"]]))  # should be "counts"

# -----------------------
# 5) Groundtruth labels
# -----------------------
tabula$groundtruth <- tabula$cell_ontology_class
combined_h5$groundtruth <- combined_h5$sample_id

# -----------------------
# 6) Save as Seurat RDS
# -----------------------
if (file.exists(combined_out)) file.remove(combined_out)
if (file.exists(tabula_out))   file.remove(tabula_out)

SaveSeuratRds(combined_h5, file = combined_out)
SaveSeuratRds(tabula,      file = tabula_out)

cat("Saved:\n", combined_out, "\n", tabula_out, "\n")
```

---

## 4) Run STEVE Module 1: Subsampling Evaluation

Example command:

```bash
Rscript scripts/master_STEVE_script.R --module subsampling \
  --userdata data/emtab9544/combined_5subsets.seurat.rds \
  --reference_annotation groundtruth \
  --ratios 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 \
  --n_repeats 1 \
  --outputdir output/subsampling_pbmc
```

### Example output

```text
Warning messages:
1: package ‘Seurat’ was built under R version 4.4.3 
2: package ‘SeuratObject’ was built under R version 4.4.3 
3: package ‘sp’ was built under R version 4.4.3 
4: package ‘sctransform’ was built under R version 4.4.3 
5: package ‘ggplot2’ was built under R version 4.4.3 
6: package ‘Rcpp’ was built under R version 4.4.3 
7: package ‘circlize’ was built under R version 4.4.3 
8: package ‘tidyr’ was built under R version 4.4.3 
[2026-03-01 20:43:54]  STEVE master starting. module=subsampling output=output/subsampling_pbmc 
Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
This message will be shown once per session
[2026-03-01 20:44:27]  Subsampling done: ratio_10_iter_01 (macro sens=0.699 spec=0.965) 
[2026-03-01 20:44:57]  Subsampling done: ratio_20_iter_01 (macro sens=0.716 spec=0.964) 
[2026-03-01 20:45:27]  Subsampling done: ratio_30_iter_01 (macro sens=0.677 spec=0.972) 
[2026-03-01 20:45:55]  Subsampling done: ratio_40_iter_01 (macro sens=0.712 spec=0.963) 
[2026-03-01 20:46:24]  Subsampling done: ratio_50_iter_01 (macro sens=0.729 spec=0.963) 
[2026-03-01 20:46:53]  Subsampling done: ratio_60_iter_01 (macro sens=0.725 spec=0.963) 
[2026-03-01 20:47:23]  Subsampling done: ratio_70_iter_01 (macro sens=0.727 spec=0.966) 
[2026-03-01 20:47:52]  Subsampling done: ratio_80_iter_01 (macro sens=0.740 spec=0.965) 
[2026-03-01 20:48:23]  Subsampling done: ratio_90_iter_01 (macro sens=0.742 spec=0.965) 
  ratio iteration macro_sensitivity macro_specificity
1   0.1         1         0.6987711         0.9648993
2   0.2         1         0.7155663         0.9641254
3   0.3         1         0.6774562         0.9720962
4   0.4         1         0.7123739         0.9634414
5   0.5         1         0.7285312         0.9631283
6   0.6         1         0.7246662         0.9626505
7   0.7         1         0.7269893         0.9655262
8   0.8         1         0.7397868         0.9650307
9   0.9         1         0.7421571         0.9647420
Warning message:
package ‘future’ was built under R version 4.4.3 
[2026-03-01 20:48:23]  OK!
```

---

## 5) Run STEVE Module 2: Novel Cell Evaluation

Example command:

```bash
Rscript scripts/master_STEVE_script.R --module novel \
  --userdata data/emtab9544/combined_5subsets.seurat.rds \
  --reference_annotation groundtruth \
  --n_repeats 1 \
  --outputdir output/novel_stewart
```

### Example output

```text
Warning messages:
1: package ‘Seurat’ was built under R version 4.4.3 
2: package ‘SeuratObject’ was built under R version 4.4.3 
3: package ‘sp’ was built under R version 4.4.3 
4: package ‘sctransform’ was built under R version 4.4.3 
5: package ‘ggplot2’ was built under R version 4.4.3 
6: package ‘Rcpp’ was built under R version 4.4.3 
7: package ‘circlize’ was built under R version 4.4.3 
8: package ‘tidyr’ was built under R version 4.4.3 
[2026-03-01 20:58:32]  STEVE master starting. module=novel output=output/novel_stewart 
Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
This message will be shown once per session
[2026-03-01 20:59:06]  Novel eval iter=1 omit=classical (sens=0.142 spec=0.928) 
[2026-03-01 20:59:33]  Novel eval iter=1 omit=dn (sens=0.299 spec=0.930) 
[2026-03-01 21:00:00]  Novel eval iter=1 omit=igmmem (sens=0.677 spec=0.841) 
[2026-03-01 21:00:28]  Novel eval iter=1 omit=naive (sens=0.130 spec=0.856) 
[2026-03-01 21:00:58]  Novel eval iter=1 omit=trans (sens=0.425 spec=0.833) 
  iteration omitted_type novel_sensitivity novel_specificity
1         1    classical         0.1421569         0.9283820
2         1           dn         0.2990256         0.9299291
3         1       igmmem         0.6766091         0.8409967
4         1        naive         0.1302578         0.8559819
5         1        trans         0.4251012         0.8332977
Warning message:
package ‘future’ was built under R version 4.4.3 
[2026-03-01 21:00:59]  OK!
```

As can be seen, it is difficult to transcriptionally separate these subtypes of B cells.

---

## 6) STEVE Module 3: Benchmarking (SingleR example)

In this module, we use any already-existing annotation tool (e.g., **SingleR**) to annotate a dataset and then compare the tool’s annotations against ground truth.

SingleR requires a reference dataset, so here we annotate **PBMC3k** using **Tabula Sapiens Blood** as the reference (we are **not** using Stewart here because it only contains B cells).

### Download PBMC3k (10x)

Download from:
[https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz](https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz)

Unpack it into `data/pbmc3k/` (or wherever you prefer).

---

### SingleR annotation script (PBMC3k using Tabula reference)

> This script produces: `pbmc3k.SingleR_annotated.seurat.rds`

```r
#!/usr/bin/env Rscript

# ============================================================
# SingleR annotation of PBMC3k (test) using Tabula Sapiens (reference)
# Memory-safe: pseudo-bulk Tabula by label + top genes + correct mapping back to Seurat
#
# IMPORTANT (for benchmarking later):
#   - groundtruth      = PBMC reference labels (cluster-derived, mapped to Tabula names)
#   - SingleR_label    = user_annotation (what you are evaluating)
#
# Inputs:
#   pbmc3k.seurat.rds        (test)
#   Blood_TSP1_30.seurat.rds (reference)
#
# Output:
#   pbmc3k.SingleR_annotated.seurat.rds
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)  # SaveSeuratRds / LoadSeuratRds
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(Matrix)
  library(scuttle)
  library(S4Vectors)
  library(SingleR)
})

# -----------------------
# Paths
# -----------------------
pbmc_in   <- "data/pbmc3k/pbmc3k.seurat.rds"
tabula_in <- "data/tabula/Blood_TSP1_30.seurat.rds"
pbmc_out  <- "data/pbmc3k/pbmc3k.SingleR_annotated.seurat.rds"

# -----------------------
# Load Seurat objects
# -----------------------
if (!file.exists(pbmc_in)) stop("PBMC RDS not found: ", pbmc_in)
if (!file.exists(tabula_in)) stop("Tabula RDS not found: ", tabula_in)

pbmc   <- LoadSeuratRds(pbmc_in)
tabula <- LoadSeuratRds(tabula_in)

DefaultAssay(pbmc)   <- "RNA"
DefaultAssay(tabula) <- "RNA"

# -----------------------
# PBMC reference labels (cluster identity) -> store as pbmc_cluster
# -----------------------
pbmc$pbmc_cluster <- as.character(Idents(pbmc))

cat("PBMC cluster identity breakdown (raw):\n")
print(sort(table(pbmc$pbmc_cluster), decreasing = TRUE))

# -----------------------
# Helpers: Seurat -> SCE (counts only, v5-safe)
# -----------------------
get_counts_layer <- function(seu, assay = "RNA") {
  lay <- tryCatch(Layers(seu[[assay]]), error = function(e) character(0))
  if (length(lay) == 0) return("counts") # fallback for older objects
  if ("counts" %in% lay) return("counts")
  cl <- lay[grepl("^counts", lay)]
  if (length(cl) > 0) return(cl[1])
  stop("No counts-like layer found. Layers: ", paste(lay, collapse = ", "))
}

seurat_to_sce_counts <- function(seu, assay = "RNA") {
  layer_use <- get_counts_layer(seu, assay = assay)
  counts <- GetAssayData(seu, assay = assay, layer = layer_use)

  if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "dgCMatrix")
  if (is.null(rownames(counts))) rownames(counts) <- rownames(seu[[assay]])
  if (is.null(colnames(counts))) colnames(counts) <- colnames(seu)

  SingleCellExperiment(
    assays  = S4Vectors::SimpleList(counts = counts),
    colData = S4Vectors::DataFrame(seu@meta.data)
  )
}

# -----------------------
# Build SCEs
# -----------------------
sce_test <- seurat_to_sce_counts(pbmc,   assay = "RNA")
sce_ref  <- seurat_to_sce_counts(tabula, assay = "RNA")

# -----------------------
# Choose Tabula label column for SingleR
# -----------------------
ref_label_col <- "cell_ontology_class"   # could swap to a broader label if you have one
if (!ref_label_col %in% colnames(colData(sce_ref))) {
  stop("Tabula reference missing '", ref_label_col, "'. Available columns:\n",
       paste(colnames(colData(sce_ref)), collapse = ", "))
}
ref_labels <- as.character(colData(sce_ref)[[ref_label_col]])
cat("Tabula unique labels:", length(unique(ref_labels)), "\n")

# -----------------------------
# 1) Pseudobulk reference by label (MEMORY FIX)
# -----------------------------
ref_pb_se <- scuttle::sumCountsAcrossCells(
  sce_ref,
  ids = ref_labels,
  exprs_values = "counts"
)

sum_assay <- SummarizedExperiment::assayNames(ref_pb_se)[1]  # usually "sum"
pb_sum <- SummarizedExperiment::assay(ref_pb_se, sum_assay)

labs <- make.unique(as.character(SummarizedExperiment::colData(ref_pb_se)$ids))
colnames(pb_sum) <- labs
rownames(pb_sum) <- rownames(ref_pb_se)

ncells <- as.numeric(SummarizedExperiment::colData(ref_pb_se)$ncells)
if (any(!is.finite(ncells)) || any(ncells <= 0)) stop("Bad ncells in pseudo-bulk reference")

# average-per-cell
pb_avg <- pb_sum %*% Matrix::Diagonal(x = 1 / ncells)
dimnames(pb_avg) <- list(rownames(pb_sum), colnames(pb_sum))

ref_pb <- SingleCellExperiment(
  assays  = S4Vectors::SimpleList(counts = pb_avg),
  colData = S4Vectors::DataFrame(label = labs, ncells = ncells)
)
colnames(ref_pb) <- labs

# -----------------------------
# 2) Intersect genes + optional top genes
# -----------------------------
common_genes <- intersect(rownames(sce_test), rownames(ref_pb))
if (length(common_genes) < 2000) {
  warning("Low overlap genes (", length(common_genes),
          "). Check gene naming consistency (symbols vs Ensembl).")
}

sce_test2 <- sce_test[common_genes, ]
ref_pb2   <- ref_pb[common_genes, ]

top_n_genes <- 5000  # try 3000–10000
ref_means <- Matrix::rowMeans(SummarizedExperiment::assay(ref_pb2, "counts"))
keep_genes <- names(sort(ref_means, decreasing = TRUE))[seq_len(min(top_n_genes, length(ref_means)))]

sce_test2 <- sce_test2[keep_genes, ]
ref_pb2   <- ref_pb2[keep_genes, ]

# -----------------------------
# 3) log-normalize (small now)
# -----------------------------
sce_test2 <- scuttle::logNormCounts(sce_test2)
ref_pb2   <- scuttle::logNormCounts(ref_pb2)

# -----------------------------
# 4) Run SingleR
# -----------------------------
pred <- SingleR(
  test  = sce_test2,
  ref   = ref_pb2,
  labels = colnames(ref_pb2),
  assay.type.test = "logcounts",
  assay.type.ref  = "logcounts"
)

# -----------------------------
# 5) Prune + NAME vectors correctly (so mapping works)
# -----------------------------
pruned <- pruneScores(pred)

final_labels <- pred$labels
final_labels[pruned] <- "unknown"

cell_ids <- rownames(pred)  # test cell IDs stored here
if (is.null(cell_ids) || length(cell_ids) != length(final_labels)) {
  stop("Unexpected: rownames(pred) missing/wrong length; cannot map labels back to PBMC.")
}

names(final_labels) <- cell_ids
names(pruned)       <- cell_ids

# -----------------------------
# 6) Attach SingleR back to PBMC
# -----------------------------
seu_cells <- colnames(pbmc)
overlap <- intersect(seu_cells, names(final_labels))

cat("PBMC cells:", length(seu_cells), "\n")
cat("SingleR cells:", length(final_labels), "\n")
cat("Overlap:", length(overlap), "\n")

if (length(overlap) == 0) {
  cat("Example PBMC cell IDs:\n"); print(head(seu_cells))
  cat("Example SingleR cell IDs:\n"); print(head(names(final_labels)))
  stop("No overlap between PBMC cell names and SingleR output cell IDs.")
}

pbmc$SingleR_label  <- NA_character_
pbmc$SingleR_pruned <- NA

pbmc$SingleR_label[overlap]  <- final_labels[overlap]
pbmc$SingleR_pruned[overlap] <- pruned[overlap]

cat("\nSingleR label counts (PBMC):\n")
print(head(sort(table(pbmc$SingleR_label, useNA = "ifany"), decreasing = TRUE), 25))

cat("\nPruned flags:\n")
print(table(pbmc$SingleR_pruned, useNA = "ifany"))

# ============================================================
# 7.5) Harmonize PBMC groundtruth to Tabula ontology + filter
#      Groundtruth = PBMC reference labels (cluster-derived, mapped to Tabula)
#      User annotation = SingleR_label
# ============================================================

# Tabula ontology labels we allow
tabula_types <- unique(as.character(tabula$cell_ontology_class))
tabula_types <- tabula_types[!is.na(tabula_types)]

# Groundtruth starts from PBMC cluster identity
pbmc$groundtruth <- as.character(pbmc$pbmc_cluster)

# *** REQUIRED mapping: PBMC cluster names -> Tabula names ***
map_pbmc_to_tabula <- c(
  "B"            = "b cell",
  "NK"           = "natural killer cell",
  "CD14+ Mono"   = "classical monocyte",
  "FCGR3A+ Mono" = "non-classical monocyte",
  "CD8 T"        = "cd8-positive, alpha-beta t cell",
  "Naive CD4 T"  = "naive thymus-derived cd4-positive, alpha-beta t cell",
  "Memory CD4 T" = "cd4-positive, alpha-beta t cell",
  "Platelet"     = "platelet",
  # Choose ONE (or comment out to drop DC cells from PBMC)
  "DC"           = "myeloid dendritic cell"
  # "DC"        = "plasmacytoid dendritic cell"
)

# Apply mapping; anything unmapped becomes NA (and will be dropped)
pbmc$groundtruth <- unname(map_pbmc_to_tabula[pbmc$groundtruth])
pbmc$groundtruth <- as.character(pbmc$groundtruth)

# Diagnostics: show unmapped PBMC clusters
unmapped <- setdiff(unique(pbmc$pbmc_cluster), names(map_pbmc_to_tabula))
if (length(unmapped) > 0) {
  cat("\nUnmapped PBMC cluster labels (will be dropped):\n")
  print(unmapped)
}

# Keep only PBMC cells where:
# - groundtruth is a Tabula label (after mapping)
# - groundtruth not NA
# - SingleR_label not NA
# - optionally drop SingleR unknowns (recommended for benchmarking)
keep_cells <- !is.na(pbmc$groundtruth) &
  pbmc$groundtruth %in% tabula_types &
  !is.na(pbmc$SingleR_label) &
  pbmc$SingleR_label != "unknown"

cat("\nKeeping", sum(keep_cells), "of", length(keep_cells),
    "PBMC cells after Tabula groundtruth harmonization.\n")

if (sum(keep_cells) == 0) {
  cat("\nDEBUG:\n")
  cat("Unique mapped groundtruth labels:\n")
  print(sort(unique(pbmc$groundtruth)))
  cat("Example Tabula labels:\n")
  print(head(sort(tabula_types), 30))
  stop("No cells survived filtering. Mapping likely mismatched.")
}

pbmc <- subset(pbmc, cells = colnames(pbmc)[keep_cells])

# Helpful: set identity to reference groundtruth
Idents(pbmc) <- pbmc$groundtruth

cat("\nPBMC groundtruth (Tabula-matched) breakdown:\n")
print(sort(table(pbmc$groundtruth), decreasing = TRUE))

cat("\nPBMC SingleR_label (user annotation) breakdown:\n")
print(head(sort(table(pbmc$SingleR_label, useNA = "ifany"), decreasing = TRUE), 25))

# -----------------------
# Save annotated PBMC object
# -----------------------

# Optional: keep a harmonized column name for STEVE benchmarking
pbmc$SingleR_label_harmonized <- pbmc$SingleR_label

if (file.exists(pbmc_out)) file.remove(pbmc_out)
SaveSeuratRds(pbmc, file = pbmc_out)

cat("\nSaved annotated PBMC to:\n", pbmc_out, "\n")

# Optional: show Tabula counts
cat("\nTabula cell_ontology_class breakdown:\n")
print(sort(table(tabula$cell_ontology_class), decreasing = TRUE))
```

---

### Run STEVE benchmarking (Module 3)

Now benchmark SingleR’s annotations (`SingleR_label_harmonized`) against your ground truth (`groundtruth`):

```bash
Rscript scripts/master_STEVE_script.R --module benchmark \
  --userdata data/pbmc3k/pbmc3k.SingleR_annotated.seurat.rds \
  --reference_annotation groundtruth \
  --user_annotation SingleR_label_harmonized \
  --n_repeats 1 \
  --outputdir output/benchmarking_pbmc_singleR
```

#### Example output

```text
Warning messages:
1: package ‘Seurat’ was built under R version 4.4.3 
2: package ‘SeuratObject’ was built under R version 4.4.3 
3: package ‘sp’ was built under R version 4.4.3 
4: package ‘sctransform’ was built under R version 4.4.3 
5: package ‘ggplot2’ was built under R version 4.4.3 
6: package ‘Rcpp’ was built under R version 4.4.3 
7: package ‘circlize’ was built under R version 4.4.3 
8: package ‘tidyr’ was built under R version 4.4.3 
[2026-03-02 01:31:10]  STEVE master starting. module=benchmark output=output/benchmarking_pbmc_singleR 
Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
This message will be shown once per session
[2026-03-02 01:31:23]  Benchmark iter=1 candidate(sens=0.627 spec=0.941) steve(sens=0.804 spec=0.994) 
  iteration candidate_macro_sensitivity candidate_macro_specificity
1         1                   0.6267394                   0.9412333
  steve_macro_sensitivity steve_macro_specificity
1               0.8039559               0.9938641
Warning message:
package ‘future’ was built under R version 4.4.3 
[2026-03-02 01:31:23]  OK!
```

---

## 7) STEVE Module 4: Reference Transfer Annotation

While not used for evaluating a pipeline, STEVE includes a final module (**reference transfer**) that can be used to annotate a dataset using another dataset. In this example, we use **Tabula Sapiens** to annotate **PBMC3k**.

```bash
Rscript scripts/master_STEVE_script.R --module transfer \
  --userdata data/pbmc3k/pbmc3k.seurat.rds \
  --groundtruthdata data/tabula/Blood_TSP1_30.seurat.rds \
  --reference_annotation groundtruth \
  --outputdir output/transfer_run1
```

---

If you want, paste your `export_h5ad_to_10x.py` usage header (help text) and I can add a short “Script arguments” section that reads like polished documentation.
