# A Novel Anti-Influenza Combined Therapy Assessed by Single Cell RNA-sequencing

This repository contains the R code (Jupyter notebooks) producing the analysis and
the figures in the manuscript. The starting point of the analysis in terms of data
are the filtered gene count matrices produced by Cell Ranger - they are stored in
another repository and are imported as a git submodule.

## Notebooks

### 00-ingest-and-annotate

Ingests the matrix files merging them into a single Seurat object, filters the cells
and genes, and annotates the cells with the following metadata:
  - Experimental condition (corresponding to the sample id)
  - t-SNE projection used for visualization
  - Viral load (% of UMIs) and viral load class (Zero, Noise, Low, Medium, High)
  - Global clusters and identified cell type (pre-computed)
  - Cell type-speciic subcluster (pre-computed)

The resulting Seurat object is saved under `out/cache/annotated.Robj` and the metadata
is additionally saved as `out/cache/annotated.metadata.Robj`.

### 01-Figures-2-and-S5

Produces the Figure 2 (clustering results and cell identification) and Figure S5 (
expression of hallmark genes across global clusters).

### 03-Figure-3

Produces the Figure 3 (viral load across conditions and cell types).

### Figures-4-and-S8

Produces the Figure 4 (proportions of cell type-specific subclusters across
conditions) and Figure S8 (heatmaps illustrating these subclusters).

### 05-host-factors-Figure-S6

Finds the host factors differentially expressed in secretory vs ciliated cells
and produces Figure S6 visualizing expression of selected host factors.

### 30-clusters-subclusters-proportions-tests

Chi-Squared contigency tests assessing the differences in the proportions of
clusters / subclusters across conditions.

### 31-viral-load-classes-proportion-tests

Chi-Squared contigency tests assessing the differences in the proportions of
clusters / subclusters across conditions.

### 32-viral-load-regressions

Analysis of viral load across cell types and conditions using Beta-Binomial
regression model as well as tests for significance of viral load differences
across cell types or conditions. Also contains alternative models that have
been considered and basic model selection metrics.

### 98-global-clustering

Reproduces the global clustering used in the paper
(sensitive to versions of R packages).

### 99-subclustering

Reproduces the subclustering within cell-types
(sensitive to versions of R packages).
