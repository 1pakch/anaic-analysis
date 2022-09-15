# A Novel Anti-Influenza Combined Therapy Assessed by Single Cell RNA-sequencing

This repository contains the R code (Jupyter notebooks) producing the analysis and
the figures in the manuscript. The starting point of the analysis in terms of data
are either raw or the filtered gene count matrices produced by Cell Ranger. 

## Usage

### Obtaining matrices

If you are cloning from a git repository either use `git clone --recurse-submodules`
command or run `git submodule init` followed by `git submodule update` commands in the
repository directory after cloning.

Alternatively, one can download either filtered matrices manually from
[anaic-filtered-matrices](https://github.com/ilya-kolpakov/anaic-filtered-matrices)
or raw matrices from GEO submission
[GSE191176](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE191176)
and then adjust their expected paths and filenames in the function `matrix_path()`
of [00-ingest-and-annotate](notebooks/ingest-and-annotate.ipynb) notebook.

When the repository is cloned one needs additionally edit the notebooks path in
`notebooks/config.R`.

### R and R packages

The code performing data ingestion, plotting and the basic statistical analysis is
not sensitive to R versions and has been run with R versions 4.0.3 and 4.2.1 having
Seurat versions 3.2.3 and 4.0.5 respectively. The clustering results in the notebooks

 - [98-global-clustering.ipynb](notebooks/98-global-clustering.ipynb)
 - [99-subclustering.ipynb](notebooks/99-subclustering.ipynb)

were observed to produce slightly different results depending on the exact version of
R and R versions (e.g. t-SNE projections might be reflected across the axes and the
resulting coordinates might be not byte-for-byte identical). The minor dependence of
the results on R package versions might be also true for the regression analysis in
[32-viral-load-regressions.ipynb](notebooks/32-viral-load-regressions.ipynb)

The singularity containers with R environments used for the analysis are available
from the authors upon request.

## Notebooks

### 00-ingest-and-annotate

Ingests the matrix files merging them into a single Seurat object, filters the cells
and genes, and annotates the cells with the following metadata:
  - Experimental condition (corresponding to the sample id)
  - t-SNE projection used for visualization
  - Viral load (% of UMIs) and viral load class (Zero, Noise, Low, Medium, High)
  - Global clusters and identified cell type (pre-computed)
  - Cell type-specific subcluster (pre-computed)

The resulting Seurat object is saved under `out/cache/annotated.Robj` and the metadata
is additionally saved as `out/cache/annotated.metadata.Robj`.

### 01-Figures-2-and-S5

Produces the Figure 2 (clustering results and cell identification) and Figure S5 (
expression of hallmark genes across global clusters).

### 03-Figure-3

Produces the Figure 3 (viral load across conditions and cell types).

### 04-Figures-4-and-S8

Produces the Figure 4 (proportions of cell type-specific subclusters across
conditions) and Figure S8 (heatmaps illustrating these subclusters).

### 05-host-factors-Figure-S6

Finds the host factors differentially expressed in secretory vs ciliated cells
and produces Figure S6 visualizing expression of selected host factors.

### 30-clusters-subclusters-proportions-tests

Chi-Squared contingency tests assessing the differences in the proportions of
clusters / subclusters across conditions.

### 31-viral-load-classes-proportion-tests

Chi-Squared contingency tests assessing the differences in the proportions of
clusters / subclusters across conditions.

### 32-viral-load-regressions

Analysis of viral load across cell types and conditions using Beta-Binomial
regression model as well as tests for significance of viral load differences
across cell types or conditions. Also contains alternative models that have
been considered and basic model selection metrics.
The byte-for-byte reproduction might require using exactly the same versions
of R packages - see the `sessionInfo()` output within the notebook.

### 98-global-clustering

Reproduces the global clustering used in the paper.
The byte-for-byte reproduction likely requires using exactly the same versions
of R packages - see the `sessionInfo()` output within the notebook.

### 99-subclustering

Reproduces the subclustering within cell-types
The byte-for-byte reproduction likely requires using exactly the same versions
of R packages - see the `sessionInfo()` output within the notebook.
