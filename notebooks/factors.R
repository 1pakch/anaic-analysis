library(tidyverse)

make_factor <- function(values, levels=NULL, ordered=TRUE) (
    factor(
        values,
        levels = if (is.null(levels)) values,
        ordered = ordered
    )
)

## Conditions linked to library/sample ids

condition_levels <- c(
    'Mock steady state',
    'Infected untreated',
    'Infected + 6\u2032SLN-CD',
    'Infected + IFN \u03BB1',
    'Infected + IFN \u03BB1 + 6\u2032SLN-CD',
    'Mock + IFN \u03BB1 + 6\u2032SLN-CD'
)

condition_text_levels <- c(
     'Mock steady state',
     'Infected untreated',
     'Infected cyclodextrin',
     'Infected interferon lambda',
     'Infected combined treatment',
     'Mock combined treatment'
)

condition_short_levels <- c(
     'M',
     'V',
     'VCy',
     'VIf',
     'VCyIf',
     'MCyIf'
)

conditions_meta <- tibble(
    library_id = c(5, 1, 2, 3, 4, 6),
    condition = make_factor(condition_levels),
    condition_text = make_factor(condition_text_levels),
    condition_short = make_factor(condition_short_levels),
)


## Viral load classes

vlc_levels <- c('Zero', 'Noise', 'Low', 'Medium', 'High')

vlc_colors <- c('#BBBBBB', '#d8c08b', '#FECC5C', '#FD8D3C', '#E31A1C')

vlc_meta = tibble(
    viral_load_class = vlc_levels,
    color_fill = vlc_colors,
    color_text = c(rep('black', 4), rep('white', 1))
)


## Global clusters and cell types

global_clusters_levels <- c('Basal', 'BdiS', 'Secretory', 'Undefined', 'Ciliated')
cell_type_levels <- c('Basal', 'BdiS', 'Secretory', 'Ciliated')
.cell_type_levels2 <- c('Basal', 'BdiS', 'Secretory', NA, 'Ciliated')

global_clusters_meta <- tibble(
    global_clusters = make_factor(global_clusters_levels),
    seurat_clusters = factor(c(0, 3, 2, 4, 1)), # original cluster labels from Seurat
    cell_type = .cell_type_levels2 %>% factor(cell_type_levels, ordered=TRUE)
)


## Sublusters metadata

make_subcluster_meta <- function(df) (
    df
    %>% arrange(order)
    %>% mutate(subcluster = labels %>% make_factor())
    %>% select(subcluster, color)
)


subcluster_metas <- list(
    Basal = tibble(
        order = c(1, 2, 5, 3, 4, 6),
        color = c('#00ccff', '#6699cc', '#99ff99', '#9999ff', '#66cc00', '#ffcccc'),
        labels = c(
            'Steady state', 'Differentiating', 'Highly proliferating',
            'Diff. KRT14/KRT16 high', 'Low proliferating', 'Inflamed'
        )
    ),
    Ciliated = tibble(
        order = c(1, 2, 3, 4),
        color = c('#00ccff', '#6699cc', '#ffcccc', '#ff0033'),
        labels = c(
            'Steady state', 'Immature', 'Inflamed IFN-\U03BB\U2212', 'Inflamed IFN-\U03BB+'
        )
    ),   
    BdiS = tibble(
        order = c(1, 2, 4, 3, 6, 5),
        color = c('#00ccff', '#6699cc', '#99ff99', '#66cc00', '#ffcccc', '#9999ff'),
        labels = c(
            'Steady state', 'SERPINB3 high', 'KRT24 high', 'CXCL1 high', 'Inflamed', 'KRT23 high'
        )
    ),
    Secretory = tibble(
        order = c(2, 1, 3),
        color = c('#6699cc', '#00ccff', '#ffcccc'),
        labels = c(
            'SERBINB3/PI3 high', 'Steady state', 'Inflamed'
        )
    )
) %>% map(make_subcluster_meta)
