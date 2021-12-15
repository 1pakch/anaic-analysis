library(tidyverse)
library(ggplot2)

vhead <- function(l) (if (length(l) > 0) l[[1]] else NULL)
vtail <- function(l) (if (length(l) > 1) l[2:length(l)] else NULL)

# Recursively split objects into nested lists using `split_one` at each level
split_rec <- function(input, by, split_one) (
    if (length(by) == 0)
        input
    else
        split_one(input, vhead(by)) %>% map(split_rec, vtail(by), split_one)
)

# Recursively reduce nested lists using `reduce_one` functor at each level
reduce_rec <- function(input, reduce_one, level=1) (
    if (!is_tibble(input))
        input %>% map(reduce_rec, reduce_one, level+1) %>% reduce_one(level)
    else
        input %>% mutate(level_change=0)
)

split_df <- function(input, by) (
    input %>% group_by_(by) %>% group_split())

group_separator <- function(level) (
    tibble(level_change=level))

repeat_tibble <- function(data, n) (
    list(data) %>% rep(n) %>% bind_rows())
    
repeat_if_separator <- function(data, level_spaces) (
    (data$level_change[[1]])
    %>% {if (.) repeat_tibble(data, level_spaces[[.]]) else data}
)
    
combine_at <- function(level) (
    function(x, y) bind_rows(x, group_separator(level), y)
)

reduce_flat <- function(values, level) (
    values %>% map(mutate) %>% reduce(combine_at(level))
)

fill_non_changing <- function(x) {
    l <- lag(x)
    f <- rev(lag(rev(x)))
    ifelse(
        !is.na(x),
        x,
        ifelse(f == l, f, NA))
}
    
prepare_axis_meta1 <- function(df, group_vars, level_spaces) (
    df
    %>% split_rec(group_vars, split_df)
    %>% reduce_rec(reduce_flat)
    %>% mutate(pos=1:nrow(.))
    %>% mutate_at(group_vars, fill_non_changing)
    %>% group_by_at(c(group_vars, 'level_change'))
    %>% group_split()
    %>% map(repeat_if_separator, level_spaces)
    %>% bind_rows()
    %>% arrange(pos)
    %>% mutate(pos=1:nrow(.))
)

restore_factor <- function(x) (
  if (!is.factor(x)) (function(y) y)
  else function(y) {
    l <- levels(x)
    factor(l[y], levels=l)
  }
)

relative_break_sizes <- function(fractions=c(0.005, 0.0025)) (
    function(n_items) (
        (n_items * fractions)
        %>% as.integer
    )
)

absolute_break_sizes <- function(sizes) (
    function(n_items) sizes
)

prepare_axis_meta <- function(df, group_vars, break_sizes) (
    prepare_axis_meta1(df, group_vars, break_sizes(nrow(df)))
    %>% mutate(across(all_of(group_vars), ~restore_factor(df[[cur_column()]])(.x)))
)

if.null <- function(x, t, f) ifelse(is.null(x), t, f)

groups.barplot0 <- function(axis.meta, pos_var, group_vars, level) (
    ggplot(data=axis.meta)
    + geom_raster(aes(x=get(pos_var), y=group_vars[[level]], fill=get(group_vars[[level]])))
    + scale_x_discrete(limits=axis.meta %>% pull(pos_var) %>% factor())
    + theme(
        legend.position='top',
        axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks=element_blank())
)
    
groups.barplot <- function(axis.meta, pos_var, group_vars, level, palette=NULL) (
    groups.barplot0(axis.meta, pos_var, group_vars, level)
    + scale_fill_brewer(
        type='qual',
        name=group_vars[[level]],
        breaks=levels(axis.meta[[group_vars[[level]]]]),
        palette=if.null(palette, level, palette),
        guide='legend',
        drop=FALSE)
)

groups.barplot2 <- function(axis.meta, pos_var, group_vars, level, name, palette=NULL) (
    groups.barplot0(axis.meta, pos_var, group_vars, level)
    + scale_fill_brewer(
        type='qual',
        name=name,
        breaks=levels(axis.meta[[group_vars[[level]]]]),
        palette=if.null(palette, level, palette),
        guide='legend',
        drop=FALSE)
)
    
level.barplot <- function(axis.meta, pos_var, level_var, color=NULL) (
    ggplot(data=axis.meta)
    + geom_col(
        aes(x=get(pos_var),
            y=get(level_var)),
        fill=color,
        size=0)
    + scale_x_discrete(
        limits=axis.meta %>% pull(pos_var) %>% factor())
    + scale_y_continuous(
        limits=NULL,
        breaks=(function(limits) mean(limits)),
        labels=level_var)
    + theme(
        axis.text.x=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank()) 
)

get.plot.cells.matrix <- function(data, genes, barcodes) (
    data@assays$SCT@scale.data[genes,]
    %>% as.matrix()
    %>% t()
    %>% {.[barcodes,]}
)

get.plot.cells.data1 <- function(data, genes, barcodes, exmax=2.5, exmin=-2.5) {
    present.genes <- intersect(genes, rownames(data))
    missing.genes <- setdiff(genes, rownames(data))
    m <- get.plot.cells.matrix(data, present.genes, barcodes)
    (
        m
        %>% tibble::as_tibble(rownames='barcode')
        %>% tidyr::pivot_longer(cols=all_of(present.genes), names_to='gene')
        %>% mutate(value = pmin(pmax(value, exmin), exmax))
        %>% full_join(crossing(barcode=barcodes, gene=missing.genes, value=-0.01))
    )
}

get.plot.cells.data <- function(data, genes.meta, cells.meta, exmax=2.5, exmin=-2.5) (
    get.plot.cells.data1(
        data,
        genes.meta %>% pull(gene) %>% na.exclude,
        cells.meta %>% pull(barcode) %>% na.exclude,
        exmax=exmax,
        exmin=exmin)
    %>% left_join(genes.meta %>% select(gene, pos_y), by='gene')
    %>% left_join(cells.meta %>% select(barcode, pos_x), by='barcode')
)
    
make.heatmap <- function(plot.cells.data, genes.meta, cells.meta, limits=c(-2.5, 2.5)) (
    ggplot(data=plot.cells.data) +
    geom_raster(aes(x=pos_x, y=pos_y, fill=value)) +
    scale_fill_distiller(palette='RdBu', limits=limits, na.value='white', name='Normalized\nexpression') +
    scale_x_discrete(limits=cells.meta$pos_x %>% as.factor()) +
    scale_y_discrete(limits=rev(genes.meta$pos_y), labels=rev(genes.meta$gene) %>% str_replace_na('')) +
    theme(axis.text.x=element_blank(), axis.title=element_blank(), axis.ticks.x=element_blank())
)

smaller_legend <- function(p, point_size=1.5, text_size=8, title_size=9) (
    p
    + guides(
        shape = guide_legend(override.aes = list(size = point_size)),
        color = guide_legend(override.aes = list(size = point_size)))
    + theme(
        legend.text = element_text(size=text_size),
        legend.title = element_text(size=title_size),
    )
)

heatmap_size <- function(genes.meta, cells.meta, stretch_x=1, stretch_y=1) {
    width <- nrow(cells.meta) / 100 * stretch_x
    height <- nrow(genes.meta) / 8 * stretch_y
    return(c(width, height))
}

title_height <- function() (0.5)
legend_height <- function() (0.7)
top_subplot_height <- function() (0.15)
top_subplot_spacer_height <- function() top_subplot_height()/3

as_unit <- function(x, unit) (
    unit(x, rep(unit, length(x)))
)

conditions_scale_fill <- scale_fill_brewer(type='qual', palette=6)
conditions_scale_color <- scale_color_brewer(type='qual', palette=6)
global_clusters_scale_fill <- scale_fill_brewer(type='qual', palette=2, drop=FALSE)
global_clusters_scale_color <- scale_color_brewer(type='qual', palette=2, drop=FALSE)
subcluster_scale_fill <- scale_fill_brewer(type='qual', palette=3)
subcluster_scale_color <- scale_color_brewer(type='qual', palette=3)

    
condition_colors <- c('#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', '#999999')
global_cluster_colors <- c('#1B9E77', '#D95F02', '#7570B3', '#E7298A', '#66A61E')
cell_type_colors <- c('#1B9E77', '#D95F02', '#7570B3', '#66A61E')


