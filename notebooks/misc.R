library(purrr)

load_object <- function(path) {
    loaded = new.env()
    load(path, env=loaded)
    stopifnot(length(loaded) == 1)
    loaded[[names(loaded)[[1]]]]
}

get_metadata <- function(obj) (
    obj@meta.data %>% as_tibble(rownames = 'barcode')
)
    
restore_subcluster <- function(obj, cell_type) UseMethod("restore_subcluster")
restore_subcluster.tbl <- function(m) {
    cell_type <- m %>% pull(cell_type) %>% unique() %>% as.character()
    stopifnot(1 == length(cell_type))
    colname = paste(cell_type, 'subcluster', sep='_')
    m$subcluster <- m[[colname]]
    return(m)
}
restore_subcluster.Seurat <- function(obj) {
    m <- get_metadata(obj) %>% restore_subcluster()
    AddMetaData(obj, m$subcluster, 'subcluster')
}