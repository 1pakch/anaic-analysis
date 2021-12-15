library(purrr)

source('config.R')

join_with <- function(sep) (
    function(...) paste(..., sep=sep))

ensure_dir_exists <- function(path) {
    dir.create(path, recursive=TRUE, showWarnings=FALSE)
    path
}
    
path_at <- function(...) {
    folder <- join_with('/')(...) %>% ensure_dir_exists()
    function(...) join_with('/')(folder, join_with('.')(...))
}
