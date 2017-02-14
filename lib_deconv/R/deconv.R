#' Deconvolve expressions. 
#' @param Y The expression matrix. Data frame or matrix. Each row contains expression measurements for a particular sample. 
#' @param pure_samples The pure sample indicies. List of vectors. The i-th element of list is a vector of indicies (rows of Y) that are pure samples of type i.
#' @param data_type A string indicating the data type. Use to choose pre-estimated gamma value. Current support for probe-level microarray as ``microarray-probe'', gene-level microarray as ``microarray-gene'' or rna-seq as ``rna-seq''.
#' @param n_choose How many markers to use. Can either be a single integer or a vector of integers, one for each cell type. If a single integer then all cell types use that number of markers. If not supplied number of markers will be adaptively chosen using \code{\link{est_n}}.
#' @param gamma Expression sensitivity parameter. If provided as a single positive number then that value will be used for gamma and over-ride the ``data_type`` option.  
#' @param full_markers User supplied markers. List of vectors. List should be same length as \code{pure_samples}. Each element of list is vector of indicies (columns of Y) that will be considered markers of that particular type. If not supplied then list of markers will be determined by \code{\link{find_markers}}.  
#' @param marker_method The method used to determine which measurements are markers.
#' @return List. Matrix of estimated mixing proportions. One row for each sample (row) in \code{Y}. One column for each cell type. Markers of each cell type used. Number of markers used for each cell type. Gamma parameter selected.  
#' @seealso \code{\link{find_markers}}
#' @export
deconv <- function(Y, pure_samples, n_choose, data_type = NULL, gamma = NULL, full_markers = NULL, 
    marker_method = "eta") {
    
    stopifnot(all(n_choose > 0))
    stopifnot(!is.null(c(data_type, gamma)))
    
    Y <- as.matrix(Y)
    K <- length(pure_samples)
    
    if (length(n_choose) == 1) 
        n_choose <- rep(n_choose, K)
    
    if (is.null(gamma)) 
        gamma <- get_gamma(data_type)
    
    if (is.null(full_markers)) {
        full_markers <- find_markers(Y, pure_samples, data_type = data_type, gamma = gamma, 
            method = marker_method)$L
    }
    markers <- lapply(1:K, function(i) {
        full_markers[[i]][1:n_choose[i]]
    })
    
    baseline <- baseline_exprs(Y, pure_samples, markers)
    phats <- est_phats(Y, markers, baseline, gamma)
    
    return(list(estimates = phats, markers = markers, n_choose = n_choose, gamma = gamma))
}
