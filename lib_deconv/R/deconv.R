#' Deconvolve expressions. 
#' @param Y The expression matrix. Data frame or matrix. Each row contains expression measurements for a particular sample. 
#' @param pure_samples The pure sample indicies. List of vectors. The i-th element of list is a vector of indicies (rows of Y) that are pure samples of type i. 
#' @param n_choose How many markers to use. Can either be a single integer or a vector of integers, one for each cell type. If a single integer then all cell types use that number of markers. If not supplied number of markers will be adaptively chosen using \code{\link{est_n}}.
#' @param full_markers User supplied markers. List of vectors. List should be same length as \code{pure_samples}. Each element of list is vector of indicies (columns of Y) that will be considered markers of that particular type. If not supplied then list of markers will be determined by \code{\link{find_markers}}.  
#' @param gamma Expression sensitivity parameter. If provided as a single positive number then that value will be used for gamma. Otherwise can pass as a string indicating data type. Current support for probe-level microarray as ``microarray-probe'', gene-level microarray as ``microarray-gene'' or rna-seq as ``rna-seq''.
#' @param n_test How many observations to generate for training sample used to optiimize \code{gamma} or \code{n_choose} when they are not specified. Single integer. 
#' @param n_test_reps How many iterations to use in optimization of \code{n_choose} when it is not speicified by user. Single integer. 
#' @param gamma_test_reps How many iterations to use in optimization of \code{gamma} when it is not speicified by user. Single integer.
#' @param marker_method The method used to determine which measurements are markers.
#' @return Matrix of estimated mixing proportions. One row for each sample (row) in \code{Y}. One column for each cell type.
#' @seealso \code{\link{find_markers}}, \code{\link{est_n}}, \code{\link{est_gamma}}
#' @export
deconv <- function(Y, pure_samples, gamma, n_choose = NULL, full_markers = NULL, 
    n_test = 500, n_test_reps = 2, marker_method = "eta", seed = NULL) {
    
    if (!is.null(seed)) {
        set.seed(seed)
    }
    
    Y <- as.matrix(Y)
    K <- length(pure_samples)
    
    if (length(n_choose) == 1) 
        n_choose <- rep(n_choose, K)
    
    if (is.null(full_markers)) {
        full_markers <- find_markers(Y, pure_samples, method = marker_method)$L
    }
    
    if (!is.numeric(gamma)) {
        if (gamma == "microarray-probe") {
            gamma <- 0.9392593
        } else if (gamma == "microarray-gene") {
            gamma <- 0.9902935
        } else if (gamma == "rng-seq") {
            gamma <- 1
        }
        
    }
    
    if (is.null(n_choose)) {
        td <- gen_test_data(Y, pure_samples, n_test)
        curr_n <- rpois(K, 100)
        n_choose <- est_n(td, pure_samples, full_markers, gamma, curr_n, n_test_reps)
    }
    
    markers <- lapply(1:K, function(i) {
        full_markers[[i]][1:n_choose[i]]
    })
    
    baseline <- baseline_exprs(Y, pure_samples, markers)
    phats <- est_phats(Y, markers, baseline, gamma)
    
    return(list(estimates = phats, markers = markers, n_choose = n_choose, gamma = gamma))
}
