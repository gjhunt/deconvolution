#' Deconvolve expressions. 
#' @param Y The expression matrix. 
#' @param pure_samples The pure sample indicies.
#' @param n_choose How many markers to use. 
#' @return Matrix of estimated mixing proportions.
#' @export
deconv <- function(Y, pure_samples, n_choose) {
    Y <- as.matrix(Y)
    full_markers <- find_markers(Y, pure_samples)
    markers <- lapply(full_markers, function(x) {
        x[1:n_choose]
    })
    
    gamma <- est_gamma()
    thetas <- est_thetas(Y, pure_samples, markers)
    phats <- est_phats(Y, markers, thetas, gamma)
    return(list(estimates = phats, markers = markers))
}

