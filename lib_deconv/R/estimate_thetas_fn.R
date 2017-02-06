#' Estimate the offset terms.
#' @param Y Matrix of expressions to deconvolve. Rows are observations. 
#' @param markers List of vectors. One list entry for each type. Each vector is columns of Y that are markers for respective type. 
#' @param pure_samples List of vectors. List element for each type. Each vector contains pure samples (rows of Y) of respective type. 
#' @return List of vectors. Each vector is estimated estimated offsets of markers for each group, resp.
#' @export
est_thetas <- function(Y, pure_samples, markers) {
    K <- length(pure_samples)
    
    thetas <- list()
    for (i in 1:K) {
        compute_theta <- function(x) {
            colMeans(Y[pure_samples[[i]], x, drop = FALSE])
        }
        thetas[[i]] <- base::sapply(markers[[i]], compute_theta)
    }
    return(thetas)
}
