#'Estimate the gene type proportions.
#' @param Y Matrix of expressions to deconvolve. Rows are observations. 
#' @param markers List of vectors. One list entry for each type. Each vector is columns of Y that are markers for respective type. 
#' @param thetas_ests List of vectors (same as markers). One list entry for each type. Vector is estimated offset for each marker of the respective type. 
#' @param gamma_est Scalar. Estimated slope. 
#' 
#' @return Estimated matrix of mixing proportions.
est_phats <- function(Y, markers, theta_ests, gamma_est) {
    
    K <- length(markers)
    
    contrib_est <- function(i) {
        Y_i <- Y[, markers[[i]], drop = FALSE]
        theta_adj <- sweep(Y_i, 2, theta_ests[[i]])
        amt <- 2^rowMeans(theta_adj/gamma_est)
        return(amt)
    }
    
    contribs <- sapply(1:K, contrib_est)
    phats <- t(apply(contribs, 1, function(x) {
        x/sum(x)
    }))
    
    return(phats)
}
