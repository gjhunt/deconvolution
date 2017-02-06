#' Generate resampled test data.
#' @param Y The known expression matrix.
#' @param pure_samples Thee pure sample indicies.
#' @param N The number of test data points to generate.
#' @return A list of the generated test data and mixing proportions. 
#' @export
gen_test_data <- function(Y, pure_samples, N = 500) {
    K <- length(pure_samples)
    
    psd_mix <- NULL
    for (j in 1:K) {
        l <- length(pure_samples[[j]])
        tmp <- array(0, c(l, K))
        tmp[, j] <- rep(1, l)
        psd_mix <- do.call(rbind, list(psd_mix, tmp))
    }
    
    t_mix <- array(rbeta(N * K, 1, 1:(N * K)), c(N, K))
    t_mix <- diag(1/(K * rowMeans(t_mix))) %*% t_mix
    t_mix[unlist(pure_samples), ] <- psd_mix
    
    mean_psd <- do.call(rbind, lapply(pure_samples, function(x) {
        colMeans(Y[x, ])
    }))
    
    t_Y <- t_mix %*% mean_psd
    return(list(t_Y = t_Y, t_mix = t_mix))
}
