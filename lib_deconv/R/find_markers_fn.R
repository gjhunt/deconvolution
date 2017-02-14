#' Find marker measurements for each group.
#' @param Y the expression data set. 
#' @param pure_samples List of pure samples for each group.
#' @param data_type The data type. 
#' @param gamma The sensitivity tuning parameter if we want to over-ride choice by ``data-type``.
#' @param data_type A string indicating the data type. Use to choose pre-estimated gamma value. Current support for probe-level microarray as ``microarray-probe'', gene-level microarray as ``microarray-gene'' or rna-seq as ``rna-seq''.
#' @param method Method to select markers. 
#' @return List with two elements. ``L'' is respective ranked markers for each group and ``V'' is the corresponding values of the ranking method (lower are better). 
#' @export
find_markers <- function(Y, pure_samples, data_type = NULL, gamma = NULL, method = "eta") {
    if (is.null(gamma)) 
        gamma <- get_gamma(data_type)
    
    K <- length(pure_samples)
    N <- dim(Y)[2]
    pure <- unlist(pure_samples)
    
    C <- array(0, c(K, N))
    colnames(C) <- colnames(Y)
    
    if (method == "eta") {
        eta_hats <- t(sapply(pure_samples, function(x) {
            colMeans(exp(Y[x, , drop = FALSE]))/gamma
        }))
        C <- t(sapply(1:K, function(i) {
            apply(eta_hats[-i, ], 2, sum)/eta_hats[i, ]
        }))
        
        pick_top <- function(x) {
            m <- which.min(x)
            return(c(m, -x[m]))
        }
        
        M <- apply(C, 2, pick_top)
    }
    
    if (method == "regression") {
        for (i in 1:K) {
            X <- as.numeric(pure %in% pure_samples[[i]])
            m <- lm(Y[pure, ] ~ 1 + X)
            cfdf <- data.frame(t(coef(m)))
            C[i, ] <- cfdf$X
        }
        
        pick_top <- function(x) {
            m <- which.max(x)
            return(c(m, x[m]))
        }
        
        M <- apply(C, 2, pick_top)
    }
    
    if (method == "diff") {
        for (i in 1:K) {
            C[i, ] <- apply(Y[pure_samples[[i]], , drop = FALSE], 2, median)
        }
        tops <- apply(C, 2, which.max)
        
        f <- function(i) {
            A <- C[order(C[, i], decreasing = TRUE)[1:2], i]
            A[1]/A[2]
        }
        res <- sapply(1:N, f)
        
        pick_top <- function(i) {
            m <- which.max(C[, i])
            return(c(m, res[i]))
        }
        
        M <- sapply(1:N, pick_top)
    }
    
    M <- data.frame(t(M))
    colnames(M) <- c("top", "value")
    M$rn <- 1:N
    
    sM <- M[order(M$top, -M$value), ]
    L <- lapply(1:K, function(i) {
        sM[sM$top == i, "rn"]
    })
    
    V <- lapply(1:K, function(i) {
        sM[sM$top == i, "value"]
    })
    
    
    return(list(L = L, V = V))
}
