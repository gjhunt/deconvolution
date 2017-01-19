#' Find marker measurements for each group.
#' @param Y the expression data set. 
#' @param pure_samples List of pure samples for each group.
#' @return List. Respective ranked markers for each group.
find_markers <- function(Y, pure_samples) {
    K <- length(pure_samples)
    N <- dim(Y)[2]
    pure <- unlist(pure_samples)
    
    C <- array(0, c(K, N))
    colnames(C) <- colnames(Y)
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
    M <- data.frame(t(M))
    colnames(M) <- c("top", "value")
    M$rn <- 1:N
    
    sM <- M[order(M$top, -M$value), ]
    L <- lapply(1:K, function(i) {
        sM[sM$top == i, "rn"]
    })
    
    return(L)
}
