sigmoid <- function(x, t) {
    return(t[1] + (t[2] - t[1])/(1 + exp(-(t[4] * (x - t[3])))))
}

logistic_fit <- function(x, y) {
    fail <- TRUE
    i <- 0
    max_try <- 10
    t <- c(7, 13, 5, 0.4)
    
    while (fail) {
        i <- i + 1
        if (i > max_try) {
            fail <- 2
            break
        }
        
        fail <- tryCatch({
            if (i != 1) {
                t <- t + rnorm(length(t), 0, sqrt(i))
            }
            init <- list(t1 = t[1], t2 = t[2], t3 = t[3], t4 = t[4])
            mod <- nls(y ~ sigmoid(x, c(t1, t2, t3, t4)), start = init, control = list(warnOnly = FALSE, 
                maxiter = 1000, minFactor = 1e-12))
            fail <- FALSE
        }, error = function(e) {
            return(TRUE)
        })
    }
    
    if (fail < 2) {
        return(list(dat = cbind(x, y), fail = FALSE, coef = coef(mod)))
    }
    
    return(list(dat = cbind(x, y), fail = TRUE, coef = rep(0, 4)))
}
