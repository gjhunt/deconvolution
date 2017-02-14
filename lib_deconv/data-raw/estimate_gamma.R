source("latin_data_read.R")
source("logistic_fit.R")

eps <- 1e-07

x <- log2(all_titration[, 1] + eps)
y <- sd[, 1]

gfit <- lapply(1:dim(sdg)[2], function(i) {
    logistic_fit(x, sdg[, i])
})
pfit <- lapply(1:dim(sd)[2], function(i) {
    logistic_fit(x, sd[, i])
})

gamma_est <- function(t) {
    0.25 * t[4] * (t[2] - t[1])
}

p_gamma <- median(sapply(pfit, function(x) {
    gamma_est(x$coef)
}))
g_gamma <- median(sapply(gfit, function(x) {
    gamma_est(x$coef)
}))
gma <- list(ma_probe = p_gamma, ma_gene = g_gamma, rna_seq = 1)

devtools::use_data(gma, internal = TRUE, overwrite = TRUE)

