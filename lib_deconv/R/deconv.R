# --> Y
# A single DNA microarray to deconvolve
# --> marker_probes
# list of list of vectors of probes. Level 1 is cell type, level two is genes.
# --> thetas
# same structure as marker_probes but contains theta estimates for probes.
# --> gamma
# scalar estimate of gamma
phats = function(Y, marker_probes, thetas, gamma){
    K = length(marker_probes)
    
    contrib_est = function(i){
        Y_i = lapply(marker_probes[[i]],function(x){Y[,x]})
        theta_adj = mapply('-',Y_i,thetas[[i]],SIMPLIFY=FALSE)
        ests = lapply(theta_adj,function(x){mean(2^(x/gamma))})
        amt = Reduce("+",ests)/length(ests)
        return(amt)
    }

    contribs = sapply(1:K,contrib_est)
    phats = contribs / sum(contribs)

    return(phats)
}
