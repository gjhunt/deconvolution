deconv = function(affy_batch, pure_samples){
    ips = indexProbes(affy_batch,which="both")
    lengths = sapply(ips,length)
    gnames_by_probe = unlist(sapply(1:length(ips),function(i){rep(names(ips)[i],lengths[i])}))
    
    d = log2(t(intensity(affy_batch)))
    d = d[,unlist(ips)]

    M = rank_marker_probes(d,pure_samples)
    M$gene = gnames_by_probe

    marker_probes = choose_marker_probes(M)
    theta_hats = estimate_thetas(d,pure_samples,marker_probes)
    gamma_hat = estimate_gamma()

    estimator = function(Y){
        return(phats(Y, marker_probes, theta_hats, gamma_hat))
    }

    return(estimator(d))
}

## Estimate slope gamma.
estimate_gamma = function(){
    return(.7)
}

## Estimate the intercept terms
# d -- the probe level data set.
# pure_samples -- list of pure samples for each group.
# marker_probes -- list of lists of marker probes by gene
estimate_thetas = function(d,pure_samples,marker_probes){
    K = length(pure_samples)

    thetas = list()
    for(i in 1:K){
        compute_theta = function(x){colMeans(d[pure_samples[[i]],x,drop=FALSE])}
        thetas[[i]] = lapply(marker_probes[[i]],compute_theta)
    }
    return(thetas)
}

## Choose the number of marker probes and organize.
# M -- the matrix output from rank_marker_probes with gene info. 
choose_marker_probes = function(M,q_top=.9985){
    K = length(unique(M$top_probe))
    marker_probes = list()
    for(i in 1:K){
        type_i = which(M$top_probe==i)
        M_i = M[type_i,]
        top_M = which(M_i$slope>=quantile(M_i$slope,q_top))
        tops = M_i[top_M,]
        tmp_list = list()
        for(j in unique(tops$gene)){
            tmp_list[[j]] = rownames(tops)[which(tops$gene==j)]
        }
        marker_probes[[i]] = tmp_list
    }
    return(marker_probes)
}

## Rank the probes as markers for each group.
# d -- the probe level data set.
# pure_samples -- list of pure samples for each group.
rank_marker_probes = function(d,pure_samples){
    K = length(pure_samples)
    N = dim(d)[2]
    pure = unlist(pure_samples)
    
    C = array(0,c(K,N))
    colnames(C) = colnames(d)
    for(i in 1:3){
        X = as.numeric(pure %in% pure_samples[[i]])
        m = lm(d[pure,] ~ 1+X)
        C[i,] = coef(m)[2,]
    }

    pick_top = function(x){
        m = which.max(x)
        return(c(m,x[m]))
    }

    ordered_C = apply(C,2,pick_top)
    #slope_ord = order(ordered_C[2,],decreasing=TRUE)
    #ordered_C = ordered_C[,slope_ord]
    ordered_C = data.frame(t(ordered_C))
    colnames(ordered_C) = c("top_probe","slope")
    return(ordered_C)
}

## Estimate the gene type proportions. 
# Y -- A single DNA microarray to deconvolve.
# marker_probes -- list of list of vectors of probes. Level 1 is cell type, level two is genes.
# thetas -- same structure as marker_probes but contains theta estimates for probes.
# gamma -- scalar estimate of gamma.
phats = function(Y, marker_probes, thetas, gamma){    
    K = length(marker_probes)
    
    contrib_est = function(i){
        Y_i = lapply(marker_probes[[i]],function(x){Y[,x,drop=FALSE]})
        theta_adj = mapply('-',Y_i,thetas[[i]],SIMPLIFY=FALSE)
        ests = lapply(theta_adj,function(x){rowMeans(2^(x/gamma))})
        amt = Reduce("+",ests)/length(ests)
        return(amt)
    }

    contribs = sapply(1:K,contrib_est)
    phats = t(apply(contribs,1,function(x){x/sum(x)}))
    return(phats)
}
