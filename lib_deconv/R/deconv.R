deconv = function(affy_batch, pure_samples,
                  marker_probes = NULL,marker_info_only=FALSE,
                  n_choose = 200,model='linear'){
    ips = indexProbes(affy_batch,which="both")
    lengths = sapply(ips,length)
    gnames_by_probe = unlist(sapply(1:length(ips),function(i){rep(names(ips)[i],lengths[i])}))
    
    d = log2(t(intensity(affy_batch)))
    d = d[,unlist(ips)]

    if(is.null(marker_probes)){
        M = rank_marker_probes(d,pure_samples)
        M$gene = gnames_by_probe
        marker_info = sort_marker_probes(M,pure_samples)
        
        marker_probes = lapply(marker_info,function(x){rownames(x[1:n_choose,])})
    }

    if(marker_info_only){
        return(list("marker_info"=marker_info))
    }
    
    theta_hats = estimate_thetas(d,pure_samples,marker_probes)
    gamma_hat = estimate_gamma()

    estimator = function(Y){
        return(phats(Y, marker_probes, theta_hats, gamma_hat,model))
    }

    ests = estimator(d)
    colnames(ests) = names(pure_samples)
    return(list("estimates"=ests,"marker_probes"=marker_probes))
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
        thetas[[i]] = sapply(marker_probes[[i]],compute_theta)
    }
    return(thetas)
}

## Choose the number of marker probes and organize.
# M -- the matrix output from rank_marker_probes with gene info. 
sort_marker_probes = function(M,pure_samples){
    K = length(unique(M$top_probe))
    marker_probes = list()
    for(i in 1:K){
        type_i = which(M$top_probe==i)
        M_i = M[type_i,]
        top_M = order(M_i$slope,decreasing=TRUE)
        tops = M_i[top_M,]
        marker_probes[[i]] = tops
    }
    names(marker_probes)=names(pure_samples)
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
    for(i in 1:K){
        X = as.numeric(pure %in% pure_samples[[i]])
        m = lm(d[pure,] ~ 1+X)
        cfdf = data.frame(t(coef(m)))
        C[i,] = cfdf$X
    }

    pick_top = function(x){
        m = which.max(x)
        return(c(m,x[m]))
    }

    ordered_C = apply(C,2,pick_top)
    ordered_C = data.frame(t(ordered_C))
    colnames(ordered_C) = c("top_probe","slope")
    return(ordered_C)
}

## Estimate the gene type proportions. 
# Y -- A single DNA microarray to deconvolve.
# marker_probes -- list of list of vectors of probes. Level 1 is cell type, level two is genes.
# thetas -- same structure as marker_probes but contains theta estimates for probes.
# gamma -- scalar estimate of gamma.
phats = function(Y, marker_probes, thetas, gamma, model){

    if(model=='linear'){
        inv_fn = function(X){
            X/gamma
        }
    } else if(model=='sigmoid'){
        param = c(-7,  7.9989337, -2.7370245,  0.4616581)
        sigmoid_inv = function(x,t){return((-1/t[4])*(log((t[2]/(x-t[1]))-1)+t[3]))}
        esigmoid_inv = function(x,t){
            smin = t[1]+1E-4
            smax = t[1]+t[2]-1E-1
            if(x>smax){
                return(sigmoid_inv(smax,t))
            } else if(x<smin){
                return(sigmoid_inv(smin,t))
            } else {
                return(sigmoid_inv(x,t))
            }
        }
        esig_inv = function(x) sapply(x,function(y) esigmoid_inv(y,param))
        inv_fn = function(X) apply(X,2,esig_inv)
    }
    
    K = length(marker_probes)
    
    contrib_est = function(i){
        Y_i = lapply(marker_probes[[i]],function(x){Y[,x,drop=FALSE]})
        theta_adj = mapply('-',Y_i,thetas[[i]],SIMPLIFY=FALSE)
        amt = 2^rowMeans(inv_fn(do.call(cbind,theta_adj)))
        return(amt)
    }

    contribs = sapply(1:K,contrib_est)
    phats = t(apply(contribs,1,function(x){x/sum(x)}))
    return(phats)
}

