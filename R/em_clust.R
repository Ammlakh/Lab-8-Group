#' Preform EM clustering on a dataset.
#'
#'
#' @usage EM(dat, num)
#'
#'
#' @param dat data that can be made into a matrix for em clustering
#'
#' @param num the number of clusters
#'
#' @return A list containing means of the clusters, clustering vector, and variance.
#'
#' @import tidyverse
#' @import mvtnorm
#'
#' @export

em_clust = function(dat, num, show_cluster_probs = F){

    priors = rep(1/num, num)
    dat = as.matrix(dat)

    if(ncol(dat) > 1){
        dat = data.frame(dat)
        means = sample_n(dat, num)
    }else{
        means = sample(dat, num)}
    dat = as.matrix(dat)
    varlist = vector('list', length = num)

    for(i in 1:num){
        if(ncol(dat) > 1){
            varlist[[i]] = cov(dat)
        }else{
            varlist[[i]] = sd(dat)}

    }

    prob = 0
    filler = 1

    while(sum((prob - filler)^2) > 0.00001){
        filler = prob
        prob = NULL
        if(ncol(dat) > 1){
            for(i in 1:num){
                initprob = dmvnorm(dat, as.numeric(means[i,]), varlist[[i]])
                prob = cbind(prob, initprob)
            }
        }else{
            for(i in 1:num){
                initprob = dnorm(dat, means[i], varlist[[i]])
                prob = cbind(prob, initprob)
            }
        }

        posterior = 0

        for(i in 1:num){
            postupdate = priors[i]*prob[,i]
            posterior = posterior + postupdate
        }

        posts = NULL

        for(i in 1:num){
            initpost = (priors[i]*prob[,i])/posterior
            posts = cbind(posts, initpost)
        }

        priors = apply(posts, 2, mean)
        means = NULL

        if(ncol(dat) > 1){
            for(i in 1:num){
                tempM = colSums(posts[,i]*dat)/sum(posts[,i])
                means = rbind(means, tempM)
            }
        }else{
            for(i in 1:num){
                tempM = sum(posts[,i]*dat)/sum(posts[,i])
                means = c(means, tempM)
            }
        }

        parameter = NULL

        for(i in 1:num){
            para_update = posts[,i]/(length(posts[,i])*priors[i])
            parameter = cbind(parameter, para_update)
        }

        varlist = vector('list', length = num)

        if(ncol(dat) > 1){
            for(i in 1:num){
                tempvar = 0
                for(l in 1:nrow(posts)){
                    varcal = parameter[l,i] * ((dat[l,] - means[i,]) %*% t(dat[l,] - means[i,]))
                    tempvar = tempvar + varcal
                }
                varlist[[i]] = tempvar
            }
        }else{
            varlist = vector('list', length = num)
            for(i in 1:num){
                tempvar = 0
                for(l in 1:nrow(posts)){
                    varcal = parameter[l,i] * ((dat[l,] - means[i]) %*% t(dat[l,] - means[i]))
                    tempvar =tempvar + varcal
                    }
                varlist[[i]] = sqrt(tempvar)
            }
        }

        if(show_cluster_probs == T){
            print(matrix(priors)) # cluster probabilities
            }
    }

    clusters = rep(NA, nrow(posts))
    for(i in 1:nrow(posts)){
        clusters[i] = which(posts[i,] == max(posts[i,]), arr.ind = T)}
    rownames(means) = NULL
    return(list('Clusters' = clusters, 'Means' = means, 'Variance' = varlist))
}
