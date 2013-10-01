require('Matrix')
require('mvtnorm')
require('ggplot2')
require('reshape2')
require('kernlab')

clif <- function(X, Y, Pi, clr.mat, nS, no.pr.val, weights.mat, cores) {
  if (!all(apply(Pi, 1, identical, Pi[1,]))) {
    stop('BBSR not implemented for biclusters. Use CallBestSubSetRegression instead')
  }
  
  perm <- Pi[1, ]
  
  # let us ignore the TF responses for now
  is.tf <- rownames(Y) %in% rownames(X)
  Y <- Y[!is.tf ,]
  clr.mat <- clr.mat[!is.tf, ]
  weights.mat <- weights.mat[!is.tf, ]
  
  # Scale and permute design and response matrix
  X <- t(scale(t(X[, perm])))
  Y <- t(scale(t(Y[, perm])))
  
  G <- nrow(Y)  # number of genes
  K.ini <- 7  # number of clusters
  
  # no clustering
  nc.beta.mat <- no.clustering(X, Y, weights.mat, no.pr.val, clr.mat, nS, cores)
  
  # spectral clustering based on response
  spec.res <- spectral.clustering(Y, K.ini)
  
  # an initial clustering
  #mem.probs <- diag(G)
  mem.probs <- random.clustering(G, K.ini)
  
  # variances of conditions
  sigma <- diag(apply(Y, 2, var))
  
  p.th <- 1e-2
  hard.clustering <- max.col(mem.probs)
  
  for (it in 1:12) {
    cat('\nIteration', it, '\n')
    
    # M-step
    ########
  
    cl.size <- apply(mem.probs, 2, sum)
    cl.rel.size <- cl.size / G
    too.small <- cl.size < p.th
    cat('These clusters are (almost) empty:', which(too.small), '\n')
    
    # create a new response
    Y.new <- (t(mem.probs) %*% Y) / cl.size
    Y.new <- Y.new[!too.small, ]
    
    K <- nrow(Y.new)  # number of clusters
    cat('Number of clusters', K, '\n')
    
    # create a new weights matrix
    weights.mat.new <- (t(mem.probs) %*% weights.mat) / cl.size
    weights.mat.new <- weights.mat.new[!too.small, ]
    
    # create a new clr matrix
    clr.mat.new <- (t(mem.probs) %*% clr.mat) / cl.size
    clr.mat.new <- clr.mat.new[!too.small, ]
    
    pp <- potential.predictors(weights.mat.new, no.pr.val, clr.mat.new, nS)
    
    beta.mat <- betas.from.bbsr(X, Y.new, pp, weights.mat.new, nS, cores)
    
    
    # E-step
    ########
    
    # cluster means
    mean.mat <- beta.mat %*% X
    
    # cluster membership log-probs, each column is a cluster, each row a gene
    logprobs <- matrix(0, G, K)
    for (i in 1:K) {
      logprobs[, i] <- dmvnorm(Y, mean.mat[i, ], sigma, log=T)
    }
    
    # complete data log likelihood
    ll <- sum(apply(t(t(logprobs) + cl.rel.size), 1, logsumexp))
    cat('Complete data log likelihood', ll, '\n')
    
    # cluster membership probabilities
    lp.sums <- apply(logprobs, 1, logsumexp)
    mem.probs.new <- exp(logprobs - lp.sums)
    #print(max.col(mem.probs))
    #print(max.col(mem.probs.new))
    if (identical(max.col(mem.probs.new), max.col(mem.probs))) {
      cat('Hard cluster assignment did not change\n')
    }
    mem.probs <- mem.probs.new
    hard.clustering <- cbind(hard.clustering, max.col(mem.probs))

    #dev.new()
    #heatmap(mem.probs, Rowv=NA, Colv=NA, scale='none')
    
  }
  vis.clustering(X, Y, max.col(mem.probs), beta.mat, mean.mat, nc.beta.mat)
  #vis.clustering(X, Y, spec.res, beta.mat, mean.mat, nc.beta.mat)
  #vis.clustering(X, Y, spectral.clustering(nc.beta.mat, K.ini), beta.mat, mean.mat, nc.beta.mat)
  print(table(spec.res, max.col(mem.probs)))
  browser()
}

betas.from.bbsr <- function(X, Y, pp, weights.mat, nS, cores) {
  out.list <- mclapply(1:nrow(Y), BBSRforOneGene, X, Y, pp, weights.mat, nS, mc.cores=cores)
  betas <- Matrix(0, nrow(Y), nrow(X))
  for (res in out.list) {
    betas[res$ind, res$pp] <- res$betas
  }
  return(betas)
}

logsumexp <- function(x) {
  x.max <- max(x)
  return(log(sum(exp(x - x.max))) + x.max)
}

potential.predictors <- function(weights.mat, no.pr.val, clr.mat, nS) {
  G <- nrow(weights.mat)
  K <- ncol(weights.mat)
  pp <- matrix(FALSE, G, K)  # predictors that will be used in the regression
  
  # keep all predictors that we have priors for
  pp[abs(weights.mat - no.pr.val) > 3 * .Machine$double.eps] <- TRUE
  
  # for each gene, add the top nS predictors of the list to possible predictors
  for (ind in 1:G) {
    clr.order <- order(clr.mat[ind, ], decreasing=TRUE)
    pp[ind, clr.order[1:min(K, nS)]] <- TRUE
  }
  return(pp)
}

random.clustering <- function(G, K) {
  ret <- matrix(0, G, K)
  for (i in 1:nrow(ret)) {
    ret[i, sample(K, 1)] <- 1
  }
  return(ret)
}

spectral.clustering <- function(Y, K) {
  spec.res <- specc(Y, centers=K, kernel='splinedot')
  return(as.numeric(as.character(spec.res)))
}

vis.clustering <- function(X, Y, lables, beta.mat, mean.mat, nc.beta.mat) {
  
  gene.order <- order(lables)
  
  dev.new()
  beta.m <- melt(as.matrix(beta.mat))
  g <- ggplot(beta.m, aes(factor(Var2), factor(Var1))) + 
    geom_tile(aes(fill = value), color='gray') + 
    scale_fill_gradient2(low='red', mid='white', high='green', midpoint=0, name='Beta') +
    #scale_fill_continues(name="Experimental\nCondition")
    xlab('TF') + ylab('Cluster') + theme(legend.position="bottom")
  print(g)
  
  dev.new()
  hlines <- melt(which(as.logical(diff(lables[gene.order]))) + 0.5)
  rownames(nc.beta.mat) <- paste(lables, 1:nrow(nc.beta.mat), sep='_')
  df <- melt(nc.beta.mat[gene.order, ])
  df$Var1 <- factor(df$Var1, ordered=TRUE, levels=as.character(unique(df$Var1)))
  g <- ggplot(df, aes(factor(Var2), Var1)) + 
    geom_tile(aes(fill = value), color='gray') + 
    scale_fill_gradient2(low='red', mid='white', high='green', midpoint=0, name='Beta') + 
    geom_hline(data=hlines, aes(yintercept=value)) + 
    xlab('TF') + ylab('Cluster_Gene') + theme(legend.position="bottom")
  print(g)
  
  dev.new()
  df <- melt(Y[gene.order, ])
  df$Var1 <- factor(df$Var1, ordered=TRUE, levels=as.character(unique(df$Var1)))
  g <- ggplot(df, aes(Var2, Var1)) + geom_tile(aes(fill = value)) + 
    scale_fill_gradient2(low='red', mid='white', high='green', midpoint=0, name='Expr') + 
    geom_hline(data=hlines, aes(yintercept=value)) + 
    xlab('Condition') + ylab('Gene') + theme(legend.position="bottom")
  print(g)
  
  return()
  dev.new()
  #pdf(width=10, height=8, useDingbats=F)
  par(mar=c(0, 0, 0, 0) + 0.1)
  par(mfrow=c(ceiling(K/4), 4))
  for (k in 1:K) {
    if (sum(lables == k) > 0) {
      Y.k <- matrix(Y[lables == k, ], sum(lables == k))
      plot(c(1, ncol(Y)), range(Y), type='n', yaxt='n', xaxt='n', ann=FALSE)
      for (i in 1:nrow(Y.k)) {
        lines(Y.k[i, ], type='p', col='cyan', lty=1, pch=16)
      }
      lines(apply(Y.k, 2, mean), type='p', col='green', lty=1, pch=1)
      #matplot(t(Y[lables == k, ]), type='l', col='red', lty=1)
      lines(mean.mat[k, ], type='p', col='blue', lty=1, pch=1)
    }
  }
  #dev.close()
}

no.clustering <- function(X, Y, weights.mat, no.pr.val, clr.mat, nS, cores) {
  pp <- potential.predictors(weights.mat, no.pr.val, clr.mat, nS)
  beta.mat <- betas.from.bbsr(X, Y, pp, weights.mat, nS, cores)
  return(as.matrix(beta.mat))
}

