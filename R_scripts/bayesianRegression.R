# TODO: Add comment
# 
# Author: Christoph
###############################################################################
require(mgcv)

BBSR <- function(X, Y, clr.mat, nS, no.pr.val, weights.mat, sc.mat, cores) {

  # Scale and permute design and response matrix
  X <- t(scale(t(X)))
  Y <- t(scale(t(Y)))
  
  G <- nrow(Y)  # number of genes
  K <- nrow(X)  # max number of possible predictors (number of TFs)
  
  pp <- matrix(FALSE, G, K)  # predictors that will be used in the regression
  
  # keep all predictors that we have priors for
  pp[weights.mat != no.pr.val] <- TRUE
  
  # for each gene, add the top nS predictors of the list to possible predictors
  clr.mat[clr.mat == 0] <- NA
  for (ind in 1:G) {
    clr.order <- order(clr.mat[ind, ], decreasing=TRUE, na.last=NA)
    pp[ind, clr.order[1:min(K, nS, length(clr.order))]] <- TRUE
  }
  diag(pp) <- FALSE
  
  out.list <- mclapply(1:G, BBSRforOneGene, X, Y, pp, weights.mat, sc.mat, nS, mc.cores=cores)
  return(out.list)
}

BBSRforOneGene <- function(ind, X, Y, pp, weights.mat, sc.mat, nS) {
  if (ind %% 100 == 0) {
    cat('Progress: BBSR for gene', ind, '\n')
  }
  
  pp.i <- pp[ind, ]
  
  if (sum(pp.i) == 0) {
    return(list(ind=ind, pp=rep(TRUE, length(pp.i)), betas=0, betas.resc=0))
  }

  # create BestSubsetRegression input  
  y <- as.vector(Y[ind, ], mode="numeric")
  x <- t(matrix(X[pp.i, ], ncol=ncol(X)))
  g <- matrix(weights.mat[ind, pp.i], ncol=sum(pp.i))
  if (is.null(sc.mat)) {
    sc <- NULL
  }
  else {
    sc <- sc.mat[ind, pp.i]
  }
  
  # experimental stuff
  spp <- ReduceNumberOfPredictors(y, x, g, nS, sc)
  pp.i[pp.i == TRUE] <- spp
  x <- t(matrix(X[pp.i, ], ncol=ncol(X)))
  g <- matrix(weights.mat[ind, pp.i], ncol=sum(pp.i))
  if (!is.null(sc.mat)) {
    sc <- sc.mat[ind, pp.i]
  }
  
  
  betas <- BestSubsetRegression(y, x, g, sc)
  betas.resc <- PredErrRed(y, x, betas, sc)
  
  return(list(ind=ind, pp=pp.i, betas=betas, betas.resc=betas.resc))
}


ReduceNumberOfPredictors <- function(y, x, g, n, sc) {
  K <- ncol(x)
  
  if (K <= n) {
    return(rep(TRUE, K))
  }

  combos <- cbind(diag(K) == 1, CombCols(diag(K)))
  bics <- ExpBICforAllCombos(y, x, g, combos, sc)
  bics.avg <- apply(t(t(combos) * bics), 1, sum)
  ret <- rep(FALSE, K)
  ret[order(bics.avg)[1:n]] <- TRUE
  
  return(ret)
}

BayesianModelAveraging <- function(y, x, g) {
  
  mprior.size <- g
  bms.out <- bms(cbind(y,x), nmodel=0, mprior='pip', 
                 mprior.size=mprior.size, user.int=F)
  #bms.out <- bms(cbind(y,x), nmodel=0, user.int=F)
  
  tmp <- coef(bms.out)
  tmp <- tmp[order(tmp[, 'Idx']), ]
  
  ret <- tmp[, 'Post Mean']
  ret[tmp[, 'PIP'] < 0.9] <- 0
  #print(sum(tmp[, 'PIP'] > 0.5))
  #return(as.numeric(tmp[, 'PIP'] > 0.9)[order(tmp[, 'Idx'])])
  return(ret)
}

BestSubsetRegressionAllWeights <- function(y, x, g){
  #Do best subset regression without any weight (if this option is chosen)
  #and with the combination of all weights that are chosen
  #return the results in a list
  
  # Q CH 11|18|2011: Is this ever going to be used?
}

BestSubsetRegression <- function(y, x, g, sc=NULL) {
  # Do best subset regression by using all possible combinations of columns of
  # x as predictors of y. Model selection criterion is BIC using results of
  # Bayesian regression with Zellner's g-prior.
  #
  # Args:
  #   y: dependent variable
  #   x: independent variable
  #   g: value for Zellner's g-prior; can be single value or vector
  #  sc: vector of sign constraints, each entry is element of (-1, 0, 1)
  #
  # Returns:
  #   Beta vector of best model

  K <- ncol(x)
  N <- nrow(x)
  ret <- c()

  combos <- AllCombinations(K)
  bics <- ExpBICforAllCombos(y, x, g, combos, sc)
  
  not.done <- TRUE
  while (not.done) {
    best <- which.min(bics)
    #cat('best is', best, '\n')
  
    # For the return value, re-compute beta ignoring g-prior.
    betas <- rep(0, K)
    if (best > 1) {
      x.tmp <- matrix(x[,combos[, best]], N)
      if (is.null(sc)) {
        tryCatch({
          bhat <- solve(crossprod(x.tmp), crossprod(x.tmp, y))
          betas[combos[, best]] <- bhat
          not.done <- FALSE
        }, error = function(e) {
          if (any(grepl('solve.default', e$call)) & grepl('singular', e$message)) {
            # error in solve - system is computationally singular
            cat(bics[best], 'at', best, 'replaced\n')
            bics[best] <<- Inf
          } else {
            stop(e)
          }
        })
      } else {
        sc.tmp <- sc[combos[, best]]
        k <- length(sc.tmp)
        M<-list(X   = x.tmp,
                y   = y,
                p   = sc.tmp,
                Ain = diag(sc.tmp, k, k),
                bin = rep(0, k),
                w   = rep(1, N),
                off = array(0,0),
                S   = list(),
                C   = matrix(0,0,0),
                sp  = array(0,0))
        bhat <- pcls(M)
        betas[combos[, best]] <- bhat
        not.done <- FALSE
      }
    }
    else {
      not.done <- FALSE
    }
  }
  return(betas)
}


AllCombinations <- function(k) {
  # Create a boolean matrix with all possible combinations of 1:k.
  # Output has k rows and 2^k columns where each column is one combination.
  # Note that the first column is all FALSE and corresponds to the null model.
  if (k < 1) {
    stop("No combinations for k < 1")
  }
  
  N <- 2^k
  out <- matrix(FALSE, k, N)
  out[1, 2] <- TRUE
  
  row <- 2
  col <- 3
  while (col < N) {
    out[row, col] <- TRUE
    
    for (i in 1:(col-2)) {
      out[, col + i] <- out[, col] | out[, i + 1]
    }
    
    row <- row + 1
    col <- col * 2 - 1
  }
  
  return(out)
}


CombCols <- function(m) {
  K <- ncol(m)
  ret <- matrix(TRUE, nrow(m), K * (K - 1) / 2)
  ret.col <- 1
  for (i in 1:(K - 1)) {
    for (j in (i + 1):K) {
      ret[, ret.col] <- m[, i] | m[, j]
      ret.col <- ret.col + 1
    }
  }
  return(ret)
}


ExpBICforAllCombos <- function(y, x, g, combos, sc) {
  # For a list of combinations of predictors do Bayesian linear regression,
  # more specifically calculate the parametrization of the inverse gamma 
  # distribution that underlies sigma squared using Zellner's g-prior method.
  # Parameter g can be a vector. The expected value of the log of sigma squared
  # is used to compute expected values of BIC.
  # Returns list of expected BIC values, one for each model.
  K <- ncol(x)
  N <- nrow(x)
  
  C <- ncol(combos)
  bics <- rep(0, C)
  
  # is the first combination the null model?
  first.combo <- 1
  if (sum(combos[, 1]) == 0) {  
    bics[1] <- N * log(var(y))
    first.combo <- 2
  }
  
  # shape parameter for the inverse gamma sigma squared would be drawn from  
  shape <- N / 2
  
  # compute digamma of shape here, so we can re-use it later
  dig.shape <- digamma(shape)
  
  # In Zellner's formulation there is a factor in the calculation of the rate 
  # parameter: 1 / (g + 1)
  # Here we replace the factor with the approriate matrix since g is a vector
  # now.
  var.mult <- matrix(sqrt(1 / (g + 1)), K, K)
  var.mult <- var.mult * t(var.mult)

  # pre-compute the crossproducts that we will need to solve for beta  
  xtx <- crossprod(x)
  xty <- crossprod(x, y)
    
  if (is.null(sc)) {
    for (i in first.combo:C){
      comb <- combos[, i]
      x.tmp <- matrix(x[, comb], N)
      k <- sum(comb)
      
      tryCatch({
        # this is faster than calling lm
        bhat <- solve(xtx[comb, comb], xty[comb])
        
        ssr <- sum((y - x.tmp %*% bhat)^2)  # sum of squares of residuals
      
        # rate parameter for the inverse gamma sigma squared would be drawn from
        # our guess on the regression vector beta is all 0 for sparse models
        rate <- (ssr + 
                (0 - t(bhat)) %*% 
                (xtx[comb, comb] * var.mult[comb, comb]) %*% 
                t(0 - t(bhat))) / 2
        
        # the expected value of the log of sigma squared based on the 
        # parametrization of the inverse gamma by rate and shape
        exp.log.sigma2 <- log(rate) - dig.shape
        
        # expected value of BIC
        bics[i] <- N * exp.log.sigma2 + k * log(N)
        
      }, error = function(e) {
        if (any(grepl('solve.default', e$call)) & grepl('singular', e$message)) {
          # error in solve - system is computationally singular
          bics[i] <<- Inf
        } else {
          stop(e)
        }
      })
      
    }
  } else {
    for (i in first.combo:C){
      comb <- combos[, i]
      x.tmp <- matrix(x[, comb], N)
      sc.tmp <- sc[comb]
      k <- sum(comb)
      
      M<-list(X=x.tmp,
              y=y,
              p=sc.tmp,
              Ain=diag(sc.tmp, k, k),
              bin=rep(0, k),
              w=rep(1, N),
              off=array(0,0),
              S=list(),
              C=matrix(0,0,0),
              sp=array(0,0))
      bhat <- pcls(M)
        
      ssr <- sum((y - x.tmp %*% bhat)^2)  # sum of squares of residuals
      
      # rate parameter for the inverse gamma sigma squared would be drawn from
      # our guess on the regression vector beta is all 0 for sparse models
      rate <- (ssr + 
              (0 - t(bhat)) %*% 
              (xtx[comb, comb] * var.mult[comb, comb]) %*% 
              t(0 - t(bhat))) / 2
      
      # the expected value of the log of sigma squared based on the 
      # parametrization of the inverse gamma by rate and shape
      exp.log.sigma2 <- log(rate) - dig.shape
      
      # expected value of BIC
      bics[i] <- N * exp.log.sigma2 + k * log(N)
      
    }
  }
  
  return(bics)
}





PredErrRed <- function(y, x, beta, sc) {
  # Calculates the error reduction (measured by variance of residuals) of each
  # predictor - compare full model to model without that predictor
  N <- nrow(x)
  K <- ncol(x)
  pred <- beta != 0
  P <- sum(pred)

  # compute sigma^2 for full model
  residuals <- y - x %*% beta
  sigma.sq.full <- var(residuals)

  # this will be the output
  err.red <- rep(0, K)

  # special case if there is only one predictor
  if (P == 1) {
    err.red[pred] <- 1 - sigma.sq.full / var(y)
    return(err.red)
  }

  k <- P - 1
  
  # one by one leave out each predictor and re-compute the model with the
  # remaining ones
  for (i in (1:K)[pred]) {
    pred.tmp <- pred
    pred.tmp[i] <- FALSE
    x.tmp <- matrix(x[,pred.tmp], N, k)
    if (is.null(sc)) {
      #bhat <- solve(t(x.tmp) %*% x.tmp) %*% t(x.tmp) %*% y
      bhat <- solve(crossprod(x.tmp), crossprod(x.tmp, y))
    } else {
      sc.tmp <- sc[pred.tmp]
      M<-list(X=x.tmp,
              y=y,
              p=sc.tmp,
              Ain=diag(sc.tmp, k, k),
              bin=rep(0, k),
              w=rep(1, N),
              off=array(0,0),
              S=list(),
              C=matrix(0,0,0),
              sp=array(0,0))
      bhat <- pcls(M)
    }
    residuals <- y - x.tmp %*% bhat
    sigma.sq <- var(residuals)

    err.red[i] <- 1 - sigma.sq.full / sigma.sq
  }
  return(err.red)
}


