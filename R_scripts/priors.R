require('Matrix')

getPriors <- function(exp.mat, tf.names, priors.mat, gs.mat, eval.on.subset, 
                      job.seed, perc.tp, perm.tp, perc.fp, perm.fp, sel.mode) {
  # priors.pars is list of prior parameters; every entry is a vector of four
  # elements: perc.tp, tp permutation number, perc.fp, fp permutation number

  priors.pars <- list()

  # Note: We require the perc.*p and perm.*p vectors to be of the same length as
  # in the following code only the same positions of the vectors will be paired.

  lengths <- sapply(list(perc.tp, perm.tp, perc.fp, perm.fp), length)
  if (length(unique(lengths)) != 1) {
    stop("Error parsing prior parameters: perc.tp, perm.tp, perc.fp, perm.fp \
          don't have the same length.")
  }

  for (pos in 1:unique(lengths)) {
    
    # there need not be the same number of permutations for tp and fp
    # but we have to make shure that the smaller one can be repeated to match 
    # the other one
    
    rmndr <- max(perm.tp[pos], perm.fp[pos]) %% min(perm.tp[pos], perm.fp[pos])
    if (rmndr != 0) {
      stop("Error parsing prior parameters: Larger number of permutations is \
            not multiple of smaller number.")
    }
    
    priors.pos <- cbind(perc.tp[pos], 1:perm.tp[pos], perc.fp[pos], 1:perm.fp[pos])

    # every row represents one combination of priors we want to test;
    # append to list of priors parameters
    priors.pars <- c(priors.pars, as.list(as.data.frame(t(priors.pos))))
    
  }
    
  priors <- list()
  for (i in 1:length(priors.pars)) {
    pp <- priors.pars[[i]]
    
    priors[[i]] <- Matrix(0, nrow(exp.mat), length(tf.names), dimnames=dimnames(priors.mat))
    names(priors)[i] <- paste('frac_tp_', pp[1], '_perm_', pp[2], '--frac_fp_', 
                              pp[3], '_perm_', pp[4], sep="")
    
    if (pp[1] > 0 | pp[3] > 0) {
      priors[[i]] <- getPriorMatrix(priors.mat, pp, gs.mat, eval.on.subset, 
                                    job.seed, sel.mode)
    }
    
  }
  return(priors)
}

getPriorMatrix <- function(priors, prior.pars, gs, from.subset, seed, sel.mode) {
  # Given the gold standard, path of the input data, and a set of prior
  # parameters, this function returns a single -1, 0, 1 matrix of random 
  # priors.
  
  perc.tp <- prior.pars[1]
  perm.tp <- prior.pars[2]
  perc.fp <- prior.pars[3]
  perm.fp <- prior.pars[4]

  if (!from.subset) {
    if (sel.mode == 'random') {
      p.mat <- makePriorMat(priors, perc.tp, perm.tp + seed) + 
               makePriorMat(priors, perc.fp, perm.fp + seed, false.priors = TRUE)
    }
    if (sel.mode == 'tf') {
      p.mat <- makePriorMatTF(priors, perc.tp, perm.tp + seed)
    }
  } else {
    p.mat <- Matrix(0, nrow(priors), ncol(priors))
    rows  <- apply(gs,1,sum) > 0
    cols <- apply(gs,2,sum) > 0
    p.mat[rows, cols] <- getPriorMatrix(priors[rows, cols], prior.pars, NULL, 
                                        FALSE, seed, sel.mode)
  }
  dimnames(p.mat) <- dimnames(priors)
  return(p.mat)
  
}

makePriorMat <- function(priors, perc, perm, false.priors = FALSE) {
  # Creates prior matrix based on:
  #
  # priors - matrix with priors
  # perc - percentage of priors to use
  # perm - the permutation of this priors-perc combination
  # false.priors - whether to use false priors
  #
  # Notes: By setting the RNG seed using perm, we make sure that, for a given
  # priors-perm combination, matrices with higher perc values include all the 
  # priors of those with lower perc values.

  # save the state of the RNG
  rng.state <- .Random.seed

  set.seed(perm, "Mersenne-Twister", "Inversion")

  n.priors <- floor(sum(priors != 0) * perc / 100)
  p.mat <- Matrix(0, nrow(priors), ncol(priors))
  
  if (n.priors > 0) {
    if (false.priors) {
      prior.order <- sample(which(priors == 0))
      p.mat[prior.order[1:n.priors]] <- sample(priors[priors != 0], n.priors, replace=TRUE)
    } else {
      prior.order <- sample(which(priors != 0))
      p.mat[prior.order[1:n.priors]] <- priors[prior.order[1:n.priors]]
    }
    
    # The above code ignores whether perc is set too high, but we should warn
    # the user.
    if (n.priors > length(prior.order)) {
      warning("Percent of priors set too high. Only the max used.", call.=TRUE)
    }
  }

  # return RNG to old state
  .Random.seed <<- rng.state
  
  return(p.mat)
}

makePriorMatTF <- function(priors, perc, perm, false.priors = FALSE) {
  # Creates prior matrix based on:

  # save the state of the RNG
  rng.state <- .Random.seed

  set.seed(perm, "Mersenne-Twister", "Inversion")

  # how many TFs do we pick
  tf.out.deg <- apply(priors != 0, 2, sum)
  n.tfs <- floor(sum(tf.out.deg > 0) * perc / 100)
  p.mat <- Matrix(0, nrow(priors), ncol(priors))
  
  if (n.tfs > 0) {
    if (false.priors) {
      stop('False priors not supported in all-or-nothing TF mode')
    } else {
      tf.order <- sample(which(tf.out.deg > 0))
      p.mat[, tf.order[1:n.tfs]] <- priors[, tf.order[1:n.tfs]]
    }
  }

  # return RNG to old state
  .Random.seed <<- rng.state
  
  return(p.mat)
}

