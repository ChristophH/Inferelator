## Bonneau lab 
## NYU - Center for Genomics and Systems Biology

# Call this script with a job config file as arguments
# Example call: Rscript inferelator.R jobs/dream4_cfg.R


require(Matrix)

rm(list=ls())

source('R_scripts/utils.R')
source('R_scripts/design_and_response.R')
source('R_scripts/priors.R')
source('R_scripts/mi_and_clr.R')
source('R_scripts/bayesianRegression.R')
source('R_scripts/men.R')
source('R_scripts/evaluate.R')
source('R_scripts/tfa.R')


date.time.str <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")


# default job parameters
PARS <- list()

PARS$input.dir <- 'input/dream4'

PARS$exp.mat.file <- 'expression.tsv'
PARS$tf.names.file <- 'tf_names.tsv'
PARS$meta.data.file <- NULL
PARS$priors.file <- NULL
PARS$gold.standard.file <- NULL
PARS$leave.out.file <- NULL
PARS$sc.file <- NULL


PARS$job.seed <- 42  # set to NULL if a random seed should be used
PARS$save.to.dir <- NULL
PARS$num.boots <- 20
PARS$max.preds <- 10
PARS$mi.bins <- 10
PARS$cores <- 8

PARS$delT.max <- 110
PARS$delT.min <- 0
PARS$tau <- 45

PARS$perc.tp <- 0
PARS$perm.tp <- 1
PARS$perc.fp <- 0
PARS$perm.fp <- 1

PARS$eval.on.subset <- FALSE

PARS$method <- 'BBSR'  # 'BBSR' or 'MEN'
PARS$prior.weight <- 1

PARS$use.tfa <- FALSE


# some of the elastic net parameters that are essentially constants;
# only override in config script if you know what you are doing
PARS$enet.sparseModels <- TRUE    # sparser models
PARS$enet.nCv <- 10               # number of cross-validations
PARS$enet.lambda <- c(0, 1, 100)  # l2 weights
PARS$enet.verbose <- FALSE        # print progress to screen
PARS$enet.plot.it <- FALSE        # generate cross-validation plots
PARS$enet.plot.file.name <- NULL  # file name for plots


# input argument is the job config script which overrides the default parameters
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1) {
  job.cfg <- args[1]
} else {
  #job.cfg <- 'jobs/ecoli_bbsr_noprior.R'
  job.cfg <- 'jobs/dream4_bbsr_noprior.R'
  job.cfg <- 'jobs/dream4_bbsr_noprior_sc.R'
}

# load job specific parameters from input config file
if (!is.null(job.cfg)) {
  source(job.cfg)
}


# read input data
IN <- read.input(PARS$input.dir, PARS$exp.mat.file, PARS$tf.names.file, 
                 PARS$meta.data.file, PARS$priors.file, PARS$gold.standard.file,
                 PARS$leave.out.file, PARS$sc.file)

# keep only TFs that are part of the expression data
IN$tf.names <- IN$tf.names[IN$tf.names %in% rownames(IN$exp.mat)]

# order genes so that TFs come before the other genes
gene.order <- rownames(IN$exp.mat)
gene.order <- c(gene.order[match(IN$tf.names, gene.order)], 
                gene.order[which(!(gene.order %in% IN$tf.names))])

IN$exp.mat <- IN$exp.mat[gene.order, ]
if (!is.null(IN$priors.mat)) {
  IN$priors.mat <- IN$priors.mat[gene.order, IN$tf.names]
}
if (!is.null(IN$sc.mat)) {
  IN$sc.mat <- IN$sc.mat[gene.order, IN$tf.names]
}
if (!is.null(IN$gs.mat)) {
  IN$gs.mat <- IN$gs.mat[gene.order, IN$tf.names]
}


# no meta data given - assume all steady state measurements
if (is.null(IN$meta.data)) {
  IN$meta.data <- trivial.meta.data(colnames(IN$exp.mat))
}


# create dummy clusterStack - a real clusterStack is only needed when inferring 
# on bi-clusters
clusterStack <- trivial.cluster.stack(IN$exp.mat)


# set the random seed
if(!is.null(PARS$job.seed)) {
  set.seed(PARS$job.seed, "Mersenne-Twister", "Inversion")
  cat("RNG seed has been set to ", PARS$job.seed, "\n")
} else {
  ignore <- runif(1)
}
SEED <- .Random.seed

if(is.null(PARS$save.to.dir)) {
  PARS$save.to.dir <- file.path(PARS$input.dir, date.time.str)
}
cat("Output dir:", PARS$save.to.dir, "\n")
if (!file.exists(PARS$save.to.dir)){
  dir.create(PARS$save.to.dir, recursive=TRUE)
}


##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# create design and response matrix
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.

# get design matrix
x <- get_usr_chosen_design(IN$meta.data, IN$exp.mat, PARS$delT.min, 
                           PARS$delT.max, time_delayed=T, all_intervals=F, 
                           use_t0_as_steady_state=F, 
                           use_delt_bigger_than_cutoff_as_steady_state=T)
IN$design_matrix_steady_state <- x[[1]]
IN$design_matrix_time_series <- x[[2]]

# get response matrix
x <- get_usr_chosen_response(IN$meta.data, IN$exp.mat, PARS$delT.min, 
                             PARS$delT.max, 'inf_1', PARS$tau, 
                             use_t0_as_steady_state=F, 
                             use_delt_bigger_than_cutoff_as_steady_state=T)
IN$response_matrix_steady_state <- x[[1]]
IN$response_matrix_time_series <- x[[2]]

# make final design/response matrices
x <- make_final_design_and_response_matrix(IN$design_matrix_steady_state,
                                           IN$design_matrix_time_series,
                                           IN$response_matrix_steady_state,
                                           IN$response_matrix_time_series,
                                           'all', clusterStack, IN$exp.mat, 
                                           IN$tf.names, make.des.red.exp=T)
IN$final_response_matrix <- x[[1]]
IN$final_design_matrix <- x[[2]]
resp.idx <- x[[3]]


if (!all(apply(resp.idx, 1, identical, resp.idx[1,]))) {
    stop('This version of the Inferelator does not support biclusters. Sorry.')
}
    

##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# set up the bootstrap permutations
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.

IN$bs.pi <- matrix(0, nrow=PARS$num.boots, ncol=ncol(resp.idx))
if (PARS$num.boots == 1) {
  IN$bs.pi[1, ] <- resp.idx[1, ]
} else {
  for (bootstrap in 1:PARS$num.boots) {
    IN$bs.pi[bootstrap, ] <- resp.idx[1, sample(ncol(resp.idx), replace=TRUE)]
  }
}


##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# parse priors parameters and set up priors list
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.

priors <- getPriors(IN$exp.mat, IN$tf.names, IN$priors.mat, IN$gs.mat, 
                    PARS$eval.on.subset, PARS$job.seed, PARS$perc.tp, 
                    PARS$perm.tp, PARS$perc.fp, PARS$perm.fp)


##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# main loop
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.

for (prior.name in names(priors)) {
  cat('Method:', PARS$method, '\nWeight:', PARS$prior.weight, '\nPriors:', prior.name, '\n')
  prior <- priors[[prior.name]]
  
  # set the prior weights matrix
  no.pr.weight <- 1
  if (sum(prior != 0) > 0) {
    if (PARS$prior.weight == no.pr.weight) {
      warning(paste('Priors present, but they will not be used, because \
                    PARS$prior.weight is set to ', no.pr.weight, '.', sep=''), 
              immediate. = TRUE)
    }
    if (PARS$method == 'BBSR') {
      no.pr.weight <- 1 / PARS$prior.weight
    }
  }
  weights.mat <- matrix(no.pr.weight, nrow(IN$exp.mat), length(IN$tf.names))
  weights.mat[prior != 0] <- PARS$prior.weight
  
  betas <- list()
  betas.resc <- list()
  for (bootstrap in 1:PARS$num.boots) {
    cat("Bootstrap", bootstrap, "of", PARS$num.boots, "\n")
    
    # set up bootstrap specific design and response
    X <- IN$final_design_matrix[, IN$bs.pi[bootstrap, ]]
    Y <- IN$final_response_matrix[, IN$bs.pi[bootstrap, ]]

    if (nrow(X) > 6000) {
      X <- X[IN$tf.names, ]  # speeds up MI calculation for large datasets
    }
    
    if(PARS$use.tfa) {
      X <- tfa(prior, Y, X, PARS$cores)
    }

    # fill mutual information matrices
    cat("Calculating MI\n")	
    Ms <- mi(t(Y), t(X), nbins=PARS$mi.bins, cpu.n=PARS$cores)
    diag(Ms) <- 0
    cat("Calculating Background MI\n")
    Ms_bg <- mi(t(X), t(X), nbins=PARS$mi.bins, cpu.n=PARS$cores)
    diag(Ms_bg) <- 0
    
    # get CLR matrix
    cat("Calculating CLR Matrix\n")
    clr.mat = mixedCLR(Ms_bg,Ms)
    dimnames(clr.mat) <- list(rownames(Y), rownames(X))
    clr.mat <- clr.mat[, IN$tf.names]
    
    # DREAM8 induced change:
    for (tf1 in IN$tf.names) {
      for (tf2 in IN$tf.names) {
        if (tf1 != tf2) {
          #if (clr.mat[tf1, tf2] > clr.mat[tf2, tf1]) {
          if (Ms[tf1, tf2] > Ms[tf2, tf1]) {
            clr.mat[tf2, tf1] <- min(clr.mat)
          } else if (Ms[tf1, tf2] < Ms[tf2, tf1]) {
            clr.mat[tf1, tf2] <- min(clr.mat)
          }
        }
      }
    }
    
    # get the sparse ODE models
    X <- X[IN$tf.names, ]
    cat('Calculating sparse ODE models\n')
    if (PARS$method == 'BBSR') {
      x <- BBSR(X, Y, clr.mat, PARS$max.preds, no.pr.weight, weights.mat, 
                IN$sc.mat, PARS$cores)
    }
    if (PARS$method == 'MEN' ) {
      x <- mclapply(1:nrow(Y), callMEN, Xs=X, Y=Y, 
                    clr.mat=clr.mat, nS=PARS$max.preds, nCv=PARS$enet.nCv,
                    lambda=PARS$enet.lambda, verbose=PARS$enet.verbose, 
                    plot.it=PARS$enet.plot.it, 
                    plot.file.name=PARS$enet.plot.file.name, 
                    weights.mat=weights.mat, no.pr.val=no.pr.weight, 
                    mc.cores=PARS$cores)
    }
    cat('\n')
    
    # our output will be a list holding two matrices: betas and betas.resc
    bs.betas <- Matrix(0, nrow(Y), nrow(X), 
                       dimnames=list(rownames(Y), rownames(X)))
    bs.betas.resc <- Matrix(0, nrow(Y), nrow(X), 
                            dimnames=list(rownames(Y), rownames(X)))
    for (res in x) {
      bs.betas[res$ind, res$pp] <- res$betas
      bs.betas.resc[res$ind, res$pp] <- res$betas.resc
    }
    betas[[bootstrap]] <- bs.betas
    betas.resc[[bootstrap]] <- bs.betas.resc
    
  }  # end bootstrap for loop
  
  res.file <- paste(PARS$save.to.dir, "/betas_", prior.name, "_", PARS$prior.weight, ".RData", sep="")
  save(betas, betas.resc, file = res.file)
  
  # rank-combine the rescaled betas (confidence scores) of the bootstraps
  confs.file <- sub('/betas_', '/combinedconf_', res.file)
  comb.confs <- Matrix(0, nrow(betas.resc[[1]]), ncol(betas.resc[[1]]), 
                       dimnames=dimnames(betas.resc[[1]]))
  for (beta.resc in betas.resc) {
    comb.confs <- comb.confs + rank(as.matrix(beta.resc), ties.method='average')
  }
  save(comb.confs, file=confs.file)
  
}  # end prior.name loop

save(PARS, IN, SEED, file = paste(PARS$save.to.dir, "/params_and_input.RData", sep=""))

if (!is.null(IN$gs.mat)) {
  cat('Using gold standard to evaluate results. Evaluate on subset is set to', PARS$eval.on.subset, '. \n')
  summarizeResults(PARS$save.to.dir, PARS$eval.on.subset)
}

