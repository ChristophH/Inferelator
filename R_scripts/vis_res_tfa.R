require(Matrix)
source('R_scripts/priors.R')

ChristophsPR <- function(ord.idx, gs) {
  prec <- cumsum(gs[ord.idx]) / cumsum(rep(1, length(ord.idx)))
  rec <- cumsum(gs[ord.idx]) / sum(gs)

  prec <- c(prec[1], prec)
  rec <- c(0, rec)

  auc <- ChristophsAUC(rec, prec)
  return(list(prec=prec, rec=rec, auc=auc))
}


ChristophsAUC <- function(x, y) {
  dx <- diff(x)
  my <- y[1:(length(y) - 1)] + diff(y) / 2
  return(sum(dx * my))
}


auprs <- function(mat, gs, prior, eval.on.subset=FALSE) {
  rows <- rep(TRUE, nrow(gs))
  cols <- rep(TRUE, ncol(gs))
  if (eval.on.subset) {
    rows <- apply(gs, 1, sum) > 0
    cols <- apply(gs, 2, sum) > 0
  }
  mat <- mat[rows, cols]
  gs <- gs[rows, cols]
  prior <- prior[rows, cols]
  
  total <- ChristophsPR(order(mat, decreasing=TRUE), gs)$auc
  trainset <- ChristophsPR(order(mat[prior!=0], decreasing=TRUE), gs[prior!=0])$auc
  testset <- ChristophsPR(order(mat[prior==0], decreasing=TRUE), gs[prior==0])$auc
  return(c(total, trainset, testset))
}


dirs <- c('/home/ch1421/Projects/Inferelator/output/bsusnew_spc_tfabs_tp50_nomd_60_15_15_BBSR_1_FALSE/',
          '/home/ch1421/Projects/Inferelator/output/bsusnew_spc_tfabs_tp50_nomd_60_15_15_BBSR_1.1_FALSE/',
          '/home/ch1421/Projects/Inferelator/output/bsusnew_spc_tfabs_tp50_nomd_60_15_15_BBSR_1.6_FALSE/',
          '/home/ch1421/Projects/Inferelator/output/bsusnew_spc_tfabs_tp50_nomd_60_15_15_BBSR_1_TRUE/',
          '/home/ch1421/Projects/Inferelator/output/bsusnew_spc_tfabs_tp50_nomd_60_15_15_BBSR_1.1_TRUE/',
          '/home/ch1421/Projects/Inferelator/output/bsusnew_spc_tfabs_tp50_nomd_60_15_15_BBSR_1.6_TRUE/')
          
res <- c()
for (dir in dirs) {
  print(dir)
  files <- list.files(dir, "combinedconf_.+\\.RData$")
  load(paste(dir, 'params_and_input.RData', sep=''))
  
  .Random.seed <- SEED
  weight <- PARS$prior.weight
  tfa <- PARS$use.tfa
  eval.on.subset <- PARS$eval.on.subset
  gs <- IN$gs.mat
  
  priors <- getPriors(IN$exp.mat, IN$tf.names, IN$priors.mat, IN$gs.mat, 
                      PARS$eval.on.subset, PARS$job.seed, PARS$perc.tp, 
                      PARS$perm.tp, PARS$perc.fp, PARS$perm.fp)
  
  prior.order <- sapply(names(priors), function(x) which(grepl(x, files)))
  files <- files[prior.order]
  
  for (i in 1:length(files)) {
    res.file <- files[i]
    load(file.path(dir, res.file))
    prior <- priors[[i]]
    aupr.vals <- auprs(comb.confs, gs, prior, eval.on.subset)
    print(aupr.vals)
    res <- rbind(res, c(weight, tfa, aupr.vals))
  } 
}

colnames(res) <- c('weight', 'tfa', 'aupr.tot', 'aupr.tr', 'aupr.te')
res <- data.frame(res)
res$aupr.tot <- as.numeric(res$aupr.tot)
res$aupr.te <- as.numeric(res$aupr.te)

library(ggplot2)
library(reshape2)

dfm = subset(melt(res, id.vars=c('weight', 'tfa')), variable!='aupr.tr')

ggplot(dfm, aes(x=as.factor(weight), y=value)) + geom_point(size = 4, position=position_dodge(width=0.5), aes(color = as.factor(tfa), shape = variable))
