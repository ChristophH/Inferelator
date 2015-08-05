library('Matrix')
source('R_scripts/evaluate.R')

get.mean.and.lh <- function(mat) {
  ret <- list()
  ret$mean <- apply(mat, 1, mean)
  ret$median <- apply(mat, 1, median)
  tmp <- apply(mat, 1, quantile, probs=c(0.05, 0.95))
  ret$low <- tmp[1, ]
  ret$high <- tmp[2, ]
  return(ret)
}

to.confidence.score <- function(mat) {
  tmp <- rank(mat)
  return(1 - (tmp - 1) / (max(tmp) - 1))
}

sum.net <- function(betas, betas.resc, comb.confs, IN, cc.file) {
  
  if (is.null(colnames(comb.confs))) {
    colnames(comb.confs) <- IN$tf.names
    rownames(comb.confs) <- rownames(IN$exp.mat)
  }

  beta.sign <- as.matrix(betas[[1]] * 0)
  beta.non.zero <- as.matrix(betas[[1]] * 0)
  for (n in 1:length(betas)) {
    beta.sign <- beta.sign + sign(betas[[n]])
    beta.non.zero <- beta.non.zero + (betas[[n]] != 0)
  }
  dimnames(beta.sign) <- dimnames(comb.confs)
  dimnames(beta.non.zero) <- dimnames(comb.confs)
  
  # we only care about interactions that are present in at least one bootstrap
  non.zero <- which(beta.non.zero != 0, arr.ind=TRUE)
  vec.ind <- which(beta.non.zero != 0)
  
  # get mean beta score
  betas.bs <- lapply(betas, function(x) x[vec.ind])
  betas.bs <- matrix(unlist(betas.bs), length(betas.bs[[1]]))
  betas.stats <- get.mean.and.lh(betas.bs)
  rm(betas.bs)

  # get mean rescaled beta
  resc.betas <- lapply(betas.resc, function(x) x[vec.ind])
  resc.betas <- matrix(unlist(resc.betas), length(resc.betas[[1]]))
  resc.betas.stats <- get.mean.and.lh(resc.betas)

  # get mean of ranks
  resc.betas.ranks <- apply(-resc.betas, 2, rank)
  resc.betas.ranks.stats <- get.mean.and.lh(resc.betas.ranks)
  rm(resc.betas, resc.betas.ranks)

  # get confidence scores
  conf.scores <- lapply(betas.resc, function(x) to.confidence.score(-as.matrix(x))[vec.ind])
  conf.scores <- matrix(unlist(conf.scores), length(conf.scores[[1]]))
  conf.scores.stats <- get.mean.and.lh(conf.scores)
  rm(conf.scores)
  
  # match the ccfile to the prior that was used
  p.ind <- which(sapply(names(IN$priors), function(n) grepl(n, cc.file)))
  p.name <- names(IN$priors)[p.ind]
  p.mat <- IN$grouped.pred[[p.name]]$prior.mat[rownames(beta.sign), colnames(beta.sign)]
  
  gp <- IN$grouped.pred[[p.name]]
  cat('A total of', length(gp$pred.has.na), 'predictors contained NA and were removed.\n')
  cat('A total of', length(gp$pred.is.const), 'predictors were constant and were removed.\n')
  cat('A total of', length(unique(unlist(gp$pred.groups))), 'predictors formed', length(gp$pred.groups), 'groups.\n')
  source('R_scripts/evaluate.R')
  cc.aupr <- aupr(comb.confs, p.mat, eval.on.subset=TRUE)
  cat('AUPR is', cc.aupr, '\n')

  # create data frame of non-zero scores 
  net.df <- data.frame(regulator=colnames(beta.sign)[non.zero[, 2]], 
                       target=rownames(beta.sign)[non.zero[, 1]], 
                       beta.sign.sum=beta.sign[vec.ind],
                       beta.non.zero=beta.non.zero[vec.ind] / length(betas),
                       beta.low=betas.stats$low,
                       beta.mean=betas.stats$mean,
                       beta.high=betas.stats$high,
                       var.exp.low=resc.betas.stats$low,
                       var.exp.median=resc.betas.stats$median,
                       var.exp.mean=resc.betas.stats$mean,
                       var.exp.high=resc.betas.stats$high,
                       var.exp.rank.low=resc.betas.ranks.stats$low,
                       var.exp.rank.mean=resc.betas.ranks.stats$mean,
                       var.exp.rank.high=resc.betas.ranks.stats$high,
                       conf.score.low=conf.scores.stats$low,
                       conf.score.mean=conf.scores.stats$mean,
                       conf.score.high=conf.scores.stats$high,
                       prior=p.mat[vec.ind])
  
  iao <- order(-net.df$conf.score.mean, net.df$var.exp.rank.mean, -net.df$var.exp.mean, -abs(net.df$beta.non.zero), net.df$regulator, net.df$target)
  net.df <- net.df[iao, ]

  net.df$precision <- cumsum((net.df$prior != 0)) / (1:length(net.df$prior))

  fn <- gsub('combinedconf', 'summary', cc.file)
  save(net.df, file=fn)
  fn <- gsub('.RData', '.tsv', fn)
  tmp.net.df <- net.df[net.df$beta.non.zero >= 0.5, ]
  for (i in 3:ncol(tmp.net.df)) {
    tmp.net.df[, i] <- round(tmp.net.df[, i], 3)
  }
  write.table(tmp.net.df, file=fn, sep='\t', row.names=FALSE, quote=FALSE)
}

sum.net.cc <- function(cc.file) {
  in.dir <- dirname(cc.file)

  load(cc.file)
  load(gsub('combinedconf', 'betas', cc.file))
  load(paste(in.dir, 'params_and_input.RData', sep='/'))
  
  tmp <- net.stab(betas.resc, N = 2290)
  dev.new()
  plot(tmp[, 1], tmp[, 2])
  browser()
  sum.net(betas, betas.resc, comb.confs, IN, cc.file)
}

# network stability over the bootstraps
net.stab <- function(betas.resc, N=5000) {
  #el <- edge.list(betas.resc, N)
  el <- edge.list.var(betas.resc, 0.2)
  nbs <- length(betas.resc)
  
  x <- 2:nbs
  #overlap <- sapply(x, function(i) length(intersect(el[[nbs+i]], el[[nbs+i-1]])))
  overlap <- sapply(x, function(i) length(intersect(el[[nbs+i]], el[[nbs+i-1]])) / length(el[[nbs+i-1]]))
  return(cbind(x, overlap))
  
  #x <- 1:nbs
  #overlap <- sapply(x, function(i) length(intersect(el[[nbs+i]], el[[nbs+nbs]])))
  #return(cbind(x, overlap))
}

edge.list <- function(betas.resc, N) {
  nbs <- length(betas.resc)
  cat('edge.list;', nbs, 'bootstraps; N =', N, '\n')
  indiv.edges <- list()
  cum.edges <- list()
  cum.ranking <- matrix(0, nrow(betas.resc[[1]]), ncol(betas.resc[[1]]))
  cum.conf.score <- matrix(0, nrow(betas.resc[[1]]), ncol(betas.resc[[1]]))
  for (bs in 1:length(betas.resc)) {
    cat('.')
    bs.mat <- as.matrix(betas.resc[[bs]])
    #cum.ranking <- cum.ranking + rank(-bs.mat)
    #cum.edges[[bs]] <- order(cum.ranking)[1:N]
    cum.conf.score <- cum.conf.score + to.confidence.score(-bs.mat)
    cum.edges[[bs]] <- order(cum.conf.score, decreasing=TRUE)[1:N]
    indiv.edges[[bs]] <- order(bs.mat, decreasing=TRUE)[1:N]
  }
  cat('done\n')
  return(c(indiv.edges, cum.edges))
}

edge.list.var <- function(betas.resc, var.exp.th = 0.2) {
  nbs <- length(betas.resc)
  cat('edge.list.var;', nbs, 'bootstraps; var.exp.th =', var.exp.th, '\n')

  beta.non.zero <- as.matrix(betas.resc[[1]] * 0)
  for (n in 1:nbs) {
    beta.non.zero <- beta.non.zero + (betas.resc[[n]] != 0)
  }
  
  # we only care about interactions that are present in at least one bootstrap
  vec.ind <- which(beta.non.zero != 0)

  resc.betas <- lapply(betas.resc, function(x) x[vec.ind])
  resc.betas <- matrix(unlist(resc.betas), length(resc.betas[[1]]))

  indiv.edges <- list()
  cum.edges <- list()
  for (bs in 1:nbs) {
    cat('.')
    indiv.edges[[bs]] <- which(resc.betas[, bs] > var.exp.th)
    cum.edges[[bs]] <- which(apply(resc.betas[, 1:bs, drop=FALSE], 1, median) > var.exp.th)
  }

  cat('done\n')
  return(c(indiv.edges, cum.edges))
}

#v13.medium <- '/home/ch1421/Projects/Rice/inferelator_output/ALL_1505_upstr_or_0_9_50_adjpth001_incl_autoreg_filtered_16_61224_518_3761_tfa_Peking_plus_lit_BBSR_1_20/combinedconf_frac_tp_100_perm_1--frac_fp_0_perm_1_1.RData'
#sum.net.cc(v13.medium)
#sum.net.cc('output/tmp_201502_final_BBSR_1.1_100_0/combinedconf_frac_tp_100_perm_1--frac_fp_0_perm_1_1.1.RData')
#sum.net.cc('output/bsubtilis_us_201506_stfa_noself3_filtered_newresp_BBSR_1.1_100_0/combinedconf_frac_tp_100_perm_1--frac_fp_0_perm_1_1.1.RData')
#sum.net.cc('output/bsubtilis_eu_201506_stfa_noself3_filtered_newresp_BBSR_1.1_100_0/combinedconf_frac_tp_100_perm_1--frac_fp_0_perm_1_1.1.RData')

#eu.allpr <- 'output/bsubtilis_eu_201506_allpr_noself3_filtered_newresp_BBSR_1.1_100_0/combinedconf_frac_tp_100_perm_1--frac_fp_0_perm_1_1.1.RData'
#us.allpr <- 'output/bsubtilis_us_201502_final_allpr_BBSR_1.1_100_0/combinedconf_frac_tp_100_perm_1--frac_fp_0_perm_1_1.1.RData'

#v21.small.pr <- '/home/ch1421/Projects/Rice/inferelator_output/150622_allmotifs_694bp_open_expr_adjpth001_incl_autoreg_merged005_filtered_32_47842_462_3364_BBSR_1.1_20/combinedconf_frac_tp_100_perm_1--frac_fp_0_perm_1_1.1.RData'
#v22.small.pr <- '/home/ch1421/Projects/Rice/inferelator_output/150622_allmotifs_694bp_open_expr_adjpth001_incl_autoreg_merged005_filtered_32_47842_462_3364_sspr_BBSR_1.1_33/combinedconf_frac_tp_100_perm_1--frac_fp_0_perm_1_1.1.RData'

#v23.small.pr <- '/home/ch1421/Projects/Rice/inferelator_output/150622_allmotifs_694bp_open_expr_adjpth001_incl_autoreg_merged005_filtered_32_47842_462_3364_462tfs_BBSR_1.1_20/combinedconf_frac_tp_100_perm_1--frac_fp_0_perm_1_1.1.RData'
#v24.small.pr <- '/home/ch1421/Projects/Rice/inferelator_output/150622_allmotifs_694bp_open_expr_adjpth001_incl_autoreg_merged005_filtered_32_47842_462_3364__462tfs_sspr_BBSR_1.1_20/combinedconf_frac_tp_100_perm_1--frac_fp_0_perm_1_1.1.RData'
#v24.small.pr.41bs <- '/home/ch1421/Projects/Rice/inferelator_output/150622_allmotifs_694bp_open_expr_adjpth001_incl_autoreg_merged005_filtered_32_47842_462_3364_lesstfs_sspr_BBSR_1.1_41/combinedconf_frac_tp_100_perm_1--frac_fp_0_perm_1_1.1.RData'

#sum.net.cc(v23.small.pr)
#sum.net.cc(v24.small.pr.41bs)

