library('Matrix')
library('ggplot2')
library('reshape2')
library('gplots')
if ('parallel' %in% installed.packages()[, 'Package']) {
  library('parallel')
} else {
  library('multicore')
}

# for heatmaps
my.hclust <- function(d) hclust(d, method="ward.D")
my.hclust.co <- function(d) hclust(d, method="complete")
my.hclust.av <- function(d) hclust(d, method="average")
dist.p.cor <- function(x) as.dist(1 - cor(t(x), method='pearson'))
dist.p.cor.abs <- function(x) as.dist(1 - abs(cor(t(x), method='pearson')))
dist.s.cor <- function(x) as.dist(1 - cor(t(x), method='spearman'))
my.cols <- colorRampPalette(c("#053061", "#f7f7f7", "#67001f"))(31)
my.cols.BrBG <- colorRampPalette(c("#a6611a", "#f5f5f5", "#018571"))(31)

get.cols <- function(val.min, val.median, val.max, col.min="#a6611a", col.median="#f5f5f5", col.max="#018571", N=32) {
  seg1 <- val.median - val.min
  seg2 <- val.max - val.median
  seg.tot <- val.max - val.min
  steps1 <- round(seg1 / seg.tot * N, 0)
  steps2 <- round(seg2 / seg.tot * N, 0)
  pal1 <- colorRampPalette(c(col.min, col.median))(steps1)
  pal2 <- colorRampPalette(c(col.median, col.max))(steps2)
  pal <- c(pal1, pal2)
  return(pal[!duplicated(pal)])
}

plot.tf <- function(des.mat, res.mat, exp.mat, tf, net, prior, tfa.mat, fn, width=18, height=9, dpi=100) {
  if (file.exists(fn)) return(NULL)
  targets <- rownames(net)
  net.tar <- targets[net[, tf] != 0]
  p.tar <- targets[prior[targets, tf] != 0]
  in.net <- targets[apply(net != 0, 1, sum) > 0]
  
  p.tar <- setdiff(p.tar, tf)
  p.no <- setdiff(p.tar, net.tar)
  novel.tar <- setdiff(net.tar, p.tar)
  p.yes <- intersect(p.tar, net.tar)
  p.no.in.net <- intersect(p.no, in.net)
  
  cat(tf, 'targets', length(net.tar), 'priors', length(p.tar), 'not pred', length(p.no), 'novel', length(novel.tar), 'pred', length(p.yes))
  cat('\n')
  
  if (!is.null(tfa.mat)) {
    des.mat <- tfa.mat
  }
  
  data.df <- melt(t(res.mat[union(net.tar, p.tar), , drop=FALSE]), varnames=c('sample', 'gene'))
  data.df$sample <- factor(data.df$sample, levels=colnames(res.mat), ordered=TRUE)
  data.df$type <- sprintf('novel targets (n=%d)', length(novel.tar))
  data.df$type[data.df$gene %in% p.no] <- sprintf('priors not in net (n=%d)', length(p.no))
  data.df$type[data.df$gene %in% p.yes] <- sprintf('priors in net (n=%d)', length(p.yes))
  
  tf.df <- data.frame(gene=tf, sample=colnames(des.mat), value=des.mat[tf,])
  
  #browser()
  g <- ggplot(data.df, aes(x=sample, y=value, group=gene)) +
    geom_line(aes(color=gene), alpha=0.2, size=2) +
    scale_color_discrete(guide = FALSE) +
    #geom_point(data=tfa.df, aes(x=sample, y=value, group=gene), color='black', size=2)
    geom_line(data=tf.df, aes(x=sample, y=value, group=gene), color='black', size=1)
  
  if (tf %in% rownames(exp.mat)) {
    tf.exp.df <- data.frame(gene=tf, sample=colnames(exp.mat), value=exp.mat[tf,])
    g <- g + geom_line(data=tf.exp.df, aes(x=sample, y=value, group=gene), color='black', size=0.5, linetype=2)
  }
    
  g <- g +
    facet_grid(type ~ .) +
    ylab('Expression') + xlab('Sample') + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = rel(0.6)))
  
  ggsave(fn, g, width=width, height=height, dpi=dpi)
  
  rm(data.df)
  gc()
  
  # for the heatmap
  fn.hm <- gsub('.png', '_heatmap.png', fn)
  goi <- union(net.tar, p.tar)
  if (length(goi) > 1) {
    hm.res.mat <- res.mat[goi, , drop=FALSE]
    cols <- rep('white', length(goi))
    cols[goi %in% novel.tar] <- '#33a02c'  #'#b2df8a'
    cols[goi %in% p.yes] <- '#1f78b4'  #'#a6cee3'
    cols[goi %in% p.no.in.net] <- 'lightgrey'

    vcols <- get.cols(min(hm.res.mat), median(hm.res.mat), max(hm.res.mat))
    png(fn.hm, width=1800, height=900)
    heatmap.2(hm.res.mat, scale='none', trace='none', hclustfun=my.hclust, col=vcols, 
      distfun=dist.p.cor, keysize = 0.66, cexRow=0.4, cexCol=0.3, RowSideColors=cols, Colv=NA,
      symbreaks=FALSE)
    dev.off()
  }
  
  return(NULL) 
  #print(g)
}


plot.TFs.and.targets <- function(net, prior, des.mat, res.mat, exp.mat, out.dir, tfa.mat, CORES) {
  tfs <- colnames(net)
  targets <- rownames(net)
  
  tf.out.deg <- apply(net != 0, 2, sum)
  tf.out.deg.p <- apply(prior[targets, tfs] != 0, 2, sum)
  
  ignore <- mclapply(which(tf.out.deg > 0), function(i) plot.tf(des.mat, res.mat, exp.mat, tfs[i], net, prior, tfa.mat, sprintf('%s/%s_%03d_%03d.png', out.dir, tfs[i], tf.out.deg[i], tf.out.deg.p[i])), mc.cores=CORES)
  return(NULL)
  for (i in which(tf.out.deg > 0)) {
    #if (tfs[i] != 'YBL005W') next
    net.tar <- targets[net[, i] != 0]
    p.tar <- targets[prior[targets, tfs[i]] != 0]

    #png(file=sprintf('%s/%s_%03d_%03d.png', out.dir, tfs[i], tf.out.deg[i], tf.out.deg.p[i]), width=1800, height=900)
    fn <- sprintf('%s/%s_%03d_%03d.png', out.dir, tfs[i], tf.out.deg[i], tf.out.deg.p[i])
    g <- plot.tf(des.mat, res.mat, tfs[i], net, prior, tfa.mat, fn)
    
  }
}


get.net <- function(betas, betas.resc, th, method) {
  if (method == 'frequency') {
    beta.sign <- as.matrix(betas[[1]] * 0)
    for (n in 1:length(betas)) {
      beta.sign <- beta.sign + sign(betas[[n]])
    }
    net <- sign(beta.sign) * (abs(beta.sign) >= length(betas)*th)
  }
  if (method == 'var.exp') {
    beta.non.zero <- as.matrix(betas.resc[[1]] * 0)
    for (n in 1:length(betas.resc)) {
      beta.non.zero <- beta.non.zero + (betas.resc[[n]] != 0)
    }
    vec.ind <- which(beta.non.zero != 0)
    resc.betas <- lapply(betas.resc, function(x) x[vec.ind])
    resc.betas <- matrix(unlist(resc.betas), length(resc.betas[[1]]))
    net <- as.matrix(betas.resc[[1]] * 0)
    net[vec.ind] <- apply(resc.betas, 1, median)
    net[net <= th] <- 0
  }
  return(net)
}


vis.tfs.and.targets <- function(ccfile, CORES=4) {
  in.dir <- dirname(ccfile)
  
  load(ccfile)
  load(gsub('combinedconf', 'betas', ccfile))
  load(paste(in.dir, 'params_and_input.RData', sep='/'))
  
  net <- get.net(betas, betas.resc, 0.2, 'var.exp')
  
  # match the ccfile to the prior that was used
  p.ind <- which(sapply(names(IN$priors), function(n) grepl(n, ccfile)))
  p.name <- names(IN$priors)[p.ind]
  
  out.dir <- file.path(in.dir, sprintf('TF_plots_%s', p.name))
  dir.create(out.dir)
  
  p.mat <- IN$grouped.pred[[p.name]]$prior.mat
  
  des.mat <- IN$grouped.pred[[p.name]]$des.mat  # if TFA was used, this will be activities
  des.mat <- t(scale(t(des.mat)))
  res.mat <- t(scale(t(IN$final_response_matrix)))

  exp.mat <- IN$exp.mat
  exp.mat <- t(scale(t(exp.mat)))
  
  tfa.mat <- NULL
  #if (PARS$use.tfa) {
  #  tfa.mat <- IN$tf.activities[[p.name]]
  #  tfa.mat <- t(scale(t(tfa.mat)))
  #}
  
  plot.TFs.and.targets(net, p.mat, des.mat, res.mat, exp.mat, out.dir, tfa.mat, CORES)
}


#ccfile <- '/home/ch1421/Projects/Kostya/inferelator_output/test_2587_conds_tfa_1_BBSR_1_20/combinedconf_frac_tp_100_perm_1--frac_fp_0_perm_1_1.RData'
#ccfile <- '/home/ch1421/Projects/Rice/inferelator_output/ALL_hughes_plus_150326_v8_incl_autoreg_filtered_32_61101_264_7371_tfa_Peking_plus_lit_BBSR_1_20/combinedconf_frac_tp_100_perm_1--frac_fp_0_perm_1_1.RData'
#ccfile <- '/home/ch1421/Projects/Rice/inferelator_output/ALL_hughes_plus_150327_v8_incl_autoreg_tfa_Peking_plus_lit_BBSR_1_20/combinedconf_frac_tp_100_perm_1--frac_fp_0_perm_1_1.RData'

#v8.small <- '/home/ch1421/Projects/Rice/inferelator_output/ALL_hughes_plus_150416_v8_incl_autoreg_filtered_32_61101_264_7371_tfa_Peking_plus_lit_BBSR_1_20/combinedconf_frac_tp_100_perm_1--frac_fp_0_perm_1_1.RData'
#v8.large <- '/home/ch1421/Projects/Rice/inferelator_output/ALL_hughes_plus_150416_v8_incl_autoreg_tfa_Peking_plus_lit_BBSR_1_20/combinedconf_frac_tp_100_perm_1--frac_fp_0_perm_1_1.RData'
#v8.medium <- '/home/ch1421/Projects/Rice/inferelator_output/ALL_hughes_plus_150416_v8_incl_autoreg_filtered_16_105631_340_13434_tfa_Peking_plus_lit_BBSR_1_20/combinedconf_frac_tp_100_perm_1--frac_fp_0_perm_1_1.RData'

#d5.small <- '/home/ch1421/Projects/Inferelator/output/dream5_net1_small/combinedconf_frac_tp_0_perm_1--frac_fp_0_perm_1_1.RData'

#v9.medium <- '/home/ch1421/Projects/Rice/inferelator_output/ALL_1505_or_0_9_50_incl_autoreg_filtered_16_279593_519_6620_tfa_Peking_plus_lit_BBSR_1_20/combinedconf_frac_tp_100_perm_1--frac_fp_0_perm_1_1.RData'

#v12.medium <- '/home/ch1421/Projects/Rice/inferelator_output/ALL_1505_upstr_or_0_9_50_incl_autoreg_filtered_16_200513_546_6320_tfa_Peking_plus_lit_BBSR_1_20/combinedconf_frac_tp_100_perm_1--frac_fp_0_perm_1_1.RData'

#v13.medium <- '/home/ch1421/Projects/Rice/inferelator_output/ALL_1505_upstr_or_0_9_50_adjpth001_incl_autoreg_filtered_16_61224_518_3761_tfa_Peking_plus_lit_BBSR_1_20/combinedconf_frac_tp_100_perm_1--frac_fp_0_perm_1_1.RData'
#v13.medium.pr <- '/home/ch1421/Projects/Rice/inferelator_output/ALL_1505_upstr_or_0_9_50_adjpth001_incl_autoreg_filtered_16_61224_518_3761_tfa_Peking_plus_lit_BBSR_1.1_20/combinedconf_frac_tp_100_perm_1--frac_fp_0_perm_1_1.1.RData'

#v16.medium.pr <- '/home/ch1421/Projects/Rice/inferelator_output/ALL_150616_upstr_or_0_9_50_adjpth001_incl_autoreg_filtered_16_49869_520_3741_tfa_Peking_plus_lit_BBSR_1.1_20/combinedconf_frac_tp_100_perm_1--frac_fp_0_perm_1_1.1.RData'
#v17.medium.pr <- '/home/ch1421/Projects/Rice/inferelator_output/150616_upstr_or_0_9_50_adjpth001_incl_autoreg_filtered_16_64303_504_4004_BBSR_1.1_20/combinedconf_frac_tp_100_perm_1--frac_fp_0_perm_1_1.1.RData'
#v22.small.pr <- '/home/ch1421/Projects/Rice/inferelator_output/150622_allmotifs_694bp_open_expr_adjpth001_incl_autoreg_merged005_filtered_32_47842_462_3364_sspr_BBSR_1.1_33/combinedconf_frac_tp_100_perm_1--frac_fp_0_perm_1_1.1.RData'

#v23.small.pr <- '/home/ch1421/Projects/Rice/inferelator_output/150622_allmotifs_694bp_open_expr_adjpth001_incl_autoreg_merged005_filtered_32_47842_462_3364_462tfs_BBSR_1.1_20/combinedconf_frac_tp_100_perm_1--frac_fp_0_perm_1_1.1.RData'
#v24.small.pr <- '/home/ch1421/Projects/Rice/inferelator_output/150622_allmotifs_694bp_open_expr_adjpth001_incl_autoreg_merged005_filtered_32_47842_462_3364__462tfs_sspr_BBSR_1.1_20/combinedconf_frac_tp_100_perm_1--frac_fp_0_perm_1_1.1.RData'

#vis.tfs.and.targets(v24.small.pr, CORES=4)

