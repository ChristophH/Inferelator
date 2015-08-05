
pretty.much.the.same <- function(x, tol=.Machine$double.eps) {
  return(all(abs(x - x[1]) < tol))
}

# get all nodes that can be reached from cur.node
get.cc <- function(cur.node, adj.mat, visited=c()) {
  visited <- union(visited, cur.node)
  new.adj.nodes <- setdiff(which(adj.mat[cur.node, ] > 0), visited)
  for (a.node in new.adj.nodes) {
    visited <- union(visited, get.cc(a.node, adj.mat, visited))
  }
  return(visited)
}

# get all connected components (with size > 1)
get.all.cc <- function(adj.mat) {
  N <- nrow(adj.mat)
  in.list <- rep(FALSE, N)
  cc.list <- list()
  for (i in 1:N) {
    if (!in.list[i]) {
      cc <- get.cc(i, adj.mat)
      if (length(cc) > 1) {
        cc.list[[length(cc.list) + 1]] <- cc
      }
      in.list[cc] <- TRUE
    }
  }
  return(cc.list)
}

group.prior <- function(groups, prior.mat, tfs.keep) {
  if (is.null(prior.mat)) {
    return(prior.mat)
  }
  prior.mat.grp <- matrix(0, nrow(prior.mat), length(groups))
  colnames(prior.mat.grp) <- names(groups)
  g.n <- 1
  for (group in groups) {
    #when consolidating the prior, check for inconsistencies
    tmp.p.mat <- prior.mat[, group]
    tmp.p.mat.incons <- any(apply(tmp.p.mat, 1, function(x) any(x > 0) & any(x < 0)))
    if (tmp.p.mat.incons) {
      stop(sprintf('inconsistency when consolidating the prior/GS in group.predictors for group %d', g.n))
    }
    prior.mat.grp[, g.n] <- sign(apply(tmp.p.mat, 1, sum))
    g.n <- g.n + 1
  }
  new.prior.mat <- prior.mat[, sort(setdiff(tfs.keep, unlist(groups)))]
  new.prior.mat <- cbind(new.prior.mat, prior.mat.grp)
  return(new.prior.mat)
}

group.expression <- function(groups, exp.mat) {
  grp.exp.mat <- matrix(0, length(groups), ncol(exp.mat), 
    dimnames=list(names(groups), colnames(exp.mat)))
  for (group in names(groups)) {
    genes <- intersect(groups[[group]], rownames(exp.mat))
    grp.exp.mat[group, ] <- apply(scale(t(exp.mat[genes, , drop=FALSE])), 1, mean)
  }
  return(grp.exp.mat)
}

# des.mat is TFs x conditions matrix
# bs.pi is bs x conditions matrix of indices of conditions per bootstrap
group.predictors <- function(des.mat, prior.mat, gs.mat, bs.pi, cor.th=0.99,
                             grp.pre='pred.group.') {
  tfs <- rownames(des.mat)
  have.na <- apply(is.na(des.mat), 1, sum) > 0
  const <- rep(FALSE, nrow(des.mat))
  const[!have.na] <- apply(des.mat[!have.na, ], 1, pretty.much.the.same)
  
  # of the TFs that have no NA and are not constant, which ones are highly correlated
  # in any of the sub-samples of the data
  keep <- !have.na & !const
  tfs.keep <- tfs[keep]
  cor.mat <- matrix(0, sum(keep), sum(keep))
  for (i in 1:nrow(bs.pi)) {
    bs.cor.mat <- cor(t(des.mat[keep, bs.pi[i, ]]))
    cor.mat <- cor.mat + (bs.cor.mat > cor.th)
  }
  cc <- get.all.cc(cor.mat)
  groups <- lapply(cc, function(x) tfs.keep[x])
  
  if (length(groups) > 0) {
    names(groups) <- sprintf('%s%d', grp.pre, 1:length(groups))
  }
  
  # apply the grouping to design matrix, exclude NA containing and constant predictors
  des.mat.grp <- matrix(0, length(groups), ncol(des.mat))
  rownames(des.mat.grp) <- names(groups)
  g.n <- 1
  for (group in groups) {
    exemplar <- group[which.max(apply(des.mat[group, ], 1, function(x) sd(x)/mean(x)))]
    des.mat.grp[g.n, ] <- des.mat[exemplar, ]
    g.n <- g.n + 1
  }
  new.des.mat <- des.mat[sort(setdiff(tfs.keep, unlist(groups))), ]
  new.des.mat <- rbind(new.des.mat, des.mat.grp)

  # also adjust the prior matrix and gold standard to the grouped predictors
  new.prior.mat <- group.prior(groups, prior.mat, tfs.keep)
  new.gs.mat <- group.prior(groups, gs.mat, tfs.keep)
  
  
  # if no TF was excluded and no TFs were grouped, return original design matrix
  # and original prior matrix
  if (identical(sort(rownames(des.mat)), sort(rownames(new.des.mat)))) {
    new.des.mat <- des.mat
    new.prior.mat <- prior.mat
    new.gs.mat <- gs.mat
  }
  
  return(list(pred.has.na=tfs[have.na], pred.is.const=tfs[const], 
              pred.groups=groups, des.mat=new.des.mat, prior.mat=new.prior.mat,
              gs.mat=new.gs.mat, tf.names=rownames(new.des.mat)))
}



