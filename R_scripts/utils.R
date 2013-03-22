# given an expression matrix, create a trivial cluster stack - no (bi)clusters
trivial.cluster.stack <- function(exp.mat) {
  clusterStack <- list()
  for (i in 1:nrow(exp.mat)) {
    clusterStack[[i]] <- list(cols=NA, ncols=NA, rows=NA, nrows=NA, resid=NA, k=NA, redExp=NA)
    clusterStack[[i]]$cols <- colnames(exp.mat)
    clusterStack[[i]]$ncols <- ncol(exp.mat)
    clusterStack[[i]]$rows <- rownames(exp.mat)[i]
    clusterStack[[i]]$nrows <- 1
    clusterStack[[i]]$k <- i
    clusterStack[[i]]$redExp <- exp.mat[i, ]
  }
  return(clusterStack)
}

# given the condition names, create meta data data frame that assumes all 
# observations are steady state measurements
trivial.meta.data <- function(cond.names) {
  meta.data <- data.frame(condName=cond.names)
  meta.data$isTs <- FALSE
  meta.data$is1stLast <- 'e'
  meta.data$prevCol <- NA
  meta.data$del.t <- NA
  return(meta.data)
}
