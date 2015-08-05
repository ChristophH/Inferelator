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

reshape.prior <- function(priors.mat, gene.names, tf.names) {
  bg.val <- as.numeric(names(which.max(table(priors.mat))))
  res <- matrix(bg.val, length(gene.names), length(tf.names), 
                dimnames=list(gene.names, tf.names))
  p.genes <- intersect(rownames(priors.mat), gene.names)
  p.tfs <- intersect(colnames(priors.mat), tf.names)
  res[p.genes, p.tfs] <- priors.mat[p.genes, p.tfs]
  return(res)
}

read.input <- function(input.dir, exp.mat.file, tf.names.file, meta.data.file, 
                       priors.file, gold.standard.file, leave.out.file, randomize.expression) {
  IN <- list()
  
  cat('Reading input: expression matrix')
  if (grepl('.RData$', exp.mat.file)) {
    IN$exp.mat <- as.matrix(local(get(load(file.path(input.dir, exp.mat.file)))))
  } else {
    IN$exp.mat <- as.matrix(read.table(file=file.path(input.dir, exp.mat.file),
                            row.names=1, header=T, sep='\t', check.names=F))
  }
  
  cat(', TF names')
  if (grepl('.RData$', tf.names.file)) {
    IN$tf.names <- local(get(load(file.path(input.dir, tf.names.file))))
  } else {
    IN$tf.names <- unique(as.vector(as.matrix(read.table(file.path(input.dir, tf.names.file)))))
  }
  IN$tf.with.expr <- IN$tf.names[IN$tf.names %in% rownames(IN$exp.mat)]

  cat(', meta data')
  IN$meta.data <- NULL
  if (!is.null(meta.data.file)) {
    if (grepl('.RData$', meta.data.file)) {
      IN$meta.data <- local(get(load(file.path(input.dir, meta.data.file))))
    } else {
      IN$meta.data <- read.table(file=file.path(input.dir, meta.data.file), 
                                 header=T, sep='\t', check.names=F)
    }
  }
  
  cat(', leave-out file')
  # if there is a leave-out file, ignore some conditions
  if (!is.null(leave.out.file)) {
    leave.out <- as.vector(as.matrix(read.table(file.path(input.dir, leave.out.file))))
    cat('Leaving out the following conditions:', leave.out, '\n')
    cat(sprintf('%d of %d are present in the expression data\n', 
      sum(leave.out %in% colnames(IN$exp.mat)), length(leave.out)))
    lo <- (IN$meta.data$prevCol %in% leave.out) | (IN$meta.data$condName %in% leave.out)
    cat('Total number of conditions being dropped:', sum(lo), '\n')
    IN$meta.data.lo <- IN$meta.data[lo, ]
    IN$meta.data <- IN$meta.data[!lo, ]
    IN$exp.mat.lo <- IN$exp.mat[, as.character(IN$meta.data.lo$condName)]
    IN$exp.mat <- IN$exp.mat[, as.character(IN$meta.data$condName)]
  }
  
  cat(', priors matrix')
  IN$priors.mat <- NULL
  if (!is.null(priors.file)) {
    if (grepl('.RData$', priors.file)) {
      IN$priors.mat <- as.matrix(local(get(load(file.path(input.dir, priors.file)))))
    } else {
      IN$priors.mat <- as.matrix(read.table(file=file.path(input.dir, priors.file),
                                 row.names=1, header=T, sep='\t', check.names=F))
    }
    IN$priors.mat <- reshape.prior(IN$priors.mat, rownames(IN$exp.mat), IN$tf.names)
  }

  cat(', gold standard matrix')
  IN$gs.mat <- NULL
  if (!is.null(gold.standard.file)) {
    if (grepl('.RData$', gold.standard.file)) {
      IN$gs.mat <- as.matrix(local(get(load(file.path(input.dir, gold.standard.file)))))
    } else {
      IN$gs.mat <- as.matrix(read.table(file=file.path(input.dir, gold.standard.file),
                                    row.names=1, header=T, sep='\t', check.names=F))
    }
    IN$gs.mat <- reshape.prior(IN$gs.mat, rownames(IN$exp.mat), IN$tf.names)
  }
  cat(' ... done.\n')
  
  if (randomize.expression) {
    cat('randomize.expression is set to TRUE; randomizing expression matrix ...')
    IN$exp.mat[sample(length(IN$exp.mat))] <- as.numeric(IN$exp.mat)
    cat(' done.\n')
  }
  return(IN)
}

