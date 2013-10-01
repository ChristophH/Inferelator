
dream8.out <- function(scores, cell.type, stimulus, 
                       team.name='Netzwerk', 
                       outdir='~/Projects/DREAM8/submission/') {

  sif.lines <- c()
  eda.lines <- 'EdgeScore'
  
  for (i in 1:ncol(scores)) {
    for (j in 1:nrow(scores)) {
      if (scores[j, i] > 0) {
        sif.lines <- c(sif.lines, paste(colnames(scores)[i], '1', rownames(scores)[j], sep='\t'))
        eda.lines <- c(eda.lines, paste(colnames(scores)[i], ' (1) ', rownames(scores)[j], ' = ', scores[j, i], sep=''))
      }
    }
  }
  degree <- apply(scores, 1, sum)
  degree[1:ncol(scores)] <- degree[1:ncol(scores)] + apply(scores, 2, sum)
  for (i in 1:length(degree)) {
    if (degree[i] == 0) {
      sif.lines <- c(sif.lines, rownames(scores)[i])
    }
  }
  sif.lines <- sort(unique(sif.lines))
  
  if (is.null(cell.type)) {
    base.name <- paste(team.name, 'Network-Insilico', sep='-')
  } else {
    base.name <- paste(team.name, cell.type, stimulus, 'Network', sep='-')
  }
  
  cat(sif.lines, file=paste(outdir, base.name, '.sif', sep=''), sep='\n')
  cat(eda.lines, file=paste(outdir, base.name, '.eda', sep=''), sep='\n')
  
}
