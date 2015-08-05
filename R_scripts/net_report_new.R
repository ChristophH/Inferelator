
net.report <- function(cc.file) {
  require('markdown')
  require('knitr')

  cc.dir <- normalizePath(dirname(cc.file))
  
  mdf <- normalizePath('R_scripts/net_report.Rmd')
  kmdf <- file.path(cc.dir, 'net_report.md')
  htf <- file.path(cc.dir, 'net_report.html')
  
  # knit
  opts_knit$set(base.dir=cc.dir, width=112)
  opts_chunk$set(fig.path=file.path(cc.dir, 'report_figures/'), fig.width=14, fig.height=8)
  knit(mdf, kmdf)
  
  # render markdown
  markdownToHTML(kmdf, htf, options=c('base64_images', 'highlight_code'))
  file.remove(kmdf)
}
