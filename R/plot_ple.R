### PLE Waterfall ###
plot_ple_waterfall <- function(x, ple.label) {
  
  if (is.null(x$mu_train)){
    stop("PLE:Waterfall plot uses PLEs: Check ple argument")
  }
  y.label <- paste("Estimates:", ple.label)
  mu_hat = x$mu_train
  mu_hat$id = 1:nrow(mu_hat)
  mu_hat = data.frame(mu_hat, Subgrps = factor(x$out.train$Subgrps) )
  res = ggplot2::ggplot(mu_hat, aes(x=reorder(id, PLE), y=PLE, fill=Subgrps)) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::ggtitle( paste("Waterfall Plot: Patient-level Estimates,", ple.label) ) +
    ggplot2::ylab(y.label) + 
    ggplot2::xlab("")+
    ggplot2::theme_bw() +
    ggplot2::theme(axis.line.x = element_blank(), axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.title.y = element_text(face="bold",angle=90))
  return(res)
}
### PLE Density ###
plot_ple_density <- function(x, ple.label) {
  
  if (is.null(x$mu_train)){
    stop("PLE:density plot uses PLEs: Check ple argument")
  }
  x.label <- paste("Estimates:", ple.label)
  mu_hat = x$mu_train
  mu_hat = data.frame(mu_hat, Subgrps = factor(x$out.train$Subgrps) )
  res = ggplot2::ggplot(mu_hat, aes(PLE, fill=Subgrps)) + 
    ggplot2::geom_density(alpha=0.30) +
    ggplot2::xlab(x.label) +
    ggplot2::ggtitle( paste("Density Plot: Patient-Level Estimates,", ple.label))+
    ggplot2::theme(plot.title=element_text(size=16,face="bold"),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank(),
                   axis.text.x=element_text(face="bold"),
                   axis.title=element_text(size=12,face="bold"))+
    ggplot2::theme_bw()
  return(res)
}