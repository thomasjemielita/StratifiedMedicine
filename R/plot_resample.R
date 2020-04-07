### Resampling ###
plot_resample <- function(x, estimand=NULL) {
  plot.dat = x$resamp.dist
  if (is.null(x$Rules)) {
    plot.dat$Rules = ifelse(plot.dat$Subgrps==0, "Overall",
                            as.character(plot.dat$Subgrps))
  }
  if (!is.null(x$Rules)) {
    rules = rbind(data.frame(Subgrps=0,Rules="Overall"), x$Rules)
    plot.dat = left_join(plot.dat, rules, by="Subgrps")
  }
  if (is.null(estimand)) {
    if (x$family %in% c("gaussian", "binomial")) {
      estimand <- "E(Y|A=1)-E(Y|A=0)"
    }
    if (x$family=="survival") {
      if (x$param=="param_cox") {
        estimand <- "HR(A=1 vs A=0)" 
      }
      if (x$param=="param_rmst") {
        estimand <- "RMST(A=1-A=0)" 
      }
    }
  }
  if (x$param=="param_cox") {
    plot.dat$est = exp(plot.dat$est)
    plot.dat$estimand = "HR(A=1 vs A=0)" 
  }
  plot.dat = plot.dat[plot.dat$estimand==estimand,]
  res = ggplot2::ggplot(plot.dat, aes(est)) + 
    ggplot2::geom_density() +
    ggplot2::xlab( paste("Bootstrap Estimates:", estimand)  ) +
    ggplot2::facet_wrap(~Rules) +
    ggplot2::ggtitle("Bootstrap Distribution of Overall/Subgroup Estimates")+
    ggplot2::theme(plot.title=element_text(size=16,face="bold"),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank(),
                   axis.text.x=element_text(face="bold"),
                   axis.title=element_text(size=12,face="bold"))+
    theme_bw()
  return(res)
}