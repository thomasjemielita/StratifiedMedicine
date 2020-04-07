### Forest Plot ###
plot_forest <- function(x) {
  
  # Combine parameter-estimates with rules ##
  parm <- x$param.dat
  if (is.null(x$Rules)){
    plot.dat <- parm
    plot.dat$Rules <- ifelse(plot.dat$Subgrps==0, "Overall",
                             as.character(plot.dat$Subgrps))
  }
  if (!is.null(x$Rules)){
    rules <- rbind(data.frame(Subgrps=0,Rules="Overall"), x$Rules)
    plot.dat <- left_join(parm, rules, by="Subgrps")
  }
  plot.dat$label = with(plot.dat, paste( sprintf("%.2f", round(est0,2)),
                                         " [",
                                         sprintf("%.2f", round(LCL0,2)), ",",
                                         sprintf("%.2f", round(UCL0,2)), "]", sep=""))
  # Plot #
  res = ggplot2::ggplot(data=plot.dat, aes(x=estimand, y=est0, ymin=LCL0, ymax=UCL0)) +
    ggplot2::geom_pointrange(aes(col=estimand)) + 
    ggplot2::geom_text(aes(label = label, col=estimand), size=3, 
                       position = position_nudge(x = 0.3)) +
    ggplot2::facet_wrap(~Rules, strip.position = "left", 
                        nrow = length(unique(plot.dat$Rules)),
                        scales = "free_y") + 
    ggplot2::xlab("Subgroup") + ylab("Estimate (95% CI)") + ggtitle("PRISM Forest Plot") +
    ggplot2::theme_bw() + 
    ggplot2::theme(plot.title=element_text(size=16,face="bold"),
                   axis.text.y=element_blank(),
                   axis.text.x=element_text(face="bold"),
                   axis.title=element_text(size=12,face="bold"),
                   strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
    ggplot2::coord_flip()
  return(res)
} 