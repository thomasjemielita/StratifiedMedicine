globalVariables(c("Rules", "est", "LCL", "UCL"))
#' plot.PRISM
#'
#' Plots PRISM results, either forest plot (estimate with CIs) or resampling distribution.
#'
#' @param x PRISM object
#' @param type Type of plot (default="CI"). Alternatively, if the PRISM object used resampling,
#' type="resample" plots the resampling distributions.
#' @param ... Additional arguments (currently ignored).
#' @return Plot (ggplot2) object
#' @method plot PRISM
#' @export
#' @importFrom stats reorder


plot.PRISM = function(x, type="CI", ...){

  if (type=="CI"){
    # Combine parameter-estimates with rules ##
    parm = x$param.dat
    if (is.null(x$Rules)){
      plot.dat = parm
      plot.dat$Rules = ifelse(plot.dat$Subgrps==0, "Overall",
                              as.character(plot.dat$Subgrps))
    }
    if (!is.null(x$Rules)){
      rules = rbind(data.frame(Subgrps=0,Rules="Overall"), x$Rules)
      plot.dat = left_join(parm, rules, by="Subgrps")
    }

    # Order subgroups in descending order ##
    plot.dat$order = c(-1, rank(plot.dat$est[plot.dat$Subgrps>0]) )
    res = ggplot(data=plot.dat, aes(x=reorder(Rules,order) , y=est, ymin=LCL, ymax=UCL)) +
      geom_pointrange() +
      coord_flip() +
      xlab("Subgroup") + ylab("Estimate (95% CI)") + ggtitle("PRISM Forest Plot") +
      theme(plot.title=element_text(size=16,face="bold"),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.text.x=element_text(face="bold"),
            axis.title=element_text(size=12,face="bold"),
            strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
      theme_bw()
  }
  if (type=="resample"){
    plot.dat = x$resamp.dist
    if (is.null(x$Rules)){
      plot.dat$Rules = ifelse(plot.dat$Subgrps==0, "Overall",
                              as.character(plot.dat$Subgrps))
    }
    if (!is.null(x$Rules)){
      rules = rbind(data.frame(Subgrps=0,Rules="Overall"), x$Rules)
      plot.dat = left_join(plot.dat, rules, by="Subgrps")
    }
    res = ggplot(plot.dat, aes(est)) + geom_density() +
      facet_wrap(~Rules) +
      ggtitle("Bootstrap Distribution of Overall/Subgroup Estimates")+
      theme(plot.title=element_text(size=16,face="bold"),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.text.x=element_text(face="bold"),
            axis.title=element_text(size=12,face="bold"))+
      theme_bw()

  }
  return(res)
}
