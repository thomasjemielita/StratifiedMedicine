globalVariables(c("Rules", "est", "LCL", "UCL", "PLE"))
#' plot.PRISM
#'
#' Plots PRISM results, either forest plot (estimate with CIs) or resampling distribution.
#'
#' @param x PRISM object
#' @param type Type of plot (default="forest", forest plot for overall and subgroups). Other
#' options include "PLE:waterfall" (waterfall plot of PLEs), "PLE:density" (density plot
#' of PLEs), "resample" (resampling distribution of parameter estimates for overall
#' and subgroups), and "heatmap" (heatmap of ple estimates/probabilities).
#' @param grid.data Input grid of values for 2-3 covariates (if 3, last variable cannot
#' be continuous). This is required for type="heatmap". Default=NULL.
#' @param grid.thres Threshold for PLE, ex: I(PLE>thres). Used to estimate P(PLE>thres) for
#' type="heatmap". Default is ">0". Direction can be reversed and can include equality
#' sign (ex: "<=").
#' @param ... Additional arguments (currently ignored).
#' @return Plot (ggplot2) object
#' @method plot PRISM
#' @export
#' @importFrom stats reorder
#' @seealso \code{\link{PRISM}}


plot.PRISM = function(x, type="forest", grid.data=NULL, grid.thres=">0", ...){

  if (type=="forest"){
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
  if (type=="PLE:waterfall"){
    if (is.null(x$mu_train)){
      stop("PLE:Waterfall plot uses PLEs: Check ple argument")
    }
    mu_hat = x$mu_train
    mu_hat$id = 1:nrow(mu_hat)
    res = ggplot(mu_hat, aes(x=reorder(id, PLE), y=PLE)) +
      geom_bar(stat="identity", width=0.7) +
      ggtitle("Waterfall Plot: Patient-level Estimates") +
      ylab("Patient Level Estimates (PLEs)") + xlab("")+
      theme_bw() +
      theme(axis.line.x = element_blank(), axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(face="bold",angle=90))

  }
  if (type=="PLE:density"){
    if (is.null(x$mu_train)){
      stop("PLE:density plot uses PLEs: Check ple argument")
    }
    mu_hat = x$mu_train
    res = ggplot(mu_hat, aes(PLE)) + geom_density() +
      xlab("Patient Level Estimates (PLEs)") +
      ggtitle("Density Plot: Patient-Level Estimates")+
      theme(plot.title=element_text(size=16,face="bold"),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.text.x=element_text(face="bold"),
            axis.title=element_text(size=12,face="bold"))+
      theme_bw()

  }
  if (type=="submod"){
    if (!( "party" %in%  class(x$submod.fit$mod)) ){
      plot(x$submod.fit$mod)
    }
    if ("party" %in% class(x$submod.fit$mod)){
      param.dat = res0$param.dat[res0$param.dat$Subgrps>0, ]
      param.dat$est.CI = with(param.dat, paste(sprintf("%.2f",round(est,2)),
                                               " [", sprintf("%.2f",round(LCL,2)),
                                               ",", sprintf("%.2f",round(UCL,2)),
                                               "]",sep="") )
      ct = x$submod.fit$mod
      plot(ct)
      for(gg in grid.ls(print=F)[[1]]) {
        if (grepl("text", gg)) {
          print(paste(gg, grid.get(gg)$label,sep=": "))
        }
      }
      ct_node <- as.list(ct$node)
      for(i in 1:nrow(param.dat)) {
        ct_node[[param.dat[i,1]]]$info$est <- sprintf("%.2f",round(param.dat$est[i],2))
        ct_node[[param.dat[i,1]]]$info$CI <- paste("[", sprintf("%.2f",round(param.dat$LCL[i],2)),
                                                   ",", sprintf("%.2f",round(param.dat$UCL[i],2)),
                                                   "]",sep="")
        ct_node[[param.dat[i,1]]]$info$prediction <- param.dat$est.CI[i]
      }
      ct$node <- as.partynode(ct_node)
      plot(ct, terminal_panel = node_terminal, tp_args = list(
        FUN = function(node) c("E(Y|A=1)-E(Y|A=0):",
                               node$est, node$CI,
                               paste("N=", node$nobs)) ) )
    }
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
  if (type=="heatmap"){
    if (is.null(x$ple.fit)){
      stop("Heatmap requires ple model fit: Check ple argument")
    }
    if (dim(grid.data)[2]>3 | dim(grid.data)[2]<2){
      stop("Heatmap only applicable for grid.data with 2 or 3 variables")
    }
    if ( dim(grid.data)[2]==3 ){
      if ( !is.factor(grid.data[,3]) | length(unique(grid.data[,3]))>6 ){
        stop("Third column in grid.data should be factor or <=6 unique values")
      }
    }
    # Extract training set covariate space #
    X.train = x$out.train[,!(colnames(x$out.train) %in% c("Y", "A", "Subgrps"))]
    name.var1 = colnames(grid.data)[1]
    name.var2 = colnames(grid.data)[2]
    name.var3 = colnames(grid.data)[3]

    # Create stacked covariate space #
    stack_grid = function(i){
      var1 = grid.data[i,1]
      var2 = grid.data[i,2]
      var3 = grid.data[i,3]
      newdata = X.train
      newdata[name.var1] = var1
      newdata[name.var2] = var2
      if (!is.null(var3)){
        newdata[name.var3] = var3
      }
      newdata$counter = i
      return(newdata)
    }
    X.grid = lapply(1:dim(grid.data)[1], stack_grid)
    X.grid = do.call(rbind, X.grid)
    counter.vec = X.grid$counter
    X.grid = X.grid[,!(colnames(X.grid) %in% "counter")]
    ## Next, predict PLEs across grid ##
    grid.ple = predict(x, newdata = X.grid, type="ple")
    grid.ple$ind.ple = eval(parse(text=paste("ifelse(grid.ple$PLE",
                                             grid.thres, ", 1, 0)")))
    avg.ple = aggregate(grid.ple$PLE ~ counter.vec, FUN="mean")
    prob.ple = aggregate(grid.ple$ind.ple ~ counter.vec, FUN="mean")

    ### Plot Heat-map ###
    est.dat = data.frame(grid.data, est = avg.ple$`grid.ple$PLE`,
                         prob = prob.ple$`grid.ple$ind.ple`)
    res.est = ggplot(data = est.dat, aes_string(x=name.var1, y=name.var2, fill="est")) +
      geom_tile() + labs(fill = "PLE") +
      scale_fill_gradient2(low="navy", mid="white", high="red")
    res.prob = ggplot(data = est.dat, aes_string(x=name.var1, y=name.var2, fill="prob")) +
      geom_tile() + labs(fill = paste("Prob(PLE", grid.thres, ")",sep="")) +
      scale_fill_gradient2(low="navy", mid="white", high="red",
                           midpoint=0.5)
    # Lastly, if there were 3 input variables, facet_wrap by third variable #
    if ( dim(grid.data)[2]==3  ){
      res.est = res.est + facet_wrap(as.formula(paste("~", name.var3)))
      res.prob = res.prob + facet_wrap(as.formula(paste("~", name.var3)))
    }
    res = list(heatmap.est=res.est, heatmap.prob=res.prob)
  }
  return(res)
}
