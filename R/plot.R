globalVariables(c("Rules", "est", "LCL", "UCL", "PLE", "label", "N", "estimand",
                  "splitvar", "p.value", "surv", "A"))
#' plot.PRISM
#'
#' Plots PRISM results, either forest plot (estimate with CIs) or resampling distribution.
#'
#' @param x PRISM object
#' @param type Type of plot (default="submod", tree plot + parameter estimates). Other 
#' optins include  "forest" ( forest plot for overall and subgroups),"PLE:waterfall" 
#' (waterfall plot of PLEs), "PLE:density" (density plot of PLEs), "resample" (resampling 
#' distribution of parameter estimates for overall and subgroups), and "heatmap" 
#' (heatmap of ple estimates/probabilities). For "submod" and "forest", CIs are based on 
#' the observed data unless bootstrap resampling. If calibrate=TRUE (the default), then
#' calibrated CIs are shown, otherse CIs based on the percentile method are shown.
#' @param estimand For "resample" plot only, must be specify which estimand to visualize.
#' Default=NULL.
#' @param grid.data Input grid of values for 2-3 covariates (if 3, last variable cannot
#' be continuous). This is required for type="heatmap". Default=NULL.
#' @param grid.thres Threshold for PLE, ex: I(PLE>thres). Used to estimate P(PLE>thres) for
#' type="heatmap". Default is ">0". Direction can be reversed and can include equality
#' sign (ex: "<=").
#' @param ... Additional arguments (currently ignored).
#' @return Plot (ggplot2) object
#' @method plot PRISM
#' @export
#' @importFrom stats reorder as.formula
#' @seealso \code{\link{PRISM}}


plot.PRISM = function(x, type="submod", estimand=NULL, grid.data=NULL, grid.thres=">0",
                      ...){

  if (type=="submod"){
    cls = class(x$submod.fit$mod)
    if ("party" %in% cls ){
      res = do.call("plot_submod", list(object=x))
    }
    if (!("party" %in% cls)){
      stop( paste("Plots for custom submod wrappers not currently supported.") )
    }
  }
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
    # Create label: Use calibrated bootstrap interval if available #
    if( !is.null(plot.dat$LCL.calib)  ){
      plot.dat$LCL = plot.dat$LCL.calib
      plot.dat$UCL = plot.dat$UCL.calib
    }
    if (is.null(plot.dat$LCL.calib) & !is.null(plot.dat$LCL.pct)){
      plot.dat$LCL = plot.dat$LCL.pct
      plot.dat$UCL = plot.dat$UCL.pct
    }
    if (x$param=="param_cox"){
      plot.dat$est = exp(plot.dat$est)
      plot.dat$LCL = exp(plot.dat$LCL)
      plot.dat$UCL = exp(plot.dat$UCL)
      plot.dat$estimand = "HR(A=1 vs A=0)"
    }
    plot.dat$label = with(plot.dat, paste( sprintf("%.2f", round(est,2)),
                                           " [",
                                           sprintf("%.2f", round(LCL,2)), ",",
                                           sprintf("%.2f", round(UCL,2)), "]", sep=""))
    # Plot #
    res = ggplot(data=plot.dat, aes(x=estimand, y=est, ymin=LCL, ymax=UCL)) +
      geom_pointrange(aes(col=estimand)) + 
      geom_text(aes(label = label, col=estimand), size=3, 
                position = position_nudge(x = 0.3)) +
      facet_wrap(~Rules, strip.position = "left", nrow = length(unique(plot.dat$Rules)),
                 scales = "free_y") + 
      xlab("Subgroup") + ylab("Estimate (95% CI)") + ggtitle("PRISM Forest Plot") +
      theme_bw() + 
      theme(plot.title=element_text(size=16,face="bold"),
            axis.text.y=element_blank(),
            axis.text.x=element_text(face="bold"),
            axis.title=element_text(size=12,face="bold"),
            strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
      coord_flip()
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
    if (is.null(estimand)){
      stop("Must provide estimand for resampling plot")
    }
    if (x$param=="param_cox"){
      plot.dat$est = exp(plot.dat$est)
      plot.dat$estimand = "HR(A=1 vs A=0)" 
    }
    plot.dat = plot.dat[plot.dat$estimand==estimand,]
    res = ggplot(plot.dat, aes(est)) + geom_density() +
      xlab( paste("Bootstrap Estimates:", estimand)  ) +
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

## PRISM submod plot ##
plot_submod = function(object){
  
  n <- dim(object$out.train)[1]
  alpha_s <- object$alpha_s
  family <- object$family
  # Extract tree fit #
  ct <- object$submod.fit$mod
  # Extract parameter estimates #
  param.dat <- object$param.dat
  if (length(unique(param.dat$Subgrps))==1){
    param.dat = param.dat
    param.dat$Subgrps=1
  }
  if (length(unique(param.dat$Subgrps))>1){
    param.dat <- param.dat[param.dat$Subgrps>0,]
  }
  Subgrps <- unique(param.dat$Subgrps)
  if (object$param=="param_cox"){
    param.dat$est = exp(param.dat$est)
    param.dat$LCL = exp(param.dat$LCL)
    param.dat$UCL = exp(param.dat$UCL)
    param.dat$estimand = "HR(A=1 vs A=0)" 
  }
  param.dat$label <- with(param.dat, paste( sprintf("%.2f", round(est,2)),
                                            " [",
                                            sprintf("%.2f", round(LCL,2)), ",",
                                            sprintf("%.2f", round(UCL,2)), "]", sep=""))
  param.dat$id <- param.dat$Subgrps
  
  # Add estimates into tree #
  ct_node <- as.list(ct$node)
  for (s in unique(param.dat$Subgrps)) {
    ct_node[[s]]$info$label <- param.dat$label[param.dat$Subgrps==s] 
    ct_node[[s]]$info$N <- param.dat$N[param.dat$Subgrps==s] 
    ct_node[[s]]$info$estimand <- param.dat$estimand[param.dat$Subgrps==s] 
    ct_node[[s]]$info$label <- param.dat$label[param.dat$Subgrps==s] 
  }
  # Round Edge Labels #
  for (ii in 1:length(ct_node)) {
    brk = ct_node[[ii]]$split$breaks
    if (!is.null(brk)){
      ct_node[[ii]]$split$breaks <- ifelse(is.numeric(brk), round(brk,2), brk)
    }
  }
  ct$node <- as.partynode(ct_node)
  
  # Create ggparty plot #
  plt <- ggparty(ct, horizontal = TRUE,
                 add_vars = list(N = "$node$info$N",
                                 est = "$node$info$label",
                                 estimand = "$node$info$estimand",
                                 label = "$node$info$label")) +
    geom_edge() +
    geom_edge_label() +
    geom_node_splitvar() + 
    geom_node_label(line_list = list(aes(label = splitvar),
                                     aes(label = pval_convert(p.value))
                    ),
                    line_gpar = list(list(size = 8, col = "black", fontface = "bold"),
                                     list(size = 8)
                    ),
                    ids = "inner")
  if (family!="survival"){
    plt <- plt +
      geom_node_label(
        line_list = list(
          aes(label = paste("Node ", id) ),
          aes(label = paste("N = ", N, " (", round(N/n*100,1), "%)", sep=""))
        ),
        line_gpar = list(list(size = 8, fontface="bold"),
                         list(size = 8) ),
        ids = "terminal") + 
      geom_node_plot(gglist =
                       list(geom_pointrange(data=param.dat,
                                            aes(x=estimand, y=est,
                                                ymin=LCL, ymax=UCL,col=estimand)),
                            geom_text(data=param.dat,
                                      aes(x=estimand, y=est,
                                          label = label, col=estimand), size=3,
                                      position = position_nudge(x = 0.6),
                                      check_overlap = TRUE),
                            scale_x_discrete(expand=c(0.4, 0.4)), 
                            xlab(""),
                            ylab( paste("Estimate (", (1-alpha_s)*100,"% CI)",sep="")),
                            theme_bw(), theme( axis.text.y = element_blank()),
                            coord_flip()),
                     nudge_x = 0.05,
                     shared_axis_labels = TRUE,
                     legend_separator = TRUE, scales = "fixed") 
  }
  if (family=="survival"){
    # Create KM survival prediction data-set #
    out.train = object$out.train
    pred.surv = NULL
    for (s in unique(out.train$Subgrps)){
      hold.dat = out.train[out.train$Subgrps==s,]
      if (is.null(out.train$A)){
        km_mod = survfit(Y ~ 1, data=hold.dat)
        pred.s = data.frame(A = "",
                            time=c(0,km_mod$time), surv=c(1, km_mod$surv))
      }
      if (!is.null(out.train$A)){
        km0 = survfit(Y ~ 1, data=hold.dat[hold.dat$A==0,])
        surv0 = data.frame(A="A=0", time=c(0,km0$time), surv=c(1, km0$surv))
        km1 = survfit(Y ~ 1, data=hold.dat[hold.dat$A==1,])
        surv1 = data.frame(A="A=1", time=c(0,km1$time), surv=c(1, km1$surv))
        pred.s = suppressWarnings( bind_rows(surv0, surv1) )
      }
      pred.s = data.frame(id=s, pred.s)
      pred.surv = suppressWarnings( bind_rows(pred.surv, pred.s) )
    }
    plt <- plt +
      geom_node_label(
        line_list = list(
          aes(label = paste("Node ", id) ),
          aes(label = paste("N = ", N, " (", round(N/n*100,1), "%)", sep="")),
          aes(label = paste(estimand, " [", (1-alpha_s)*100,"% CI]", sep="")),
          aes(label = label)
        ),
        line_gpar = list(list(size = 8, fontface="bold"),
                         list(size = 8),
                         list(size = 8, fontface="bold"),
                         list(size = 8) ),
        ids = "terminal") +
      geom_node_plot(gglist = list(geom_line(data=pred.surv,
                                             aes(x=time, y=surv,
                                                 col=A)),
                                   xlab("Time"),
                                   ylab( "Survival Probability" ),
                                   theme_bw(),
                                   theme(legend.title = element_blank())),
                     nudge_x = 0.10,
                     shared_axis_labels = TRUE,
                     legend_separator = TRUE, scales = "fixed") 
    
  }
  plt <- plt + ggtitle("PRISM Tree Plot")
  return(plt)
}
