globalVariables(c("Rules", "est", "LCL", "UCL", "PLE", "label", "N", "estimand",
                  "splitvar", "p.value", "surv", "A", "Subgrps",
                  "LCL0", "UCL0", "est0", "events", "prob.est",
                  "x", "y"))
#' plot.PRISM
#'
#' Plots PRISM results. Options include "tree", "forest", "resample", and "PLE:waterfall".
#'
#' @param x PRISM object
#' @param type Type of plot (default="tree", tree plot + parameter estimates/outcome and 
#' or probability density plots). Other options include  "forest" 
#' (forest plot for overall and subgroups),"PLE:waterfall" 
#' (waterfall plot of PLEs), "PLE:density" (density plot of PLEs), "resample" (resampling 
#' distribution of parameter estimates for overall and subgroups), and "heatmap" 
#' (heatmap of ple estimates/probabilities). For "tree" and "forest", CIs are based on 
#' the observed data unless resampling is used. For bootstrap resampling, if 
#' calibrate=TRUE, then calibrated CIs along are shown, otherse CIs based on the 
#' percentile method are shown.
#' @param estimand For "resample" plot only, must be specify which estimand to visualize.
#' Default=NULL.
#' @param grid.data Input grid of values for 2-3 covariates (if 3, last variable cannot
#' be continuous). This is required for type="heatmap". Default=NULL.
#' @param grid.thres Threshold for PLE, ex: I(PLE>thres). Used to estimate P(PLE>thres) for
#' type="heatmap". Default is ">0". Direction can be reversed and can include equality
#' sign (ex: "<=").
#' @param tree.plots Type of plots to include in the "tree" plot. Default="outcome"
#' (boxplots of treatment-specific outcomes, or counterfactual estimates if PLE!=NULL).
#' For "density", the estimated probability density of the treatment effects is shown 
#' (normal approximation, unless resampling is used). "both" combines both plots.
#' @param tree.thres Probability threshold, ex: P(Mean(A=1 vs A=0)>c. Default=NULL, 
#' which defaults to using ">0", unless param="param_cox", which  "P(HR(A=1 vs A=0))<1". 
#' If a density plot is included, setting tree.thres=">c" will use green colors 
#' for values above c, and red colors for values below c. If tree.thres="<c", the 
#' reverse color scheme is used.
#' @param est.resamp Should plot present resampling based estimates? Default=TRUE if 
#' bootstrap or CV  based resampling is used. Only applicable for type="submod". 
#' If bootstrap calibration is used, calibrated CIs are presented. If no calibration,
#' then percentile Cis are presented with the smoothed bootstrap point-estimates.
#' @param nudge_out Nudge tree outcome plot (see ggparty for details)
#' @param nudge_dens Nudge tree density plot
#' @param width_out Width of tree outcome plot (see ggparty for details)
#' @param width_dens Width of density tree outcome plot
#' @param ... Additional arguments (currently ignored).
#' @return Plot (ggplot2) object
#' @method plot PRISM
#' @export
#' @importFrom stats reorder as.formula density
#' @seealso \code{\link{PRISM}}


plot.PRISM = function(x, type="tree", estimand=NULL, grid.data=NULL, grid.thres=">0",
                      tree.thres=NULL,
                      est.resamp=TRUE, tree.plots="outcome",
                      nudge_out=0.1, width_out=0.5,
                      nudge_dens=ifelse(tree.plots=="both", 0.3, 0.1),
                      width_dens=0.5, ...) {
  # Default Setup #
  param.dat <- x$param.dat
  param.dat$est0 <- param.dat$est
  param.dat$SE0 <- param.dat$SE
  param.dat$LCL0 <- param.dat$LCL
  param.dat$UCL0 <- param.dat$UCL
  param.dat$prob.est <- param.dat$`Prob(>0)`
  resample <- ifelse(is.null(x$resample), "None", x$resample)
  bayes <- ifelse(is.null(x$bayes.fun), FALSE, TRUE)
  label.param <- ""
  if (est.resamp & (resample %in% c("Bootstrap", "CV"))) {
    if (resample=="Bootstrap" & is.null(param.dat$LCL.calib)) {
      param.dat$est0 <- param.dat$est_resamp
      param.dat$LCL0 <- param.dat$LCL.pct
      param.dat$UCL0 <- param.dat$UCL.pct
      label.param <- "(Boot,Pct)"
    }
    if (resample=="Bootstrap" & !is.null(param.dat$LCL.calib)) {
      param.dat$LCL0 <- param.dat$LCL.calib
      param.dat$UCL0 <- param.dat$UCL.calib
      label.param <- "(Boot,Calib)"
    }
    if (resample=="CV"){
      param.dat$est0 <- param.dat$est_resamp
      param.dat$LCL0 <- param.dat$LCL.CV
      param.dat$UCL0 = param.dat$UCL.CV
      label.param <- "(CV)"
    }
  }
  if (bayes) {
    param.dat$est0 <- param.dat$est.bayes
    param.dat$SE0 <- param.dat$SE.bayes
    param.dat$LCL0 <- param.dat$LCL.bayes
    param.dat$UCL0 <- param.dat$UCL.bayes
  }
  if (x$family=="survival") {
    if (x$param=="param_cox") {
      param.dat$est0 = exp(param.dat$est0)
      param.dat$LCL0 = exp(param.dat$LCL0)
      param.dat$UCL0 = exp(param.dat$UCL0)
      param.dat$estimand = gsub("logHR", "HR", param.dat$estimand)
      param.dat$prob.est = 1-param.dat$`Prob(>0)`
    }
  }
  x$param.dat <- param.dat
  
  if (type=="submod" & length(unique(x$out.train$Subgrps))==1) {
    message("No Subgroups found: Forest Plot is for Overall Population Only")
    type = "forest"
  }
  if (type=="tree"){
    if (!is.null(tree.thres)) {
      thres.name <- paste("Prob(",tree.thres, ")", sep="")
      x2 <- prob_calculator(x, thres=tree.thres)
      colnames(x2$param.dat)[which(colnames(x2$param.dat)==thres.name)] <- "prob.est"
    }
    if (is.null(tree.thres)) {
      x2 <- x
      tree.thres <- ifelse(x2$param=="param_cox", "<1", ">0")
    }
    cls <- class(x2$submod.fit$mod)
    if ("party" %in% cls) {
      res <- do.call("plot_tree", list(object=x2, plots=tree.plots,
                                       prob.thres = tree.thres,
                                       nudge_out=nudge_out, width_out=width_out,
                                       nudge_dens=nudge_dens, width_dens=width_dens))
    }
    if (!("party" %in% cls)) {
      stop( paste("Plots for non partykit models not currently supported.") )
    }
  }
  if (type=="forest"){
    res <- plot_forest(x)
  }
  # PLE Label #
  ple.label <- ""
  if (x$family %in% c("gaussian", "binomial")) {
    ple.label <- "E(Y|A=1,X)-E(Y|A=0,X)"
  }
  if (x$family %in% "survival") {
    if (x$ple %in% c("ple_ranger", "ple_rfsrc")){
      ple.label <- "RMST(A=1 vs A=0)"
    }
    if (x$ple %in% c("ple_bart")) {
      ple.label <- "ETR(A=1 vs A=0)"
    }
    if (x$ple %in% c("ple_glmnet")) {
      ple.label <- "logHR(A=1 vs A=0)"
    }
  }
  if (type=="PLE:waterfall") {
    res <- plot_ple_waterfall(x=x, ple.label=ple.label)
  }
  if (type=="PLE:density") {
    res <- plot_ple_density(x=x, ple.label=ple.label)
  }
  if (type=="resample") {
    res <- plot_resample(x=x, estimand=estimand)
  }
  if (type=="heatmap"){
    res <- plot_heatmap(x=x, grid.data=grid.data, grid.thres=grid.thres)
  }
  return(res)
}
## PRISM submod plot (with resampling) ##
plot_tree = function(object, plots, prob.thres, width_out, nudge_out,
                     width_dens, nudge_dens) {
  # Extract #
  n <- dim(object$out.train)[1]
  family <- object$family
  prob.label <- paste("Prob(",prob.thres, ")", sep="")
  if (family=="survival"){
    n.events <- object$param.dat$events[object$param.dat$Subgrps==0] 
  }
  alpha_s <- object$alpha_s
  resample <- ifelse(is.null(object$resample), "None", object$resample)
  ct <- object$submod.fit$mod
  # Parameter Estimates #
  param.dat <- object$param.dat
  Subgrps <- unique(param.dat$Subgrps[param.dat$Subgrps>0])
  
  param.dat$label <- with(param.dat, paste( sprintf("%.2f", round(est0,2)),
                                            " [",
                                            sprintf("%.2f", round(LCL0,2)), ",",
                                            sprintf("%.2f", round(UCL0,2)), "]", sep=""))
  param.dat$prob.est <- with(param.dat, sprintf("%.2f", round(prob.est,2)))
  param.dat$id <- param.dat$Subgrps
  param.ovrl <- param.dat[param.dat$Subgrps==0,]
  param.subs <- param.dat[param.dat$Subgrps>0,]
  estimand <- unique(param.dat$estimand)
  
  # A levels (if not NULL) #
  if (!is.null(object$out.train$A)) {
    A_lvls <- with(object$out.train, unique(A)[order(unique(A))])
  }

  if (family!="survival" & !is.null(object$out.train$A)) {
    E_A0 <- paste("E(Y|A=", A_lvls[1], ")", sep="")
    E_A1 <- paste("E(Y|A=", A_lvls[2], ")", sep="")
    E_diff <- paste(E_A1, "-", E_A0, sep="")
    mu_A0 <- paste("mu", A_lvls[1], sep="_")
    mu_A1 <- paste("mu", A_lvls[2], sep="_")
    param.subs <- param.dat[param.subs$estimand==E_diff,]
    param.subs$estimand <- as.character(param.subs$estimand)
    estimand <- E_diff
    if (object$ple!="None") {
      mu_hat <- object$mu_train
      plot.dat <- rbind( data.frame(estimand=E_A0, est=mu_hat[,mu_A0], 
                                    Subgrps = object$out.train$Subgrps),
                         data.frame(estimand=E_A1, est=mu_hat[,mu_A1],
                                    Subgrps = object$out.train$Subgrps))
      plot.dat$id <- plot.dat$Subgrps
    }
    if (object$ple=="None") {
      plot.dat <- object$out.train[,colnames(object$out.train) %in% 
                                     c("Y", "A", "Subgrps")]
      colnames(plot.dat)[1] <- c("est")
      plot.dat$estimand <- with(plot.dat, paste("E(Y|A=",A,")",sep=""))
      plot.dat$estimand <- factor(plot.dat$estimand, levels=c(E_A0, E_A1))
      plot.dat$id <- plot.dat$Subgrps
    }
  }
  # If no Treatment #
  if (is.null(object$out.train$A)) {
    plot.dat <- object$out.train[,colnames(object$out.train) %in% 
                                   c("Y", "Subgrps")]
    colnames(plot.dat)[1] <- c("est")
    plot.dat$estimand <- "E(Y)"
    plot.dat$id <- plot.dat$Subgrps
  }
  
  # Add estimates into tree #
  ct_node <- as.list(ct$node)
  for (s in Subgrps) {
    ct_node[[s]]$info$label <- param.subs$label[param.subs$Subgrps==s] 
    ct_node[[s]]$info$N <- param.subs$N[param.subs$Subgrps==s] 
    ct_node[[s]]$info$estimand <- param.subs$estimand[param.subs$Subgrps==s]
    ct_node[[s]]$info$prob.est <- param.subs$prob.est[param.subs$Subgrps==s]
    if (family=="survival") {
      ct_node[[s]]$info$events <- param.subs$events[param.subs$Subgrps==s] 
    }
  }
  # Round Edge Labels #
  for (ii in 1:length(ct_node)) {
    brk = ct_node[[ii]]$split$breaks
    if (!is.null(brk)){
      ct_node[[ii]]$split$breaks <- ifelse(is.numeric(brk), round(brk,2), brk)
    }
  }
  ct$node <- as.partynode(ct_node)
  
  ## Create density data ##
  post.prob <- object$resamp.dist
  if (is.null(post.prob)) {
    post.prob <- NULL
    if (object$param=="param_cox") {
      param.subs$est0 <- log(param.subs$est0)
    }
    for (s in unique(Subgrps)) {
      dat.s <- rnorm(10000, mean = param.subs$est0[param.subs$Subgrps==s],
                     sd = param.subs$SE0[param.subs$Subgrps==s])
      dat.s <- data.frame(Subgrps=s, est=dat.s)
      post.prob <- rbind(post.prob, dat.s)
    }
  }
  post.prob <- post.prob[post.prob$Subgrps>0,]
  if (family=="survival" & object$param=="param_cox") {
    post.prob$est <- exp(post.prob$est)
  }
  dat.dens <- NULL
  for (s in unique(post.prob$Subgrps)){
    hold.s = post.prob %>% filter(Subgrps==s)
    dat.s = with(density(hold.s$est), data.frame(x, y))
    dat.s = data.frame(id = s, dat.s)
    dat.dens = rbind(dat.dens, dat.s)
  }
  max.dens <- max(dat.dens$y)*1.05
  
  ## Plot setup ##
  add_vars_list <- list(N = "$node$info$N", est = "$node$info$label", 
                        estimand = "$node$info$estimand", 
                        prob.est = "$node$info$prob.est",
                        label = "$node$info$label")
  if (family=="survival"){
    add_vars_list[["events"]] = "$node$info$events"
  }
  
  plt <- ggparty(ct, horizontal = TRUE, add_vars = add_vars_list) +
    geom_edge() +
    geom_edge_label() +
    geom_node_label(line_list = list(aes(label = splitvar),
                                     aes(label = pval_convert(p.value))),
                    line_gpar = list(list(size = 8, col = "black", fontface = "bold"),
                                     list(size = 8)),ids = "inner")
  # Set up prob cutoff params #
  dir <- substr(prob.thres, 1, 1)
  thres <- substr(prob.thres, 2, nchar(prob.thres))
  if (suppressWarnings(is.na(as.numeric(thres)))) {
    thres <- substr(prob.thres, 3, nchar(prob.thres))
  }
  thres <- as.numeric(thres)
  fill_L <- ifelse(dir=="<", "green", "red")
  fill_R <- ifelse(dir==">", "green", "red")
  ## Continuous/Binary Plots ##
  if (family!="survival") {
    # Add Node Info #
    plt <- plt + geom_node_label(
      line_list = list(
        aes(label = paste("Node ", id) ),
        aes(label = paste("N = ", N, " (", round(N/n*100,1), "%)", sep="")),
        aes(label = paste(estimand, " [", (1-alpha_s)*100,"% CI]", sep="")),
        aes(label = label),
        aes(label = paste(prob.label, "=", prob.est))
      ),
      line_gpar = list(list(size = 8, fontface="bold"),
                       list(size = 8),
                       list(size = 8, fontface="bold"),
                       list(size = 8),
                       list(size = 8)),
      ids = "terminal")
    # Outcome Plot #
    out_plot <- geom_node_plot(gglist =
                                 list(geom_boxplot(data=plot.dat,
                                                   aes(x=estimand, y=est, fill=estimand)),
                                      scale_x_discrete(expand=c(0.4, 0.4)), 
                                      xlab(""),
                                      ylab("Estimate"),
                                      theme_bw(), theme( axis.text.y = element_blank()),
                                      coord_flip()),
                               nudge_x = nudge_out, width = width_out,
                               shared_axis_labels = TRUE,
                               legend_separator = TRUE, scales = "fixed")
    # Density Plot #
    dens_plot <- geom_node_plot(
      gglist = list( geom_area(data=dat.dens, 
                               aes(x = ifelse(x < thres , x, thres),y=y), fill = fill_L),
                     geom_area(data=dat.dens, aes(x = ifelse(x > thres , x, thres),
                                                  y=y), fill = fill_R),
                     geom_line(data=dat.dens, aes(x=x,y=y)),
                     ylim(0,max.dens),
                     theme_bw(),
                     ggtitle(estimand),
                     theme(plot.title = element_text(size = 8, face = "bold")),
                     xlab("Density"), ylab("") ), shared_axis_labels = TRUE,
      width=width_dens, nudge_x = nudge_dens)
    # Outcome Only #
    if (plots=="outcome") {
      plt.out <- plt + out_plot 
    }
    if (plots=="density") {
      plt.out <- plt + dens_plot
    }
    if (plots=="both") {
      plt.out <- plt + out_plot + dens_plot
    }
  }
  ### Survival Plot ###
  if (family=="survival") {
    # Add Node Info #
    plt <- plt +
      geom_node_label(
        line_list = list(
          aes(label = paste("Node ", id) ),
          aes(label = paste("N = ", N, " (", round(N/n*100,1), "%)", sep="")),
          aes(label = paste("Events = ", events,
                            " (", round(events/n.events*100,1), "%)", sep="")),
          aes(label = paste(estimand, " [", (1-alpha_s)*100,"% CI]", sep="")),
          aes(label = label),
          aes(label = paste(prob.label, "=", prob.est))
        ),
        line_gpar = list(list(size = 8, fontface="bold"),
                         list(size = 8),
                         list(size = 8),
                         list(size = 8, fontface="bold"),
                         list(size = 8), 
                         list(size = 8)),
        ids = "terminal")
    # Create KM survival prediction data-set #
    out.train <- object$out.train
    pred.surv <- NULL
    for (s in unique(out.train$Subgrps)) {
      hold.dat <- out.train[out.train$Subgrps==s,]
      if (is.null(out.train$A)){
        km_mod <- survfit(Y ~ 1, data=hold.dat)
        pred.s <- data.frame(A = "",
                            time=c(0,km_mod$time), surv=c(1, km_mod$surv))
      }
      if (!is.null(out.train$A)) {
        km0 <- survfit(Y ~ 1, data=hold.dat[hold.dat$A==A_lvls[1],])
        surv0 <- data.frame(A=A_lvls[1], time=c(0,km0$time), surv=c(1, km0$surv))
        km1 <- survfit(Y ~ 1, data=hold.dat[hold.dat$A==A_lvls[2],])
        surv1 <- data.frame(A=A_lvls[2], time=c(0,km1$time), surv=c(1, km1$surv))
        pred.s <- suppressWarnings(bind_rows(surv0, surv1))
      }
      pred.s <- data.frame(id=s, pred.s)
      pred.surv <- suppressWarnings( bind_rows(pred.surv, pred.s) )
    }
    out_plot <- geom_node_plot(gglist = list(geom_line(data=pred.surv,
                                                       aes(x=time, y=surv,
                                                           col=A), size=1.5),
                                             xlab("Time"),
                                             ylab( "Survival Probability" ),
                                             theme_bw(),
                                             theme(legend.title = element_blank())),
                               shared_axis_labels = TRUE,
                               legend_separator = TRUE, scales = "fixed", 
                               width=width_out, nudge_x=nudge_out)
    # Density #
    dens_plot <- geom_node_plot(
      gglist = list( geom_area(data=dat.dens, 
                               aes(x = ifelse(x < thres , x, thres),y=y), fill = fill_L),
                     geom_area(data=dat.dens, aes(x = ifelse(x > thres , x, thres),
                                                  y=y), fill = fill_R),
                     geom_line(data=dat.dens, aes(x=x,y=y)),
                     ylim(0,max.dens), theme_bw(), xlab("Density"), ylab(""),
                     ggtitle(estimand),
                     theme(plot.title = element_text(size = 8, face = "bold"))),
      shared_axis_labels = TRUE, width = width_dens, nudge_x = nudge_dens)
    if (plots=="outcome") {
      plt.out <- plt + out_plot
    }
    if (plots=="density") {
      plt.out <- plt + dens_plot
    }
    if (plots=="both") {
      plt.out <- plt + out_plot + dens_plot
    }
  }
  plt.out
  return(plt.out)
}

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
### Heat Map ###
plot_heatmap <- function(x, grid.data, grid.thres) {
  
  .Deprecated("plot_dependence")
  
  if (is.null(x$ple.fit)) {
    stop("Heatmap requires ple model fit: Check ple argument")
  }
  if (dim(grid.data)[2]>3 | dim(grid.data)[2]<2) {
    stop("Heatmap only applicable for grid.data with 2 or 3 variables")
  }
  if (dim(grid.data)[2]==3) {
    if (!is.factor(grid.data[,3]) | length(unique(grid.data[,3]))>6 ) {
      stop("Third column in grid.data should be factor or <=6 unique values")
    }
  }
  # Extract training set covariate space #
  X.train = x$out.train[,!(colnames(x$out.train) %in% c("Y", "A", "Subgrps"))]
  name.var1 = colnames(grid.data)[1]
  name.var2 = colnames(grid.data)[2]
  name.var3 = colnames(grid.data)[3]
  
  # Create stacked covariate space #
  stack_grid = function(i) {
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
  res.est = ggplot2::ggplot(data = est.dat,
                            aes_string(x=name.var1, y=name.var2, fill="est")) +
    ggplot2::geom_tile() + 
    ggplot2::labs(fill = "PLE") +
    ggplot2::scale_fill_gradient2(low="navy", mid="white", high="red")
  res.prob = ggplot2::ggplot(data = est.dat, 
                             aes_string(x=name.var1, y=name.var2, fill="prob")) +
    ggplot2::geom_tile() + 
    ggplot2::labs(fill = paste("Prob(PLE", grid.thres, ")",sep="")) +
    ggplot2::scale_fill_gradient2(low="navy", mid="white", high="red",
                         midpoint=0.5)
  # Lastly, if there were 3 input variables, facet_wrap by third variable #
  if ( dim(grid.data)[2]==3  ){
    res.est = res.est + ggplot2::facet_wrap(as.formula(paste("~", name.var3)))
    res.prob = res.prob + ggplot2::facet_wrap(as.formula(paste("~", name.var3)))
  }
  res = list(heatmap.est=res.est, heatmap.prob=res.prob)
  return(res)
}
