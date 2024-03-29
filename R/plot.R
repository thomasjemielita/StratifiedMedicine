globalVariables(c("Rules", "est", "LCL", "UCL", "PLE", "label", "N", "estimand",
                  "splitvar", "p.value", "surv", "A", "Subgrps",
                  "LCL0", "UCL0", "est0", "events", "prob.est",
                  "x", "y"))
#' plot.PRISM
#'
#' Plots PRISM results. Options include "tree", "forest", "resample", and "PLE:waterfall".
#'
#' @param x PRISM object
#' @param type Type of plot (default="tree", \code{ggparty} based plot with parameter 
#' estimates, along with options for including outcome or probability based plots). 
#' Other options include  "forest" (forest plot for overall and subgroups),"PLE:waterfall" 
#' (waterfall plot of PLEs), "PLE:density" (density plot of PLEs), "resample" (resampling 
#' distribution of parameter estimates for overall and subgroups), and "heatmap" 
#' (heatmap of ple estimates/probabilities). For "tree" and "forest", CIs are based on 
#' the observed data unless resampling is used. For bootstrap resampling, if 
#' calibrate=TRUE, then calibrated CIs along are shown, otherse CIs based on the 
#' percentile method are shown.
#' @param target For "resample" plot only, must be specify which estimand to visualize.
#' Default=NULL.
#' @param grid.data Input grid of values for 2-3 covariates (if 3, last variable cannot
#' be continuous). This is required for type="heatmap". Default=NULL.
#' @param grid.thres Threshold for PLE, ex: I(PLE>thres). Used to estimate P(PLE>thres) for
#' type="heatmap". Default is ">0". Direction can be reversed and can include equality
#' sign (ex: "<=").
#' @param tree.plots Type of plots to include in each node of the "tree" plot. Default="outcome". 
#' For non-survival data, if the fitted PRISM object (x) does not include patient-level 
#' estimates (ple="None"), or if param="lm", this will plot the observed outcomes (Y) by the 
#' treatment assignment (A). If the fitted PRISM object includes patient-level estimates 
#' (ex: ple="ranger"), this includes box-plots of the model-based (if param="ple") or double-robust
#' based (if param="dr") counter-factual estimates of E(Y|X,A=a) for continuous outcomes or 
#' Prob(Y=1|X,A=a) for binary outcomes (truncated to 0,1). For survival data, Kaplan-Meier 
#' based survival estimates are plotted by treatment group. For "density", the estimated 
#' probability density of the treatment effects is shown (normal approximation, unless resampling is used). 
#' "both" include the "outcome" and "density" plots. If tree.plots = "none", then only the 
#' tree structure is shown.
#' @param prob.thres Probability threshold, ex: P(Mean(A=1 vs A=0)>c. Default=NULL, 
#' which defaults to using ">0", unless param="cox", which  "P(HR(A=1 vs A=0))<1". 
#' If a density plot is included, setting prob.thres=">c" will use green colors 
#' for values above c, and red colors for values below c. If prob.thres="<c", the 
#' reverse color scheme is used.
#' @param nudge_out Nudge tree outcome plot (see ggparty for details)
#' @param nudge_dens Nudge tree density plot
#' @param width_out Width of tree outcome plot (see ggparty for details)
#' @param width_dens Width of density tree outcome plot
#' @param ... Additional arguments (currently ignored).
#' @return Plot (ggplot2) object
#' @method plot PRISM
#' @export
#' @importFrom stats reorder as.formula density
#' @importFrom ggplot2 geom_pointrange geom_text xlab theme theme_bw coord_flip
#' @importFrom ggplot2 position_nudge ylab element_text element_blank
#' @importFrom ggplot2 geom_density geom_point geom_line
#' @import ggparty
#' @importFrom partykit nodeapply nodeids split_node data_party as.partynode
#' @seealso \code{\link{PRISM}}


plot.PRISM = function(x, type="tree", target=NULL, grid.data=NULL, grid.thres=">0",
                      prob.thres=NULL,
                      tree.plots="outcome",
                      nudge_out=0.1, width_out=0.5,
                      nudge_dens=ifelse(tree.plots=="both", 0.3, 0.1),
                      width_dens=0.5, ...) {
  
  if (type=="PLE:waterfall") {
    ple.fit <- list(ple = x$ple, mu_train = x$mu_train, 
                    family = x$family, 
                    treetype = x$ple.fit$mod$fit0[[1]]$mod$mod$treetype)
    ple.fit$mu_train$Subgrps <- factor(x$out.train$Subgrps)
    res <- plot_ple(object=ple.fit, type="waterfall")
    return(res)
  }
  if (type=="PLE:density") {
    ple.fit <- list(ple = x$ple, mu_train = x$mu_train, 
                    family = x$family,
                    treetype = x$ple.fit$mod$fit0[[1]]$mod$mod$treetype)
    ple.fit$mu_train$Subgrps <- factor(x$out.train$Subgrps)
    res <- plot_ple(object=ple.fit, type="density")
    return(res)
  }
  # Default Setup #
  x <- default_trt_plots(obj=x)
  param.dat <- x$param.dat
  resample <- ifelse(is.null(x$resample), "None", x$resample)
  bayes <- ifelse(is.null(x$bayes.fun), FALSE, TRUE)
 
  if (type=="tree" & length(unique(x$out.train$Subgrps))==1) {
    type = "forest"
  }
  if (type=="tree"){
    if (!is.null(prob.thres)) {
      thres.name <- paste("Prob(",prob.thres, ")", sep="")
      x2 <- prob_calculator(x, thres=prob.thres)
    }
    if (is.null(prob.thres)) {
      x2 <- x
      prob.thres <- ifelse(x2$param=="cox", "<1", ">0")
    }
    cls <- class(x2$submod.fit$mod)
    if ("party" %in% cls) {
      # print(x2$param.dat)
      res <- do.call("plot_tree", list(object=x2, plots=tree.plots,
                                       prob.thres = prob.thres,
                                       nudge_out=nudge_out, width_out=width_out,
                                       nudge_dens=nudge_dens, width_dens=width_dens))
    }
    if (!("party" %in% cls)) {
      stop(paste("Tree Plots for non partykit models not currently supported."))
    }
  }
  if (type=="forest"){
    res <- plot_forest(x)
  }
  if (type=="resample") {
    res <- plot_resample(x=x, target=target)
  }
  if (type=="heatmap"){
    res <- plot_heatmap(x=x, grid.data=grid.data, grid.thres=grid.thres)
  }
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
