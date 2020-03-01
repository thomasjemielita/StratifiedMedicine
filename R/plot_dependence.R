#' Partial dependence plots: Single Variable (marginal effect) or heat map (2 to 
#' 3 variables). 
#'
#' @param object Fitted PRISM object
#' @param vars Variables to visualize (ex: c("var1", "var2", "var3)). If no grid.data 
#' provided, defaults to using seq(min(var), max(var)) for each continuous variables. 
#' For categorical, uses all categories. 
#' @param grid.data Input grid of values for 2-3 covariates (if 3, last variable cannot
#' be continuous). This is required for type="heatmap". Default=NULL.
#' @param grid.thres Threshold for PLE, ex: I(PLE>thres). Used to estimate P(PLE>thres) for
#' type="heatmap". Default is ">0". Direction can be reversed and can include equality
#' sign (ex: "<=").
#' @param estimand Estimand for which to generate dependendence or heat map plots.
#' @param ... Additional arguments (currently ignored).
#' @return Plot (ggplot2) object
#' @references 
#' \itemize{
#' \item Friedman, J. Greedy function approximation: A gradient boosting machine.
#'  Annals of statistics (2001): 1189-1232
#' \item Zhao, Qingyuan, and Trevor Hastie. Causal interpretations of black-box models.
#'  Journal of Business & Economic Statistics, to appear. (2017).
#' }
#' @export

plot_dependence <- function(object, vars, grid.data=NULL, grid.thres=">0", 
                            estimand=NULL, ...) {
  
  if (is.null(grid.data)) {
    numb_vars <- length(vars)
    var1_type <- ifelse(is.numeric(object$out.train[,vars[1]]), "ctns", "fact")
    if (var1_type=="ctns") {
      min.v <- min(object$out.train[vars[1]])
      max.v <- max(object$out.train[vars[1]])
      by.v <- (max.v-min.v)/20
      vals.1 <- seq(min.v, max.v, by=by.v)
    }
    if (var1_type=="fact") {
      vals.1 <- unique(object$out.train[vars[1]])
    }
    if (numb_vars>1) {
      var2_type <- ifelse(is.numeric(object$out.train[,vars[2]]), "ctns", "fact")
      if (var2_type=="ctns") {
        min.v <- min(object$out.train[vars[2]])
        max.v <- max(object$out.train[vars[2]])
        by.v <- (max.v-min.v)/20
        vals.2 <- seq(min.v, max.v, by=by.v)
      }
      if (var2_type=="fact") {
        vals.2 <- unique(object$out.train[vars[2]])
      }
      if (numb_vars==3) {
        var3_type <- ifelse(is.numeric(object$out.train[,vars[3]]), "ctns", "fact")
        vals.3 <- unique(object$out.train[vars[3]])
      }
    }
    if (numb_vars==1) {
      grid.data <- data.frame(var = vals.1)
      colnames(grid.data) <- vars
    }
    if (numb_vars>1) {
      grid.data = expand.grid(var1 = vals.1, var2 =vals.2)
      if (numb_vars==3) {
        for (v in vals.3) {
          grid.data <- data.frame(grid.data, var3 = v)
        }
      }
    }
    colnames(grid.data) <- vars
  }
  if (!is.null(grid.data)) {
    numb_vars <- dim(grid.data)[2] 
  }
  if (is.null(object$ple.fit)) {
    stop("Heatmap requires ple model fit: Check ple argument")
  }
  if (numb_vars>3) {
    stop("Heatmap only applicable for grid.data up to 3 variables")
  }
  if (numb_vars==3) {
    if (!is.factor(grid.data[,3]) | length(unique(grid.data[,3]))>6 ) {
      stop("Third column in grid.data should be factor or <=6 unique values")
    }
  }
  # Extract training set covariate space #
  X.train = object$out.train[,!(colnames(object$out.train) %in% c("Y", "A", "Subgrps"))]
  name.var1 = colnames(grid.data)[1]
  name.var2 = colnames(grid.data)[2]
  name.var3 = colnames(grid.data)[3]
  
  if (is.null(estimand)) {
    A_lvls <- as.character(with(object$out.train, unique(A)[order(unique(A))]))
    if (object$family %in% c("gaussian", "binomial")) {
      E_A0 <- paste("E(Y|A=", A_lvls[1], ")", sep="")
      E_A1 <- paste("E(Y|A=", A_lvls[2], ")", sep="")
      E_diff <- paste(E_A1, "-", E_A0, sep="")
      e.name <- E_diff
    }
    if (object$family %in% "survival") {
      if (object$ple %in% c("ple_ranger", "ple_rfsrc")) {
        if (is.null(object$out.train$A)) {
          e.name <- "rmst"
        }
        if (!is.null(object$out.train$A)) {
          e.name <- paste("rmst(A=", A_lvls[2], " vs ", "A=", A_lvls[1], ")", sep="")
        }
      }
      if (object$ple=="ple_glmnet") {
        e.name <- paste("logHR(A=", A_lvls[2], " vs ", "A=", A_lvls[1], ")", sep="")
      }
    }
  }
  
  # Create stacked covariate space #
  stack_grid = function(i) {
    var1 = grid.data[i,1]
    var2 = grid.data[i,2]
    var3 = grid.data[i,3]
    newdata = X.train
    newdata[name.var1] = var1
    if (!is.null(var2)){
      newdata[name.var2] = var2
    }
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
  
  # Univariate (Marginal Effect) #
  if (numb_vars==1) {
    grid.ple = predict(object, newdata = X.grid, type="ple")
    grid.ple$ind.ple = eval(parse(text=paste("ifelse(grid.ple$PLE",
                                             grid.thres, ", 1, 0)")))
    avg.ple = aggregate(grid.ple$PLE ~ counter.vec, FUN="mean")
    prob.ple = aggregate(grid.ple$ind.ple ~ counter.vec, FUN="mean")
    est.dat = data.frame(grid.data, est = avg.ple$`grid.ple$PLE`,
                         prob = prob.ple$`grid.ple$ind.ple`)
    res.est = ggplot(data = est.dat, aes_string(x=name.var1, y="est")) + 
      geom_point() + geom_smooth(se=FALSE) + xlab(name.var1) + ylab(e.name) +
      geom_rug(data=object$out.train, aes_string(name.var1), sides="b", inherit.aes = F)
    res <- list(res.est=res.est)
  }
  # Heat Map (2 or 3 variables) #
  if (numb_vars>1) {
    ## predict across grid ##
    grid.ple = predict(object, newdata = X.grid, type="ple")
    grid.ple$ind.ple = eval(parse(text=paste("ifelse(grid.ple$PLE",
                                             grid.thres, ", 1, 0)")))
    avg.ple = aggregate(grid.ple$PLE ~ counter.vec, FUN="mean")
    prob.ple = aggregate(grid.ple$ind.ple ~ counter.vec, FUN="mean")
    ### Plot Heat-map ###
    est.dat = data.frame(grid.data, est = avg.ple$`grid.ple$PLE`,
                         prob = prob.ple$`grid.ple$ind.ple`)
    res.est = ggplot(data = est.dat, aes_string(x=name.var1, y=name.var2, fill="est")) +
      geom_tile() + labs(fill = "est") + ggtitle(paste("Heat Map:", e.name))+
      geom_rug(data=object$out.train, aes_string(name.var1, name.var2), inherit.aes = F) + 
      scale_fill_gradient2(low="navy", mid="white", high="red")
    res.prob = ggplot(data = est.dat, aes_string(x=name.var1, y=name.var2, fill="prob")) +
      geom_tile() + labs(fill = paste("Prob(", grid.thres, ")",sep="")) +
      geom_rug(data=object$out.train, aes_string(name.var1, name.var2), inherit.aes = F) + 
      ggtitle(paste("Heat Map:", e.name))+
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