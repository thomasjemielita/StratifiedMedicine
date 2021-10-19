#' Partial dependence plots: Single Variable (marginal effect) or heat map (2 to 
#' 3 variables). 
#'
#' @param object Fitted \code{ple_train} or \code{PRISM} object
#' @param X input covariate space. Default=NULL. 
#' @param target Which patient-level estimate to target for PDP based plots. Default=NULL, 
#' which uses the estimated treatment difference. 
#' @param vars Variables to visualize (ex: c("var1", "var2", "var3)). If no grid.data 
#' provided, defaults to using seq(min(var), max(var)) for each continuous variables. 
#' For categorical, uses all categories. 
#' @param grid.data Input grid of values for 2-3 covariates (if 3, last variable cannot
#' be continuous). This is required for type="heatmap". Default=NULL.
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
#' @importFrom ggplot2 ggplot aes_string aes geom_tile ggtitle geom_rug 
#' @importFrom ggplot2 scale_fill_gradient2 facet_wrap
#' 
#' @examples
#' \donttest{
#' library(StratifiedMedicine)
#' ## Continuous ##
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#'
#'
#' # Fit through ple_train wrapper #
#' mod = ple_train(Y=Y, A=A, X=X, Xtest=X, ple="ranger", meta="X-learner")
#' plot_dependence(mod, X=X, vars="X1")
#' }
#' 
plot_dependence <- function(object, X=NULL, target=NULL, vars, 
                            grid.data=NULL, ...) {
  
  if (class(object)=="PRISM") {
    out.train <- object$out.train
  }
  if (class(object)=="ple_train") {
    if (is.null(X)) {
      stop("Must supply covariate data (X)")
    }
    out.train <- data.frame(X)
  }
  family <- object$family
  ## Target? ##
  if (is.null(target)) {
    ple_name <- colnames(object$mu_train)[grepl("diff", colnames(object$mu_train))]
    ple_name <- ple_name[1]
  }
  if (!is.null(target)) { ple_name <- target }
  # Set up Labels #
  pieces <- unlist(strsplit(ple_name, "_"))
  if (pieces[1] %in% c("diff")) {
    if (family=="gaussian") {
      ple.label <- paste("E(Y|X,A=", pieces[2], ")-", "E(Y|X,A=", pieces[3], ")", sep="")
    }
    if (family=="binomial") {
      ple.label <- paste("P(Y=1|X,A=", pieces[2], ")-", "P(Y=1|X,A=", pieces[3], ")", sep="")
    }
    if (family=="survival") {
      if (object$ple=="ranger") {
        ple.label <- paste("RMST(X,A=", pieces[2], ")-", "RMST(X,A=", pieces[3], ")", sep="")
      }
      if (object$ple %in% c("linear", "glmnet")) {
        ple.label <- paste("logHR(X,A=", pieces[2], ")-", "logHR(X,A=", pieces[3], ")", sep="")
      }
    } 
  }
  if (pieces[1] != "diff") {
    ple.label <- ple_name
  }
  # Set up data for estimation #
  if (is.null(grid.data)) {
    numb_vars <- length(vars)
    var1_type <- ifelse(is.numeric(out.train[,vars[1]]), "ctns", "fact")
    if (var1_type=="ctns") {
      min.v <- min(out.train[vars[1]])
      max.v <- max(out.train[vars[1]])
      by.v <- (max.v-min.v)/20
      vals.1 <- seq(min.v, max.v, by=by.v)
    }
    if (var1_type=="fact") {
      vals.1 <- unique(out.train[vars[1]])
    }
    if (numb_vars>1) {
      var2_type <- ifelse(is.numeric(out.train[,vars[2]]), "ctns", "fact")
      if (var2_type=="ctns") {
        min.v <- min(out.train[vars[2]])
        max.v <- max(out.train[vars[2]])
        by.v <- (max.v-min.v)/20
        vals.2 <- seq(min.v, max.v, by=by.v)
      }
      if (var2_type=="fact") {
        vals.2 <- unique(out.train[vars[2]])
      }
      if (numb_vars==3) {
        var3_type <- ifelse(is.numeric(out.train[,vars[3]]), "ctns", "fact")
        vals.3 <- unique(out.train[vars[3]])
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
  if (numb_vars>3) {
    stop("Heatmap only applicable for grid.data up to 3 variables")
  }
  if (numb_vars==3) {
    if (!is.factor(grid.data[,3]) | length(unique(grid.data[,3]))>6 ) {
      stop("Third column in grid.data should be factor or <=6 unique values")
    }
  }
  # Extract training set covariate space #
  X.train = out.train[,!(colnames(out.train) %in% c("Y", "A", "Subgrps"))]
  name.var1 = colnames(grid.data)[1]
  name.var2 = colnames(grid.data)[2]
  name.var3 = colnames(grid.data)[3]
  
  # Create stacked covariate space #
  stack_grid = function(i, X.train) {
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
  X.grid = lapply(1:dim(grid.data)[1], stack_grid, X.train=X.train)
  X.grid = do.call(rbind, X.grid)
  counter.vec = X.grid$counter
  X.grid = X.grid[,!(colnames(X.grid) %in% "counter")]

  # Univariate (Marginal Effect) #
  if (numb_vars==1) {
    if (class(object)=="PRISM") {
      grid.ple = predict(object, newdata = X.grid, type="ple")
    }
    if (class(object)=="ple_train") {
      grid.ple = predict(object, newdata = X.grid) 
    }
    y.label <- paste("Estimates:", ple.label)
    grid.ple$PLE <- grid.ple[[ple_name]]
    avg.ple = aggregate(grid.ple$PLE ~ counter.vec, FUN="mean")
    est.dat = data.frame(grid.data, est = avg.ple$`grid.ple$PLE`)
    res.est = ggplot2::ggplot(data = est.dat, 
                              ggplot2::aes_string(x=name.var1, y="est")) + 
      ggplot2::geom_point() + ggplot2::geom_smooth(se=FALSE) + 
      ggplot2::xlab(name.var1) +
      ggplot2::ylab(y.label)
      ggplot2::geom_rug(data=out.train, 
                        ggplot2::aes_string(name.var1), sides="b", inherit.aes = F)
    res <- res.est
  }
  # Heat Map (2 or 3 variables) #
  if (numb_vars>1) {
    ## predict across grid ##
    grid.ple = predict(object, newdata = X.grid, type="ple")
    grid.ple$PLE <- grid.ple[[ple_name]]
    avg.ple = aggregate(grid.ple$PLE ~ counter.vec, FUN="mean")
    ### Plot Heat-map ###
    est.dat = data.frame(grid.data, est = avg.ple$`grid.ple$PLE`)
    res.est = ggplot2::ggplot(data = est.dat, 
                              ggplot2::aes_string(x=name.var1, y=name.var2, fill="est")) +
      ggplot2::geom_tile() + ggplot2::labs(fill = "est") + 
      ggplot2::ggtitle(paste("Heat Map Estimates:", ple.label))+
      ggplot2::geom_rug(data=out.train, 
                        ggplot2::aes_string(name.var1, name.var2), inherit.aes = F) + 
      ggplot2::scale_fill_gradient2(low="navy", mid="white", high="red")
    # Lastly, if there were 3 input variables, facet_wrap by third variable #
    if (dim(grid.data)[2]==3){
      res.est = res.est + ggplot2::facet_wrap(as.formula(paste("~", name.var3)))
    }
    res = res.est
  }
  return(res)
}