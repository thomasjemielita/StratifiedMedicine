#' Patient-level Estimates: randomForestSRC
#'
#' Uses treatment-specific (or with explicit X*A interactions) random forest models 
#' (randomForestSRC package) to obtain patient-level estimates. Used for continuous, 
#' binary, or survival outcomes.
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate space.
#' @param Xtest Test set
#' @param ntree Number of trees (default=1000)
#' @param byTrt If TRUE, fit treatment-specific rfsrc models. If FALSE, fit a single
#' rfsrc model with covariate space (X, A, X*A).
#' @param upweight Whether to upweight the probability that the treatment variable is
#' included as a splitting variable (through rfsrc's xvar.wt argument). Default=100 
#' (other variables receive weight of 1). Only applicable for single rfsrc model 
#' (byTrt=FALSE).
#' @param min.node.pct Minimum sample size in forest nodes (n*min.node.pct)
#' @param family Outcome type ("gaussian", "binomial", "survival"), 
#' default is "gaussian".
#' @param ... Any additional parameters, not currently passed through.
#'
#'
#' @return Trained random forest (rfsrc) model(s).
#'  \itemize{
#'   \item mod - trained model(s)
#'   \item A - treatment variable (training set)
#'   \item X - covariate space (training set)
#' }
#' @references 
#' \itemize{
#' \item Breiman L. (2001). Random forests, Machine Learning, 45:5-32.
#' \item Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S. (2008). Random 
#' survival forests, Ann. App. Statist., 2:841-860.
#' }
#' @examples
#' \donttest{
#' library(StratifiedMedicine)
#' ## Continuous ##
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#'
#' # Counter-factual Random Forest (treatment-specific ranger models) #
#' mod1 = ple_rfsrc(Y, A, X, Xtest=X)
#'
#' }
#'
#' @export
#'
#### Counterfactual Forest: rfsrc ####
ple_rfsrc = function(Y, A, X, Xtest, ntree=1000, byTrt=TRUE, upweight=100,
                     min.node.pct=0.10, family="gaussian", ...){
  
  if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
    stop("Package randomForestSRC needed for ple_rfsrc. Please install.")
  }
  
  miss = sum(is.na(X))
  if (miss==0){ na.action = "na.omit" }
  if (miss==1){ na.action = "na.impute" }
  # Convert to model matrix #
  X = model.matrix(~., data = X )[,-1]
  if (is.Surv(Y)){
    train.init <- data.frame(time=Y[,1],status=Y[,2], X)
    form <- as.formula("Surv(time,status)~.")
  }
  if (!is.Surv(Y)){
    train.init <- data.frame(Y, X)
    form <- as.formula("Y~.")
  }
  if (is.null(A)){
    mod <- randomForestSRC::rfsrc(form, data=train.init, ntree=ntree,
                 nodesize = min.node.pct*dim(X)[1], na.action=na.action )
    mod = list(mod=mod)
    # Prediction Function #
    pred.fun = function(mod, X, tau=NULL){
      X = model.matrix(~., data = X )
      X = X[,-1]
      family = mod[[1]]$family
      if (family!="surv"){
          mu_hat = data.frame(PLE = predict(mod$mod, X)$predicted )
      }
      if (family=="surv"){
        preds = predict( mod$mod, data.frame(X) )
        if (is.null(tau)){
          tau.t <- max(preds$time.interest)
        }
        looper_rmst <- function(i, surv, time){
          est.rmst <- rmst_calc(surv = surv[i,],
                                time = time,
                                tau=tau.t)
          return(est.rmst)
        }
        rmst <- lapply(1:dim(X)[1], looper_rmst, surv=preds$survival,
                        time=preds$time.interest)
        mu_hat <- do.call(rbind, rmst)
        mu_hat = data.frame(PLE=mu_hat)
      }
      return(mu_hat)
    }
  }
  if (!is.null(A)) {
    A_lvls <- unique(A)[order(unique(A))]
    ## Random Forest models for each Treatment ##
    if (byTrt){
      looper <- function(a) {
        train_a <- train.init[A==a,]
        mod_a <- randomForestSRC::rfsrc(form, data = train_a,
                                        nodesize = min.node.pct*dim(train_a)[1], 
                                        ntree=ntree, na.action=na.action)
        mod_a$A_lvl <- a
        return(mod_a)
      }
      mod <- lapply(A_lvls, looper)
      names(mod) <- paste("mod", unique(A)[order(unique(A))], sep="_")
      # Prediction Function #
      pred.fun = function(mod, X, tau=NULL){
        X = model.matrix(~., data = X)[,-1]
        family = mod[[1]]$family
        if (family!="surv"){
          looper <- function(a) {
            mu_a <- data.frame(predict(mod[[a]], data.frame(X))$predicted)
            colnames(mu_a) <- paste("mu", mod[[a]]$A_lvl, sep="_")
            return(mu_a)
          }
          mu_hat <- lapply(1:length(mod), looper)
          mu_hat <- do.call(cbind, mu_hat)
          mu_hat$PLE = as.numeric(mu_hat[,2]-mu_hat[,1])
        }
        if (family=="surv"){
          looper <- function(a) {
            pred_a <- predict(mod[[a]], data.frame(X))
            return(pred_a)
          }
          preds <- lapply(1:length(mod), looper)
          if (is.null(tau)) {
            max.times <- NULL
            for (a in 1:length(preds)) {
              max.times <- c(max.times, max(preds[[a]]$time.interest))
            }
            tau.t <- min(max.times)
          }
          looper_rmst <- function(i, surv, time){
            est.rmst <- rmst_calc(surv = surv[i,],
                                  time = time,
                                  tau=tau.t)
            return(est.rmst)
          }
          # Across treatments #
          rmst_preds <- function(a) {
            est <- lapply(1:dim(X)[1], looper_rmst, surv=preds[[a]]$survival,
                          time=preds[[a]]$time.interest)
            est <- data.frame(do.call(rbind, est))
            colnames(est) <- paste("mu", mod[[a]]$A_lvl, sep="_")
            return(est)
          }
          mu_hat <- lapply(1:length(preds), rmst_preds)
          mu_hat <- data.frame(do.call(cbind, mu_hat))
          mu_hat$PLE <- mu_hat[,2]-mu_hat[,1]
        }
        return(mu_hat)
      }
    }
    ## Single Random Forest Model: Generate A*X interactions manually ##
    if (!byTrt){
      A.mat <- model.matrix(~., data=data.frame(A))[,-1]
      len.A <- ifelse(is.vector(A.mat), 1, dim(A.mat))
      X_inter <- NULL
      x.weights <- NULL
      for (a in 1:len.A) {
        hold_inter = X*data.frame(A.mat)[,a]
        colnames(hold_inter) = paste(colnames(X), A_lvls[a+1], sep="_")
        X_inter <- cbind(X_inter, hold_inter)
      }
      train.inter = data.frame(train.init, A=A.mat, X_inter)
      x.weights <- c(rep(1, dim(X)[2]), upweight,
                     rep(rep(upweight, dim(X_inter)[2]), length(A_lvls[-1])))
      ## Fit RF ##
      mod.inter <- randomForestSRC::rfsrc(form, data = train.inter, xvar.wt = x.weights,
                    nodesize = min.node.pct*dim(train.inter)[1], ntree=ntree,
                    na.action=na.action)
      mod.inter$A_lvls <- A_lvls
      mod = list(mod.inter=mod.inter)
      # Prediction Function #
      pred.fun = function(mod, X, tau=NULL){
        X = model.matrix(~., data = X)[,-1]
        mod.inter = mod$mod.inter
        family = mod$mod.inter$family
        A_lvls <- mod$mod.inter$A_lvls
        X0 = data.frame(0, X, X*0)
        colnames(X0) = c( "A", colnames(X), paste(colnames(X), A_lvls[2], sep="_"))
        X1 = data.frame(1, X, X*1)
        colnames(X1) = c( "A", colnames(X), paste(colnames(X), A_lvls[2], sep="_"))
        if (family!="surv"){
          mu_hat = data.frame(mu0 = predict( mod$mod.inter, data.frame(X0) )$predicted,
                              mu1 = predict( mod$mod.inter, data.frame(X1) )$predicted)
          mu_hat$PLE = with(mu_hat, mu1-mu0)
        }
        if (family=="surv"){
          pred0 = predict( mod$mod.inter, data.frame(X0) )
          pred1 = predict( mod$mod.inter, data.frame(X1) )
          if (is.null(tau)){
            tau.t <- min( max(pred0$time.interest), max(pred1$time.interest)  )
          }
          looper_rmst <- function(i, surv, time){
            est.rmst <- rmst_calc(surv = surv[i,],
                                  time = time,
                                  tau=tau.t)
            return(est.rmst)
          }
          # A = 0 #
          rmst0 <- lapply(1:dim(X)[1], looper_rmst, surv=pred0$survival,
                          time=pred0$time.interest)
          rmst0 <- do.call(rbind, rmst0)
          # A = 1 #
          rmst1 <- lapply(1:dim(X)[1], looper_rmst, surv=pred1$survival,
                          time=pred1$time.interest)
          rmst1 <- do.call(rbind, rmst1)
          mu_hat <- data.frame(mu0 = rmst0, mu1 = rmst1)
        }
        mu_hat$PLE <- with(mu_hat, mu1-mu0)
        colnames(mu_hat)[1:2] <- paste("mu", A_lvls, sep="_")
        return(mu_hat)
      }
    }
  }
  res = list(mod=mod, pred.fun=pred.fun)
  class(res) = "ple_rfsrc"
  ## Return Results ##
  return( res )
}