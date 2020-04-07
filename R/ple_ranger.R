#' Patient-level Estimates: Ranger
#'
#' Uses treatment-specific (or with explicit X*A interactions) random forest models (ranger)
#' to obtain patient-level estimates. Used for continuous, binary, or survival outcomes.
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate space.
#' @param Xtest Test set
#' @param byTrt If TRUE, fit treatment-specific ranger models. If FALSE, fit a single ranger
#' model with covariate space (X, A, X*A). For "gaussian" or "binomial", default is TRUE.
#' For "survival", default is FALSE. 
#' @param min.node.pct Minimum sample size in forest nodes (n*min.node.pct)
#' @param family Outcome type ("gaussian", "binomial"), default is "gaussian"
#' @param ... Any additional parameters, not currently passed through.
#'
#' @import ranger
#'
#' @return Trained random forest (ranger) model(s).
#'  \itemize{
#'   \item mod - trained model(s)
#'   \item A - treatment variable (training set)
#'   \item X - covariate space (training set)
#' }
#' @references Wright, M. N. & Ziegler, A. (2017). ranger: A fast implementation of 
#' random forests for high dimensional data in C++ and R. J Stat Softw 77:1-17. 
#' \url{https://doi.org/10.18637/jss.v077.i01}.
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
#' # Counter-factual Random Forest (treatment-specific ranger models) #
#' mod1 = ple_ranger(Y, A, X, Xtest=X)
#'
#' }
#'
#' @export
#' @seealso \code{\link{PRISM}}, \code{\link{ranger}}
#'
#### Counterfactual Forest: Ranger ####
ple_ranger = function(Y, A, X, Xtest, family="gaussian", 
                      byTrt=ifelse(family=="survival", FALSE, TRUE),
                      min.node.pct=0.10, ...) {

  if (is.null(A)){
    mod <- ranger(Y ~ ., data = data.frame(Y, X), 
                   min.node.size = min.node.pct*dim(X)[2] )
    mod = list(mod=mod)
    # Prediction Function #
    pred.fun = function(mod, X, tau=NULL){
      treetype = mod[[1]]$treetype
      if (treetype!="Survival"){
          mu_hat = data.frame(PLE = predict( mod$mod, X )$predictions )
      }
      if (treetype=="Survival"){
        preds = predict( mod$mod, X )
        if (is.null(tau)){
          tau.t <- max(preds$unique.death.times)
        }
        looper_rmst <- function(i, surv, time){
          est.rmst <- rmst_calc(surv = surv[i,],
                                time = time,
                                tau=tau.t)
          return(est.rmst)
        }
        rmst <- lapply(1:dim(X)[1], looper_rmst, surv=preds$survival,
                       time=preds$unique.death.times)
        mu_hat <- do.call(rbind, rmst)
        mu_hat <- data.frame(PLE = mu_hat)
      }
      return(mu_hat)
    }
  }
  if (!is.null(A)){
    A_lvls <- unique(A)[order(unique(A))]
    ## Random Forest models for each Treatment ##
    if (byTrt){
      looper <- function(a) {
        train_a <- data.frame(Y=Y[A==a], X[A==a, , drop=FALSE])
        mod_a <- ranger(Y ~ 1+., data = train_a, 
                        min.node.size = min.node.pct*dim(train_a)[1])
        mod_a$A_lvl <- a
        return(mod_a)
      }
      mod <- lapply(A_lvls, looper)
      names(mod) <- paste("mod", unique(A)[order(unique(A))], sep="_")
      # Prediction Function #
      pred.fun = function(mod, X, tau=NULL){
        treetype = mod[[1]]$treetype
        if (treetype!="Survival"){
          looper <- function(a) {
            mu_a <- data.frame(predict(mod[[a]], X)$predictions)
            colnames(mu_a) <- paste("mu", mod[[a]]$A_lvl, sep="_")
            return(mu_a)
          }
          mu_hat <- lapply(1:length(mod), looper)
          mu_hat <- do.call(cbind, mu_hat)
          mu_hat$PLE = as.numeric(mu_hat[,2]-mu_hat[,1])
        }
        if (treetype=="Survival") {
          looper <- function(a) {
            pred_a <- predict(mod[[a]], X)
            return(pred_a)
          }
          preds <- lapply(1:length(mod), looper)
          if (is.null(tau)) {
            max.times <- NULL
            for (a in 1:length(preds)) {
              max.times <- c(max.times, max(preds[[a]]$unique.death.times))
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
                          time=preds[[a]]$unique.death.times)
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
      X = model.matrix(~., data = X )[,-1]
      A.mat <- model.matrix(~., data=data.frame(A))[,-1]
      len.A <- ifelse(is.vector(A.mat), 1, dim(A.mat))
      X_inter <- NULL
      for (a in 1:len.A) {
        hold_inter = X*data.frame(A.mat)[,a]
        colnames(hold_inter) = paste(colnames(X), A_lvls[a+1], sep="_")
        X_inter <- cbind(X_inter, hold_inter)
      }
      train.inter = data.frame(Y, A=A.mat, X, X_inter)
      ## Fit RF ##
      mod.inter <- ranger(Y ~ ., data = train.inter,
                          min.node.size = min.node.pct*dim(train.inter)[1])
      mod.inter$A_lvls <- A_lvls
      mod = list(mod.inter=mod.inter)
      # Prediction Function #
      pred.fun = function(mod, X, tau=NULL) {
        X = model.matrix(~., data = X )[,-1]
        mod.inter <- mod$mod.inter
        treetype <- mod$mod.inter$treetype
        A_lvls <- mod$mod.inter$A_lvls
        # Initial Interaction #
        X0 = data.frame(0, X, X*0)
        colnames(X0) = c( "A", colnames(X), paste(colnames(X), A_lvls[2], sep="_"))
        X1 = data.frame(1, X, X*1)
        colnames(X1) = c( "A", colnames(X), paste(colnames(X), A_lvls[2], sep="_"))
        if (treetype!="Survival") {
          mu_hat = data.frame(mu0 = predict(mod$mod.inter, X0)$predictions,
                              mu1 = predict(mod$mod.inter, X1)$predictions)
        }
        if (treetype=="Survival") {
          pred1 = predict( mod$mod.inter, X1 ) 
          pred0 = predict( mod$mod.inter, X0 )
          if (is.null(tau)){
            tau.t <- min( max(pred0$unique.death.times), 
                          max(pred1$unique.death.times)  )
          }
          looper_rmst <- function(i, surv, time){
            est.rmst <- rmst_calc(surv = surv[i,],
                                  time = time,
                                  tau=tau.t)
            return(est.rmst)
          }
          # A = 0 #
          rmst0 <- lapply(1:dim(X)[1], looper_rmst, surv=pred0$survival,
                          time=pred0$unique.death.times)
          rmst0 <- do.call(rbind, rmst0)
          # A = 1 #
          rmst1 <- lapply(1:dim(X)[1], looper_rmst, surv=pred1$survival,
                          time=pred1$unique.death.times)
          rmst1 <- do.call(rbind, rmst1)
          mu_hat <- data.frame(mu0 = rmst0, mu1 = rmst1)
        }
        mu_hat$PLE = with(mu_hat, mu1-mu0)
        colnames(mu_hat)[1:2] <- paste("mu", A_lvls, sep="_")
        return(mu_hat)
      }
    }
  }
  res = list(mod=mod, pred.fun=pred.fun)
  class(res) = "ple_ranger"
  ## Return Results ##
  return( res )
}