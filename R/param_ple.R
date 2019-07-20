#' Parameter Estimation: Patient-Level Estimates
#'
#' Parameter estimation and inference through patient-level estimates. Usable for continuous
#' and binary outcomes (possibly survival, needs further evaluation).
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate space.
#' @param mu_hat Patient-level estimates (See PLE_models)
#' @param Subgrps Identified subgroups (can be the overall population)
#' @param alpha_ovrl Two-sided alpha level for overall population
#' @param alpha_s Two-sided alpha level at subgroup
#' @param ... Any additional parameters, not currently passed through.
#'
#' @return Data-set with parameter estimates and corresponding variability metrics,
#' for overall and subgroups. Subgrps=0 corresponds to the overall population by default.
#'  \itemize{
#'   \item param.dat - Parameter estimates and variability metrics (est, SE,
#'   LCL/UCL = lower/upper confidence limits, pval = p-value).
#'   }
#' @export
#' @examples
#' library(StratifiedMedicine)
#'
#' ## Continuous ##
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#' train = data.frame(Y, A, X)
#'
#' ## Estimate PLEs (ranger) ##
#' res_ranger = ple_train(Y, A, X, Xtest=X, ple = "ple_ranger")
#'
#' ## Identify Subgroups: MOB (lmtree) ##
#' res_lmtree = submod_train(Y, A, X, Xtest=X, submod="submod_lmtree")
#'
#' ## Parameter-estimation ##
#' params = param_ple(Y, A, X, mu_hat = res_ranger$mu_train,
#'                   Subgrps = res_lmtree$Subgrps.train, alpha_ovrl=0.05,
#'                   alpha_s=0.05)
#' params
#'
#'
### PLE Param: Plug-in estimator using PLE estimates, Use EIF for SEs
param_ple = function(Y, A, X, mu_hat, Subgrps, alpha_ovrl, alpha_s, ...){

  indata = data.frame(Y=Y, A=A, X)
  ## Estimate overall and within subgroups ##
  looper = function(s){
    Y.s = indata$Y[Subgrps %in% s]
    A.s = indata$A[Subgrps %in% s]
    n.s = length(Y.s)
    probA = mean(A.s)
    mu.s = mu_hat[Subgrps %in% s,]
    # Average PLEs for point-estimate #
    est = mean(mu.s$PLE, na.rm=TRUE)
    # EIF for variance estimate #
    eif = ( A.s*Y.s - (A.s-probA)*mu.s$mu1 )/ probA -
      ( (1-A.s)*Y.s + (A.s-probA)*mu.s$mu0 ) / (1-probA) - est
    SE = sqrt( n.s^(-2) * sum( eif^2 )  )
    LCL = est-qt( (1-alpha_s/2), df=n.s-1 )*SE
    UCL = est+qt( (1-alpha_s/2), df=n.s-1 )*SE
    pval = 2*pt(-abs(est/SE), df=n.s-1)
    return( c(n.s, est, SE, LCL, UCL, pval) )
  }
  ## Across subgroups ##
  S_levels = as.numeric( names(table(Subgrps)) )
  S_N = as.numeric( table(Subgrps) )
  param.dat = lapply(S_levels, looper)
  param.dat = do.call(rbind, param.dat)
  param.dat = data.frame( S = S_levels, param.dat)
  # Overall #
  param.dat0 = looper(s = S_levels)
  param.dat0 = c(0, param.dat0)
  param.dat = rbind(param.dat0, param.dat)
  colnames(param.dat) = c("Subgrps", "N", "est", "SE", "LCL", "UCL", "pval")
  # return results #
  return( param.dat )
}
