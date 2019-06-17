#' Parameter Estimation: Patient-Level Estimates
#'
#'Parameter estimation and inference through patient-level estimates.
#'Usable for continuous and binary outcomes.
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate matrix. Must be numeric.
#' @param mu_hat Patient-level estimates (See PLE_models)
#' @param Subgrps Identified subgroups (can be the overall population)
#' @param alpha_ovrl Two-sided alpha level for overall population
#' @param alpha_s Two-sided alpha level at subgroup
#' @param ... Any additional parameters, not currently passed through.
#'
#' @return Data-set with parameter estimates and corresponding
#' variability metrics, for overall and subgroups.
#'  \itemize{
#'   \item param.dat - Parameter estimates and variability metrics.
#'   By convention, Subgrps=0 corresponds to overall population.
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
#' res_ranger = ple_ranger(Y, A, X, Xtest=X)
#'
#' ## Identify Subgroups: MOB (lmtree) ##
#' res_lmtree = submod_lmtree(Y, A, X, Xtest=X)
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
  ## Overall Estimate ##
  Y.p = indata$Y
  A.p = indata$A
  n = dim(indata)[1]
  probA = mean(A.p)
  mu = mu_hat
  est0 = mean(mu$PLE, na.rm=TRUE)
  # EIF for variance estimate #
  eif = ( A.p*Y.p - (A.p-probA)*mu$mu1 )/ probA -
    ( (1-A.p)*Y.p + (A.p-probA)*mu$mu0 ) / (1-probA) - est0
  SE0 = sqrt( n^(-2) * sum( eif^2 )  )
  LCL0 = est0-qt( (1-alpha_ovrl/2), df=n-1 )*SE0
  UCL0 = est0+qt( (1-alpha_ovrl/2), df=n-1 )*SE0
  param.dat0 = data.frame(Subgrps=0, N=n, est=est0, SE=SE0, LCL=LCL0, UCL=UCL0)
  param.dat0$pval = with(param.dat0, 2*pt(-abs(est/SE), df=N-1) )

  ## Subgroup-Specific Estimate ##
  looper = function(s){
    Y.s = indata$Y[Subgrps==s]
    A.s = indata$A[Subgrps==s]
    n.s = length(Y.s)
    probA = mean(A.s)
    mu.s = mu_hat[Subgrps==s,]
    # Average PLEs for point-estimate #
    est = mean(mu.s$PLE, na.rm=TRUE)
    # EIF for variance estimate #
    eif = ( A.s*Y.s - (A.s-probA)*mu.s$mu1 )/ probA -
      ( (1-A.s)*Y.s + (A.s-probA)*mu.s$mu0 ) / (1-probA) - est
    SE = sqrt( n.s^(-2) * sum( eif^2 )  )
    LCL = est-qt( (1-alpha_s/2), df=n.s-1 )*SE
    UCL = est+qt( (1-alpha_s/2), df=n.s-1 )*SE
    pval = 2*pt(-abs(est/SE), df=n.s-1)
    return( c(est, SE, LCL, UCL, pval) )
  }
  S_levels = as.numeric( names(table(Subgrps)) )
  S_N = as.numeric( table(Subgrps) )
  param.dat = lapply(S_levels, looper)
  param.dat = do.call(rbind, param.dat)
  param.dat = data.frame( S = S_levels, N=S_N, param.dat)
  colnames(param.dat) = c("Subgrps", "N", "est", "SE", "LCL", "UCL", "pval")
  param.dat = rbind( param.dat0, param.dat )
  return( param.dat )
}
