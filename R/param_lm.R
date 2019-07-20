#' Parameter Estimation: Linear Regression
#'
#' For each identified subgroup, fit separate linear regression models. Point-estimates and
#' variability metrics in the overall population are obtained by aggregating subgroup
#' specific results (adaptive weighting or sample size weighting).
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate space.
#' @param mu_hat Patient-level estimates (See PLE_models)
#' @param Subgrps Identified subgroups (can be the overall population)
#' @param alpha_ovrl Two-sided alpha level for overall population
#' @param alpha_s Two-sided alpha level at subgroup
#' @param combine For overall population, method of combining subgroup-specific results.
#' Default is "adaptive", "SS" corresponds to sample size weighting.
#' @param ... Any additional parameters, not currently passed through.
#'
#' @return Data-set with parameter estimates (average treatment effect) and corresponding
#' variability metrics, for overall and subgroups. Subgrps=0 corresponds to the overall
#' population by default.
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
#'
#' ## Identify Subgroups: MOB (lmtree) ##
#' res_lmtree = submod_train(Y, A, X,  Xtest=X, submod="submod_lmtree")
#'
#' ## Parameter-estimation ##
#' params = param_lm(Y, A, X, Subgrps = res_lmtree$Subgrps.train, alpha_ovrl=0.05,
#'                   alpha_s=0.05)
#' params
#'
#' @seealso \code{\link{param_combine}}
#'
### Linear Regression: Estimate E(Y|A=1) - E(Y|A=0) ###
param_lm = function(Y, A, X, mu_hat, Subgrps, alpha_ovrl, alpha_s, combine="adaptive",
                    ...){

  indata = data.frame(Y=Y,A=A, X)
  ## Subgroup Specific Estimate ##
  looper = function(s){
    lm.mod = tryCatch( lm(Y ~ A , data=indata[Subgrps==s,]),
                       error = function(e) "param error" )
    if (is.character(lm.mod)){
      est = NA; SE = NA; pval = NA; LCL = NA; UCL = NA;
    }
    if (is.list(lm.mod)){
      est = summary(lm.mod)$coefficients[2,1]
      SE = summary(lm.mod)$coefficients[2,2]
      LCL =  confint(lm.mod, level=1-alpha_s)[2,1]
      UCL =  confint(lm.mod, level=1-alpha_s)[2,2]
      pval = summary(lm.mod)$coefficients[2,4]
    }
    return( c(est, SE, LCL, UCL, pval) )
  }
  S_levels = as.numeric( names(table(Subgrps)) )
  S_N = as.numeric( table(Subgrps) )
  param.dat = lapply(S_levels, looper)
  param.dat = do.call(rbind, param.dat)
  param.dat = data.frame( S = S_levels, N=S_N, param.dat)
  colnames(param.dat) = c("Subgrps", "N", "est", "SE", "LCL", "UCL", "pval")
  # Combine results and estimate effect in overall population #
  param.dat0 = param_combine(param.dat = param.dat, alpha_ovrl=alpha_ovrl,
                             combine=combine)
  param.dat = rbind(param.dat0, param.dat)
  return( param.dat )
}
