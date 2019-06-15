#' Parameter Estimation: RMST
#'
#'For each identified subgroups and in the overall population, estimate the restricted
#'mean survival time (RMST). IN PROGRESS
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
#' @import survRM2
#'
#' @return Data-set with parameter estimates (RMST) and corresponding
#' variability metrics, for overall and subgroups.
#'  \itemize{
#'   \item param.dat - Parameter estimates and variability metrics
#'   }
#' @export
#' @examples
#'
#' \donttest{
#' library(StratifiedMedicine)
#' ## Survival Data ##
#` # (uncomment to run PRISM with resampling for survival) ##
#` # ## Load TH.data (no treatment; generate treatment randomly to simulate null effect) ##
#` # data("GBSG2", package = "TH.data")
#` # surv.dat = GBSG2
#` # ## Design Matrices ###
#` # Y = with(surv.dat, Surv(time, cens))
#` # X = surv.dat[,!(colnames(surv.dat) %in% c("time", "cens")) ]
#` # A = rbinom( n = dim(X)[1], size=1, prob=0.5  )
#'
#' # MOB-Weibull Subgroup Model ##
#' res_weibull = SubMod_weibull(Y, A, X, Xtest=X, family="survival")
#' plot(res_weibull$mod)
#'
#' ## Parameter-Estimation ##
#' params = Param_RMST(Y, A, X, Subgrps = res_weibull$Subgrps.train, alpha_ovrl=0.05,
#'                    alpha_s=0.05)
#' params
#' }
#'
#'
### RMST (unadjusted) ###
Param_RMST = function(Y, A, X, mu_hat, alpha_ovrl, alpha_s, Subgrps, ...){
  indata = data.frame(Y=Y, A=A, X)
  ### Loop through subgroups ##
  looper = function(s){
    time = indata$Y[Subgrps==s,1]
    status = indata$Y[Subgrps==s,2]
    arm = indata$A[Subgrps==s]
    obj = tryCatch( rmst2(time, status, arm),
                    error = function(e) "param error" )
    if (is.character(obj)){
      est = NA; pval = NA; LCL = NA; UCL = NA;
    }
    if (is.list(obj)){
      est = obj$unadjusted.result[1,1]
      pval = obj$unadjusted.result[1,4]
      LCL =  obj$unadjusted.result[1,2]
      UCL =  obj$unadjusted.result[1,3]
    }
    return( c(est, pval, LCL, UCL) )
  }
  S_levels = as.numeric( names(table(Subgrps)) )
  S_N = as.numeric( table(Subgrps) )
  param.dat = lapply(S_levels, looper)
  param.dat = do.call(rbind, param.dat)
  param.dat = data.frame( S = S_levels, N=S_N, param.dat)
  colnames(param.dat) = c("Subgrps", "N", "rmst (A=1)-(A=0)", "pval", "LCL", "UCL")
  return( param.dat )
}
