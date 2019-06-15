#' Parameter Estimation: Cox Regression
#'
#'For each identified subgroups and in the overall population, fit separate cox regression models.
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
#' @return Data-set with parameter estimates (hazard ratio) and corresponding
#' variability metrics, for overall and subgroups.
#'  \itemize{
#'   \item param.dat - Parameter estimates and variability metrics
#'   }
#' @export
#' @examples
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
#' params = Param_cox(Y, A, X, Subgrps = res_weibull$Subgrps.train, alpha_ovrl=0.05,
#'                    alpha_s=0.05)
#' params
#' }
#'


### Cox Regression: Hazard Ratios ###
Param_cox = function(Y, A, X, mu_hat, Subgrps, alpha_ovrl, alpha_s, ...){

  indata = data.frame(Y=Y, A=A, X)
  ## Overall estimate ##
  mod.ovrl = coxph(Y ~ A , data=indata)
  param.dat0 = data.frame( Subgrps=0, N = dim(indata)[1],
                           est = exp( summary(mod.ovrl)$coefficients[1] ),
                           SE = summary(mod.ovrl)$coefficients[3],
                           LCL = exp( confint(mod.ovrl, level=1-alpha_ovrl)[1] ),
                           UCL = exp( confint(mod.ovrl, level=1-alpha_ovrl)[2] ),
                           pval = summary(mod.ovrl)$coefficients[5] )

  ### Loop through subgroups ##
  looper = function(s){
    ## Extract HR, SE, 95% CI, and p-value for Subgroup Specific Treatment Effect ##
    cox.mod = tryCatch( coxph(Y ~ A , data=indata[Subgrps==s,]),
                        error = function(e) "param error" )
    if (is.character(cox.mod)){
      est = NA; SE = NA; pval = NA; LCL = NA; UCL = NA;
    }
    if (is.list(cox.mod)){
      est = exp( summary(cox.mod)$coefficients[1] )
      SE = summary(cox.mod)$coefficients[3]
      LCL = exp( confint(cox.mod, level=1-alpha_s)[1] )
      UCL = exp( confint(cox.mod, level=1-alpha_s)[2] )
      pval = summary(cox.mod)$coefficients[5]
    }
    return( c(est, SE, LCL, UCL, pval) )
  }
  S_levels = as.numeric( names(table(Subgrps)) )
  S_N = as.numeric( table(Subgrps) )
  param.dat = lapply(S_levels, looper)
  param.dat = do.call(rbind, param.dat)
  param.dat = data.frame( S = S_levels, N=S_N, param.dat)
  colnames(param.dat) = c("Subgrps", "N", "est", "SE", "LCL", "UCL", "pval")
  param.dat = rbind( param.dat0, param.dat  )
  return( param.dat )
}
