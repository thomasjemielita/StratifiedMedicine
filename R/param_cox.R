#' Parameter Estimation: Cox Regression
#'
#' For each identified subgroup, fit separate cox regression models. Point-estimates and
#' variability metrics in the overall population are obtained by aggregating subgroup
#' specific results (adaptive weighting or sample size weighting).
#'
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
#' @return Data-set with parameter estimates (hazard ratio) and corresponding
#' variability metrics, for overall and subgroups. Subgrps=0 corresponds to the overall
#' population by default.
#'  \itemize{
#'   \item param.dat - Parameter estimates and variability metrics (est=HR, SE=SE(logHR),
#'   LCL/UCL = lower/upper confidence limit on HR scale, pval = p-value).
#'   }
#' @export
#' @examples
#' \donttest{
#' library(StratifiedMedicine)
#' # Survival Data #
#' require(TH.data); require(coin)
#` data("GBSG2", package = "TH.data")
#` surv.dat = GBSG2
#` # Design Matrices ###
#` Y = with(surv.dat, Surv(time, cens))
#` X = surv.dat[,!(colnames(surv.dat) %in% c("time", "cens")) ]
#` A = rbinom( n = dim(X)[1], size=1, prob=0.5  ) ## simulate null treatment
#'
#' # MOB-Weibull Subgroup Model ##
#' res_weibull = submod_train(Y, A, X, Xtest=X, family="survival",
#'                             submod="submod_weibull")
#' plot(res_weibull$mod)
#'
#' ## Parameter-Estimation ##
#' params = param_cox(Y, A, X, Subgrps = res_weibull$Subgrps.train, alpha_ovrl=0.05,
#'                    alpha_s=0.05)
#' params
#' }
#'
#' @seealso \code{\link{param_combine}}

### Cox Regression: Hazard Ratios ###
param_cox = function(Y, A, X, mu_hat, Subgrps, alpha_ovrl, alpha_s, combine="adaptive",
                     ...){

  indata = data.frame(Y=Y, A=A, X)
  ### Loop through subgroups ##
  looper = function(s){
    ## Extract HR, SE, 95% CI, and p-value for Subgroup Specific Treatment Effect ##
    cox.mod = tryCatch( coxph(Y ~ A , data=indata[Subgrps==s,]),
                        error = function(e) "param error",
                        warning = function(w) "convergence issues")
    if (is.character(cox.mod)){
      est = NA; SE = NA; pval = NA; LCL = NA; UCL = NA;
    }
    if (is.list(cox.mod)){
      est = summary(cox.mod)$coefficients[1]
      SE = summary(cox.mod)$coefficients[3]
      LCL = confint(cox.mod, level=1-alpha_s)[1]
      UCL = confint(cox.mod, level=1-alpha_s)[2]
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
  # Combine results and estimate effect in overall population #
  if ( sum(is.na(param.dat$est))>0 | length(unique(Subgrps))==1  ){
    cox.mod = coxph(Y ~ A , data=indata)
    param.dat0 = data.frame(Subgrps = 0, N = dim(X)[2],
                            est = summary(cox.mod)$coefficients[1],
                            SE = summary(cox.mod)$coefficients[3],
                            LCL = confint(cox.mod, level=1-alpha_ovrl)[1],
                            UCL = confint(cox.mod, level=1-alpha_ovrl)[2],
                            pval = summary(cox.mod)$coefficients[5] )
  }
  if ( sum(is.na(param.dat$est))==0){
    param.dat0 = param_combine(param.dat = param.dat, alpha_ovrl=alpha_ovrl,
                               combine=combine)
  }
  param.dat = rbind(param.dat0, param.dat)
  # convert log(HR) to HR #
  param.dat$est = exp(param.dat$est)
  param.dat$LCL = exp(param.dat$LCL)
  param.dat$UCL = exp(param.dat$UCL)
  param.dat$estimand = "HR(A=1 vs A=0)"
  param.dat = param.dat[,c("Subgrps", "N", "estimand", "est", "SE", "LCL", "UCL", "pval")]
  return( param.dat )
}
