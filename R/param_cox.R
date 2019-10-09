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
#' @return Data-set with parameter estimates (log hazard ratio) and corresponding
#' variability metrics, for overall and subgroups. Subgrps=0 corresponds to the overall
#' population by default.
#'  \itemize{
#'   \item param.dat - Parameter estimates and variability metrics (est=logHR, 
#'   SE=SE(logHR), LCL/UCL = lower/upper confidence limit on logHR scale, pval = p-value).
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
  if (is.null(A)){
    stop("param_cox not applicable for no treatment (A=NULL)")
  }
  indata = data.frame(Y=Y, A=A, X)
  ### Loop through subgroups ##
  looper = function(s, alpha){
    n.s = dim(indata[Subgrps %in% s,])[1]
    ## Extract HR, SE, 95% CI, and p-value for Subgroup Specific Treatment Effect ##
    cox.mod = tryCatch( coxph(Y ~ A , data=indata[Subgrps %in% s,]),
                        error = function(e) "fit error",
                        warning = function(w) "convergence issues")
    if (is.character(cox.mod)){
      summ = data.frame(Subgrps = ifelse(n.s==dim(indata)[1], 0, s), N=n.s,
                        est=NA, SE=NA, LCL=NA, UCL=NA, 
                        pval=NA)
    }
    if (is.list(cox.mod)){
      est = summary(cox.mod)$coefficients[1]
      SE = summary(cox.mod)$coefficients[3]
      LCL = confint(cox.mod, level=1-alpha)[1]
      UCL = confint(cox.mod, level=1-alpha)[2]
      pval = summary(cox.mod)$coefficients[5]
      summ = data.frame( Subgrps = ifelse(n.s==dim(indata)[1], 0, s),
                         N = n.s, est, SE, LCL, UCL, pval)
    }
    return( summ )
  }
  # Across Subgroups #
  S_levels <- as.numeric( names(table(Subgrps)) )
  param.dat <- lapply(S_levels, looper, 
                      alpha = ifelse( length(unique(Subgrps))==1, alpha_ovrl, alpha_s))
  param.dat <- do.call(rbind, param.dat)
  param.dat <- data.frame( param.dat )
  # Combine results and estimate effect in overall population #
  if ( sum(is.na(param.dat$est))>0){
    if (length(unique(Subgrps))>1){
      param.dat0 <- looper(s=S_levels, alpha=alpha_ovrl)
      param.dat = rbind(param.dat0, param.dat)
    }
  }
  if ( sum(is.na(param.dat$est))==0){
    if (length(unique(Subgrps))>1){
      param.dat0 <- param_combine(param.dat = param.dat, alpha_ovrl=alpha_ovrl,
                                  combine=combine)
      param.dat <- rbind(param.dat0, param.dat) 
    }
  }
  param.dat$estimand = "logHR(A=1 vs A=0)"
  param.dat = param.dat[,c("Subgrps", "N", "estimand", "est", "SE", "LCL", "UCL", "pval")]
  return( param.dat )
}
