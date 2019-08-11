#' Parameter Estimation: Restricted Mean Survival Time (RMST)
#'
#' For each identified subgroup, estimate the restricted mean survival time (RMST), based
#' on the method described in the R package "survRM2". Point-estimates and
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
#' # Survival Data #
#' require(TH.data); require(coin)
#' data("GBSG2", package = "TH.data")
#' surv.dat = GBSG2
#' # Design Matrices ###
#' Y = with(surv.dat, Surv(time, cens))
#' X = surv.dat[,!(colnames(surv.dat) %in% c("time", "cens")) ]
#' A = rbinom( n = dim(X)[1], size=1, prob=0.5  ) ## simulate null treatment
#'
#' # MOB-Weibull Subgroup Model ##
#' res_weibull = submod_train(Y, A, X, Xtest=X, family="survival",
#'                            submod = "submod_weibull")
#' plot(res_weibull$mod)
#'
#' # Parameter-Estimation ##
#' require(survRM2)
#' params = param_rmst(Y, A, X, Subgrps = res_weibull$Subgrps.train, alpha_ovrl=0.05,
#'                    alpha_s=0.05)
#' params
#' }
#'
#' @seealso \code{\link{param_combine}}
#'
### RMST ###
param_rmst = function(Y, A, X, mu_hat, Subgrps, alpha_ovrl, alpha_s, combine="adaptive",
                      ...){

  if (!requireNamespace("survRM2", quietly = TRUE)) {
    stop("Package survRM2 needed for param_rmst. Please install.")
  }

  indata = data.frame(Y=Y, A=A, X)
  ### Loop through subgroups ##
  looper = function(s){
    time = indata$Y[Subgrps %in% s,1]
    status = indata$Y[Subgrps %in% s,2]
    arm = indata$A[Subgrps %in% s]
    n.s = length(arm)
    obj = tryCatch( survRM2::rmst2(time, status, arm),
                    error = function(e) "param error" )
    if (is.character(obj)){
      est = NA; SE = NA; pval = NA; LCL = NA; UCL = NA;
    }
    if (is.list(obj)){
      est = obj$unadjusted.result[1,1]
      SE = sqrt( obj$RMST.arm1$rmst.var + obj$RMST.arm0$rmst.var )
      LCL =  est - qnorm(1-alpha_s/2)*SE
      UCL =  est + qnorm(1-alpha_s/2)*SE
      pval = obj$unadjusted.result[1,4]
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
    param.dat0 = looper(s=S_levels)
    param.dat0 = data.frame( Subgrps=0, N=dim(X)[2],
                             param.dat0 )
  }
  if ( sum(is.na(param.dat$est))==0){
    param.dat0 = param_combine(param.dat = param.dat, alpha_ovrl=alpha_ovrl,
                               combine=combine)
  }
  param.dat = rbind(param.dat0, param.dat)
  param.dat$estimand = "RMST(1-0)"
  param.dat = param.dat[,c("Subgrps", "N", "estimand", "est", "SE", "LCL", "UCL", "pval")]
  return( param.dat )
}
