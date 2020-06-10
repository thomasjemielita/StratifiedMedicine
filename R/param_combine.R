#' Overall Population Estimate: Aggregating Subgroup-Specific Parameter Estimates
#'
#' Function that combines subgroup-specific estimates to obtain an overall population
#' estimate. Options including sample size weighting and adaptive weighting (default; as
#' described in Marceau-West and Mehrotra (to appear)).
#'
#' @param param.dat Parameter data-set with subgroup-specific point estimates, SEs, and
#' sample sizes.
#' @param combine Method to combine subgroup-specific estimates. Default is "adaptive".
#' combine="SS" uses sample size weighting.
#' @param alpha_ovrl Two-sided alpha level for overall population. Default=0.05
#' @param ... Any additional parameters, not currently passed through.
#'
#' @importFrom stats cor
#' @return Data-frame with overall population point estimate, SE, and CI
#' @export

param_combine = function(param.dat, combine="SS", alpha_ovrl=0.05, ...){

  # Set weights #
  if (combine=="adaptive"){
    est.corr = cor(param.dat$est, param.dat$SE^2, method="spearman")
    if (est.corr<0){ w.vec = param.dat$N }
    if (est.corr>=0){ w.vec = param.dat$N / param.dat$SE }
    w.vec = w.vec / sum(w.vec)
  }
  if (combine=="SS"){
    w.vec = param.dat$N / sum(param.dat$N)
  }
  # Combine based on weights #
  beta0 = as.numeric( w.vec %*% param.dat$est )
  SE.beta0 = sqrt(  sum(w.vec^2*param.dat$SE^2)    )
  beta0.LCL = beta0 - qt(1-alpha_ovrl/2, df=sum(param.dat$N)-1) * SE.beta0
  beta0.UCL = beta0 + qt(1-alpha_ovrl/2, df=sum(param.dat$N)-1) * SE.beta0
  # Return results #
  param.dat0 = data.frame(Subgrps="ovrl", N = sum(param.dat$N), 
                          est = beta0, SE = SE.beta0,
                          LCL = beta0.LCL, UCL = beta0.UCL)
  param.dat0$pval = with(param.dat0, 2*pt(-abs(est/SE), df=N-1) )
  return( param.dat0 )
}
