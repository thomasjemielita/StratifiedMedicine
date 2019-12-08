#' Parameter Estimation: Double-robust estimator
#'
#' For each identified subgroup and in the overall population, use the double robust
#' estimator (Funk et al 2011). For continuous and binary outcomes, this outputs 
#' estimates for E(Y|A=1), E(Y|A=0), and E(Y|A=1)-E(Y|A=0). For survival, point 
#' estimates correspond to E(logT|A=1), E(logT|A=0), and E(logT|A=1)-E(logT|A=0). 
#' Here, patient-level estimates (mu_hat) must correspond to the accelerated failure 
#' time (AFT) survival model framework (default uses AFT-based bart, "abart").
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
#' for overall and subgroups. Subgrps=0 corresponds to the overall population by 
#' default.
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
#' ## Estimate PLEs (ranger) ##
#' res_ranger = ple_train(Y, A, X, Xtest=X, ple="ple_ranger")
#'
#' ## Identify Subgroups: MOB (lmtree) ##
#' res_lmtree = submod_train(Y, A, X, Xtest=X, submod="submod_lmtree")
#'
#' ## Parameter-estimation ##
#' params = param_dr(Y, A, X, mu_hat = res_ranger$mu_train,
#'                   Subgrps = res_lmtree$Subgrps.train, alpha_ovrl=0.05,
#'                   alpha_s=0.05)
#' params
#'

### AIPTW (Double-Robust) Param ###
param_dr = function(Y, A, X, mu_hat, Subgrps, alpha_ovrl, alpha_s, ...){
  
  if (is.null(A)){
    stop("param_dr not applicable for no treatment (A=NULL)") 
  }
  if (is.Surv(Y)){
    estimands <- c("E(logT|A=0)", "E(logT|A=1)", "E(logT|A=1)-E(logT|A=0)")
  }
  if (!is.Surv(Y)){
    estimands <- c("E(Y|A=0)", "E(Y|A=1)", "E(Y|A=1)-E(Y|A=0)")
  }
  indata <- data.frame(Y=Y, A=A, X)
  # Subgroup and overall estimates #
  looper <- function(s, alpha){
    Y.s <- indata$Y[Subgrps %in% s]
    if (is.Surv(Y)){
      Y.s = log(Y.s[,1])
    }
    A.s <- indata$A[Subgrps %in% s]
    n.s <- length(Y.s)
    probA <- mean(A.s)
    mu.s <- mu_hat[Subgrps %in% s,]
    # EIF ###
    eif.0 = ( (1-A.s)*Y.s + (A.s-probA)*mu.s$mu0 ) / (1-probA)
    eif.1 = ( A.s*Y.s - (A.s-probA)*mu.s$mu1 )/ probA
    eif = eif.1 - eif.0
    # Double robust estimator: Average eifs #
    est = c( mean(eif.0, na.rm=TRUE), 
             mean(eif.1, na.rm=TRUE),
             mean(eif, na.rm=TRUE) )
    # EIF for variance estimate #
    SE = sqrt( n.s^(-2) * c( sum( (eif.0-est[1])^2 ),
                             sum( (eif.1-est[2])^2 ),
                             sum( (eif-est[3])^2 ) ) )
    LCL = est-qt( (1-alpha/2), df=n.s-1 )*SE
    UCL = est+qt( (1-alpha/2), df=n.s-1 )*SE
    pval = 2*pt(-abs(est/SE), df=n.s-1)
    summ = data.frame( Subgrps = ifelse( mean(S_levels %in% s)==1, 0, s),
                       N = n.s, 
                       estimand = estimands,
                       est, SE, LCL, UCL, pval)
    return( summ )
  }
  S_levels = as.numeric( names(table(Subgrps)) )
  ## Fit Overall and across subgroups ##
  param.dat0 = looper(s = S_levels, alpha = alpha_ovrl)
  if (length(unique(S_levels))>1){
    param.dat = lapply(S_levels, looper, alpha = alpha_s)
    param.dat = do.call(rbind, param.dat)
    param.dat = data.frame( param.dat )
    param.dat = rbind(param.dat0, param.dat)
  }
  if (length(unique(S_levels))==1){
    param.dat = param.dat0
  }
  return( param.dat )
}
