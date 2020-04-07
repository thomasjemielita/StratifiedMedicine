#' Parameter Estimation: Patient-Level Estimates
#'
#' For each identified subgroup and in the overall population, average the patient-level
#' estimates of E(Y|A=1), E(Y|A=0), and E(Y|A=1)-E(Y|A=0). Pseudo-outcomes are 
#' used for variance estimates (Jemielita and Mehrotra 2019).
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
#'
#' @references Jemielita T, Mehrotra D. PRISM: Patient Response Identifiers for 
#' Stratified Medicine. \url{https://arxiv.org/abs/1912.03337}
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
param_ple = function(Y, A, X, mu_hat, Subgrps, alpha_ovrl, alpha_s, ...) {

  noA = FALSE
  if (is.null(A)) {
    A_num = rep(1, length(Y))
    noA = TRUE
  }
  if (!noA) {
    A_lvls <- unique(A)[order(unique(A))]
    E_A0 <- paste("E(Y|A=", A_lvls[1], ")", sep="")
    E_A1 <- paste("E(Y|A=", A_lvls[2], ")", sep="")
    E_diff <- paste(E_A1, "-", E_A0, sep="")
    estimands <- c(E_A0, E_A1, E_diff)
    mu_A0 <- paste("mu", A_lvls[1], sep="_")
    mu_A1 <- paste("mu", A_lvls[2], sep="_")
    A_num <- model.matrix(~., data=data.frame(A))[,-1] 
  }
  
  indata = data.frame(Y=Y, A=A_num, X)
  ## Estimate overall and within subgroups ##
  looper = function(s, alpha){
    Y.s = indata$Y[Subgrps %in% s]
    A.s = indata$A[Subgrps %in% s]
    n.s = length(Y.s)
    probA = mean(A.s)
    mu.s = mu_hat[Subgrps %in% s,]
    if (noA){
      est = mean(mu.s)
      SE = sqrt( n.s^(-2)*sum((est-Y)^2) ) 
      LCL = est-qt( (1-alpha/2), df=n.s-1 )*SE
      UCL = est+qt( (1-alpha/2), df=n.s-1 )*SE
      pval = 2*pt(-abs(est/SE), df=n.s-1)
      summ = data.frame( Subgrps = ifelse(n.s==length(Y), 0, s),
                         N = n.s, estimand = "E(Y)", est, SE, LCL, UCL, pval)
    }
    if (!noA){
      # PLEs for point-estimates #)
      est = c( mean(mu.s[,mu_A0], na.rm=TRUE), 
               mean(mu.s[,mu_A1], na.rm=TRUE),
               mean(mu.s$PLE, na.rm=TRUE))
      # EIF for variance estimate #
      eif.0 = ( (1-A.s)*Y.s + (A.s-probA)*mu.s[,mu_A0] ) / (1-probA)
      eif.1 = ( A.s*Y.s - (A.s-probA)*mu.s[,mu_A1] )/ probA
      eif = eif.1 - eif.0
      SE = sqrt( n.s^(-2) * c( sum( (eif.0-est[1])^2), 
                               sum( (eif.1-est[2])^2),
                               sum( (eif-est[3])^2 )   ) )
      LCL = est-qt( (1-alpha/2), df=n.s-1 )*SE
      UCL = est+qt( (1-alpha/2), df=n.s-1 )*SE
      pval = 2*pt(-abs(est/SE), df=n.s-1)
      summ = data.frame( Subgrps = ifelse(n.s==length(Y), 0, s),
                         N = n.s, 
                         estimand = estimands,
                         est, SE, LCL, UCL, pval)
    }
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
