#' Generate Subgroup Data-sets
#'
#' Simulation/real data-sets; useful for testing new models and PRISM configurations.
#'
#' @param family Outcome type ("gaussian", "binomial", "survival")
#' @param ... Any additional parameters, not currently passed through.
#'
#' @importFrom stats rbinom rnorm time
#' @importFrom utils data
#' @import mvtnorm
#' @import TH.data
#' @return Simulation data set (Y=outcome, A=treatment, X=covariates)
#' @export
#'

generate_subgrp_data = function(family, ...){

  # Sample size #
  n = 800
  ### Covariate Space (P=50) ####
  jd = function(n,m) matrix(c(1),nrow=n,ncol=m)
  Sigma_X = diag(50)+jd(50,50)*0.10 - diag(50)*0.10
  X = data.frame(rmvnorm(n=n, sigma=Sigma_X) )
  #### Treatment Variable: Binary #####
  A = rbinom(n=n, size=1, prob=0.5)
  ## Create subgroup based on two continuous covariates ##
  X1_cut = with(X, ifelse(X1>qnorm(0.80),1,0) )
  X2_cut = with(X, ifelse( X2<qnorm(0.42), 1, 0) )
  ## Subgroups ##
  subgrps = NA
  subgrps = ifelse( X1_cut==0 & X2_cut==0, "[X1- X2-]", subgrps)
  subgrps = ifelse( X1_cut==0 & X2_cut==1, "[X1- X2+]", subgrps)
  subgrps = ifelse( X1_cut==1 & X2_cut==0, "[X1+ X2-]", subgrps)
  subgrps = ifelse( X1_cut==1 & X2_cut==1, "[X1+ X2+]", subgrps)
  if (family=="gaussian"){
    ## Mean Function ##
    mu.e = with(X, 1.5 +   A*( 0.40*ifelse(subgrps=="[X1+ X2+]", 1, 0)+
                                 0.35*ifelse(subgrps=="[X1+ X2-]", 1, 0)+
                                 0.30*ifelse(subgrps=="[X1- X2+]", 1, 0) )+
                  0.10*X1 - 0.20*(X2-mean(X2) )/sd(X2) + 0.15*(X3-mean(X3))/sd(X3)+
                  0.08*(X5-mean(X5) ) / sd(X5)  + 0.09*(X7-mean(X7))/ sd(X7) )
    ## Outcome: Normal ##
    Y = rnorm(n = length(mu.e), mean = mu.e, sd = 1)
  }
  if (family=="binomial"){
    ## TBD ##
  }
  if (family=="survival"){
    ## TBD ##
  }
  return( list(Y=Y, A=A, X=X) )
}
