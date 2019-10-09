#' Generate Subgroup Data-sets
#'
#' Simulation/real data-sets; useful for testing new models and PRISM configurations.
#'
#' @param n sample size (default=800)
#' @param seed seed number (default=513413)
#' @param family Outcome type ("gaussian", "binomial", "survival")
#' @param null Simulate null hypothesis of no treatment effect and no subgruops. Default
#' is FALSE.
#' @param ... Any additional parameters, not currently passed through.
#'
#' @importFrom stats rbinom rnorm time
#' @importFrom utils data
#' @import mvtnorm
#' @return Simulation data set (Y=outcome, A=treatment, X=covariates)
#' @export
#'

generate_subgrp_data = function(n=800, seed=513413, family, null=FALSE, ...){

  # Sample size #
  set.seed(seed)
  ### Covariate Space (P=50) ####
  jd = function(n,m) matrix(c(1),nrow=n,ncol=m)
  Sigma_X = diag(50)+jd(50,50)*0.10 - diag(50)*0.10
  X = data.frame(rmvnorm(n=n, sigma=Sigma_X) )
  #### Treatment Variable: Binary #####
  A = rbinom(n=n, size=1, prob=0.5)
  ## Create subgroup based on two continuous covariates ##
  X1_cut = with(X, ifelse(X1>qnorm(0.50),1,0) )
  X2_cut = with(X, ifelse( X2>qnorm(0.50), 1, 0) )
  ## Subgroups ##
  subgrps = NA
  subgrps = ifelse( X1_cut==0 & X2_cut==0, "[X1- X2-]", subgrps)
  subgrps = ifelse( X1_cut==0 & X2_cut==1, "[X1- X2+]", subgrps)
  subgrps = ifelse( X1_cut==1 & X2_cut==0, "[X1+ X2-]", subgrps)
  subgrps = ifelse( X1_cut==1 & X2_cut==1, "[X1+ X2+]", subgrps)
  if (family=="gaussian"){
    if (!null){
      trt.effect = A*( 0.65*ifelse(subgrps=="[X1+ X2+]", 1, 0)+
          0.45*ifelse(subgrps=="[X1+ X2-]", 1, 0)+
          0.25*ifelse(subgrps=="[X1- X2+]", 1, 0) )
    }
    if (null){
      trt.effect = 0
    }
    ## Mean Function ##
    mu.e = with(X, 1.5 + trt.effect +
                  0.15*(X1-mean(X1))/sd(X1) + 0.20*(X2-mean(X2) )/sd(X2) +
                  0.10*(X3-mean(X3))/sd(X3)+
                  0.10*(X5-mean(X5) ) / sd(X5)  + 0.10*(X7-mean(X7))/ sd(X7) )
    ## Outcome: Normal ##
    Y = rnorm(n = length(mu.e), mean = mu.e, sd = 1)
  }
  if (family=="binomial"){
    ## TBD ##
  }
  if (family=="survival"){
    ## TBD ##
  }
  return( list(Y=Y, A=A, X=X, subgrps = subgrps) )
}
