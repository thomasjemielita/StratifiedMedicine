#' Filter: Random Forest Variable Importance
#'
#' Filtering through Random Forest Variable Importance with p-values.
#' Default is to remove variables with p-values >= 0.10.
#' Used for continuous, binary, or survival outcomes
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate matrix. Must be numeric.
#' @param b Subsample size (n^b)
#' @param K Number of samples (default=200)
#' @param DF2 2-DF test statistic (default=FALSE)
#' @param FDR FDR correction for p-values (default=FALSE)
#' @param pval.thres p-value threshold for filtering (default=0.10)
#' @param family Outcome type ("gaussian", "binomial", "survival"), default is "gaussian"
#' @param ... Any additional parameters, not currently passed through.
#'
#' @return Filter model and variables that remain after filtering.
#'  \itemize{
#'   \item mod - Filtering model
#'   \item filter.vars - Variables that remain after filtering (could be all)
#' }
#'
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
#' \donttest{
#' mod1 = filter_RF(Y, A, X, K=200) # Same as default #
#' mod1$filter.vars
#' mod1$mod # summary of variable importance outputs
#' }
#'
#'
#' ## Survival (TBD) ##
#'
#' @seealso \code{\link{PRISM}}, \code{\link{ranger}}
#'
##### RF Variable Importance ######
filter_RF = function(Y, A, X, b=0.66, K=200, DF2=FALSE, FDR=FALSE, pval.thres=0.10,
                     family="gaussian", ...){

  if (DF2==TRUE){ #Generate the interaction covariates #
    X_inter = X*A
    colnames(X_inter) = paste(colnames(X), "_A", sep="")
    W = cbind(X, A, X_inter)
  }
  if (DF2==FALSE){
    W = cbind(X)
  }

  ## Calculcate observed variable importance ##
  set.seed(15314)
  mod0 <- ranger(Y ~ ., data = data.frame(Y,W), importance = "permutation")
  VI0 = as.numeric(mod0$variable.importance)

  ### Subsample function ###
  b1 = ( dim(X)[1] )^(b)
  sub_VI = function(s){
    ## Subsample data ##
    hold = data.frame(Y, W)
    seedz = s+3141
    set.seed(seedz)
    hold = hold[sample(nrow(hold), size=b1),]
    mod.s <- ranger(Y ~ ., data = hold, importance = "permutation")
    VI = c( as.numeric(mod.s$variable.importance) )
    return(VI)
  }
  res = lapply(1:K, sub_VI)
  res = do.call(rbind, res)
  VI.var = NULL
  for (j in 1:dim(res)[2]){
    VI.var = c(VI.var, b1/((length(Y)-b1)*K) * sum( (res[,j] - VI0[j])^2, na.rm = TRUE ) )
  }
  #### Initial Output of VI and SE(VI) for each variable ###
  out = data.frame(Variables = colnames(W), est = VI0, SE = sqrt(VI.var) )
  if (DF2==FALSE){
    out.F = out
    ### T-Statistics and one-sided p-values ###
    out.F$Tstat = with(out.F, est / SE )
    out.F$pval = with(out.F, 1-pnorm(Tstat) )
    out.F$pval.fdr = p.adjust(out.F$pval, "fdr")
  }
  if (DF2==TRUE){
    out.F = data.frame(Variables=colnames(X), est.2DF=NA, SE.2DF=NA)
    for (var in colnames(X)){
      #print(var)
      #var = "X187"
      inter.name = paste(var, "_A", sep="")
      ## Extract VI estimates and SEs ##
      est0 = out$est[out$Variables==var]; SE0 = out$SE[out$Variable==var]
      estI = out$est[out$Variables==inter.name]; SEI = out$SE[out$Variable==inter.name]
      ## Calculate empirical covariance between main/interaction VIs ##
      indices = match(c(var, inter.name), colnames(W))
      cov.est = b1/((length(Y)-b1)*K) *
        sum( (res[,indices[1]]-est0)*(res[,indices[2]]-estI), na.rm=TRUE)
      ## Calculcate "2 DF" VI Estimate and SE ##
      est.F = estI+est0
      SE.F = sqrt( SE0^2 + SEI^2 + 2*cov.est )
      ## Store ##
      out.F$est.2DF[out.F$Variables==var] =  est.F
      out.F$SE.2DF[out.F$Variables==var] =  SE.F
    }
    ### T-Statistics and One-sided Pval ###
    out.F$Tstat = with(out.F, est.2DF / SE.2DF)
    out.F$pval = with(out.F, 1-pnorm(Tstat) )
    out.F$pval.fdr = p.adjust(out.F$pval, "fdr")
  }

  ## FDR Adjustment? ##
  if (FDR){
    out.F$filter_vars = with(out.F, ifelse(pval.fdr<pval.thres, 1, 0) )
  }
  if (!FDR){
    out.F$filter_vars = with(out.F, ifelse(pval<pval.thres, 1, 0) )
  }
  # ggplot(out.F[out.F$pval<0.05,], aes(x=reorder(Variables,-pval), y=est)) +
  #   geom_errorbar(aes(ymin=est-1.96*SE, ymax=est+1.96*SE), width=.1)+geom_point()+
  #   ylab("Variable Importance")+
  #   xlab("Variables")+
  #   ggtitle("Variable Importance with 95% CI")+coord_flip()
  ## Return VI information and variables that we filter ##
  mod = out.F
  filter.vars = as.character( mod$Variables[mod$filter_vars==1] )
  return( list(mod=mod, filter.vars=filter.vars) )
}
