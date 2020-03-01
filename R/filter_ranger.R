#' Filter: Random Forest (ranger) Variable Importance
#'
#' Filtering through Random Forest Variable Importance with p-values. P-values are obtained
#' through subsampling based T-statistics, as described in Ishwaran and Lu 2017. Default is
#' to remove variables if one-sided (VI>0) p-values >= 0.10. Used for continuous, binary, 
#' or survival outcomes.
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate space.
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
#' @references Wright, M. N. & Ziegler, A. (2017). ranger: A fast implementation of 
#' random forests for high dimensional data in C++ and R. J Stat Softw 77:1-17. 
#' \url{https://doi.org/10.18637/jss.v077.i01}.
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
#' mod1 = filter_ranger(Y, A, X, K=200) # Same as default #
#' mod1$filter.vars
#' mod1$mod # summary of variable importance outputs
#' }
#'
#'
##### RF Variable Importance ######
filter_ranger = function(Y, A, X, b=0.66, K=200, DF2=FALSE, FDR=FALSE, pval.thres=0.10,
                     family="gaussian", ...) {

  if (DF2==TRUE) { #Generate the interaction covariates #
    A_lvls <- unique(A)[order(unique(A))]
    ## Generate the interaction covariates ###
    A.mat <- model.matrix(~., data=data.frame(A))[,-1]
    X_inter = X*A.mat
    colnames(X_inter) = paste(colnames(X), A_lvls[2], sep="_")
    W = cbind(A=A.mat, X, X_inter)
    # X_inter = X*A
    # colnames(X_inter) = paste(colnames(X), "_A", sep="")
    # W = cbind(X, A, X_inter)
  }
  if (DF2==FALSE) {
    W = cbind(X)
  }

  ## Calculcate observed variable importance ##
  mod0 <- ranger(Y ~ ., data = data.frame(Y,W), importance = "permutation")
  VI0 = as.numeric(mod0$variable.importance)

  ### Subsample function ###
  b1 = ( dim(X)[1] )^(b)
  sub_VI = function(s) {
    ## Subsample data ##
    hold = data.frame(Y, W)
    hold = hold[sample(nrow(hold), size=b1),]
    mod.s <- ranger(Y ~ ., data = hold, importance = "permutation")
    VI = c( as.numeric(mod.s$variable.importance) )
    return(VI)
  }
  res = lapply(1:K, sub_VI)
  res = do.call(rbind, res)
  VI.var = NULL
  for (j in 1:dim(res)[2]){
    VI.var = c(VI.var, b1/((length(Y)-b1)*K) * sum((res[,j] - VI0[j])^2, na.rm = TRUE ))
  }
  #### Initial Output of VI and SE(VI) for each variable ###
  out = data.frame(Variables = colnames(W), est = VI0, SE = sqrt(VI.var) )
  if (DF2==FALSE) {
    out.F = out
    ### T-Statistics and one-sided p-values ###
    out.F$Tstat = with(out.F, est / SE)
    out.F$pval = with(out.F, 1-pnorm(Tstat))
    out.F$pval.fdr = p.adjust(out.F$pval, "fdr")
  }
  if (DF2==TRUE) {
    out.F = data.frame(Variables=colnames(X), est.2DF=NA, SE.2DF=NA)
    for (var in colnames(X)) {
      inter.name = paste(var, A_lvls[2], sep="_")
      # Extract VI estimates and SEs #
      est0 = out$est[out$Variables==var]; SE0 = out$SE[out$Variable==var]
      estI = out$est[out$Variables==inter.name]; SEI = out$SE[out$Variable==inter.name]
      # Calculate empirical covariance between main/interaction VIs #
      indices = match(c(var, inter.name), colnames(W))
      cov.est = b1/((length(Y)-b1)*K) *
        sum((res[,indices[1]]-est0)*(res[,indices[2]]-estI), na.rm=TRUE)
      # Calculcate "2 DF" VI Estimate and SE #
      est.F = estI+est0
      SE.F = sqrt( SE0^2 + SEI^2 + 2*cov.est )
      # Store #
      out.F$est.2DF[out.F$Variables==var] =  est.F
      out.F$SE.2DF[out.F$Variables==var] =  SE.F
    }
    # T-Statistics and One-sided Pval #
    out.F$Tstat = with(out.F, est.2DF/SE.2DF)
    out.F$pval = with(out.F, 1-pnorm(Tstat))
    out.F$pval.fdr = p.adjust(out.F$pval, "fdr")
    out.F$est <- out.F$est.2DF
    out.F$SE <- out.F$SE.2DF
  }

  # FDR Adjustment? #
  if (FDR) {
    out.F$filter_vars = with(out.F, ifelse(pval.fdr<pval.thres, 1, 0))
  }
  if (!FDR) {
    out.F$filter_vars = with(out.F, ifelse(pval<pval.thres, 1, 0))
  }
  out.F$LCL <- with(out.F, est-qnorm(1-pval.thres)*SE)
  out.F$UCL <- with(out.F, est+qnorm(1-pval.thres)*SE)
  filter.vars = with(out.F, as.character(Variables[filter_vars==1]))
  out.F$rank <- with(out.F, rank(pval))
  out.F$pval.thres <- pval.thres
  VI.dat <- out.F
  VI.dat$family <- family
  mod = list(VI.dat=VI.dat)
  # Return VI information and variables that we filter #
  return(list(mod=mod, filter.vars=filter.vars) )
}
# VIMP Plot #
plot_vimp_ranger <- function(mod, top_n=NULL) {
  VI.dat <- mod$VI.dat
  family <- unique(VI.dat$family)
  pval.thres <- unique(VI.dat$pval.thres)
  plot.title <- paste("Random Forest (", family,
                      ") Importance Plot", sep="")
  y.label <- paste("Variable Importance ", 
                   "(", (1-pval.thres)*100, "% CI)",sep="")
  if (!is.null(top_n)) {
    VI.dat <- VI.dat[VI.dat$rank<=top_n,]
    plot.title <- paste(plot.title, " [Top ", top_n, "]", sep="")
  }
  vimp.plt <- ggplot2::ggplot(VI.dat, aes(x=reorder(.data$Variables,-.data$pval), 
                                 y=.data$est, ymin=.data$LCL, ymax=.data$UCL)) +
    ggplot2::geom_errorbar(width=.1)+geom_point()+
    ggplot2::geom_hline(yintercept=0, color="red")+
    ggplot2::ylab(y.label)+
    ggplot2::xlab("Variable")+
    ggplot2::ggtitle(plot.title)+
    ggplot2::coord_flip()
  return(vimp.plt)
}
