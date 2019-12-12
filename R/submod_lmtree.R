#' Subgroup Identification: Model-based partitioning (lmtree)
#'
#' Uses the lmtree (model-based partitioning, OLS) algorithm to identify subgroups
#' (Zeileis, Hothorn, Hornik 2008). Usable for continuous and binary outcomes.
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate space.
#' @param Xtest Test set
#' @param mu_train Patient-level estimates (See PLE_models)
#' @param alpha Significance level for variable selection (default=0.05)
#' @param minsize Minimum number of observations in a tree node.
#' Default = floor( dim(train)[1]*0.05  )
#' @param maxdepth Maximum depth of any node in the tree (default=4)
#' @param parm Model parameters included in parameter instability tests 
#' (default=NULL, all parameters)
#' @param ... Any additional parameters, not currently passed through.
#'
#' @import partykit
#'
#' @return Trained lmtree model.
#'  \itemize{
#'   \item mod - lmtree model object
#' }
#'
#' @export
#' @examples
#' 
#' \donttest{
#' library(StratifiedMedicine)
#'
#' ## Continuous ##
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#' train = data.frame(Y, A, X)
#' # Outcome/treatment must be labeled as Y/A #
#'
#' res_lmtree1 = submod_lmtree(Y, A, X, Xtest=X)
#' res_lmtree2 = submod_lmtree(Y, A, X, Xtest=X, maxdepth=2, minsize=100)
#' plot(res_lmtree1$mod)
#' plot(res_lmtree2$mod)
#'}
#'
#### lmtree (MOB) ###
submod_lmtree = function(Y, A, X, Xtest, mu_train, alpha=0.05,
                         minsize = floor( dim(X)[1]*0.10  ),
                         maxdepth = 4, parm=NULL, ...){

  ## Fit Model ##
  mod <- lmtree(Y~A | ., data = X, alpha=alpha, maxdepth = maxdepth, 
                minsize=minsize, parm=parm)
  
  # Prediction Function #
  pred.fun <- function(mod, X=NULL, type="subgrp"){
    Subgrps <- NULL; pred <- NULL
    Subgrps <- as.numeric( predict(mod, type="node", newdata = X) )
    if (type=="all"){
      L.mat <- rbind( c(1,0), c(1,1) )
      pred <- data.frame(Subgrps=Subgrps, mu0 = NA, mu1 = NA)
      for (s in unique(Subgrps)){
        hold <- suppressWarnings(  summary(mod)[[as.character(s)]] )
        hold.est <- as.numeric( L.mat %*% coef(hold)  ) 
        pred$mu0[pred$Subgrps==s] <- hold.est[1]
        pred$mu1[pred$Subgrps==s] <- hold.est[2]
      }
      pred$PLE <- with(pred, mu1-mu0) 
    }
    return(list(Subgrps=Subgrps, pred=pred))
  }
  res <- list(mod=mod, pred.fun=pred.fun)
  class(res) <- "submod_lmtree"
  ## Return Results ##
  return(  res )
}
