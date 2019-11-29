#' Subgroup Identification: Model-based partitioning (Weibull)
#'
#' Uses the MOB (with weibull loss function) algorithm to identify subgroups
#' (Zeileis, Hothorn, Hornik 2008; Seibold, Zeileis, Hothorn 2016). Usable for
#' survival outcomes.
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
#' @import survival
#'
#' @return Trained MOB (Weibull) model.
#'  \itemize{
#'   \item mod - MOB (Weibull) model object
#' }
#'
#' @export
#' @examples
#' library(StratifiedMedicine)
#'
#'
#' \donttest{
#' ## Load TH.data (no treatment; generate treatment randomly to simulate null effect) ##
#' data("GBSG2", package = "TH.data", envir = e <- new.env() )
#' surv.dat = e$GBSG2
#' ## Design Matrices ###
#' Y = with(surv.dat, Surv(time, cens))
#' X = surv.dat[,!(colnames(surv.dat) %in% c("time", "cens")) ]
#' A = rbinom( n = dim(X)[1], size=1, prob=0.5  )
#' res_weibull = submod_weibull(Y, A, X, Xtest=X, family="survival")
#' plot(res_weibull$mod)
#' }
#'
#'
## MOB: Weibull ##
submod_weibull = function(Y, A, X, Xtest, mu_train, alpha=0.05,
                          minsize = floor( dim(X)[1]*0.10  ),
                          maxdepth = 4, parm=NULL, ...){

  ## Fit Model ##
  mod <- mob(Y ~ A | ., data = X,
             fit = wbreg, control = mob_control(alpha=alpha, minsize=minsize,
                                                maxdepth=maxdepth, parm=parm))
  # Prediction Function #
  pred.fun <- function(mod, X=NULL, type="subgrp"){
    pred <- NULL
    Subgrps <- as.numeric( predict(mod, type="node", newdata = X) )
    if (type=="all"){
      L.mat <- rbind( c(1,0), c(1,1) )
      pred <- data.frame(Subgrps=Subgrps, mu0 = NA, mu1 = NA)
      for (s in unique(Subgrps)){
        hold <- summary(mod)[[as.character(s)]]
        hold.est <- exp( as.numeric( L.mat %*% coef(hold)  ) )
        pred$mu0[pred$Subgrps==s] <- hold.est[1]
        pred$mu1[pred$Subgrps==s] <- hold.est[2]
      } 
    }
    return(list(Subgrps=Subgrps, pred=pred))
  }
  res <- list(mod=mod, pred.fun=pred.fun)
  class(res) <- "submod_weibull"
  ## Return Results ##
  return(  res  )
}