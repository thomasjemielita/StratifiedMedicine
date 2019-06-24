#' Subgroup Identification: Model-based partitioning (Weibull)
#'
#' Uses the MOB (with weibull loss function) algorithm to identify subgroups.
#' Usable for survival outcomes.
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate matrix. Must be numeric.
#' @param Xtest Test set
#' @param mu_train Patient-level estimates (See PLE_models)
#' @param minsize Minimum number of observations in a tree node.
#' Default = floor( dim(train)[1]*0.05  )
#' @param maxdepth Maximum depth of any node in the tree (default=4)
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
#' @seealso \code{\link{PRISM}}, \code{\link{mob}}
#'
## MOB: Weibull ##
submod_weibull = function(Y, A, X, Xtest, mu_train, minsize = floor( dim(X)[1]*0.05  ),
                          maxdepth = 4, ...){

  #### Weibull Mob Functions ###
  wbreg <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...){
    survreg(y ~ 0 + x, weights = weights, dist = "weibull", ...)
  }
  logLik.survreg <- function(object, ...){
    structure(object$loglik[2], df = sum(object$df), class = "logLik")
  }

  ## Fit Model ##
  mod <- mob(Y ~ A | ., data = X,
             fit = wbreg, control = mob_control(parm=1:3, minsize=minsize,
                                                maxdepth=maxdepth))

  res = list(mod=mod)
  class(res) = "submod_weibull"
  ## Return Results ##
  return(  res  )
}

#' Predict submod: Model-based partitioning (Weibull)
#'
#' Predict subgroups and obtain subgroup-specific point-estimates (in pprogress).
#'
#' @param object Trained MOB (Weibull) model.
#' @param newdata Data-set to make predictions at.
#' @param ... Any additional parameters, not currently passed through.
#'
#' @import partykit
#'
#' @return Identified subgroups with subgroup-specific predictions.
#' \itemize{
#'   \item Subgrps - Identified subgroups
#'   \item pred - Predictions, by subgroup.
#'}
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
#' out = predict(res_weibull, newdata=X)
#' plot(res_weibull$mod)
#' }
#'
#'
#' @method predict submod_weibull
#' @export
#'
predict.submod_weibull = function(object, newdata, ...){

  # Extract mod #
  mod = object$mod
  ##  Predict Subgroups ##
  Subgrps = as.numeric( predict(mod, type="node", newdata = newdata) )
  ## What should we predict? TBD in progress ##
  pred =  NA

  ## Return Results ##
  return(  list(Subgrps=Subgrps, pred=pred) )
}

