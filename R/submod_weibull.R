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
#' @return MOB (Weibull) model, predictions, identified subgroups, and subgroup rules/definitions.
#'  \itemize{
#'   \item mod - MOB (Weibull) model object
#'   \item Subgrps.train - Identified subgroups (training set)
#'   \item Subgrps.test - Identified subgroups (test set)
#'   \item pred.train - Predictions (training set)
#'   \item pred.test - Predictions (test set)
#'   \item Rules - Subgroups rules/definitions
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
  ##  Predict Subgroups for Train/Test ##
  Subgrps.train = as.numeric( predict(mod, type="node") )
  Subgrps.test = as.numeric( predict(mod, type="node", newdata = Xtest) )
  Rules = list_rules(mod)
  if (length(unique(Subgrps.train))==1){
    Rules = data.frame(Subgrps = unique(Subgrps.train), Rules = "All" )
  }
  if (length(unique(Subgrps.train))>1){
    Rules = data.frame(Subgrps = as.numeric(names(Rules)), Rules = as.character(Rules) )
  }
  # ## Predict Hazard Ratio across subgroups ##
  # for (sub in unique(Subgrps.train)){
  #   # Extract Model #
  #   mod.s = summary(mod)[]
  #
  # }
  ## Predict E(Y|X=x, A=1)-E(Y|X=x,A=0) ##
  pred.train = NA
  pred.test =  NA
  ## Return Results ##
  return(  list(mod=mod, Subgrps.train=Subgrps.train, Subgrps.test=Subgrps.test,
                pred.train=pred.train, pred.test=pred.test, Rules=Rules) )
}
