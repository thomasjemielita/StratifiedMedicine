#' Subgroup Identification: CART (rpart)
#'
#' Uses the CART algorithm (rpart) to identify subgroups. Usable for continuous and binary
#' outcomes. Option to use the observed outcome or PLEs for subgroup identification.
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate space.
#' @param Xtest Test set
#' @param mu_train Patient-level estimates (See PLE_models)
#' @param minbucket Minimum number of observations in a tree node.
#' Default = floor( dim(train)[1]*0.05  )
#' @param maxdepth Maximum depth of any node in the tree (default=4)
#' @param outcome_PLE If TRUE, use PLE as outcome (mu_train must contain PLEs).
#' Else use observed outcome Y
#' @param family Outcome type ("gaussian", "binomial"), default is "gaussian"
#' @param ... Any additional parameters, not currently passed through.
#'
#' @return Trained rpart (CART).
#'  \itemize{
#'   \item mod - rpart model as partykit object
#' }
#'
#' @export
#' @references Breiman L, Friedman JH, Olshen RA, and Stone CJ. (1984) Classification 
#' and Regression Trees. Wadsworth
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
#'
#' require(rpart)
#' res_rpart1 = submod_rpart(Y, A, X, Xtest=X)
#' res_rpart2 = submod_rpart(Y, A, X, Xtest=X, maxdepth=2, minbucket=100)
#' plot(res_rpart1$mod)
#' plot(res_rpart2$mod)
#' }
#'
#'
## CART(rpart) ###
submod_rpart = function(Y, A, X, Xtest, mu_train, minbucket = floor( dim(X)[1]*0.10  ),
                       maxdepth = 4, outcome_PLE=FALSE, family="gaussian", ...){

  if (!requireNamespace("rpart", quietly = TRUE)) {
    stop("Package rpart needed for submod_rpart. Please install.")
  }

  ## Use PLE as outcome? #
  if (outcome_PLE==TRUE){
    Y <- mu_train$PLE
  }
  ## Fit Model ##
  mod <- rpart::rpart(Y ~ ., data = X,
               control = rpart::rpart.control(minbucket=minbucket, maxdepth=maxdepth))
  mod <- as.party(mod)
  # Prediction Function #
  pred.fun <- function(mod, X=NULL, type="subgrp"){
    pred <- NULL
    Subgrps <- as.numeric( predict(mod, type="node", newdata = X) )
    if (type=="all"){
      ## Response Predictions ##
      pred <- data.frame(Subgrps=Subgrps,
                         mu = as.numeric( predict(mod, newdata = X, type="response")) )
    }
    return( list(Subgrps=Subgrps, pred=pred) )
  }
  ## Return Results ##
  res <- list(mod=mod, pred.fun=pred.fun)
  class(res) <- "submod_rpart"
  return(  res )
}
