#' Filter: Elastic Net (glmnet)
#'
#' Filter variables through elastic net (Zou and Hastie 2005). Default is to regress
#' Y~X (search for prognostic variables). Variables with estimated coefficients of zero
#' (depends on lambda choice; default is lambda.min) are filtered. Usable for continuous,
#' binary, and survival outcomes.
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate space.
#' @param lambda Lambda for elastic net model (default="lambda.min"). Other options include
#' "lambda.1se" and fixed values
#' @param family Outcome type ("gaussian", "binomial", "survival"), default is "gaussian"
#' @param interaction Regress Y~X+A+A*X (interaction between covariates and treatment)?
#' Default is FALSE. If TRUE, variables with zero coefficients (both X and X*A terms)
#' are filtered.
#' @param ... Any additional parameters, not currently passed through.
#'
#' @return Filter model and variables that remain after filtering.
#'  \itemize{
#'   \item mod - Filtering model
#'   \item filter.vars - Variables that remain after filtering (could be all)
#' }
#' @export
#' @references Friedman, J., Hastie, T. and Tibshirani, R. (2008) Regularization Paths for
#'  Generalized Linear Models via Coordinate Descent,
#'  \url{https://web.stanford.edu/~hastie/Papers/glmnet.pdf} Journal of Statistical 
#'  Software, Vol. 33(1), 1-22 Feb 2010 Vol. 33(1), 1-22 Feb 2010.

#' @examples
#'
#' library(StratifiedMedicine)
#'
#' ## Continuous ##
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#'
#' # Default: Regress Y~X (search for prognostic factors) #
#' mod1 = filter_glmnet(Y, A, X)
#' mod2 = filter_glmnet(Y, A, X, lambda = "lambda.min") # same as default
#' mod3 = filter_glmnet(Y, A, X, lambda = "lambda.1se")
#' mod1$filter.vars
#' mod2$filter.vars
#' mod3$filter.vars
#'
#' # Interaction=TRUE; Regress Y~X+A+X*A (search for prognostic and/or predictive) #
#' mod4 = filter_glmnet(Y, A, X, interaction=TRUE)
#' mod4$filter.vars

##### Elastic net (glmnet): Y~X ######
filter_glmnet = function(Y, A, X, lambda="lambda.min", family="gaussian",
                         interaction=FALSE, ...){

  ## Model Matrix #
  fact.vars <- sapply(X, is.factor)
  X.mat <- X
  colnames(X.mat)[fact.vars] <- paste(colnames(X.mat)[fact.vars], "_lvl_", sep="")
  X.mat <- model.matrix(~., data = X.mat)
  X.mat <- X.mat[, colnames(X.mat) != "(Intercept)"]
  W <- X.mat
  intercept <- TRUE

  if (interaction){
    A.mat <- model.matrix(~., data=data.frame(A))[,-1]
    X_inter = X.mat*A.mat
    colnames(X_inter) = paste(colnames(X.mat), "_trtA", sep="")
    W = cbind(X.mat, A, X_inter)
  }
  # Center and Scale #
  n <- dim(W)[1]
  W_centered <- apply(W, 2, function(x) x - mean(x))
  Ws <- apply(W_centered, 2, function(x) x / sqrt(sum(x^2) / n))

  # Fit Elastic Net #
  if (family=="survival") { family = "cox" }
  mod <- suppressWarnings( 
    cv.glmnet(x = Ws, y = Y, nlambda = 100, alpha=0.5, family=family,
                   intercept=intercept) )

  ### Extract filtered variable based on lambda ###
  VI <- coef(mod, s = lambda)[,1]
  VI <- VI[ names(VI) != "(Intercept)" ]
  # Extract variables that pass the filter ##
  filter.vars <- names(VI[VI!=0])
  filter.vars <- unique( gsub("_lvl_.*","",filter.vars) )
  if (interaction) {
    filter.vars <- unique(sub("_trtA","",filter.vars))
    filter.vars <- filter.vars[filter.vars!="A"]
  }
  # Store selected lambda in model #
  mod$sel.lambda <- lambda
  mod$sel.family <- family
  # Return model fit, filter.vars #
  return(list(mod=mod, filter.vars=filter.vars))
}
# VI Plot #
plot_vimp_glmnet <- function(mod) {
  # Extract filtered variable based on lambda #
  VI <- coef(mod, s = mod$sel.lambda)[,1]
  VI = VI[ names(VI) != "(Intercept)" ]
  
  # Importance Plot #
  plot.title <- paste("Elastic Net (", mod$sel.family, ") Importance Plot", sep="")
  VI.dat <- data.frame(covariate = names(VI),
                       est = as.numeric(VI))
  VI.dat <- VI.dat[VI.dat$est!=0,]
  VI.dat$rank <- 1
  vimp.plt <- ggplot2::ggplot(VI.dat, aes(x=reorder(.data$covariate, abs(.data$est)), 
                                 y=.data$est)) + 
    ggplot2::geom_bar(stat="identity")+
    ggplot2::ylab("Beta (standardized)")+
    ggplot2::xlab("Variable")+
    ggplot2::ggtitle(plot.title)+
    ggplot2::coord_flip()
  return(vimp.plt)
}