#' Subgroup Identification: Model-based partitioning (glmtree)
#'
#' Uses the glmtree (model-based partitioning, glm; through partykit R package) algorithm 
#' to identify subgroups (Zeileis, Hothorn, Hornik 2008). Usable for continuous and 
#' binary outcomes.
#'
#' @param Y The outcome variable. Must be numeric
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate space.
#' @param Xtest Test set
#' @param mu_train Patient-level estimates (See PLE_models)
#' @param alpha Significance level for variable selection (default=0.05)
#' @param glm.fam Family for GLM; default=binomial
#' @param link Link function for GLM; default="identity"
#' @param minsize Minimum number of observations in a tree node.
#' Default = floor( dim(train)[1]*0.05  )
#' @param maxdepth Maximum depth of any node in the tree (default=4)
#' @param parm Model parameters included in parameter instability tests 
#' (default=NULL, all parameters)
#' @param ... Any additional parameters, not currently passed through.
#'
#' @import partykit
#' @importFrom stats binomial
#'
#' @return Trained lmtree model.
#'  \itemize{
#'   \item mod - lmtree model object
#' }
#'
#' @references Zeileis A, Hothorn T, Hornik K (2008). Model-Based Recursive Partitioning. 
#' Journal of Computational and Graphical Statistics, 17(2), 492–514.
#' @export
#' @examples
#' 
#' \donttest{
#' library(StratifiedMedicine)
#'
#' ## Binomial ##
#' dat_bin = generate_subgrp_data(family="binomial")
#' Y = dat_bin$Y
#' X = dat_bin$X
#' A = dat_bin$A
#'
#' 
#' res_glmtree1 = submod_glmtree(Y, A, X, Xtest=X)
#' res_glmtree2 = submod_glmtree(Y, A, X, Xtest=X, link="logit")
#' plot(res_glmtree1$mod)
#' plot(res_glmtree2$mod)
#' }
#'
#'
#### glmtree (MOB) ###
submod_glmtree = function(Y, A, X, Xtest, mu_train, 
                         glm.fam = binomial, link="identity", alpha=0.05,
                         minsize = floor( dim(X)[1]*0.10  ),
                         maxdepth = 4, parm=NULL, ...){

  ## Fit Model ##
  mod <- glmtree(Y~A | ., data = X, family= glm.fam(link=link), 
                 alpha=alpha, maxdepth = maxdepth, 
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
  class(res) <- "submod_glmtree"
  ## Return Results ##
  return(  res )
}
