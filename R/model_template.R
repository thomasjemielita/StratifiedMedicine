#' Model Template: Generates template functions for PRISM algorithm
#'
#' Helper function that generates template functions for the individual components of 
#' PRISM (filter, ple, submod, param). Useful for creating user-specific models, for 
#' example a different ple model for estimating counterfactuals.
#'
#'
#' @param type Type of model template, accepts "filter", "ple", "submod", and "param".
#'
#' @return Function template for chosen model type. 
#' @export
#' 
#' @examples
#' model_template(type="filter")
#' model_template(type="ple")
#' model_template(type="submod")
#' model_template(type="param")

model_template = function(type){
  
  if (type=="filter"){
    filter_template = function(Y, A, X, ...){
      
      # Step 1: Fit Filter Model #
      # for example:
      # mod = cv.glmnet(Y, X)
      
      # Step 2: Extract variables that pass the filter #
      # for example:
      # VI <- coef(mod, s = "lambda.1se" )[,1]
      # VI = VI[ names(VI) != "(Intercept)" ]
      # filter.vars = names(VI[VI!=0])
      
      # # Return model fit and filtered variables #
      # res = list(mod=mod, filter.vars=filter.vars)
      # return( res )
    }
    return(filter_template)
  }
  if (type=="ple"){
    ple_template <- function(Y, A, X, Xtest, ...){
      
      # # Step 1: Fit PLE Model #
      # # for example: Estimate E(Y|A=1,X), E(Y|A=0,X), E(Y|A=1,X)-E(Y|A=0,X)
      # # Split data by treatment #
      # train0 <-  data.frame(Y=Y[A==0], X[A==0,])
      # train1 <-  data.frame(Y=Y[A==1], X[A==1,])
      # # RandomForest (A=1 patients) #
      # mod0 <- ranger(Y ~ ., data = train0)
      # # RandomForest (A=0 patients) #
      # mod1 <- ranger(Y ~ ., data = train1)
      # mod = list(mod0=mod0, mod1=mod1)
      # 
      # # Step 2: Predictions
      # # Option 1: Create a Prediction Function #
      # pred.fun <- function(mod, X){
      #   mu1_hat <- predict( mod$mod1, X )$predictions
      #   mu0_hat <- predict( mod$mod0, X )$predictions
      #   mu_hat <- data.frame(mu1 = mu1_hat, mu0 = mu0_hat, PLE = mu1_hat-mu0_hat)
      #   return(mu_hat)
      # }
      # # Option 2: Directly Output Predictions (here, we still use pred.fun) #
      # mu_train <- pred.fun(mod, X)
      # mu_test <- pred.fun(mod, Xtest)
      
      # # Return model fits and pred.fun (or just mu_train/mu_test) #
      # res <- list(mod=mod, pred.fun=pred.fun, mu_train=mu_train, mu_test=mu_test)
      # return( res )
    }
    return(ple_template)
  }
  if (type=="submod"){
    submod_template <- function(Y, A, X, Xtest, mu_train, ...){
      # # Step 1: Fit subgroup model #
      # # for example: 
      # mod <- ctree(Y ~ ., data = X)
      # 
      # # Step 2: Predictions #
      # # Option 1: Create Prediction Function #
      # pred.fun <- function(mod, X=NULL){
      #   Subgrps <- as.numeric( predict(mod, type="node", newdata = X) )
      #   return( list(Subgrps=Subgrps) )
      # }
      # # Option 2: Output Subgroups for train/test (here we use pred.fun)
      # Subgrps.train = pred.fun(mod, X)
      # Subgrps.test = pred.fun(mod, X)
      # Return fit and pred.fun (or just Subgrps.train/Subgrps.test)
      # res <- list(mod=mod, pred.fun=pred.fun, Subgrps.train=Subgrps.train,
      #             Subgrps.test=Subgrps.test)
      # return(res)
    }
    return(submod_template)
  }
  if (type=="param"){
    param_template <- function(Y, A, X, mu_hat, Subgrps, alpha_ovrl, alpha_s,...){
      # Key Outputs: Subgroup specific and overall parameter estimates
      # Simple Approach: Fit LM within each subgroup and overall #
      # indata = data.frame(Y=Y,A=A, X)
      # ## Subgroup Specific Estimate ##
      # looper = function(s, alpha){
      #   hold = indata[Subgrps %in% s,]
      #   lm.mod = lm(Y ~ A , data=hold )
      #   n.s = dim(hold)[1]
      #   est = summary(lm.mod)$coefficients[2,1]
      #   SE = summary(lm.mod)$coefficients[2,2]
      #   LCL = est - qt(1-alpha/2, df=n.s-1)*SE
      #   UCL = est + qt(1-alpha/2, df=n.s-1)*SE
      #   pval = 2*pt(-abs(est/SE), df=n.s-1)
      #   summ = data.frame( Subgrps = ifelse(n.s==dim(indata)[1], 0, s), N = n.s, 
      #                      estimand = "E(Y|A=1)-E(Y|A=0)", est, SE, LCL, UCL, pval)
      #   return( summ )
      # }
      # # Across Subgroups #
      # S_levels = as.numeric( names(table(Subgrps)) )
      # param.dat = lapply(S_levels, looper, alpha_s)
      # param.dat = do.call(rbind, param.dat)
      # param.dat = data.frame( param.dat )
      # ## Overall ##
      # param.dat0 = looper(S_levels, alpha_ovrl)
      # # Combine and return ##
      # param.dat = rbind(param.dat0, param.dat)
      # return( param.dat )
    }
    return(param_template)
  }
}
