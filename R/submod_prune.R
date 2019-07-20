#' Subgroup Identification: Pruning
#'
#' Function to prune a trained subgroup model. Can be used within PRISM. This is currently
#' experimental and includes two approaches: (1) Pruning initial subgroups based on optimal
#' treatment regime loss function, and (2) "Double" Subgroup model, which uses the initial
#' discovered subgroups as the input covariate in the subgroup model.
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate space.
#' @param Xtest Test set
#' @param type Type of pruning. Options include "OTR" and "submod2X".
#' @param submod0 Initial subgroup model fit.
#' @param submod Subgroup identification (submod) function. Maps the observed data and/or
#' PLEs to subgroups.
#' @param hyper Hyper-parameters (must be list). Default is NULL.
#' @param mu_train Patient-level estimates (See PLE_models). Default=NULL
#' @param family Outcome type ("gaussian", "binomial", "survival"). Default="gaussian".
#' @param ... Any additional parameters, not currently passed through.
#'
#'
#' @return Trained subgroup model and subgroup predictions/estimates for train/test sets.
#'
#'  \itemize{
#'   \item mod - trained subgroup model
#'   \item Subgrps.train - Identified subgroups (training set)
#'   \item Subgrps.test - Identified subgroups (test set)
#'   \item pred.train - Predictions (training set)
#'   \item pred.test - Predictions (test set)
#'   \item Rules - Definitions for subgroups, if provided in submod output.
#' }

submod_prune = function(Y, A, X, Xtest, type, submod0, submod=NULL,
                        hyper = NULL, mu_train=NULL, family="gaussian", ...){

  if (type=="OTR"){

    if (is.null(mu_train)){
      stop("OTR submod pruning: PLE estimates (mu_train) are missing.")
    }
    # Fit OTR model with subgroups (as factor) to prune #
    submod0 = step3
    Subgrps.train = data.frame(Subgrps=as.factor(step3$Subgrps.train))
    Subgrps.test = data.frame(Subgrps=factor(step3$Subgrps.test,
                                             levels = levels(Subgrps.train)) )

    mod.prune = submod_train(Y=Y, A=A, X=Subgrps.train, Xtest=Subgrps.test,
                             mu_train=mu_train, family = family, submod="submod_otr",
                             hyper = hyper)
  }
  if (type=="submod2X"){
    Rules = step3$Rules
    temp = aggregate( step3$pred.train~step3$Subgrps.train, FUN = "mean")
    # temp = aggregate(mu_train$PLE~step3$Subgrps.train, FUN="mean")
    colnames(temp) = c("Subgrps.temp", "Subgrps.est")
    temp$Rank = rank(temp$Subgrps.est)
    temp = temp[,-2]
    levels_rank = temp[order(temp$Rank),]
    levels_rank = levels_rank$Rank
    # Merge into new design matrices #
    X_new = data.frame(Subgrps.temp = step3$Subgrps.train)
    X_new_test = data.frame(Subgrps.temp = step3$Subgrps.test)
    X_new = left_join(X_new, temp, by="Subgrps.temp")
    X_new_test = left_join(X_new_test, temp, by="Subgrps.temp")
    X_new$Rank = factor(X_new$Rank, ordered = TRUE, levels = levels_rank )
    X_new_test$Rank = factor(X_new_test$Rank, ordered = TRUE, levels = levels_rank )
    X.star.temp = data.frame(Rank=X_new[,"Rank"])
    Xtest.star.temp = data.frame(Rank=X_new_test[,"Rank"])
    ## Re-Fit Subgroup Model ##
    step3X = submod_train(Y=Y, A=A, X=X.star.temp, Xtest=Xtest.star.temp,
                          mu_train=mu_train,family = family, submod=submod,
                          hyper = hyper)
    ## Obtain the "Merged" Rules ##
    if (is.null(Rules)){ #if no definitions, use the numeric predictions #
      temp.rules = data.frame(Subgrps.temp=temp$Subgrps.temp,
                              Rules = temp$Subgrps.temp)
    }
    if (!is.null(Rules)){
      temp.rules = Rules
      colnames(temp.rules) = c("Subgrps.temp", "Rules")
    }
    rules.dat = data.frame(Subgrps.temp = step3$Subgrps.train,
                           Subgrps.new = step3X$Subgrps.train)
    rules.dat = left_join(rules.dat, temp.rules, by="Subgrps.temp")
    rules.dat = unique(rules.dat[,c("Subgrps.new", "Rules")])
    rules.dat <- aggregate(Rules ~ Subgrps.new, data = rules.dat, paste,
                           collapse = " , ")
    colnames(rules.dat) = c("Subgrps", "Rules")
    step3X$Rules = rules.dat
    step3 = step3X
  }

}












