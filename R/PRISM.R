#' PRISM: Patient Response Identifier for Stratified Medicine
#'
#' PRISM algorithm. Given a data-set of (Y, A, X) (Outcome, treatment, covariates),
#' the \code{PRISM} identifies potential subgroup along with point and variability metrics.
#' This four step procedure (filter, ple, submod, param) is flexible and accepts user-inputs
#' at each step.
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (ex: a=1,...,A or a="control","new")
#' @param X Covariate space. Variables types (ex: numeric, factor, ordinal) should be set
#' to align with subgroup model (submod argument). For example, for lmtree, binary variables
#' coded as numeric (ex: 0, 1) are treated differently than the corresponding factor version
#' (ex: "A", "B"). Filter and PLE models provided in the StratifiedMedicine package can
#' accomodate all variable types.
#' @param Xtest Test set. Default is NULL which uses X (training set). Variable types should
#' match X.
#' @param family Outcome type. Options include "gaussion" (default), "binomial", and "survival".
#' @param filter Maps (Y,A,X) => (Y,A,X.star) where X.star has potentially less
#' covariates than X. Default is "Filter_ENET", NULL uses no filter.
#' @param ple PLE (Patient-Level Estimate) function. Maps the observed data to PLEs.
#' (Y,A,X) ==> PLE(X). Default for "gaussian"/"binomial" is "ple_ranger"
#' (treatment-specific random forest models). The default for "survival" is
#' "ple_glmnet" (elastic net (glmnet) cox regression).
#' @param submod Subgroup identification model function. Maps the observed data and/or PLEs
#' to subgroups. Default of "gaussian"/"binomial" is "submod_lmtree" (MOB with OLS loss).
#' Default for "survival" is "submod_weibull" (MOB with weibull loss)
#' @param param Parameter estimation and inference function. Based on the discovered subgroups,
#' perform inference through the input function (by name). Default for "gaussian"/"binomial" is
#' "param_PLE", default for "survival" is "param_cox".
#' @param alpha_ovrl Two-sided alpha level for overall population. Default=0.05
#' @param alpha_s Two-sided alpha level at subgroup level. Default=0.05
#' @param filter.hyper Hyper-parameters for the Filter function (must be list). Default is NULL.
#' @param ple.hyper Hyper-parameters for the PLE function (must be list). Default is NULL.
#' @param submod.hyper Hyper-parameters for the SubMod function (must be list). Default is NULL.
#' @param submod.prune Pruning option for subgroup model. Current options include "OTR"
#' (feed initial discovered subgroups into submod_otr) and "2X" which feeds the initial
#' subgroups into submod (based on 5-STAR, Marceau-West and Mehrotra 2019 in progress).
#' Default=NULL (no pruning).
#' @param param.hyper Hyper-parameters for the Param function (must be list). Default is NULL.
#' @param prefilter_resamp Option to filter the covariate space (based on filter model) prior
#' to resampling. Default=FALSE.
#' @param resample Resampling method for resample-based estimates and variability metrics.
#' Options include "Boostrap" and "Permutation." Default=NULL (No resampling).
#' @param stratify Stratified resampling (Default=TRUE)
#' @param R Number of resamples (default=100)
#' @param filter.resamp Filter function during resampling, default=NULL (use original Filter)
#' @param ple.resamp PLE function during resampling, default=NULL (use original PLE)
#' @param submod.resamp SubMod function for resampling, default=NULL (use original SubMod)
#' @param verbose Detail progress of PRISM? Default=TRUE
#' @param verbose.resamp Output iterations during resampling? Default=FALSE
#'
#' @return Trained PRISM object. Includes filter, ple, submod, and param outputs.
#'  \itemize{
#'   \item filter.mod - Filter model
#'   \item filter.vars - Variables remaining after filtering
#'   \item ple.fit - Fitted ple model (model fit, other fit outputs)
#'   \item mu_train - Patient-level estimates (train)
#'   \item mu_test - Patient-level estimates (test)
#'   \item submod.fit - Fitted submod model (model fit, other fit outputs)
#'   \item out.train - Training data-set with identified subgroups
#'   \item out.test - Test data-set with identified subgroups
#'   \item Rules - Subgroup rules / definitions
#'   \item param.dat - Parameter estimates and variablity metrics (depends on param)
#'   \item resamp.dist - Resampling distributions (NULL if no resampling is done)
#' }
#' @export
#' @importFrom stats aggregate coef lm model.matrix p.adjust pnorm confint
#' @importFrom stats predict pt qnorm qt quantile sd weighted.mean
#' @import dplyr
#' @import ggplot2
#' @import survival
#'
#' @examples
#' ## Load library ##
#' library(StratifiedMedicine)
#'
#' ##### Examples: Continuous Outcome ###########
#'
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#'
#' # Run Default: filter_glmnet, ple_ranger, submod_lmtree, param_ple #
#' res0 = PRISM(Y=Y, A=A, X=X)
#' res0$filter.vars # variables that pass the filter
#' plot(res0, type="PLE:density") # distribution of PLEs
#' plot(res0, type="PLE:waterfall") # PLE waterfall plot
#' plot(res0$submod.fit$mod) # Plot of subgroup model
#' res0$param.dat # overall/subgroup specific parameter estimates/inference
#' plot(res0) # Forest plot: overall/subgroup specific parameter estimates (CIs)
#'
#' # Without filtering #
#' res1 = PRISM(Y=Y, A=A, X=X, filter=NULL)
#' plot(res1$submod.fit$mod)
#' plot(res1)
#'
#' ## With bootstrap (No filtering) ##
#' \donttest{
#'   res_boot = PRISM(Y=Y, A=A, X=X, resample = "Bootstrap", R=50, verbose.resamp = TRUE)
#'   # Plot of distributions and P(est>0) #
#'   plot(res_boot, type="resample")+geom_vline(xintercept = 0)
#'   aggregate(I(est>0)~Subgrps, data=res_boot$resamp.dist, FUN="mean")
#' }
#'
#' # Survival Data ##
#' \donttest{
#'   library(survival)
#'   require(TH.data); require(coin)
#'   data("GBSG2", package = "TH.data")
#'   surv.dat = GBSG2
#'   # Design Matrices ###
#'   Y = with(surv.dat, Surv(time, cens))
#'   X = surv.dat[,!(colnames(surv.dat) %in% c("time", "cens")) ]
#'   set.seed(513)
#'   A = rbinom( n = dim(X)[1], size=1, prob=0.5  )
#'
#'   # Default: PRISM: glmnet ==> MOB (Weibull) ==> Cox; bootstrapping posterior prob/inference #
#'   res_weibull1 = PRISM(Y=Y, A=A, X=X, ple=NULL, resample="Bootstrap", R=100,
#'                        verbose.resamp = TRUE)
#'   plot(res_weibull1$submod.fit$mod)
#'   plot(res_weibull1)
#'   plot(res_weibull1, type="resample")+geom_vline(xintercept = 1)
#'   aggregate(I(est<1)~Subgrps, data=res_weibull1$resamp.dist, FUN="mean")
#'
#'   # PRISM: ENET ==> CTREE ==> Cox; bootstrapping for posterior prob/inference #
#'   res_ctree1 = PRISM(Y=Y, A=A, X=X, ple=NULL, submod = "submod_ctree",
#'                      resample="Bootstrap", R=100, verbose.resamp = TRUE)
#'   plot(res_ctree1$submod.fit$submod.fit$mod)
#'   plot(res_ctree1)
#'   plot(res_ctree1, type="resample")+geom_vline(xintercept = 1)
#'   aggregate(I(est<1)~Subgrps, data=res_ctree1$resamp.dist, FUN="mean")
#' }
#'
#' @references Jemielita and Mehrotra (2019 in progress)

##### PRISM: Patient Responder Identifiers for Stratified Medicine ########
PRISM = function(Y, A, X, Xtest=NULL, family="gaussian",
                 filter="filter_glmnet", ple=NULL, submod=NULL, param=NULL,
                 alpha_ovrl=0.05, alpha_s = 0.05,
                 filter.hyper=NULL, ple.hyper=NULL, submod.hyper = NULL,
                 submod.prune=NULL, param.hyper = NULL, prefilter_resamp=FALSE,
                 resample = NULL, stratify=TRUE,
                 R = 100, filter.resamp = NULL, ple.resamp = NULL,
                 submod.resamp = NULL, verbose=TRUE,
                 verbose.resamp = FALSE){

  ## "Test" Set ##
  if (is.null(Xtest)){ Xtest = X   }

  ## Is the Outcome Survival? ##
  if (is.Surv(Y) & family!="survival"){ family="survival"  }

  # Missing data? #
  if ( sum(is.na(Y))>0 | sum(is.na(A))>0 | sum(is.na(X))>0 ){
    message("Missing Data (in outcome, treatment, or covariates)")
  }

  ### Defaults: By Family (gaussian, binomial (Risk Difference), survival ) ##
  if (family=="gaussian" | family=="binomial"){
    if (is.null(ple) ){ ple = "ple_ranger" }
    if (is.null(submod) ){ submod = "submod_lmtree" }
    if (is.null(param) ){ param = "param_ple" }
  }
  if (family=="survival"){
    if (is.null(ple) ){ ple = "ple_glmnet" }
    if (is.null(submod) ){ submod = "submod_weibull" }
    if (is.null(param) ){ param = "param_cox" }
  }

  ### Main Model Pipeline: Feeds into CV/Bootstrap/resamples procedures ##
  main = function(Y, A, X, Xtest, ple, filter, submod, verbose){
    #### Step 1: Variable Filtering #####
    if ( !is.null(filter) ){
      if (verbose) message( paste("Filtering:", filter, sep=" ") )
      step1 = do.call( filter, append(list(Y=Y, A=A, X=X, family=family), filter.hyper)  )
      filter.mod = step1$mod
      filter.vars = step1$filter.vars
    }
    if ( is.null(filter) ){
      filter.mod = NULL; filter.vars = NULL;
    }
    # Drop variables depending on filter #
    if ( is.null(filter) ){ X.star = X; Xtest.star = Xtest }
    if ( !is.null(filter) ){
      if (length(filter.vars)==0){ X.star = X; Xtest.star = Xtest  }
      if (length(filter.vars) >0){
        X.star = X[, colnames(X) %in% c(filter.vars, "A", "Y"), drop=FALSE]
        Xtest.star = Xtest[, colnames(Xtest) %in% c(filter.vars, "A", "Y"), drop=FALSE] }
    }
    #### Step 2: PLE Estimation #####
    if ( !is.null(ple) ){
      if (verbose) message( paste("PLE:", ple, sep=" " ) )
      step2 = ple_train(Y=Y,A=A,X=X.star,Xtest=Xtest.star,family=family, ple=ple,
                        hyper = ple.hyper)
      ple.fit = step2$fit
      mu_train = step2$mu_train
      mu_test = step2$mu_test
    }
    if ( is.null(ple) ){
      ple.mods = NULL
      mu_train = NULL
      mu_test = NULL
    }
    ### Step 3: Subgroup Identification ###
    if ( !is.null(submod) ){
      if (verbose) message( paste("Subgroup Identification:",
                                submod, sep=" "))
      step3 = submod_train(Y=Y, A=A, X=X.star, Xtest=Xtest.star, mu_train=mu_train,
                           family = family, submod=submod, hyper = submod.hyper)

      # Option for subgroup pruning: experimental ##
      if (!is.null(submod.prune)){

        step3X = submod_prune(Y=Y, A=A, X=X.star, Xtest=Xtest.star, submod0 = step3,
                              submod=submod, mu_train=mu_train,family = family,
                              hyper = submod.hyper)
        step3 = step3X
      }
      submod.fit = step3$fit
      Rules=step3$Rules;
      Subgrps.train = step3$Subgrps.train; Subgrps.test = step3$Subgrps.test
    }
    if ( is.null(submod) ){
      submod.fit = NULL; Rules=NULL;
      Subgrps.train = NULL; Subgrps.test = NULL;
    }
    ### Step 4: Parameter Estimation and Inference ###
    if (verbose){ message(paste("Parameter Estimation:", param, sep=" ")) }
    param.dat = do.call( param, list(Y=Y, A=A, X=X.star, mu_hat = mu_train,
                                     Subgrps=Subgrps.train,
                                     alpha_ovrl=alpha_ovrl,
                                     alpha_s=alpha_s)  )
    ### Return Outputs ###
    return( list( mu_train = mu_train, mu_test = mu_test, filter.mod = filter.mod,
                  filter.vars = filter.vars, ple.fit = ple.fit, submod.fit = submod.fit,
                  Subgrps.train = Subgrps.train, Subgrps.test=Subgrps.test,
                  Rules = Rules, param.dat=param.dat) )
  }
  ## Run on Observed Data ##
  if (verbose){ message( "Observed Data" )   }
  res0 = main(Y=Y, A=A, X=X, Xtest=Xtest, ple=ple, filter=filter,
              submod = submod, verbose=verbose)
  Subgrps = res0$Subgrps.train
  mu_train = res0$mu_train
  param.dat = res0$param.dat

  #################################################################
  #### Resampling (Bootstrapping or Permutation)    ###############
  #################################################################

  resamp_param = NULL ## Set to null (needed if no resampling)

  if (!is.null(resample)){

    ## Resampling: Auto-checks ##
    if (is.null(filter.resamp)) {  filter.resamp = filter  }
    if (is.null(ple.resamp)) {  ple.resamp = ple  }
    if (is.null(submod.resamp)) {  submod.resamp = submod  }

    obs.data = data.frame(id = 1:nrow(X), Y=Y,A=A, X, Subgrps)
    if (length(res0$filter.vars)>0 & prefilter_resamp == TRUE){
      obs.data = obs.data[, colnames(obs.data) %in%
                            c("id", "Y", "A", res0$filter.vars, "Subgrps")]
    }
    # Xtest.star = Xtest[,colnames(Xtest) %in% res0$filter.vars]
    if (verbose){ message( paste(resample, R, "resamples"))   }

    subject.counter = data.frame(id = obs.data$id )

    ### Resampling Wrapper ###
    fetty_wop = function(R, stratify, obs.data, ple, filter, submod,
                         calibrate, verbose){

      if (verbose) message( paste(resample, "Sample", R) )
      ### Permutation resampling (shuffle treatment assignment) ###
      if (resample=="Permutation"){
        seedz = 5341+R
        set.seed(seedz)
        A_resamp = sample(obs.data$A, replace=FALSE)
        resamp.data = obs.data
        resamp.data$A = A_resamp
      }
      if (resample=="Bootstrap"){
        ### Re-sample the data ##
        if (stratify){
          resamp.data = NULL
          for (s in unique(Subgrps)){
            hold.s = obs.data[obs.data$Subgrps==s,]
            seedz = 5341+R
            set.seed(seedz)
            indices.R = sample(nrow(hold.s), size = dim(hold.s)[1], replace=TRUE)
            resamp.data = rbind(resamp.data, hold.s[indices.R,])
          }
        }
        if (!stratify){
          seedz = 5341+R
          set.seed(seedz)
          indices.R = sample(nrow(obs.data), size = dim(obs.data)[1], replace=TRUE)
          resamp.data = obs.data[indices.R,]
        }

      }
      #### PRISM on resampled data: Test set = Observed Data ####
      Y.R = resamp.data$Y
      A.R = resamp.data$A
      X.R = resamp.data[!(colnames(resamp.data) %in% c("Subgrps", "id", "Y", "A"))]
      if (family=="survival"){ Y.R = Surv(Y.R[,1], Y.R[,2])  }
      res.R = main(Y=Y.R, A=A.R, X=X.R, Xtest=Xtest, ple=ple, filter=filter,
                   submod = submod, verbose=FALSE)
      Subgrps.R = res.R$Subgrps.train
      param.R = res.R$param.dat
      if (calibrate){
        ### Calibrate alpha (need to implement) ###
        ## NEED WORKS!!! ##
      }
      ## Original Subgroups: Based on bootstrap parameter estimates ##
      hold = param.R %>% filter(param.R$Subgrps>0)
      hold = hold[, colnames(hold) %in% c("Subgrps", "est") ]
      est.resamp = left_join( data.frame(Subgrps=Subgrps.R),
                              hold, by="Subgrps")
      est.resamp = data.frame(id=resamp.data$id, Subgrps = resamp.data$Subgrps, est.resamp)
      colnames(est.resamp) = c("id", "Subgrps", "Subgrps.R", "est")
      param.resamp = aggregate(est ~ Subgrps, data=est.resamp, FUN="mean")
      param.resamp = rbind( param.R[param.R$Subgrps==0,c("Subgrps", "est")],
                            param.resamp )
      param.resamp = data.frame(R=R, param.resamp)
      ## Counter for each subject (how many times did they appear in the bootstrap sample)##
      cnt.table = table(resamp.data$id)
      counter = suppressWarnings( left_join(subject.counter,
                                            data.frame(id=as.numeric(names(cnt.table)),
                                                       count=as.numeric(cnt.table)),
                                            by = "id") )
      counter = data.frame(R=R, counter)
      ## Output counter and estimated statistics ##
      return( list(param.resamp=param.resamp, counter=counter$count) )
    }

    ### Run Resampling ##
    resamp.obj = lapply(1:R, fetty_wop, stratify=stratify, obs.data=obs.data,
                        ple=ple.resamp, filter=filter.resamp,
                        submod=submod.resamp, calibrate=FALSE,
                        verbose = verbose.resamp)
    ## Extract Resampling parameter estimates and subject-counters ##
    hold = do.call(rbind, resamp.obj)
    resamp_param = do.call(rbind, hold[,1])
    resamp_param = resamp_param %>% arrange(Subgrps)
    resamp_counter = do.call(cbind, hold[,2])

    ### Calculate Resample Metrics
    ## Smoothed estimate (average across resamples), SE, pval
    ## Bootstrap specific: Covariances, Acceleraation, PCT / BCa CIs ###
    ## Standard errors, covariances, acceleration, CIs
    resamp_metrics = function(param.dat){
      final_ests = param.dat
      final_ests$est_resamp = NA
      final_ests$SE_resamp = NA
      if (resample=="Permutation"){
        final_ests$pval_perm = NA
      }
      if (resample=="Bootstrap"){
        final_ests$SE_bootS = NA
        final_ests$bias = NA
        final_ests$accel = NA
        final_ests$LCL.pct = NA
        final_ests$UCL.pct = NA
        final_ests$LCL.BCa = NA
        final_ests$UCL.BCa = NA
      }
      for (sub in unique(final_ests$Subgrps)){
        ### Smoothed Bootstrap estimate, boot SD ###
        est0 = final_ests$est[final_ests$Subgrps==sub]
        est.vec = resamp_param$est[resamp_param$Subgrps==sub]
        final_ests$est_resamp[final_ests$Subgrps==sub] = mean(est.vec)
        final_ests$SE_resamp[final_ests$Subgrps==sub] = sd( est.vec )
        ## Permutation p-value ##
        if (resample=="Permutation"){
          final_ests$pval_perm[final_ests$Subgrps==sub] =
            (sum(abs(est.vec)>abs(est0)) + 1 ) / (length(est.vec)+1)
        }
        ## Bootstrap Covariance/acceleration/bias/smoothed SE ##
        if (resample=="Bootstrap"){
          if (sub %in% c(-1, 0)){
            sub_subjs = as.matrix( resamp_counter  )
            alpha = alpha_ovrl
          }
          if (sub>0){
            sub_subjs = as.matrix( resamp_counter[obs.data$Subgrps==sub,]  )
            alpha = alpha_s
          }
          sub_subjs = apply(sub_subjs, 2, function(x) ifelse(is.na(x)==TRUE,0,x) )
          count_mean = as.numeric( rowMeans(sub_subjs, na.rm=TRUE) )
          cov.b = suppressWarnings(
            rowMeans( (sub_subjs - count_mean) * (est.vec - mean(est.vec) ) ) )
          bias.b =  suppressWarnings(
            rowMeans( (sub_subjs - 1)^2 * (est.vec - mean(est.vec) ) ) )
          bca.accel = (1/6)*sum(cov.b^3) / (sum(cov.b^2))^(2/3)
          bca.bias = mean(bias.b)
          final_ests$accel[final_ests$Subgrps==sub] = bca.accel
          final_ests$bias[final_ests$Subgrps==sub] = bca.bias
          final_ests$SE_bootS[final_ests$Subgrps==sub] = sqrt( sum(cov.b^2) )
          ### Confidence Intervals (Pct and BCa) ###
          BCa.QL = pnorm(bca.bias +(bca.bias+qnorm(alpha/2))/
                           (1-bca.accel*(bca.bias+qnorm(alpha/2))) )
          BCa.QU = pnorm(bca.bias +(bca.bias+qnorm(1-alpha/2))/
                           (1-bca.accel*(bca.bias+qnorm(1-alpha/2))) )
          quants = as.numeric(
            quantile(est.vec, probs=c(alpha/2, (1-alpha/2), BCa.QL, BCa.QU)) )
          final_ests$LCL.pct[final_ests$Subgrps==sub] = quants[1]
          final_ests$UCL.pct[final_ests$Subgrps==sub] = quants[2]
          final_ests$LCL.BCa[final_ests$Subgrps==sub] = quants[3]
          final_ests$UCL.BCa[final_ests$Subgrps==sub] = quants[4]
        }
      }
      return(final_ests)
    }
    param.dat = resamp_metrics(param.dat=param.dat)

  }

  ### Return Results ##
  res = list( filter.mod = res0$filter.mod, filter.vars = res0$filter.vars,
              ple.fit = res0$ple.fit, mu_train=res0$mu_train, mu_test=res0$mu_test,
              submod.fit = res0$submod.fit,
              out.train = data.frame(Y, A, X, Subgrps=res0$Subgrps.train),
              out.test = data.frame(Xtest, Subgrps=res0$Subgrps.test),
              Rules=res0$Rules,
              param.dat = param.dat, resamp.dist = resamp_param,
              family = family,
              filter = filter, ple = ple, submod=submod, param=param)
  class(res) <- c("PRISM")
  return(res)
}

#' PRISM: Patient Response Identifier for Stratified Medicine (Predictions)
#'
#' Predictions for PRISM algorithm. Given the training set (Y,A,X) or new test set (Xtest),
#' output ple predictions and identified subgroups with correspond parameter estimates.
#'
#' @param object Trained PRISM model.
#' @param newdata Data-set to make predictions at (Default=NULL, predictions correspond
#' to training data).
#' @param type Type of prediction. Default is "all" (ple, submod, and param predictions).
#' Other options include "ple" (ple predictions), "submod" (submod predictions with
#' associated parameter estimates).
#' @param ... Any additional parameters, not currently passed through.
#'
#' @return Data-frame with predictions (ple, submod, or both).
#' @export
#'
#' @examples
#' ## Load library ##
#' library(StratifiedMedicine)
#'
#' ##### Examples: Continuous Outcome ###########
#'
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#'
#' # Run Default: filter_glmnet, ple_ranger, submod_lmtree, param_ple #
#' res0 = PRISM(Y=Y, A=A, X=X)
#' summary( predict(res0) ) # all #
#' summary( predict(res0, type="ple") )
#' summary( predict(res0, type="submod") )
#'
#'
#' @method predict PRISM
#' @export
#'
predict.PRISM = function(object, newdata=NULL, type="all", ...){

  if (type=="all"){
    mu_hat = predict(object$ple.fit, newdata=newdata)
    Subgrps = predict(object$submod.fit, newdata=newdata)
    res = data.frame(mu_hat, Subgrps=Subgrps$Subgrps)
    params = object$param.dat
    res = left_join(res, params[,c("Subgrps", "est")], by="Subgrps")
    res = data.frame(res)
  }
  if (type=="ple"){
    mu_hat = predict(object$ple.fit, newdata=newdata)
    res = data.frame(mu_hat)
  }
  if (type=="submod"){
    Subgrps = predict(object$submod.fit, newdata=newdata)
    res = data.frame(Subgrps=Subgrps$Subgrps)
    params = object$param.dat
    res = left_join(res, params[,c("Subgrps", "est")], by="Subgrps")
  }
  return(res)
}

