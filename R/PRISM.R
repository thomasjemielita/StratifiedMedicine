#' PRISM: Patient Responder Identifier for Stratified Medicine
#'
#' PRISM algorithm. Given a data-set of (Y, X, A) (Outcome, covariates, treatment),
#' the \code{PRISM} identifies potential subgroup along with point and variability metrics.
#' This four step procedure (Filter, PLE, SubMod, Param) is flexible and accepts user-inputs
#' at each step.
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate matrix. Must be numeric.
#' @param Xtest Test set. Default is NULL which uses X (training set)
#' @param family Outcome type. Options include "gaussion" (default), "binomial", and "survival".
#' @param Filter Maps (Y,A,X) => (Y,A,X.star) where X.star has potentially less
#' covariates than X. Default is "Filter_ENET", NULL uses no filter.
#' @param PLE PLE (Patient-Level Estimate) function. Maps the observed data to PLEs.
#' (Y,A,X) ==> PLE(X). Default for "gaussian"/"binomial" is "ranger_PLE"
#' (treatment-specific random forest models). The default for "survival" is
#' "cox_PLE" (elastic net cox regression).
#' @param SubMod Subgroup identification model function. Maps the observed data and/or PLEs
#' to subgroups. Default of "gaussian"/"binomial" is "SubMod_lmtree" (MOB with OLS loss).
#' Default for "survival" is "SubMod_weibull" (MOB with weibull loss)
#' @param Param Parameter estimation and inference function. Based on the discovered subgroups,
#' perform inference through the input function (by name). Default for "gaussian"/"binomial" is
#' "Param_PLE." Default for "survival" is "Param_cox."
#' @param alpha_ovrl Two-sided alpha level for overall population. Default=0.05
#' @param alpha_s Two-sided alpha level at subgroup level. Default=0.05
#' @param Filter.hyper Hyper-parameters for the Filter function (must be list). Default is NULL.
#' @param PLE.hyper Hyper-parameters for the PLE function (must be list). Default is NULL.
#' @param SubMod.hyper Hyper-parameters for the SubMod function (must be list). Default is NULL.
#' @param SubMod2X Option to perform "double" subgroup identification; run SubMod once,
#' rank them into ordinal groups, re-run Submod (BETA). Default=FALSE
#' @param Param.hyper Hyper-parameters for the Param function (must be list). Default is NULL.
#' @param preFilter_resamp Option to filter the covariate space (based on Filter model) prior
#' to resampling. Default=FALSE.
#' @param Resample Resampling method for resample-based estimates and variability metrics.
#' Options include "Boostrap" and "Permutation." Default=NULL (No resampling).
#' @param stratify Stratified resampling (Default=TRUE)
#' @param R Number of resamples (default=100)
#' @param Filter.resamp Filter function during resampling, default=NULL (use original Filter)
#' @param PLE.resamp PLE function during resampling, default=NULL (use original PLE)
#' @param SubMod.resamp SubMod function for resampling, default=NULL (use original SubMod)
#' @param verbose Detail progress of PRISM? Default=TRUE
#' @param verbose.resamp Output iterations during resampling? Default=FALSE
#'
#' @return Return Filter, PLE, SubMod, and Param outputs.
#'  \itemize{
#'   \item mu_train - Patient-level estimates (train)
#'   \item mu_test - Patient-level estimates (test)
#'   \item filter.mod - Filter model
#'   \item filter.vars - Variables remaining after filtering
#'   \item Sub.mod - Subgroup model
#'   \item out.train - Training data-set with identified subgroups
#'   \item out.test - Test data-set with identified subgroups
#'   \item Subpred.train - Training predictions (based on SubMod)
#'   \item Subpred.test - Test predictions (based on SubMod)
#'   \item Rules - Subgroup rules / definitions
#'   \item param.dat - Parameter estimates and variablity metrics
#'   \item resamp.dist - Resampling distributions (NULL if no resampling is done)
#'   \item Rules - Subgroups rules/definitions
#' }
#' @export
#' @importFrom stats aggregate coef lm model.matrix p.adjust pnorm
#' @importFrom stats predict pt qnorm qt quantile sd weighted.mean
#' @import dplyr
#' @import ggplot2
#' @import survival
#' @import coin
#'
#' @examples
#' ## Load library ##
#' library(StratifiedMedicine)
#' library(ggplot2)
#'
#' ##### Examples: Continuous Outcome ###########
#'
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#'
#' # Run Default: Filter_ENET, PLE_ranger, SubMod_lmtree, Param_PLE #
#' res0 = PRISM(Y=Y, A=A, X=X)
#' hist(res0$mu_train$PLE) # distribution of PLEs
#' plot(res0$Sub.mod) # Plot of subgroup model
#' res0$param.dat # overall/subgroup specific parameter estimates/inference
#' plot(res0) # Forest plot: overall/subgroup specific parameter estimates (CIs)
#'
#' ## With bootstrap (No filtering) ##
#' \donttest{
#' res_boot = PRISM(Y=Y, A=A, X=X, Resample = "Bootstrap", R=50, verbose.resamp = TRUE)
#' ggplot(res_boot$resamp.dist, aes(est)) + geom_density() +
#' facet_wrap(~Subgrps) + ggtitle("Bootstrap Distribution of Overall/Subgroup Estimates")
#' }
#'
#' # Survival Data ##
#' # Load TH.data (no treatment; generate treatment randomly to simulate null effect) ##
#' \donttest{
#' data("GBSG2", package = "TH.data")
#' surv.dat = GBSG2
#' # Design Matrices ###
#' Y = with(surv.dat, Surv(time, cens))
#' X = surv.dat[,!(colnames(surv.dat) %in% c("time", "cens")) ]
#' A = rbinom( n = dim(X)[1], size=1, prob=0.5  )
#'
#' # Default: PRISM: ENET ==> MOB (Weibull) ==> Cox; bootstrapping posterior prob/inference #
#' res_weibull1 = PRISM(Y=Y, A=A, X=X, PLE=NULL, Resample="Bootstrap", R=100)
#' plot(res_weibull1$Sub.mod)
#' res_weibull1$param.dat
#'
#' # PRISM: ENET ==> CTREE ==> Cox; bootstrapping for posterior prob/inference #
#' res_CTREE1 = PRISM(Y=Y, A=A, X=X, PLE=NULL, SubMod = "SubMod_CTREE",
#'                    Resample="Bootstrap", R=100)
#' plot(res_weibull1$Sub.mod)
#' res_CTREE1$param.dat
#' }
#'
#' @references Jemielita and Mehrotra (2019 in progress)

##### PRISM: Patient Responder Identifiers for Stratified Medicine ########
PRISM = function(Y, A, X, Xtest=NULL, family="gaussian",
                 Filter="Filter_ENET", PLE=NULL, SubMod=NULL, Param=NULL,
                 alpha_ovrl=0.05, alpha_s = 0.05,
                 Filter.hyper=NULL, PLE.hyper=NULL, SubMod.hyper = NULL, SubMod2X=FALSE,
                 Param.hyper = NULL, preFilter_resamp=FALSE, Resample = NULL, stratify=TRUE,
                 R = 100, Filter.resamp = NULL, PLE.resamp = NULL,
                 SubMod.resamp = NULL, verbose=TRUE,
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
    if (is.null(PLE) ){ PLE = "PLE_ranger" }
    if (is.null(SubMod) ){ SubMod = "SubMod_lmtree" }
    if (is.null(Param) ){ Param = "Param_PLE" }
  }
  if (family=="survival"){
    if (is.null(PLE)){ PLE = "PLE_ENET" }
    if (is.null(SubMod) ){ SubMod = "SubMod_weibull" }
    if (is.null(Param) ){ Param = "Param_cox" }
  }

  ### Main Model Pipeline: Feeds into CV/Bootstrap/resamples procedures ##
  main = function(Y, A, X, Xtest, PLE, Filter, SubMod, verbose){
    #### Step 1: Variable Filtering #####
    if ( !is.null(Filter) ){
      if (verbose) message( paste("Filtering:", Filter, sep=" ") )
      step1 = do.call( Filter, append(list(Y=Y, A=A, X=X, family=family), Filter.hyper)  )
      filter.mod = step1$mod
      filter.vars = step1$filter.vars
    }
    if ( is.null(Filter) ){
      filter.mod = NULL; filter.vars = NULL;
    }
    # Drop variables depending on filter #
    if ( is.null(Filter) ){ X.star = X; Xtest.star = Xtest }
    if ( !is.null(Filter) ){
      if (length(filter.vars)==0){ X.star = X; Xtest.star = Xtest  }
      if (length(filter.vars) >0){
        X.star = X[, colnames(X) %in% c(filter.vars, "A", "Y"), drop=FALSE]
        Xtest.star = Xtest[, colnames(Xtest) %in% c(filter.vars, "A", "Y"), drop=FALSE] }
    }
    #### Step 2: PLE Estimation #####
    if ( !is.null(PLE) ){
      if (verbose) message( paste("PLE:", PLE, sep=" " ) )
      step2 = do.call( PLE, append(list(Y=Y, A=A, X=X.star, Xtest=Xtest.star,
                                        family=family), PLE.hyper) )
      mu_train = step2$mu_train
      mu_test = step2$mu_test
    }
    if ( is.null(PLE) ){
      mu_train = NULL; mu_test = NULL;
    }
    ### Step 3: Subgroup Identification ###
    if ( !is.null(SubMod) ){
      if (verbose) message( paste("Subgroup Identification:",
                                SubMod, sep=" "))
      step3 = do.call( SubMod, append( list(Y=Y, A=A, X=X.star, Xtest=Xtest.star,
                                    mu_train = mu_train, family=family), SubMod.hyper) )
      # Option for "double" subgroup ==> Rank order subgroups by plug-in PLE #
      if (SubMod2X){
        Rules = step3$Rules
        temp = aggregate( step3$pred.train~step3$Subgrps.train, FUN = "mean")
        # temp = aggregate(mu_train$PLE~step3$Subgrps.train, FUN="mean")
        colnames(temp) = c("Subgrps.temp", "Subgrps.est")
        temp$Rank = rank(temp$Subgrps.est)
        temp = temp[,-2]
        levels_rank = temp %>% arrange(temp$Rank) %>% select(temp$Rank)
        levels_rank = levels_rank[,1]
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
        step3X = do.call( SubMod, append(list(Y=Y, A=A, X=X.star.temp, Xtest=Xtest.star.temp,
                                              mu_train = mu_train, family=family),
                                              SubMod.hyper) )
        ## Obtain the "Merged" Rules ##
        temp.rules = Rules
        colnames(temp.rules) = c("Subgrps.temp", "Rules")
        rules.dat = data.frame(Subgrps.temp = step3$Subgrps.train,
                               Subgrps.new = step3X$Subgrps.train)
        rules.dat = left_join(rules.dat, temp.rules, by="Subgrps.temp")
        rules.dat = unique(rules.dat[,c("Subgrps.new", "Rules")])
        rules.dat <- aggregate(Rules ~ Subgrps.new, data = rules.dat, paste, collapse = " , ")
        colnames(rules.dat) = c("Subgrps", "Rules")
        step3X$Rules = rules.dat
        step3 = step3X
      }
      Sub.mod = step3$mod; Rules=step3$Rules;
      Subgrps.train = step3$Subgrps.train; Subgrps.test = step3$Subgrps.test
      Subpred.train = step3$pred.train; Subpred.test = step3$pred.test
    }
    if ( is.null(SubMod) ){
      Sub.mod = NULL; Rules=NULL;
      Subgrps.train = NULL; Subgrps.test = NULL;
      Subpred.train = NULL; Subpred.test = NULL
    }
    ### Step 4: Parameter Estimation and Inference ###
    if (verbose){ message(paste("Parameter Estimation:", Param, sep=" ")) }
    param.dat = do.call( Param, list(Y=Y, A=A, X=X.star, mu_hat = mu_train,
                                     Subgrps=Subgrps.train,
                                     alpha_ovrl=alpha_ovrl,
                                     alpha_s=alpha_s)  )
    ### Return Outputs ###
    return( list( mu_train = mu_train, mu_test = mu_test, filter.mod = filter.mod,
                  filter.vars = filter.vars, Sub.mod=Sub.mod,
                  Subgrps.train = Subgrps.train, Subgrps.test=Subgrps.test,
                  Subpred.train = Subpred.train, Subpred.test = Subpred.test,
                  Rules = Rules, param.dat=param.dat) )
  }
  ## Run on Observed Data ##
  if (verbose){ message( "Observed Data" )   }
  res0 = main(Y=Y, A=A, X=X, Xtest=Xtest, PLE=PLE, Filter=Filter,
              SubMod = SubMod, verbose=verbose)
  Subgrps = res0$Subgrps.train
  mu_train = res0$mu_train
  param.dat = res0$param.dat

  #################################################################
  #### Resampling (Bootstrapping or Permutation)    ###############
  #################################################################

  resamp_param = NULL ## Set to null (needed if no resampling)

  if (!is.null(Resample)){

    ## Resampling: Auto-checks ##
    if (is.null(Filter.resamp)) {  Filter.resamp = Filter  }
    if (is.null(PLE.resamp)) {  PLE.resamp = PLE  }
    if (is.null(SubMod.resamp)) {  SubMod.resamp = SubMod  }

    obs.data = data.frame(id = 1:nrow(X), Y=Y,A=A, X, Subgrps)
    if (length(res0$filter.vars)>0 & preFilter_resamp == TRUE){
      obs.data = obs.data[, colnames(obs.data) %in%
                            c("id", "Y", "A", res0$filter.vars, "Subgrps")]
    }
    # Xtest.star = Xtest[,colnames(Xtest) %in% res0$filter.vars]
    if (verbose){ message( paste(Resample, R, "resamples"))   }

    subject.counter = data.frame(id = obs.data$id )

    ### Resampling Wrapper ###
    fetty_wop = function(R, stratify, obs.data, PLE, Filter, SubMod,
                         calibrate, verbose){

      if (verbose) message( paste(Resample, "Sample", R) )
      ### Permutation resampling (shuffle treatment assignment) ###
      if (Resample=="Permutation"){
        seedz = 5341+R
        set.seed(seedz)
        A_resamp = sample(obs.data$A, replace=FALSE)
        resamp.data = obs.data
        resamp.data$A = A_resamp
      }
      if (Resample=="Bootstrap"){
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
      if (family=="survival"){ Y = Surv(Y.R[,1], Y.R[,2])  }
      res.R = main(Y=Y.R, A=A.R, X=X.R, Xtest=Xtest, PLE=PLE, Filter=Filter,
                   SubMod = SubMod, verbose=FALSE)
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
      param.resamp = rbind( param.R[1:2,c("Subgrps", "est")], param.resamp )
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
                        PLE=PLE.resamp, Filter=Filter.resamp,
                        SubMod=SubMod.resamp, calibrate=FALSE,
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
      if (Resample=="Permutation"){
        final_ests$pval_perm = NA
      }
      if (Resample=="Bootstrap"){
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
        if (Resample=="Permutation"){
          final_ests$pval_perm[final_ests$Subgrps==sub] =
            (sum(abs(est.vec)>abs(est0)) + 1 ) / (length(est.vec)+1)
        }
        ## Bootstrap Covariance/acceleration/bias/smoothed SE ##
        if (Resample=="Bootstrap"){
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
  res = list( mu_train=res0$mu_train, mu_test=res0$mu_test,
              filter.mod = res0$filter.mod, filter.vars = res0$filter.vars,
              Sub.mod = res0$Sub.mod,
              out.train = data.frame(Y, A, X, Subgrps=res0$Subgrps.train),
              out.test = data.frame(Xtest, Subgrps=res0$Subgrps.test),
              Subpred.train = res0$Subpred.train,
              Subpred.test = res0$Subpred.test, Rules=res0$Rules,
              param.dat = param.dat, resamp.dist = resamp_param)
  class(res) <- c("PRISM")
  return(res)
}
