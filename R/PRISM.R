#' PRISM: Patient Response Identifier for Stratified Medicine
#'
#' PRISM algorithm. Given a data-set of (Y, A, X) (Outcome, treatment, covariates),
#' the \code{PRISM} identifies potential subgroups along with point-estimate and variability
#' metrics; with and without resampling (bootstrap or cross-validation based). This four 
#' step procedure (filter, ple, submod, param) is flexible and accepts user-inputs at each
#' step.
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (Default supports binary treatment, either numeric or 
#' factor). "ple_train" accomodates >2 along with binary treatments.
#' @param X Covariate space. 
#' @param Xtest Test set. Default is NULL (no test predictions). Variable types should match X.
#' @param mu_train Patient-level estimates in training set (see \code{ple_train}). 
#' Default=NULL
#' @param family Outcome type. Options include "gaussion" (default), "binomial", and "survival".
#' @param filter Filter model to determine variables that are likely associated with the 
#' outcome and/or treatment. Outputs a potential reduce list of varia where X.star 
#' has potentially less variables than X. Default is "glmnet" (elastic net). Other 
#' options include "ranger" (random forest based variable importance with p-values).
#' See \code{filter_train} for more details. "None" uses no filter. 
#' @param ple Base-learner used to estimate patient-level equantities, such as the 
#' conditional average treatment effect (CATE), E(Y|A=1,X)-E(Y|A=0, X) = CATE(X). 
#' Default is random based based through "ranger". "None" uses no ple. See below for 
#' details on estimating the treatment contrasts.
#' @param meta Using the ple model as a base learner, meta-learners can be used for 
#' estimating patient-level treatment differences. Options include "T-learner" (treatment
#' specific models), "S-learner" (single model), and "X-learner". For family="gaussian" &
#' "binomial", the default is "X-learner", which uses a two-stage regression 
#' approach (See Kunzel et al 2019). For "survival", the default is "T-learner". "X-learner" 
#' is currently not supported for survival outcomes.
#' @param submod Subgroup identification model function. Options include tree-methods that 
#' target the treatment by variable interaction directly ("lmtree", "glmtree", "mob_weib"),
#' regress the CATE ("rpart_cate", "ctree_cate"), and target prognostic variables ("rpart", "ctree").
#' Default for family="gaussian" is "lmtree" (MOB with OLS loss). For "binomial" the default is 
#' "glmtree" (MOB with binomial loss). Default for "survival" is "lmtree" (log-rank transformation 
#' on survival outcomes and then fit MOB-OLS). "None" uses no submod. Currently only available for 
#' binary treatments or A=NULL.
#' @param param Parameter estimation and inference function. Based on the discovered 
#' subgroups, estimate parameter estimates and correspond variability metrics. Options
#' include "lm" (unadjusted linear regression), "dr" (doubly-robust estimator),
#' "gcomp" (G-computation, average the patient-level estimates), "cox" (cox regression),
#' and "rmst" (RMST based estimates as in survRMST package). Default for "gaussian",
#' "binomial" is "dr", while default for "survival" is "cox". Currently only available 
#' for binary treatments or A=NULL. 
#' @param pool Whether to pool the initial identified subgroups (ex: tree nodes).
#' Default = "no". Other options include "trteff" or "trteff_boot" (check if 
#' naive or bootstrap treatment estimate is beyond clinical meaningful 
#' threshold delta, ex: trteff_boot > 0), and optimal treatment regime (OTR) pooling, 
#' "otr:logistic", "otr:rf". "otr:logistic" fits weighted logistic regression 
#' with I(mu_1-mu_0>delta) as the outcome, the candidate subgroups as covariates, 
#' and weights=abs((mu_1-mu_0) - delta). "otr:rf" follows the same approach but 
#' with weighted random forest, and also includes X in the regression. Regardless 
#' of the pooling approach, the key output is "trt_assign", a data-frame with the 
#' initial subgroups and the pooled subgroups (ex: dopt=1, patient should receive
#'  A=1, vs dopt=0, patient should receive A=0).
#' @param delta Threshold for defining benefit vs non-benefitting patients. 
#' Only applicable for submod="otr", and if pooling is used (see "pool"). 
#' Default=">0".
#' @param propensity Propensity score estimation, P(A=a|X). Default=FALSE which 
#' use the marginal estimates, P(A=a) (applicable for RCT data). If TRUE, will 
#' use the "ple" base learner to estimate P(A=a|X). 
#' @param combine Method of combining group-specific point-estimates. Options 
#' include "SS" (sample size weighting), and "maxZ" 
#' (see: Mehrotra and Marceau-West). This is used for pooling (ex: within dopt=1
#' groups, aggregate group-specific treatment estimates), and for calculating 
#' the overall population treatment effect estimate.
#' @param alpha_ovrl Two-sided alpha level for overall population. Default=0.05
#' @param alpha_s Two-sided alpha level at subgroup level. Default=0.05
#' @param filter.hyper Hyper-parameters for the filter function (must be list). 
#' Default is NULL.
#' @param ple.hyper Hyper-parameters for the PLE function (must be list). 
#' Default is NULL.
#' @param submod.hyper Hyper-parameters for the submod function (must be list). 
#' Default is NULL.
#' @param resample Resampling method for resample-based treatment effect estimates 
#' and variability metrics. Options include "Bootstrap" and 
#' "CV" (cross-validation). Default=NULL (No resampling).
#' @param stratify Stratified resampling? Default="trt" (stratify by A). Other
#' options include "sub" (stratify by the identified subgroups), "trt_sub" 
#' (stratify by A and the identified subgroups), and "no" (no stratification).
#' @param R Number of resamples (default=NULL; R=100 for Permutation/Bootstrap 
#' and R=5 for CV). This resamples the entire PRISM procedure. 
#' @param resample_submod For submod only, resampling method for treatment effect estimates.
#' Options include "Bootstrap" or NULL (no resampling).
#' @param R_submod Number of resamples for resample_submod
#' @param resample_pool For submod only, resampling method for pooling step. 
#' nly applicable if resample_submod="Bootstrap" and/or pool="trteff_boot". 
#' @param R_pool Number of resamples for resample_pool
#' @param calibrate Bootstrap calibration for nominal alpha (Loh et al 2016).
#' Default=FALSE. For TRUE, outputs the calibrated alpha level and calibrated 
#' CIs for the overall population and subgroups. Not applicable for permutation 
#' or CV resampling.
#' @param alpha.mat Grid of alpha values for calibration. Default=NULL, which 
#' uses seq(alpha/1000,alpha,by=0.005) for alpha_ovrl/alpha_s. 
#' @param filter.resamp Filter function during re-sampling. Default=NULL 
#' (uses "filter"). If "None", the "filter" model is not trained in each resample,
#' and instead use filtered variables from the observed data "filter" step. 
#' (less computationally expensive).
#' @param ple.resamp Ple function during re-sampling. Default=NULL 
#' (uses "ple"). If "None", the "ple" model is not training in each resample, 
#' and instead the original model estimates are resampled (less computationally
#'  expensive).
#' @param verbose Detail progress of PRISM? Default=TRUE
#' @param verbose.resamp Output iterations during resampling? Default=FALSE
#' @param seed Seed for PRISM run (Default=777) 
#' @param efficient If TRUE (default for PRISM), then models (filter, ple, submod) will 
#' store reduced set of outputs for faster speed. 
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
#'   \item resamp_dist - Resampling distributions (NULL if no resampling is done)
#' }
#' @export
#' @importFrom stats aggregate coef lm glm model.matrix p.adjust pnorm confint
#' @importFrom stats predict pt qnorm qt quantile sd weighted.mean vcov na.omit
#' @importFrom survival is.Surv survfit survreg Surv
#' @import dplyr
#' 
#' @details PRISM is a general framework with five key steps:
#' 
#' 0. Estimand: Determine the question of interest (ex: mean treatment difference)
#' 
#' 1. Filter (filter): Reduce covariate space by removing noise covariates. Options include 
#' elastic net ("glmnet") and random forest variable importance ("ranger").
#' 
#' 2. Patient-Level Estimates (ple): Estimate counterfactual patient-level quantities, 
#' for example, the conditional average treatment effect (CATE), E(Y|A=1,X)-E(Y|A=0,X). 
#' This calls the "ple_train" function, and follows the framework of Kunzel et al 2019. 
#' Base-learners include random forest ("ranger"), BART ("bart"), elastic net ("glmnet"), 
#' and linear models (LM, GLM, or Cox regression). Meta-learners include the "S-Learner" 
#' (single model), "T-learner" (treatment specific models), and "X-learner" (2-stage approach).
#' 
#' 3. Subgroup Model (submod): Currently uses tree-based methods to identify predictive 
#' and/or prognostic subgroups. Options include MOB OLS ("lmtree"), MOB GLM ("glmtree"), 
#' MOB Weibull ("mob_weib"), conditional inference trees ("ctree", Y~ctree(X); "ctree_cate", 
#' CATE~ctree(X)), and recursive partitioning and regression trees ("rpart", Y~rpart(X); "rpart_cate", 
#' CATE~rpart(X)), and optimal treatment regimes ("otr").
#' 
#' 4. Treatment Effect Estimation (param): For the overall population and the discovered 
#' subgroups (if any), obtain treatment effect point-estimates and variability metrics. 
#' Options include: cox regression ("cox"), double robust estimator ("dr"), linear regression 
#' ("lm"), average of patient-level estimates ("gcomp"), and restricted mean survival 
#' time ("rmst").
#' 
#' Steps 1-4 also support user-specific models. If treatment is provided (A!=NULL), 
#' the default settings are as follows:
#' 
#' Y is continuous (family="gaussian"): 
#' Elastic Net Filter ==> X-learner with random forest ==> MOB (OLS) ==> 
#' Double Robust estimator
#' 
#' Y is binary (family="binomial"): 
#' Elastic Net Filter ==> X-learner with random forest ==> MOB (GLM) ==> 
#' Double Robust estimator
#' 
#' Y is right-censored (family="survival"):
#' Elastic Net Filter ==> T-learner with random forest ==> MOB (Weibull) ==> 
#' Cox regression
#' 
#' 
#' If treatment is not provided (A=NULL), the default settings are as follows:
#' 
#' Y is continuous (family="gaussian"): 
#' Elastic Net Filter ==> Random Forest ==> ctree ==> linear regression
#' 
#' Y is binary (family="binomial"): 
#' Elastic Net Filter ==> Random Forest ==> ctree ==> linear regression
#' 
#' Y is right-censored (family="survival"):
#' Elastic Net Filter ==> Survival Random Forest ==> ctree ==> RMST
#' 
#' 
#' @references Friedman, J., Hastie, T. and Tibshirani, R. (2008) Regularization Paths for
#' Generalized Linear Models via Coordinate Descent, 
#' \url{https://web.stanford.edu/~hastie/Papers/glmnet.pdf} Journal of Statistical 
#' Software, Vol. 33(1), 1-22 Feb 2010 Vol. 33(1), 1-22 Feb 2010.
#' @references Jemielita T, Mehrotra D. PRISM: Patient Response Identifiers for 
#' Stratified Medicine. \url{https://arxiv.org/abs/1912.03337}
#' @references Hothorn T, Hornik K, Zeileis A (2006). Unbiased Recursive Partitioning: 
#' A Conditional Inference Framework. Journal of Computational and Graphical Statistics,
#' 15(3), 651–674.
#' @references Wright, M. N. & Ziegler, A. (2017). ranger: A fast implementation of 
#' random forests for high dimensional data in C++ and R. J Stat Softw 77:1-17. 
#' \doi{10.18637/jss.v077.i01}.
#' @references Zeileis A, Hothorn T, Hornik K (2008). Model-Based Recursive Partitioning. 
#' Journal of Computational and Graphical Statistics, 17(2), 492–514.
#' @examples
#' \donttest{
#' ## Load library ##
#' library(StratifiedMedicine)
#'
#' ## Examples: Continuous Outcome ##
#'
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#'
#' # Run Default: glmnet, ranger (X-learner), lmtree, dr #
#' res0 = PRISM(Y=Y, A=A, X=X)
#' summary(res0)
#' plot(res0)
#' 
#' res1 = PRISM(Y=Y, A=A, X=X, filter="None")
#' summary(res1)
#' plot(res1)
#' }
#' 
#' # Search for Prognostic Only (omit A from function) #
#' \donttest{
#' res3 = PRISM(Y=Y, X=X)
#' summary(res3)
#' plot(res3)
#' }
#'
#' ## With bootstrap (No filtering) ##
#' \donttest{
#' library(ggplot2)
#'   res_boot = PRISM(Y=Y, A=A, X=X, resample = "Bootstrap", R=50, verbose.resamp = TRUE)
#'   # Plot of distributions and P(est>0) #
#'   plot(res_boot, type="resample", estimand = "E(Y|A=1)-E(Y|A=0)")+
#'   geom_vline(xintercept = 0)
#'   aggregate(I(est>0)~Subgrps, data=res_boot$resamp_dist, FUN="mean")
#' }
#' 
#' ## Examples: Binary Outcome ##
#' \donttest{
#' dat_bin = generate_subgrp_data(family="binomial")
#' Y = dat_bin$Y
#' X = dat_bin$X
#' A = dat_bin$A
#'
#' # Run Default: glmnet, ranger, glmtree, dr #
#' res0 = PRISM(Y=Y, A=A, X=X)
#' 
#' plot(res0)
#' }
#'
#' # Survival Data ##
#' \donttest{
#'   library(survival)
#'   library(ggplot2)
#'   require(TH.data); require(coin)
#'   data("GBSG2", package = "TH.data")
#'   surv.dat = GBSG2
#'   # Design Matrices ###
#'   Y = with(surv.dat, Surv(time, cens))
#'   X = surv.dat[,!(colnames(surv.dat) %in% c("time", "cens")) ]
#'   set.seed(513)
#'   A = rbinom( n = dim(X)[1], size=1, prob=0.5  )
#'
#'   # PRISM: glmnet ==> Random Forest to estimate Treatment-Specific RMST
#'   # ==> MOB (Weibull) ==> Cox for HRs#
#'   res_weib = PRISM(Y=Y, A=A, X=X)
#'   plot(res_weib, type="PLE:waterfall")
#'   plot(res_weib)
#'   
#'   #PRISM: glmnet ==> Random Forest to estimate Treatment-Specific RMST
#'   #RPART_CATE: Regress RMST on RPART for subgroups #
#'   res_cate = PRISM(Y=Y, A=A, X=X, submod="rpart_cate")
#'   plot(res_cate)
#'
#'   # PRISM: ENET ==> CTREE ==> Cox; with bootstrap #
#'   res_ctree1 = PRISM(Y=Y, A=A, X=X, ple="None", submod = "ctree",
#'                      resample="Bootstrap", R=50, verbose.resamp = TRUE)
#'   plot(res_ctree1)
#'   plot(res_ctree1, type="resample", estimand="HR(A=1 vs A=0)")+geom_vline(xintercept = 1)
#'   aggregate(I(est<0)~Subgrps, data=res_ctree1$resamp_dist, FUN="mean")
#' }
#'

##### PRISM: Patient Responder Identifiers for Stratified Medicine ########
PRISM <- function(Y, A=NULL, X, Xtest=NULL, mu_train=NULL, family="gaussian",
                 filter="glmnet", ple="ranger", submod=NULL, param=NULL,
                 meta = ifelse(family=="survival", "T-learner", "X-learner"),
                 pool="no", delta = ">0", 
                 propensity = FALSE,
                 combine = "SS",
                 alpha_ovrl=0.05, alpha_s = 0.05,
                 filter.hyper=NULL, ple.hyper=NULL, submod.hyper = NULL,
                 resample = NULL, stratify=ifelse(!is.null(A), "trt", "no"),
                 R = NULL, resample_submod = NULL, R_submod=NULL, 
                 resample_pool = NULL, R_pool=NULL,
                 calibrate=FALSE, alpha.mat=NULL,
                 filter.resamp = NULL, ple.resamp = NULL,
                 verbose=TRUE,
                 verbose.resamp = FALSE, seed=777, 
                 efficient = TRUE){

  func_call <- match.call()
  
  if (is.null(A)){
    message("No Treatment Variable (A) Provided: Searching for Prognostic Effects")
    if (!is.null(submod)){
      if (submod %in% c("lmtree", "glmtree", "otr", "mob_weib")){
      message( paste(submod, 
              "not usable without (A): Using ctree") )
      }
    }
    meta <- "none"
  }

  ## Is the Outcome Survival? ##
  if (is.Surv(Y)) {
    if (family!="survival"){ family <- "survival" }
  }
  ## Is the outcome binary? #
  if (!is.Surv(Y)) {
    if ( mean( unique(Y) %in% c(0,1) )==1 ){ family <- "binomial" }
  }
  # Missing data? #
  if (sum(is.na(Y))>0 | sum(is.na(A))>0 | sum(is.na(X))>0) {
    stop("Missing Data (in outcome, treatment, or covariates)")
  }
  ### Defaults: By Family (gaussian, binomial (Risk Difference), survival ) ##
  if (family=="gaussian" | family=="binomial") {
    if (is.null(ple)) { ple <- "ranger" }
    if (is.null(submod)) { 
      if (family=="gaussian"){ submod <- "lmtree"}
      if (family=="binomial"){ submod <- "glmtree"}
      if (is.null(A)) { submod <- "ctree" }
    }
    if (is.null(param)) { 
      if (is.null(A)){ param <- "lm" }
      else { param <- "dr" }
    }
  }
  if (family=="survival") {
    if (is.null(ple)) {ple <- "ranger"}
    if (is.null(submod)){ 
      if (is.null(A)) {submod <- "ctree"}
      else {submod <- "lmtree"}
      }
    if (is.null(param)) { 
      if (is.null(A)) {param <- "rmst"}
      else {param <- "cox"}
    }
    # Use T-learner if X-learner is specified #
    if (meta=="X-learner") {
      meta <- "T-learner"
    }
  }
  
  list_args <- list(family = family, filter = filter, 
                    ple = ple, submod = submod, param = param,
                    meta = meta, pool = pool,
                    delta = delta, propensity = propensity,
                    combine = combine,
                    resample_submod = resample_submod, R_submod=R_submod,
                    resample_pool = resample_pool, R_pool=R_pool,
                    alpha_ovrl = alpha_ovrl, alpha_s = alpha_s, 
                    filter.hyper = filter.hyper, ple.hyper = ple.hyper,
                    submod.hyper = submod.hyper, verbose = verbose, 
                    efficient = efficient)
  
  ## Train PRISM on Observed Data (Y,A,X) ##
  if (verbose) {message("Observed Data")}
  set.seed(seed)
  fit0 <- do.call("PRISM_train", 
                 append(list(Y=Y, A=A, X=X, Xtest = Xtest), list_args))
  Subgrps <- fit0$Subgrps.train
  mu_train <- fit0$mu_train
  param.dat <- fit0$param.dat
  param.dat <- param.dat[order(param.dat$Subgrps, param.dat$estimand),]
  resamp_dist <- NULL
  resamp_subgrps <- NULL 
  resamp_calib <- NULL
  
  ## Resampling ##
  if (is.null(R) & !is.null(resample)){
    if (resample %in% c("Permutation", "Bootstrap")){R <- 50}
    if (resample == "CV" ) {R <- 5}
  }
  if ( !is.null(resample)) {
    if (verbose){ 
      if (resample=="CV") {
        message( paste("Cross Validation", R, "folds") )
      }
      if (resample!="CV") {
        message( paste(resample, R, "resamples"))   }
    }
    ## Resampling: Auto-checks ##
    if (is.null(filter.resamp)) {  filter.resamp <- filter  }
    if (is.null(ple.resamp)) {  ple.resamp <- ple  }
    list_args_R <- list_args
    list_args_R$filter <- filter.resamp
    list_args_R$ple <- ple.resamp
    
    fitR <- resampler_general(Y=Y, A=A, X=X, fit=fit0, 
                              wrapper="PRISM_train",
                              list_args=list_args_R, 
                              resample = resample, R=R, 
                              stratify=stratify, 
                              fixed_subs = "false", 
                              calibrate = calibrate, 
                              alpha.mat = alpha.mat, 
                              verbose=verbose.resamp)
    param.dat <- fitR$trt_eff_resamp
    resamp_dist <- fitR$resamp_dist
    resamp_subgrps <- fitR$resamp_subgrps
  }
  # Add Details to Trt Estimates #
  param.dat <- .add_trt_details(Y=Y, Subgrps=Subgrps, trt_dat = param.dat, family=family)
  
  if (is.null(A)) {
    out.train <- data.frame(Y, X, Subgrps=fit0$Subgrps.train)
  }
  if (!is.null(A)) {
    out.train <- data.frame(Y, A, X, Subgrps=fit0$Subgrps.train)
  }
  if (is.null(Xtest)) { 
    out.test <- NULL
  } 
  else { 
    out.test <- data.frame(Xtest, Subgrps=fit0$Subgrps.test) 
  }
    
  ### Return Results ##
  res = list(call = func_call, filter.mod = fit0$filter.mod, filter.vars = fit0$filter.vars,
             ple.fit = fit0$ple.fit, mu_train=fit0$mu_train, mu_test=fit0$mu_test,
             submod.fit = fit0$submod.fit,
             out.train = out.train,
             out.test = out.test,
             Rules=fit0$Rules,
             param.dat = param.dat, 
             pool=pool, trt_assign=fit0$trt_assign,
             trt_eff = fit0$trt_eff, trt_eff_obs = fit0$trt_eff_obs,
             trt_eff_pool = fit0$trt_eff_pool, trt_eff_dopt = fit0$trt_eff_dopt, 
             trt_eff_resamp = fit0$trt_eff_resamp, 
             resample=resample, resamp_dist = resamp_dist, 
             resamp_subgrps = resamp_subgrps, 
             resamp_calib = resamp_calib, 
             submod_rdist = fit0$submod_rdist,
             resamp_subgrps_submod = fit0$resamp_subgrps,
             family = family,
             filter = filter, ple = ple, submod=submod, param=param,
             meta = meta, combine = combine,
             alpha_ovrl = alpha_ovrl, alpha_s = alpha_s,
             filter.hyper = filter.hyper, ple.hyper = ple.hyper, 
             submod.hyper = submod.hyper)
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
#' \donttest{
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
#' # Run Default: filter_glmnet, ple_ranger, lmtree, param_ple #
#' res0 = PRISM(Y=Y, A=A, X=X)
#' summary(predict(res0, X)) # all #
#' summary(predict(res0, X, type="ple"))
#' summary(predict(res0, X, type="submod"))
#' }
#'
#'
#' @method predict PRISM
#' @export
#'
predict.PRISM = function(object, newdata=NULL, type="all", ...){

  if (type=="all"){
    mu_hat = object$ple.fit$pred.fun(object$ple.fit$mod, newdata)
    Subgrps = object$submod.fit$pred.fun(object$submod.fit$mod, newdata)$Subgrps
    res = data.frame(Subgrps = Subgrps, mu_hat)
  }
  if (type=="ple"){
    mu_hat = object$ple.fit$pred.fun(object$ple.fit$mod, newdata)
    res = data.frame(mu_hat)
  }
  if (type=="submod"){
    Subgrps = object$submod.fit$pred.fun(object$submod.fit$mod, newdata)$Subgrps
    res = data.frame(Subgrps=Subgrps)
  }
  return(res)
}

#' PRISM: Patient Response Identifier for Stratified Medicine (Summary)
#'
#' Summary for PRISM algorithm results. Outputs configuration, which variables pass the filter (if used),
#' subgroup summaries, and treatment effect estimates.
#'
#' @param object Trained PRISM model.
#' @param round_est Rounding for trt ests (default=4)
#' @param round_SE Rounding for trt SEs (default=4)
#' @param round_CI Rounding for trt CIs (default=4)
#' @param ... Any additional parameters, not currently passed through.
#'
#' @return List of key PRISM outputs: (1) Configuration, (2) Variables that pass filter 
#' (if filter is used), (3) Number of Identified Subgroups, and (4) Parameter Estimates, 
#' SEs, and CIs for each subgroup/estimand
#' 
#' @method summary PRISM
#' @export
#' 
summary.PRISM = function(object, round_est=4,
                         round_SE=4, round_CI=4,...){

  out <- NULL
  alpha_ovrl <- object$alpha_ovrl
  alpha_s <- object$alpha_s
  # Call #
  out$call <- object$call
  # Configuration #
  out$`PRISM Configuration` <- with(object, paste("[Filter]", filter, "=>",
                                                  "[PLE]", ple, "=>",
                                                  "[Subgroups]", submod, "=>",
                                                  "[Param]", param, sep=" "))
  if (!is.null(object$resample)) {
    out$`PRISM Configuration` <- paste(out$`PRISM Configuration`, object$resample, sep="=> ")
  }
  # Filter #
  if (object$filter!="None"){
    out$`Variables that Pass Filter` <- object$filter.vars
  }
  # Subgroup summary #
  out_submod <- summary_output_submod(object, round_est=round_est,
                                      round_SE=round_SE, round_CI=round_CI)
  out <- append(out, out_submod)
  
  # # Final Estimates Summary #
  # if (!is.null(object$resample)) {
  #   
  #   param.dat <- object$param.dat
  #   trt_final <- param.dat
  #   trt_final <- trt_final[,c("Subgrps", "N", "estimand", "est_resamp", "SE_resamp")]
  #   
  # }
  
  if (!is.null(object$resample)) {
    
    param.dat <- trt_data_cleaner(object$param.dat, round_est=round_est,
                                    round_SE=round_SE, round_CI=round_CI)
    if (object$resample=="Bootstrap") {
      out$`Treatment Effect Estimates (bootstrap)` <- param.dat
    }
    if (object$resample=="CV") {
      out$`Treatment Effect Estimates (CV)` <- param.dat
    }
    if (object$resample=="Permutation") {
      out$`Treatment Effect Estimates (Permutation)` <- param.dat
    }
    
  }
  class(out) <- "summary.PRISM"
  out
}