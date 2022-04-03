# PRISM (Train): Patient Response Identifier for Stratified Medicine
PRISM_train <- function(Y, A, X, Xtest=NULL, family="gaussian",
                       mu_train=NULL, filter="glmnet", 
                       ple="ranger", submod="lmtree", 
                       param="dr", meta = "X-learner", pool="no", 
                       delta=">0", propensity = FALSE, combine = "SS", 
                       resample_submod = NULL, R_submod=NULL, 
                       resample_pool = NULL, R_pool=NULL, 
                       alpha_ovrl=0.05, alpha_s = 0.05, 
                       filter.hyper=NULL, ple.hyper=NULL, 
                       submod.hyper = NULL, verbose=TRUE, mu_hat=NULL, 
                       efficient=TRUE) {
  
  if (is.null(A)) {
    pool <- "no"
  }
  ### Step 1: Variable Filtering ###
  if (!(filter=="None")) {
    if (verbose) message(paste("Filtering:", filter, sep=" "))
    step1 <- tryCatch(filter_train(Y=Y, A=A, X=X, family = family, 
                                   filter=filter, hyper = filter.hyper),
                       error = function(e) paste("filter error", e ) )
    if (is.character(step1)) {
      if (verbose) message(step1)
      if (verbose) message("Using all variables")
      filter.mod <- NULL;
      filter.vars <- NULL
    }
    if (is.list(step1)) {
      filter.mod <- step1$mod
      filter.vars <- step1$filter.vars 
    }
  }
  if (filter=="None") {
    filter.mod <- NULL; filter.vars <- NULL;
  }
  # Drop variables depending on filter #
  if (filter=="None") {X.star <- X; Xtest.star <- Xtest}
  if (!(filter=="None")) {
    if (length(filter.vars)==0) {
      X.star <- X 
      Xtest.star <- Xtest
    }
    if (length(filter.vars)>0) {
      X.star <- X[,colnames(X) %in% c(filter.vars, "A", "Y"), drop=FALSE]
      if (is.null(Xtest)) {
        Xtest.star <- Xtest
      }
      if (!is.null(Xtest)) {
        Xtest.star <- Xtest[,colnames(Xtest) %in% 
                              c(filter.vars, "A", "Y"), drop=FALSE] 
      }
    }
  }
  ### Step 2: PLE ###
  if (!(ple=="None")) {
    if (verbose) message( paste("Counterfactual Estimation:", ple, 
                                paste("(", meta, ")", sep=""), sep=" "))
    step2 <- tryCatch(
      ple_train(Y=Y,A=A,X=X.star,Xtest=Xtest.star,family=family,
                ple=ple, meta=meta, propensity=propensity, 
                hyper = ple.hyper), 
      error = function(e) paste("ple error", e ) )
    
    if (is.character(step2)) {
      if (verbose) message(step2)
      ple.fit <- NULL
      mu_train <- NULL
      mu_test <- NULL
    }
    if (is.list(step2)) {
      ple.fit <- step2$fit
      mu_train <- step2$mu_train
      mu_test <- step2$mu_test 
    }
  }
  if (ple=="None") {
    ple.fit <- NULL
    mu_train <- mu_train
    mu_test <- NULL
  }

  ### Step 3: Subgroup Identification ###
  if (!(submod=="None")) {
    if (verbose) {
      message(paste("Subgroup Identification:", submod, sep=" "))
      if (pool!="no") {
        message(paste("Pooling:", pool, sep = " "))
      }
      message(paste("Treatment Effect Estimation:", param, sep=" "))
    }
    step3 <- tryCatch(submod_train(Y=Y, A=A, X=X.star, Xtest=Xtest.star, 
                                   mu_train=mu_train, family = family, 
                                   submod = submod, hyper = submod.hyper,
                                   ple = ple, ple.hyper=ple.hyper,
                                   meta = meta,
                                   propensity = propensity, 
                                   pool=pool, delta = delta, param = param,
                                   combine = combine,
                                   resample = resample_submod,
                                   R = R_submod,
                                   resample_pool = resample_pool, 
                                   R_pool = R_pool,
                                   verbose.resamp = FALSE,
                                   alpha_ovrl = alpha_ovrl, alpha_s = alpha_s, 
                                   efficient = efficient),
                      error = function(e) paste("submod error:", e) )
    # print(step3)
    
    if (is.character(step3)) {
      if (verbose) message(step3)
      if (verbose) message("Settings Subgrps=1")
      submod.fit <- NULL
      Rules <- NULL
      Subgrps.train <- rep(1, dim(X)[1])
      Subgrps.test <- rep(1, dim(Xtest)[1])
      trt_assign <- NULL; trt_eff <- NULL
      trt_eff_obs <- NULL; trt_eff_pool <- NULL
      trt_eff_dopt <- NULL; trt_eff_resamp <- NULL
      submod_rdist <- NULL; resamp_subgrps <- NULL
    }
    if (is.list(step3)) {
      submod.fit <- step3$fit
      Rules <- step3$Rules;
      Subgrps.train <- as.character(step3$Subgrps.train)
      Subgrps.test <- as.character(step3$Subgrps.test)
      trt_assign <- step3$trt_assign
      trt_eff <- step3$trt_eff
      trt_eff_obs <- step3$trt_eff_obs
      trt_eff_pool <- step3$trt_eff_pool
      trt_eff_dopt <- step3$trt_eff_dopt
      trt_eff_resamp <- step3$trt_eff_resamp
      submod_rdist <- step3$resamp_dist
      resamp_subgrps <- step3$resamp_subgrps
    }
  }
  if (submod=="None") {
    submod.fit <- NULL; Rules <- NULL;
    Subgrps.train <- NULL; Subgrps.test <- NULL;
    trt_assign <- NULL; trt_eff <- NULL; 
    trt_eff_obs <- NULL; trt_eff_resamp <- NULL;
    trt_eff_pool <- NULL; trt_eff_dopt <- NULL;
    submod_rdist <- NULL; resamp_subgrps <- NULL
  }
  ### Step 4: Parameter Estimation and Inference ###
  param.dat <- trt_eff
  param.dat$Subgrps <- as.character(param.dat$Subgrps)
  
  if (!is.character(param.dat)) {
    if (is.null(param.dat$estimand)) {
      param.dat$estimand <- "est"
    }
  }
  
  output <- list(mu_train = mu_train, mu_test = mu_test, filter.mod = filter.mod,
                 filter.vars = filter.vars, ple.fit = ple.fit, submod.fit = submod.fit,
                 Subgrps = Subgrps.train, 
                 Subgrps.train = Subgrps.train, Subgrps.test=Subgrps.test,
                 Rules = Rules, param.dat=param.dat, trt_assign=trt_assign,
                 trt_eff = trt_eff, trt_eff_obs = trt_eff_obs,
                 trt_eff_pool = trt_eff_pool, trt_eff_dopt = trt_eff_dopt,
                 trt_eff_resamp = trt_eff_resamp, submod_rdist=submod_rdist, 
                 resamp_subgrps = resamp_subgrps)
  class(output) <- "PRISM_train"
  return(output)
}