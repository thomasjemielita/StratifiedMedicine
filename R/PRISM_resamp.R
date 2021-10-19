# PRISM (Resample): Patient Response Identifier for Stratified Medicine
PRISM_resamp <- function(PRISM.fit, Y, A, X, Xtest=NULL, family="gaussian",
                       filter="glmnet", ple="ranger", submod=NULL, param=NULL,
                       meta = "X-learner",
                       pool = "no",
                       delta = ">0", propensity = FALSE,
                       combine = "SS",
                       resample_submod = NULL, R_submod=NULL,
                       alpha_ovrl=0.05, alpha_s = 0.05,
                       filter.hyper=NULL, ple.hyper=NULL, submod.hyper = NULL,
                       verbose=TRUE,
                       resample="Bootstrap", 
                       R=NULL, stratify = "trt", calibrate=TRUE,
                       alpha.mat=NULL) {
  
  if (is.null(alpha.mat)) {
    alpha.mat = data.frame(cbind( unique(c(seq(alpha_ovrl/1000, alpha_ovrl,
                                               by=0.005),alpha_ovrl) ),
                                  unique(c(seq(alpha_s/1000, alpha_s,
                                               by=0.005),alpha_s))))
    colnames(alpha.mat) = c("alpha_ovrl", "alpha_s")
  }
  # Set up "observed" data with fitted subgroups #
  param.dat <- PRISM.fit$param.dat
  Subgrps <- PRISM.fit$Subgrps.train
  obs.data <- data.frame(id = 1:nrow(X), Y=Y,A=A, X, Subgrps = Subgrps)
  Xtest.R <- obs.data[,!(colnames(obs.data) %in% c("id", "Y", "A", "Subgrps"))]
  mu_obs <- NULL
  # mu_obs <- PRISM.fit$mu_train
  if (length(PRISM.fit$filter.vars)>0 & filter == "None") {
    obs.data <- obs.data[, colnames(obs.data) %in%
                          c("id", "Y", "A", PRISM.fit$filter.vars, "Subgrps")]
  }
  if (ple=="None") {
    mu_obs <- PRISM.fit$mu_train
  }
  
  # Generate resampling indices #
  n <- dim(X)[1]
  if (stratify=="trt" & !is.null(A)) {
    strata <- A
  }
  if (stratify=="sub") {
    strata <- Subgrps
  }
  if (stratify=="trt_sub" & !is.null(A)) {
    strata <- paste(Subgrps, ": A=", A, sep="")
  }
  if (stratify=="no") {
    strata <- NULL
  }
  if (resample=="Bootstrap") {
    indices <- bootstrap_indices(n=n, R=R, strata=strata)
  }
  if (resample=="Permutation") {
    indices <- permute_indices(n=n, R=R, strata=strata)
    calibrate <- FALSE
  }
  if (resample=="CV") {
    folds <- CV_folds(n=n, R=R, strata=strata)
    calibrate <- FALSE
  }
  
  # Resampling Wrapper #
  fetty_wop <- function(r, obs.data, mu_obs, Xtest.R, family, 
                        filter, ple, submod, param, meta, pool, 
                        delta, propensity, combine,
                        resample_submod, R_submod,
                        alpha_ovrl, alpha_s, filter.hyper,
                        ple.hyper, submod.hyper, resample,
                        calibrate, alpha.mat, verbose) {
    
    if (verbose) message( paste(resample, "Sample", r) )
    
    mu_resamp <- NULL
    if (resample=="Permutation") {
      resamp.data <- obs.data
      A_resamp <- resamp.data$A[ indices[r,] ]
      resamp.data$A <- A_resamp
      Subgrps0 <- resamp.data$Subgrps
    }
    if (resample=="Bootstrap") {
      resamp.data <- obs.data[ indices[r,], ]
      Subgrps0 <- resamp.data$Subgrps
      if (!is.null(mu_obs)) {
        mu_resamp <- mu_obs[ indices[r,], ]
      }
    }
    if (resample=="CV") {
      resamp.data <- obs.data[folds!=r,]
      test <- obs.data[folds==r,]
      Xtest.R <- test[,!(colnames(test) %in% c("Subgrps", "id", "Y", "A"))]
      Subgrps0 <- test$Subgrps
    }
    
    Y.R <- resamp.data$Y
    A.R <- resamp.data$A
    X.R <- resamp.data[!(colnames(resamp.data) %in% c("Subgrps", "id", "Y", "A"))]
    if (family=="survival") { Y.R = Surv(Y.R[,1], Y.R[,2])  }
    res.R <- tryCatch(PRISM_train(Y=Y.R, A=A.R, X=X.R, Xtest=Xtest.R, 
                                  family = family, filter=filter, ple=ple, 
                                  submod = submod, param = param, meta=meta,
                                  pool = pool,
                                  delta = delta, 
                                  propensity = propensity,
                                  combine = combine,
                                  resample_submod = resample_submod, R_submod=R_submod,
                                  alpha_ovrl = alpha_ovrl, alpha_s = alpha_s,
                                  filter.hyper = filter.hyper, 
                                  ple.hyper = ple.hyper, 
                                  submod.hyper = submod.hyper, 
                                  mu_hat = mu_resamp,
                                  verbose = FALSE), 
                      error = function(e) paste("PRISM error:", e))
      
    # Return NULL if param error #
    if (is.character(res.R)) {
      if (verbose) { message("PRISM_train error; ignoring resample")   }
      if (verbose) { message(res.R)   }
      return(NULL)
    }
    # Return NULL if submod.fit is NULL (subgroup model failed) #
    if (is.null(res.R$submod.fit)) {
      if (verbose) { message("submod fit error; ignoring resample")   }
      return(NULL)
    }
    Subgrps.R <- res.R$Subgrps.train
    param.R <- res.R$param.dat
    if (resample!="CV") {
      param.obs.R <- tryCatch(param_est(Y=obs.data$Y, A=obs.data$A, 
                                        X=Xtest.R, param=param,
                                      mu_hat=PRISM.fit$mu_train, 
                                      Subgrps=res.R$Subgrps.test,
                                      alpha_ovrl=alpha_ovrl, alpha_s=alpha_s),
                            error = function(e) "param error")
      bias.R <- param.obs.R %>% rename(est.obs=est)
      bias.R <- bias.R[,colnames(bias.R) %in% c("Subgrps", "estimand", "est.obs")]
      param.R <- suppressWarnings(
        left_join(param.R, bias.R, by=c("Subgrps", "estimand")))
      param.R$bias <- with(param.R, est-est.obs)
    }
    if (resample=="CV") {
      Subgrps.R <- res.R$Subgrps.test
      param.R <- tryCatch(param_est(Y=test$Y, A=test$A, 
                                        X=Xtest.R, param=param,
                                        mu_hat=res.R$mu_test, 
                                        Subgrps=Subgrps.R,
                                        alpha_ovrl=alpha_ovrl, alpha_s=alpha_s),
                              error = function(e) "param error")
      param.R$bias <- NA
    }
    # Return NULL if param error #
    if (is.character(param.R)) {
      if (verbose) { 
        message("param error; ignoring resample")
      }
      return( NULL )
    }
    if (sum(is.na(param.R$est))>0) {
      if (verbose) { 
        message("param error (NAs in estimates); ignoring resample")
      }
      return( NULL )
    }
    calib.dat <- NULL 
    if (calibrate) { 
      calib.dat <- coverage_counter(param.dat=param.R, alpha.mat = alpha.mat)
      calib.dat$ovrl_ind <- ifelse(calib.dat$Subgrps==0, 1, 0)
      calib.dat <- aggregate(ind~ovrl_ind*alpha*estimand, data=calib.dat, FUN="mean")
      calib.dat <- data.frame(R=r, calib.dat)
    }
    
    numb_subs <- length(unique(res.R$Subgrps.train))
    # No Subgroups #
    if (numb_subs==1) {
      hold <- param.R[param.R$Subgrps=="ovrl",
                     c("estimand", "est", "SE", "bias")]
      param.resamp <- param.dat[,c("Subgrps", "N", "estimand")]
      param.resamp <- left_join(param.resamp, hold, by = "estimand")
      param.resamp <- data.frame(R=r, param.resamp)
    }
    # >1 Subgroups #
    if (numb_subs>1) {
      hold <- param.R[param.R$Subgrps!="ovrl",]
      hold <- hold[, colnames(hold) %in% c("Subgrps", "estimand", "est", 
                                          "SE", "bias")]
     
      param.resamp <- NULL
      for (e in unique(hold$estimand)) {
        
        hold.e <- hold[hold$estimand==e,]
        est.resamp <- suppressWarnings(left_join( data.frame(Subgrps=Subgrps.R),
                                hold.e, by="Subgrps"))
        colnames(est.resamp)[1] <- "Subgrps.R"
        est.resamp <- data.frame(Subgrps = Subgrps0, est.resamp)
        est.resamp <- est.resamp %>% group_by(Subgrps, Subgrps.R) %>% mutate( N = n() )
        est.resamp <- unique(est.resamp)
        
        bias <- est.resamp %>% group_by(Subgrps) %>%
                summarise( bias = weighted.mean(bias, w=N) )
        
        param.hold <- NULL
        for (s in unique(est.resamp$Subgrps)) {
          hold.s <- param_combine(est.resamp[est.resamp$Subgrps==s,], 
                                  combine=combine)
          hold.s <- data.frame(Subgrps=s, hold.s[,colnames(hold.s) %in% 
                                                  c("N", "est", "SE")])
          param.hold <- rbind(param.hold, hold.s)
        }
        param.hold <- data.frame(R=r, param.hold, estimand=e)
        param.hold <- suppressWarnings(left_join(param.hold, bias, by="Subgrps"))
        param.resamp <- rbind(param.resamp, param.hold)
      }
      # Add in estimates from overall population #
      hold.ovrl <- param.R[param.R$Subgrps=="ovrl", 
                          colnames(param.R) %in% colnames(param.resamp)]
      param.resamp.ovrl <- data.frame(R=r, hold.ovrl)
      param.resamp <- rbind(param.resamp.ovrl, param.resamp)
    }
    param.resamp$numb_subs <- numb_subs
   
    ## Output counter and estimated statistics ##
    return( list(param.resamp=param.resamp, calib.dat = calib.dat) )
  }

  ### Run Resampling ##
  resamp.obj <- lapply(1:R, fetty_wop, obs.data=obs.data,
                      mu_obs = mu_obs,
                      Xtest.R = Xtest.R, family = family,
                      filter = filter, ple = ple, submod = submod,
                      param = param, meta = meta, pool = pool,
                      delta = delta, propensity = propensity,
                      combine = combine, 
                      resample_submod = resample_submod, R_submod=R_submod,
                      alpha_ovrl = alpha_ovrl, alpha_s = alpha_s,
                      filter.hyper = filter.hyper, ple.hyper = ple.hyper, 
                      submod.hyper = submod.hyper, resample = resample,
                      calibrate=calibrate, alpha.mat = alpha.mat,
                      verbose = verbose)
  
  
  ## Extract Resampling parameter estimates, subject-counters, and calib.dat ##
  hold <- do.call(rbind, resamp.obj)
  resamp_param <- suppressWarnings(do.call(bind_rows, hold[,1]))
  resamp_param <- resamp_param[order(resamp_param$Subgrps, resamp_param$estimand),]
  resamp_calib <- do.call(rbind,hold[,2])
  
  ## Resampling Metrics ##
  # Smoothed estimate (average across resamples), SE, pval #
  # Bootstrap: CIs (percentile, calibrated interval (if applicable))
  # Permutation: p-value
  param.dat <- resamp_metrics(param.dat=param.dat, resamp_param=resamp_param,
                              resamp_calib=resamp_calib, resample=resample)
  
  ### Return Outputs ###
  return( list(param.dat=param.dat, resamp.dist = resamp_param, 
               resamp.calib=resamp_calib) )
}