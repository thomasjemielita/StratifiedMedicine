### General Resampler: Subgroups ###
# Inputs #
# @Y: Outcome (if survival, should be Surv(time, event))
# @A: Binary Treatment (A=0 vs A=1)
# @X: Covariate Space
# @fit: Fitted model (from wrapper using observed data)
# @wrapper: Wrapper name
# @trteff_method: Approach for estimating treatment effects (ex: cox, ols, riskdiff)
# @R: Number of resamples
# @stratify: Should resampling be stratified by treatment?
# @fixed_subs: Are the subgroups fixed (ex: dopt=1 vs dopt=0). 
# If trees, use fixed_subs=="false"
# @verbose: Whether to print output
# @alpha: Alpha for CIs
resampler_general <- function(Y, A, X, fit, wrapper, list_args,
                           R=20, resample = "Bootstrap", 
                           stratify=ifelse(!is.null(A), "trt", "no"),
                           fixed_subs = "false",
                           calibrate=FALSE,
                           alpha.mat=NULL,
                           verbose=FALSE, ...) {
  
  wrapper_fn <- get(wrapper, envir = parent.frame())
  
  alpha_ovrl <- list_args$alpha_ovrl
  alpha_s <- list_args$alpha_s
  
  # Calibration alpha matrix? #
  if (is.null(alpha.mat)) {
    alpha.mat = data.frame(cbind( unique(c(seq(alpha_ovrl/1000, alpha_ovrl,
                                               by=0.005),alpha_ovrl) ),
                                  unique(c(seq(alpha_s/1000, alpha_s,
                                               by=0.005),alpha_s))))
    colnames(alpha.mat) = c("alpha_ovrl", "alpha_s")
  }
  
  # Generate resampling indices #
  n <- dim(X)[1]
  strata <- NULL
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
  
  # Initial Subgroups / mu_train #
  Subgrps <- fit$Subgrps
  mu_train <- fit$mu_train
  mu_ind <- !is.null(mu_train)
  mu_vec <- is.vector(mu_train)
  
  # Original Data #
  trt_eff <- fit$trt_eff
  id.data <- data.frame(id = 1:nrow(X))
  if (is.null(A)) {
    obs.data <- data.frame(id.data, Y=Y, X, Subgrps = Subgrps)
  }
  if (!is.null(A)) {
    obs.data <- data.frame(id.data, Y=Y, A=A, X, Subgrps = Subgrps)
  }
  Xtest.R <- obs.data[,!(colnames(obs.data) %in% c("id", "Y", "A", "Subgrps"))]
  fam_check <- ifelse(survival::is.Surv(Y), "survival", "non-survival")
  
  # Filter? #
  if (!is.null(list_args$filter)) {
    if (length(fit$filter.vars)>0 & list_args$filter == "None") {
    obs.data <- obs.data[, colnames(obs.data) %in%
                           c("id", "Y", "A", fit$filter.vars, "Subgrps")]
    }
  }
  
  
  # How to combine trt estimates #
  if (is.null(list_args$combine)) {
    combine_use <- "SS"
  }
  if (!is.null(list_args$combine)) {
    combine_use <- list_args$combine
  }
  list_args$verbose <- FALSE
  
  # Loop through resamples #
  looper <- function(r, verbose=FALSE, resample, fixed_subs, combine,
                     Xtest.R, 
                     wrapper, list_args) {
    
    
    if (verbose) message( paste(resample, "Sample", r) )
    mu_train_r <- NULL
    if (resample=="Permutation") {
      resamp.data <- obs.data
      A_resamp <- resamp.data$A[ indices[r,] ]
      resamp.data$A <- A_resamp
      Subgrps0 <- resamp.data$Subgrps
    }
    if (resample=="Bootstrap") {
      resamp.data <- obs.data[ indices[r,], ]
      Subgrps0 <- resamp.data$Subgrps
      if (mu_ind) {
        if (mu_vec) {
          mu_train_r <- mu_train[indices[r,]]
        }
        if (!mu_vec) {
          mu_train_r <- mu_train[indices[r,], ]
        }
      }
    }
    if (resample=="CV") {
      resamp.data <- obs.data[folds!=r,]
      test <- obs.data[folds==r,]
      Xtest.R <- test[,!(colnames(test) %in% c("Subgrps", "id", "Y", "A"))]
      Subgrps0 <- test$Subgrps
    }
    
    # Set up data 
    Y.R <- resamp.data$Y
    A.R <- resamp.data$A
    X.R <- resamp.data[!(colnames(resamp.data) %in% c("Subgrps", "id", "Y", "A"))]
    if (fam_check=="survival") { Y.R = Surv(Y.R[,1], Y.R[,2])  }
    
    # Fit the wrapper #
    fit.r <- tryCatch(do.call(wrapper, 
                              append(list(Y=Y.R, A=A.R, X=X.R, 
                                          Xtest=Xtest.R, mu_train=mu_train_r),
                                     list_args)),
                      error = function(e) as.character(e))
    # print(paste("resample", r))
    # print(fit.r)

    if (is.character(fit.r)) {
      if (verbose) { message("fit error; ignoring resample")   }
      return(NULL)
    }
    Subgrps.R <- fit.r$Subgrps
    trt_eff.R <- fit.r$trt_eff
    
    sub_dat <- unique(data.frame(id = resamp.data$id, Subgrps = Subgrps.R))
    sub_dat <- dplyr::left_join(id.data, sub_dat, by="id")
    Subgrp_vec <- sub_dat$Subgrps
    
    if (is.character(trt_eff.R)) {
      if (verbose) { message("trt estimate error; ignoring resample")   }
      return(NULL)
    }
    

    if (inherits(fit.r, "PRISM_train")) {
      if (is.null(fit.r$submod.fit)) {
        if (verbose) { message("submod fit error; ignoring resample")   }
        return(NULL)
      }
    }
    
    # Bias estimates #
    if (resample!="CV") {
      trt_eff_obs.R <- tryCatch(param_est(Y=obs.data$Y, A=obs.data$A, 
                                        X=Xtest.R, param=list_args$param,
                                        mu_hat=fit$mu_train, 
                                        Subgrps=fit.r$Subgrps.test,
                                        alpha_ovrl=list_args$alpha_ovrl, 
                                        alpha_s=list_args$alpha_s),
                              error = function(e) "param error")
      if (is.character(trt_eff_obs.R) | is.null(trt_eff_obs.R)) {
        trt_eff.R$bias <- NA
        trt_eff.R$est.obs <- NA
      }
      else {
        bias.R <- trt_eff_obs.R %>% rename(est.obs=est)
        bias.R <- bias.R[,colnames(bias.R) %in% c("Subgrps", "estimand", "est.obs")]
        trt_eff.R <- suppressWarnings(
          dplyr::left_join(trt_eff.R, bias.R, by=c("Subgrps", "estimand")))
        trt_eff.R$bias <- with(trt_eff.R, est-est.obs) 
      }
    }  
    if (resample=="CV") {
      Subgrps.R <- fit.r$Subgrps.test
      trt_eff.R <- tryCatch(param_est(Y=test$Y, A=test$A, 
                                    X=Xtest.R, param=list_args$param,
                                    mu_hat=fit.r$mu_test, 
                                    Subgrps=Subgrps.R,
                                    alpha_ovrl=list_args$alpha_ovrl, 
                                    alpha_s=list_args$alpha_s),
                          error = function(e) "param error")
      trt_eff.R$bias <- NA
    }
    if (sum(is.na(trt_eff.R$est))>0) {
      if (verbose) { 
        message("param error (NAs in estimates); ignoring resample")
      }
      return( NULL )
    }
    calib.dat <- NULL 
    if (calibrate) { 
      calib.dat <- coverage_counter(param.dat=trt_eff.R, alpha.mat = alpha.mat)
      calib.dat$ovrl_ind <- ifelse(calib.dat$Subgrps==0, 1, 0)
      calib.dat <- aggregate(ind~ovrl_ind*alpha*estimand, data=calib.dat, FUN="mean")
      calib.dat <- data.frame(R=r, calib.dat)
    }
    
    
    # Option 1: Original Subgroups are the "same" as bootstrap #
    # For example: groups are "Trt with A" vs "Trt with B" #
    if (fixed_subs=="true") {
      est_resamp <- data.frame(R=r, trt_eff.R)
      return(list(Subgrp_vec=Subgrp_vec, est_resamp=est_resamp, calib.dat = calib.dat))
    }
    # Option 2 (default): Subgroups can vary by bootstrap => map back to original #
    # GUIDE approach #
    if (fixed_subs=="false") {
      
      numb_subs = length(unique(trt_eff.R$Subgrps[trt_eff.R$Subgrps!="ovrl"]))
      
      # No Subgroups #
      if (numb_subs==1) {
        hold <- trt_eff.R[,c("estimand", "est", "SE", "bias")]
        est_resamp <- trt_eff[,c("Subgrps", "N", "estimand")]
        est_resamp <- dplyr::left_join(est_resamp, hold, by = "estimand")
        est_resamp <- data.frame(R=r, est_resamp)
      }
      # >1 Subgroups #
      if (numb_subs>1) {
        hold <- trt_eff.R[trt_eff.R$Subgrps!="ovrl",]
        hold <- hold[, colnames(hold) %in% c("Subgrps", "estimand", "est", 
                                             "SE", "bias")]
        
        est_resamp <- NULL
        for (e in unique(hold$estimand)) {
          
          hold.e <- hold[hold$estimand==e,]
          est.resamp <- suppressWarnings(dplyr::left_join( data.frame(Subgrps=Subgrps.R),
                                                           hold.e, by="Subgrps"))
          colnames(est.resamp)[1] <- "Subgrps.R"
          est.resamp <- data.frame(Subgrps = Subgrps0, est.resamp)
          est.resamp <- est.resamp %>% dplyr::group_by(Subgrps, Subgrps.R) %>% mutate( N = n() )
          est.resamp <- unique(est.resamp)
          
          bias <- est.resamp %>% group_by(Subgrps) %>%
            summarise( bias = weighted.mean(bias, w=N), .groups = 'drop')
          
          param.hold <- NULL
          for (s in unique(est.resamp$Subgrps)) {
            hold.s <- param_combine(est.resamp[est.resamp$Subgrps==s,], 
                                    combine=combine)
            hold.s <- data.frame(Subgrps=s, hold.s[,colnames(hold.s) %in% 
                                                     c("N", "est", "SE")])
            param.hold <- rbind(param.hold, hold.s)
          }
          param.hold <- data.frame(R=r, param.hold, estimand=e)
          param.hold <- suppressWarnings(dplyr::left_join(param.hold, bias, by="Subgrps"))
          est_resamp <- rbind(est_resamp, param.hold)
        }
        # Add in estimates from overall population #
        hold.ovrl <- trt_eff.R[trt_eff.R$Subgrps=="ovrl", 
                               colnames(trt_eff.R) %in% colnames(est_resamp)]
        est_resamp.ovrl <- data.frame(R=r, hold.ovrl)
        est_resamp <- rbind(est_resamp.ovrl, est_resamp)
      }
      est_resamp$numb_subs = numb_subs
      return(list(Subgrp_vec=Subgrp_vec, est_resamp=est_resamp, calib.dat = calib.dat))
    }
  }
  # Loop through resamples #
  resamp.obj = lapply(1:R, looper, wrapper=wrapper_fn, resample=resample,
                      list_args=list_args,
                      Xtest.R = Xtest.R, 
                      verbose = verbose, fixed_subs=fixed_subs,
                      combine=combine_use)
  
  if (any(sapply(resamp.obj, is.null))) {
    resamp.obj <- resamp.obj[-which(sapply(resamp.obj, is.null))] # remove NULL elements
  }
  hold <- do.call(rbind, resamp.obj)
  resamp_subgrps <- data.frame(id = obs.data$id, do.call(cbind, hold[,1]))
  resamp_dist <- data.frame(do.call(rbind, hold[,2]))
  resamp_calib <- do.call(rbind, hold[,3])
  
  ## Resampling Metrics ##
  # Smoothed estimate (average across resamples), SE, pval #
  # Bootstrap: CIs (percentile, calibrated interval (if applicable))
  # Permutation: p-value
  trt_eff_resamp <- resamp_metrics(trt_eff=trt_eff, resamp_dist=resamp_dist,
                              resamp_calib=resamp_calib, resample=resample)
  
  # Merge back with original treatment group #
  trt_eff_resamp <- dplyr::left_join(trt_eff[,c("Subgrps", "N", "estimand")],
                                 trt_eff_resamp, by=c("Subgrps", "estimand"))
  return(list(resamp_dist = resamp_dist, trt_eff_resamp=trt_eff_resamp, 
              resamp_subgrps = resamp_subgrps))
}