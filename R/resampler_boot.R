### Bootstrap Resampler: Subgroups ###
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
resampler_boot <- function(Y, A, X, fit, wrapper, list_args,
                           R=20, stratify="trt",
                           fixed_subs = "false", 
                           verbose=FALSE, ...) {
  
  wrapper_fn <- get(wrapper, envir = parent.frame())
  
  # Generate resampling indices #
  n <- dim(X)[1]
  strata <- NULL
  if (stratify=="trt") {
    strata <- A
  }
  indices <- bootstrap_indices(n=n, R=R, strata=strata)
  
  # Initial Subgroups #
  Subgrps <- fit$Subgrps
  
  # Original Data #
  trt_eff <- fit$trt_eff
  obs.data <- data.frame(id = 1:nrow(X), Y=Y, A=A, X, Subgrps = Subgrps)
  Xtest.R <- obs.data[,!(colnames(obs.data) %in% c("id", "Y", "A", "Subgrps"))]
  fam_check <- ifelse(is.Surv(Y), "survival", "non-survival")
  
  # How to combine trt estimates #
  if (is.null(list_args$combine)) {
    combine_use <- "SS"
  }
  if (!is.null(list_args$combine)) {
    combine_use <- list_args$combine
  }
  
  print(list_args)
  # Loop through resamples #
  looper <- function(r, verbose=FALSE, fixed_subs, combine,
                     wrapper, list_args) {
    
    if (verbose) message( paste("Resample", r) )
    resamp.data <- obs.data[ indices[r,], ]
    Subgrps0 <- resamp.data$Subgrps
    
    # Set up data 
    Y.R <- resamp.data$Y
    A.R <- resamp.data$A
    X.R <- resamp.data[!(colnames(resamp.data) %in% c("Subgrps", "id", "Y", "A"))]
    if (fam_check=="survival") { Y.R = Surv(Y.R[,1], Y.R[,2])  }
    
    # Fit the wrapper #
    fit.r <- tryCatch(do.call(wrapper, 
                              append(list(Y=Y.R, A=A.R, X=X.R, mu_train=NULL),
                                     list_args)),
                      error = function(e) as.character(e))

    if (is.character(fit.r)) {
      if (verbose) { message("fit error; ignoring resample")   }
      return(NULL)
    }
    Subgrps.R <- fit.r$Subgrps
    trt_eff.R <- fit.r$trt_eff
   
    # Return NULL if param error #
    if (is.character(trt_eff.R)) {
      if (verbose) { message("trt estimate error; ignoring resample")   }
      return(NULL)
    }
    # Option 1: Original Subgroups are the "same" as bootstrap #
    # For example: groups are "Trt with A" vs "Trt with B" #
    if (fixed_subs=="true") {
      boot_trt_eff <- data.frame(R=r, trt_eff.R)
      return(boot_trt_eff)
    }
    # Option 2 (default): Subgroups can vary by bootstrap => map back to original #
    # GUIDE approach #
    if (fixed_subs=="false") {
      
      numb_subs = length(unique(trt_eff.R$Subgrps))
      
      # No Subgroups #
      if (numb_subs==1) {
        hold <- trt_eff.R[,c("estimand", "est", "SE")]
        boot_trt_eff <- trt_eff[,c("Subgrps", "N", "estimand")]
        boot_trt_eff <- dplyr::left_join(boot_trt_eff, hold, by = "estimand")
        boot_trt_eff <- data.frame(R=r, boot_trt_eff)
      }
      # >1 Subgroups #
      if (numb_subs>1) {
        hold <- trt_eff.R[trt_eff.R$Subgrps!="ovrl",]
        hold <- hold[, colnames(hold) %in% c("Subgrps", "estimand", "est", "SE")]
        
        param.resamp <- NULL
        for (e in unique(hold$estimand)) {
          
          hold.e <- hold[hold$estimand==e,]
          est.resamp <- suppressWarnings(left_join( data.frame(Subgrps=Subgrps.R),
                                                    hold.e, by="Subgrps"))
          colnames(est.resamp)[1] <- "Subgrps.R"
          est.resamp <- data.frame(Subgrps = Subgrps0, est.resamp)
          est.resamp <- est.resamp %>% group_by(Subgrps, Subgrps.R) %>% mutate( N = n() )
          est.resamp <- unique(est.resamp)
          
          param.hold <- NULL
          for (s in unique(est.resamp$Subgrps)) {
            hold.s <- param_combine(est.resamp[est.resamp$Subgrps==s,], 
                                    combine=combine)
            hold.s <- data.frame(Subgrps=s, hold.s[,colnames(hold.s) %in% 
                                                     c("N", "est", "SE")])
            param.hold <- rbind(param.hold, hold.s)
          }
          param.hold <- data.frame(R=r, param.hold, estimand=e)
          # param.hold <- suppressWarnings(left_join(param.hold, bias, by="Subgrps"))
          param.resamp <- rbind(param.resamp, param.hold)
        }
        # Add in estimates from overall population #
        hold.ovrl <- trt_eff.R[trt_eff.R$Subgrps=="ovrl", 
                             colnames(trt_eff.R) %in% colnames(param.resamp)]
        param.resamp.ovrl <- data.frame(R=r, hold.ovrl)
        param.resamp <- rbind(param.resamp.ovrl, param.resamp)
        boot_trt_eff <- param.resamp
      }
      boot_trt_eff$numb_subs = numb_subs
      return(boot_trt_eff)
    }
  }
  # Loop through resamples #
  resamp.obj = lapply(1:R, looper, wrapper=wrapper_fn, list_args=list_args,
                      verbose = TRUE, fixed_subs=fixed_subs,
                      combine=combine_use)
  
  if (any(sapply(resamp.obj, is.null))) {
    resamp.obj <- resamp.obj[-which(sapply(resamp.obj, is.null))] # remove NULL elements
  }
  boot_dist <- data.frame(do.call(rbind, resamp.obj))
  
  # Calculate bootstrap-based treatment estimates #
  boot_trt_eff <- NULL
  for (e in unique(boot_dist$estimand)) {
    for (sub in unique(boot_dist$Subgrps)) {
      
      hold <- boot_dist[boot_dist$Subgrps==sub & boot_dist$estimand==e,]
      sub_trt <- trt_eff[trt_eff$Subgrps==sub & trt_eff$estimand==e,] 
      alpha_use <- sub_trt$alpha
      boot_est <- mean(hold$est, na.rm=TRUE)
      boot_SE <- sd(hold$est, na.rm=TRUE)
      quants = as.numeric(
        quantile(hold$est, probs=c(alpha_use/2, (1-alpha_use/2)), na.rm = TRUE) )
      pct_LCL <- quants[1]
      pct_UCL <- quants[2]
      summ <- data.frame(Subgrps = sub, estimand=e, est=boot_est, SE = boot_SE,
                         LCL = pct_LCL, UCL = pct_UCL, alpha=alpha_use)
      boot_trt_eff <- rbind(boot_trt_eff, summ)
    }
  }
  # Merge back with original treatment group #
  boot_trt_eff <- left_join(trt_eff[,c("Subgrps", "N", "estimand")],
                            boot_trt_eff, by=c("Subgrps", "estimand"))
  return(list(boot_dist=boot_dist, boot_trt_eff=boot_trt_eff))
}