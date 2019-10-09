#' PRISM (Resample): Patient Response Identifier for Stratified Medicine
#'
#' Based on initial PRISM fit (\code{PRISM_train}), run resampling (Boostrap, Permutation,
#' or cross-validation). Used directly in \code{PRISM}.
#'
#' @param PRISM.fit Fitted PRISM model
#' @inheritParams PRISM
#' 
#' @return Trained PRISM object. Includes filter, ple, submod, and param outputs.
#'  \itemize{
#'   \item param.dat - Parameter estimates and variablity metrics (depends on param)
#'   \item resamp.dist - - Resampling distributions
#' }
#' 
#' @export
#'   
##### PRISM: Patient Responder Identifiers for Stratified Medicine ########
PRISM_resamp <- function(PRISM.fit, Y, A, X, Xtest=NULL, family="gaussian",
                       filter="filter_glmnet", ple=NULL, submod=NULL, param=NULL,
                       alpha_ovrl=0.05, alpha_s = 0.05,
                       filter.hyper=NULL, ple.hyper=NULL, submod.hyper = NULL,
                       param.hyper = NULL, verbose=TRUE,
                       prefilter_resamp=FALSE, resample="Bootstrap", 
                       R=NULL, stratify = TRUE, calibrate=TRUE,
                       alpha.mat=NULL){
  if (is.null(alpha.mat)){
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
  if (length(PRISM.fit$filter.vars)>0 & prefilter_resamp == TRUE){
    obs.data <- obs.data[, colnames(obs.data) %in%
                          c("id", "Y", "A", PRISM.fit$filter.vars, "Subgrps")]
  }
  
  # Generate resampling indices #
  n <- dim(X)[1]
  if (stratify){
    strata <- Subgrps
  }
  if (!stratify){
    strata <- NULL
  }
  if (resample=="Bootstrap"){
    indices <- bootstrap_indices(n=n, R=R, strata=strata)
  }
  if (resample=="Permutation"){
    indices <- permute_indices(n=n, R=R, strata=strata)
    calibrate=FALSE
  }
  if (resample=="CV"){
    folds <- CV_folds(n=n, R=R, strata=strata)
    calibrate=FALSE
  }
  
  ### Resampling Wrapper ###
  fetty_wop <- function(r, stratify, obs.data, Xtest.R, ple, filter, submod,
                       calibrate, alpha.mat, verbose){
    
    if (verbose) message( paste(resample, "Sample", r) )
    ### Permutation resampling (shuffle treatment assignment) ###
    if (resample=="Permutation"){
      resamp.data <- obs.data
      A_resamp <- resamp.data$A[ indices[r,] ]
      resamp.data$A <- A_resamp
      Subgrps0 <- resamp.data$Subgrps
    }
    if (resample=="Bootstrap"){
      resamp.data <- obs.data[ indices[r,], ]
      Subgrps0 <- resamp.data$Subgrps
    }
    if (resample=="CV"){
      resamp.data <- obs.data[folds!=r,]
      test <- obs.data[folds==r,]
      Xtest.R <- test[,!(colnames(test) %in% c("Subgrps", "id", "Y", "A"))]
      Subgrps0 <- test$Subgrps
    }
    ### PRISM on resampled data ###
    Y.R <- resamp.data$Y
    A.R <- resamp.data$A
    X.R <- resamp.data[!(colnames(resamp.data) %in% c("Subgrps", "id", "Y", "A"))]
    if (family=="survival"){ Y.R = Surv(Y.R[,1], Y.R[,2])  }
    res.R <- PRISM_train(Y=Y.R, A=A.R, X=X.R, Xtest=Xtest.R, family=family, 
                        filter=filter, ple=ple, submod = submod, param=param,
                        alpha_ovrl = alpha_ovrl, alpha_s = alpha_s,
                        filter.hyper = filter.hyper, ple.hyper = ple.hyper,
                        submod.hyper = submod.hyper, param.hyper = param.hyper,
                        verbose = FALSE)
    Subgrps.R <- res.R$Subgrps.train
    param.R <- res.R$param.dat
    if (resample!="CV"){
      param.obs.R <- tryCatch( do.call( param, list(Y=obs.data$Y, A=obs.data$A, 
                                               X=Xtest.R, 
                                               mu_hat = PRISM.fit$mu_train,
                                               Subgrps=res.R$Subgrps.test,
                                               alpha_ovrl=alpha_ovrl,
                                               alpha_s=alpha_s)  ),
                          error = function(e) "param error" )
      bias.R <- param.obs.R %>% rename(est.obs=est)
      bias.R <- bias.R[,colnames(bias.R) %in% c("Subgrps", "estimand", "est.obs")]
      param.R <- left_join(param.R, bias.R, by=c("Subgrps", "estimand"))
      param.R$bias <- with(param.R, est-est.obs)
    }
    if (resample=="CV"){
      Subgrps.R = res.R$Subgrps.test
      param.R = tryCatch( do.call( param, list(Y=test$Y, A=test$A, X=Xtest.R, 
                                               mu_hat = res.R$mu_test,
                                               Subgrps=Subgrps.R,
                                               alpha_ovrl=alpha_ovrl,
                                               alpha_s=alpha_s)  ),
                          error = function(e) "param error" )
      param.R$bias = NA
    }
    # Return NULL if param error #
    if (is.character(param.R)){
      if (verbose){ message("param error; ignoring resample")   }
      return( NULL   )
    }
    calib.dat = NULL # Initiatilize
    if (calibrate){ # Calibration of alpha #
      # Loop across alpha values and check if est.obs in [LCL,UCL] (by estimand) #
      calib.dat = coverage_counter(param.dat=param.R, alpha.mat = alpha.mat)
      calib.dat$ovrl_ind = ifelse(calib.dat$Subgrps==0, 1, 0)
      calib.dat = aggregate(ind~ovrl_ind*alpha*estimand, data=calib.dat, FUN="mean")
      calib.dat = data.frame(R=r, calib.dat)
    }
    # Bootstrap parameter estimates (for original Subgrps) #
    numb_subs = length(unique(param.R$Subgrps))
    # No Subgroups #
    if (numb_subs==1){
      hold = param.R[,c("estimand", "est", "SE", "bias")]
      param.resamp = param.dat[,c("Subgrps", "N", "estimand")]
      param.resamp = left_join(param.resamp, hold, by = "estimand")
      param.resamp = data.frame(R=r, param.resamp)
    }
    # >1 Subgroups #
    if (numb_subs>1){
      hold = param.R %>% filter(param.R$Subgrps>0)
      hold = hold[, colnames(hold) %in% c("Subgrps", "estimand", "est", "SE", "bias") ]
      # Loop through estimands #
      param.resamp <- NULL
      for (e in unique(hold$estimand)){
        hold.e <- hold[hold$estimand==e,]
        est.resamp <- left_join( data.frame(Subgrps=Subgrps.R),
                                hold.e, by="Subgrps")
        est.resamp <- data.frame(Subgrps = Subgrps0, est.resamp)
        colnames(est.resamp) <- c("Subgrps", "Subgrps.R", "estimand", "est", "SE",
                                 "bias")
        est.resamp <- est.resamp %>% group_by(Subgrps, Subgrps.R) %>% mutate( N = n() )
        est.resamp <- unique(est.resamp)
        # Bias #
        bias <- est.resamp %>% group_by(Subgrps) %>%
                summarise( bias = weighted.mean(bias, w=N) )
        ## Obtain point-estimate / standard errors ##
        S_levels <- as.numeric( names(table(est.resamp$Subgrps)) )
        param.hold <- NULL
        for (s in S_levels){
          hold.s <- param_combine(est.resamp[est.resamp$Subgrps==s,], combine="SS")
          hold.s <- data.frame(Subgrps=s,hold.s[,colnames(hold.s) %in% c("N", "est", "SE")])
          param.hold <- rbind(param.hold, hold.s)
        }
        param.hold <- data.frame(R=r, param.hold, estimand=e)
        param.hold <- left_join(param.hold, bias, by="Subgrps")
        param.resamp <- rbind(param.resamp, param.hold)
      }
      # Add in estimates from overall population #
      hold.ovrl = param.R[param.R$Subgrps==0, 
                          colnames(param.R) %in% colnames(param.resamp)]
      param.resamp.ovrl = data.frame(R=r, hold.ovrl)
      param.resamp = rbind( param.resamp.ovrl, param.resamp )
    }
    param.resamp$numb_subs = numb_subs
   
    ## Output counter and estimated statistics ##
    return( list(param.resamp=param.resamp, calib.dat = calib.dat) )
  }

  ### Run Resampling ##
  resamp.obj = lapply(1:R, fetty_wop, stratify=stratify, obs.data=obs.data,
                      Xtest.R = Xtest.R,
                      ple=ple, filter=filter, submod=submod, calibrate=calibrate,
                      alpha.mat = alpha.mat,
                      verbose = verbose)
  ## Extract Resampling parameter estimates, subject-counters, and calib.dat ##
  hold = do.call(rbind, resamp.obj)
  resamp_param = do.call(bind_rows, hold[,1])
  resamp_param = resamp_param[order(resamp_param$Subgrps, resamp_param$estimand),]
  resamp_calib = do.call(rbind,hold[,2])
  
  ## Resampling Metrics ##
  # Smoothed estimate (average across resamples), SE, pval #
  # Bootstrap: CIs (percentile, calibrated interval (if applicable))
  # Permutation: p-value
  param.dat = resamp_metrics(param.dat=param.dat, resamp_param=resamp_param,
                              resamp_calib=resamp_calib, resample=resample)
  
  ### Return Outputs ###
  return( list(param.dat=param.dat, resamp.dist = resamp_param, 
               resamp.calib=resamp_calib) )
}