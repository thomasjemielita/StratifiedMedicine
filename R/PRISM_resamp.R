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
PRISM_resamp = function(PRISM.fit, Y, A, X, Xtest=NULL, family="gaussian",
                       filter="filter_glmnet", ple=NULL, submod=NULL, param=NULL,
                       alpha_ovrl=0.05, alpha_s = 0.05,
                       filter.hyper=NULL, ple.hyper=NULL, submod.hyper = NULL,
                       param.hyper = NULL, verbose=TRUE,
                       prefilter_resamp=FALSE, resample=NULL, R=100, stratify = TRUE){
  
  # Set up "observed" data with fitted subgroups #
  param.dat = PRISM.fit$param.dat
  Subgrps = PRISM.fit$Subgrps.train
  obs.data = data.frame(id = 1:nrow(X), Y=Y,A=A, X, Subgrps = Subgrps)
  if (length(PRISM.fit$filter.vars)>0 & prefilter_resamp == TRUE){
    obs.data = obs.data[, colnames(obs.data) %in%
                          c("id", "Y", "A", PRISM.fit$filter.vars, "Subgrps")]
  }
  subject.counter = data.frame(id = obs.data$id )
  
  # Generate resampling indices #
  n = dim(X)[1]
  if (resample=="Bootstrap" | resample=="Permutation"){
    if (!stratify){
      indices <- sample.int(n, n * R, replace = TRUE)
      dim(indices) <- c(R, n)
    }
    if (stratify){
      indices <- matrix(NA, R, n)
      for (s in unique(Subgrps)) {
        sub_i <- seq_len(n)[Subgrps == s]
        ns = length(sub_i)
        indices[, sub_i] <- sub_i[sample.int(ns, ns*R, replace = TRUE)]
      }
    }
  }
  if (resample=="CV"){
    if (!stratify){
      folds <- sample(rep(1:R,ceiling(n/R))[1:n])
    }
    if (stratify){
      folds = NULL
      for (s in unique(Subgrps)){
        sub_i <- seq_len(n)[Subgrps == s]
        ns = length(sub_i)
        folds.s <- sample(rep(1:R,ceiling(ns/R))[1:ns])
        folds[sub_i] = folds.s
      }
    }
  }
  
  ### Resampling Wrapper ###
  fetty_wop = function(R, stratify, obs.data, ple, filter, submod,
                       calibrate, verbose){
    
    if (verbose) message( paste(resample, "Sample", R) )
    resamp.data = obs.data
    ### Permutation resampling (shuffle treatment assignment) ###
    if (resample=="Permutation"){
      A_resamp = resamp.data$A[ indices[R,] ]
      resamp.data$A = A_resamp
    }
    if (resample=="Bootstrap"){
      resamp.data = obs.data[ indices[R,], ]
      Subgrps0 = resamp.data$Subgrps
    }
    if (resample=="CV"){
      resamp.data = obs.data[folds!=R,]
      test = obs.data[folds==R,]
      Xtest = test[,!(colnames(test) %in% c("Subgrps", "id", "Y", "A"))]
      Subgrps0 = test$Subgrps
    }
    ### PRISM on resampled data ###
    Y.R = resamp.data$Y
    A.R = resamp.data$A
    X.R = resamp.data[!(colnames(resamp.data) %in% c("Subgrps", "id", "Y", "A"))]
    if (family=="survival"){ Y.R = Surv(Y.R[,1], Y.R[,2])  }
    res.R = PRISM_train(Y=Y.R, A=A.R, X=X.R, Xtest=Xtest, family=family, 
                        filter=filter, ple=ple, submod = submod, param=param,
                        alpha_ovrl = alpha_ovrl, alpha_s = alpha_s,
                        filter.hyper = filter.hyper, ple.hyper = ple.hyper,
                        submod.hyper = submod.hyper, param.hyper = param.hyper,
                        verbose = FALSE)
    Subgrps.R = res.R$Subgrps.train
    param.R = res.R$param.dat
    if (resample=="CV"){
      Subgrps.R = res.R$Subgrps.test
      param.R = tryCatch( do.call( param, list(Y=test$Y, A=test$A, X=Xtest, 
                                               mu_hat = res.R$mu_test,
                                               Subgrps=Subgrps.R,
                                               alpha_ovrl=alpha_ovrl,
                                               alpha_s=alpha_s)  ),
                          error = function(e) "param error" )
    }
    # Return NULL if param error #
    if (is.character(param.R)){
      if (verbose){ message("param error; ignoring resample")   }
      return( NULL   )
    }
    if (calibrate){
      ### Calibrate alpha (need to implement) ###
      ## NEED WORKS!!! ##
    }
    # Bootstrap parameter estimates (for original Subgrps) #
    numb_subs = length(unique(param.R$Subgrps))
    # No Subgroups #
    if (numb_subs==1){
      hold = param.R[,c("estimand", "est", "SE")]
      param.resamp = param.dat[,c("Subgrps", "N", "estimand")]
      param.resamp = left_join(param.resamp, hold, by = "estimand")
      param.resamp = data.frame(R=R, param.resamp)
    }
    # >1 Subgroups #
    if (numb_subs>1){
      hold = param.R %>% filter(param.R$Subgrps>0)
      hold = hold[, colnames(hold) %in% c("Subgrps", "estimand", "est", "SE") ]
      # Loop through estimands #
      param.resamp = NULL
      for (e in unique(hold$estimand)){
        hold.e = hold[hold$estimand==e,]
        est.resamp = left_join( data.frame(Subgrps=Subgrps.R),
                                hold.e, by="Subgrps")
        est.resamp = data.frame(Subgrps = Subgrps0, est.resamp)
        colnames(est.resamp) = c("Subgrps", "Subgrps.R", "estimand", "est", "SE")
        est.resamp = est.resamp %>% group_by(Subgrps, Subgrps.R) %>% mutate( N = n() )
        est.resamp = unique(est.resamp)
        ## Obtain point-estimate / standard errors ##
        S_levels = as.numeric( names(table(est.resamp$Subgrps)) )
        param.hold = NULL
        for (s in S_levels){
          hold.s = param_combine(est.resamp[est.resamp$Subgrps==s,], combine="SS")
          hold.s = data.frame(Subgrps=s,hold.s[,colnames(hold.s) %in% c("N", "est", "SE")])
          param.hold = rbind(param.hold, hold.s)
        }
        param.hold = data.frame(R=R, param.hold, estimand=e)
        param.resamp = rbind(param.resamp, param.hold)
      }
      # Add in estimates from overall population #
      hold.ovrl = param.R[param.R$Subgrps==0, 
                          colnames(param.R) %in% colnames(param.resamp)]
      param.resamp.ovrl = data.frame(R=R, hold.ovrl)
      param.resamp = rbind( param.resamp.ovrl, param.resamp )
    }
    param.resamp$numb_subs = numb_subs
   
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
                      ple=ple, filter=filter, submod=submod, calibrate=FALSE,
                      verbose = verbose)
  ## Extract Resampling parameter estimates and subject-counters ##
  hold = do.call(rbind, resamp.obj)
  resamp_param = do.call(bind_rows, hold[,1])
  resamp_param = resamp_param[order(resamp_param$Subgrps, resamp_param$estimand),]
  resamp_counter = do.call(cbind, hold[,2])
  
  ## Resampling Metrics ##
  # Smoothed estimate (average across resamples), SE, pval #
  # Bootstrap: CIs (percentile, BCa)
  # Permutation: p-value
  resamp_metrics = function(param.dat, e){
    final_ests = param.dat[param.dat$estimand==e,]
    final_ests$est_resamp = NA
    final_ests$SE_resamp = NA
    if (resample=="Permutation"){
      final_ests$pval_perm = NA
    }
    if (resample=="Bootstrap"){
      final_ests = final_ests %>% mutate( SE_smooth = NA, SE_smooth_adj = NA,
                                          bias_smooth = NA, bca.z0 = NA,
                                          accel = NA, LCL.pct = NA, UCL.pct= NA,
                                          LCL.BCa = NA, UCL.BCa = NA)  
    }
    if (resample=="CV"){
      final_ests = final_ests %>% mutate( LCL.CV = NA, UCL.CV = NA, pval.CV = NA)  
    }
    for (sub in unique(final_ests$Subgrps)){
      if (sub<=0){ alpha = alpha_ovrl }
      if (sub>0 ){ alpha = alpha_s }
      hold = final_ests[final_ests$Subgrps==sub,]
      hold.R = resamp_param[resamp_param$Subgrps==sub & resamp_param$estimand==e,]
      ### Smoothed Bootstrap estimate, boot SD ###
      est0 = hold$est
      est.vec = hold.R$est
      est.R = mean(hold.R$est)
      final_ests$est_resamp[final_ests$Subgrps==sub] = est.R
      # if (resample=="CV"){ ## TBD, need to sort our SEs/variances ##
      #   final_ests$SE_resamp[final_ests$Subgrps==sub] = SE.R
      #   final_ests$LCL.CV[final_ests$Subgrps==sub] = est.R - qnorm(1-alpha/2)*SE.R
      #   final_ests$UCL.CV[final_ests$Subgrps==sub] = est.R + qnorm(1-alpha/2)*SE.R
      #   final_ests$pval.CV[final_ests$Subgrps==sub] = 2*pnorm(-abs(est.R/SE.R))
      # }
      if (resample!="CV"){
        final_ests$SE_resamp[final_ests$Subgrps==sub] = sd( est.vec ) 
      }
      ## Permutation p-value ##
      if (resample=="Permutation"){
        final_ests$pval_perm[final_ests$Subgrps==sub] =
          (sum(abs(est.vec)>abs(est0)) + 1 ) / (length(est.vec)+1)
      }
      ## Bootstrap Covariance/acceleration/bias/smoothed SE ##
      if (resample=="Bootstrap"){
        if (sub<=0){ sub_subjs = as.matrix( resamp_counter  ) }
        if (sub>0 ){ sub_subjs = as.matrix( resamp_counter[obs.data$Subgrps==sub,]) }
        # Bootstrap covariance #
        sub_subjs = apply(sub_subjs, 2, function(x) ifelse(is.na(x)==TRUE,0,x) )
        cov.b = suppressWarnings(
          rowMeans( (sub_subjs - 1) * (est.vec - mean(est.vec) ) ) )
        # Bias (smoothed estimate) #
        bias_smooth =  mean( suppressWarnings(
          rowMeans( (sub_subjs - 1)^2 * (est.vec - mean(est.vec) ) ) ) )
        # Acceleration #
        accel = (1/6)*sum(cov.b^3) / ( (sum(cov.b^2))^(3/2) )
        bca.z0 = qnorm( sum( est.vec < est0, na.rm=TRUE) / length(na.omit(est.vec))  )
        final_ests$accel[final_ests$Subgrps==sub] = accel
        final_ests$bca.z0[final_ests$Subgrps==sub] = bca.z0
        final_ests$bias_smooth[final_ests$Subgrps==sub] = bias_smooth
        final_ests$SE_smooth[final_ests$Subgrps==sub] = sqrt( sum(cov.b^2) )
        final_ests$SE_smooth_adj[final_ests$Subgrps==sub] = suppressWarnings(
          sqrt( sum(cov.b^2) - (length(Y) / (R^2))*sum( (est.vec-mean(est.vec))^2) ) 
        )
        ### Confidence Intervals (Pct and BCa) ###
        BCa.QL = pnorm(bca.z0 +(bca.z0+qnorm(alpha/2))/
                         (1-accel*(bca.z0+qnorm(alpha/2))) )
        BCa.QU = pnorm(bca.z0 +(bca.z0+qnorm(1-alpha/2))/
                         (1-accel*(bca.z0+qnorm(1-alpha/2))) )
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
  holder = NULL
  estimands = unique(param.dat$estimand)
  for (e in estimands){
    holder = rbind( holder, resamp_metrics(param.dat=param.dat, e=e)   )
  }
  param.dat = holder
  
  ### Return Outputs ###
  return( list(param.dat=param.dat, resamp.dist = resamp_param) )
}