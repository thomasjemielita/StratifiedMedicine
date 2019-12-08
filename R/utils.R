### Extract Predictors Used in Party/MOB Tree ###
getUsefulPredictors <- function(x) {
  varid <- nodeapply(x, ids = nodeids(x),
                     FUN = function(n) split_node(n)$varid)
  varid <- unique(unlist(varid))
  names(data_party(x))[varid]
}


## Extract summary statistics (linear regression: Y~A within subgroup) ##
lm_stats = function(summ.fit, Subgrps, s, alpha, noA){
  # Sample size #
  n.s <- length(summ.fit$residuals)
  coeff <- summ.fit$coefficients
  if (dim(coeff)==1){
    est <- coeff[1,1]
    SE <- coeff[1,2]
    LCL <- est - qt(1-alpha/2, df=n.s-1)*SE
    UCL <- est + qt(1-alpha/2, df=n.s-1)*SE
    pval <- 2*pt(-abs(est/SE), df=n.s-1)
    summ <- data.frame( Subgrps = s, N = n.s, 
                       estimand = "E(Y)", est, SE, LCL, UCL, pval)
  }
  if (dim(coeff)[1]>1){
    L.mat = rbind( c(1,0), c(1,1), c(0,1) )
    est = L.mat %*% coeff[,1]
    SE = sqrt(  diag( L.mat %*% vcov(summ.fit) %*% t(L.mat) ) )
    LCL = est - qt(1-alpha/2, df=n.s-1)*SE
    UCL = est + qt(1-alpha/2, df=n.s-1)*SE
    pval = 2*pt(-abs(est/SE), df=n.s-1)
    summ = data.frame( Subgrps = s, N = n.s, 
                       estimand = c("E(Y|A=0)", "E(Y|A=1)", "E(Y|A=1)-E(Y|A=0)"),
                       est, SE, LCL, UCL, pval) 
  }
  return(summ)
}

## Bayesian: Normal prior (at overall for subgroups), normal posterior ##
norm_norm <- function(PRISM.fit, alpha_ovrl, alpha_s,
                      scale = length(PRISM.fit$Subgrps.train), ...) {
  param.dat <- PRISM.fit$param.dat
  Subgrps <- PRISM.fit$Subgrps.train
  get_posterior_ests <- function(param.dat, e) {
    # Set prior to overall #
    param.e <- param.dat[param.dat$estimand==e,]
    prior.mu <- param.e$est[param.e$Subgrps==0]
    prior.var <- scale * ( param.e$SE[param.e$Subgrps==0] )^2
    looper <- function(s){
      ns <- param.e$N[param.e$Subgrps==s]
      obs.mu <- param.e$est[param.e$Subgrps==s]
      obs.var <- ( param.e$SE[param.e$Subgrps==s] )^2
      ### Mean/Var ##
      if (s==0){
        mu_B <- obs.mu
        var_B <- obs.var
        alpha <- alpha_ovrl
      }
      if (s>0){
        var_B <- ( 1/prior.var + 1/(obs.var) )^(-1)
        mu_B <- var_B * ( prior.mu/prior.var + obs.mu / (obs.var)    ) 
        alpha <- alpha_s
      }
      # Bayesian Intervals (q25, q75) #
      LCL <- qnorm( alpha/2, mean = mu_B, sd = sqrt(var_B) )
      UCL <- qnorm( 1-alpha/2, mean = mu_B, sd = sqrt(var_B) )
      summ <- data.frame(Subgrps=s, estimand=e,
                        est.bayes = mu_B, SE.bayes = sqrt(var_B), 
                        LCL.bayes=LCL, UCL.bayes = UCL)
      return(summ)
    }
    param.B = lapply( c(unique(Subgrps),0), looper)
    param.B = do.call(rbind, param.B)
    return(param.B)
  }
  bayes_param = NULL
  for (e in unique(param.dat$estimand)){
    param.B = get_posterior_ests(param.dat=param.dat, e = e)
    bayes_param = suppressWarnings( bind_rows(bayes_param, param.B) )
  }
  param.dat <- left_join(param.dat, bayes_param, by=c("Subgrps", "estimand"))
  bayes.sim = function(mu, sd, n=100000){
    return( rnorm(n=n, mean = mu, sd = sd)  )
  }
  return( list(param.dat=param.dat, bayes.sim=bayes.sim) )
}

## Generate bootstrap resampling indices ##
bootstrap_indices = function(n, R, strata=NULL){
  if (is.null(strata)){
    indices <- sample.int(n, n * R, replace = TRUE)
    dim(indices) <- c(R, n)
  }
  if (!is.null(strata)){
    indices <- matrix(NA, R, n)
    for (s in unique(strata)) {
      sub_i <- seq_len(n)[strata == s]
      ns = length(sub_i)
      indices[, sub_i] <- sub_i[sample.int(length(sub_i), ns*R, replace = TRUE)]
    }
  }
  return(indices)
}
## Generate permutation resampling indices ##
permute_indices = function(n, R, strata=NULL){
  if (is.null(strata)){
    indices <- matrix(NA, R, n)
    for (r in 1:R){
      indices[r,] <- sample(n, size=n, replace=FALSE)
    }
  }
  if (!is.null(strata)){
    indices <- matrix(NA, R, n)
    for (s in unique(strata)) {
      sub_i <- seq_len(n)[strata == s]
      ns = length(sub_i)
      for (r in 1:R){
        indices[r, sub_i] <- sample(sub_i, replace=FALSE)
      }
    }
  }
  return(indices)
}
## Generate CV sampling folds ###
CV_folds = function(n, R, strata=NULL){
  if (is.null(strata)){
    folds <- sample(rep(1:R,ceiling(n/R))[1:n])
  }
  if (!is.null(strata)){
    folds = NULL
    for (s in unique(strata)){
      sub_i <- seq_len(n)[strata == s]
      ns = length(sub_i)
      folds.s <- sample(rep(1:R,ceiling(ns/R))[1:ns])
      folds[sub_i] = folds.s
    }
  }
  return(folds)
}

## Coverage counter(for bootstrap calibration) ##
coverage_counter = function(param.dat, alpha.mat){
  
  counts = NULL
  counter = function(i){
    alpha_s = alpha.mat$alpha_s[i]
    alpha_ovrl = alpha.mat$alpha_ovrl[i]
    param.dat$alpha = with(param.dat, ifelse(Subgrps==0, alpha_ovrl, alpha_s))
    param.dat$LCL.a = with(param.dat, est - qt(1-alpha/2, df=N-1)*SE)
    param.dat$UCL.a = with(param.dat, est + qt(1-alpha/2, df=N-1)*SE)
    param.dat$ind = with(param.dat, ifelse(est.obs<=UCL.a & LCL.a<=est.obs, 1, 0) )
    cnt = param.dat[,c("Subgrps", "estimand", "est", "alpha", "ind")]
    return(cnt)
  }
  counts = lapply(1:dim(alpha.mat)[1], counter)
  counts = data.frame( do.call(rbind,counts) )
  return(counts)
}
calibrate_alpha = function(resamp_calib, alpha_ovrl, alpha_s){
  
  # Average across overall/subgroups and within alphas #
  calib.dat = aggregate(ind ~ ovrl_ind*alpha, data=resamp_calib, FUN="mean")
  # Linear calibration #
  fn = function(alpha, ind){
    (sum( alpha*(1-ind)) / sum(alpha^2))^(-1)
  }
  alpha.ovrl_c = with(calib.dat[calib.dat$ovrl_ind==1,], fn(alpha, ind) )
  alpha.s_c = with(calib.dat[calib.dat$ovrl_ind==0,], fn(alpha, ind) )
  out = data.frame(ovrl_ind = c(1,0), alpha_nom = c(alpha_ovrl, alpha_s),
                   alpha_calib = c(alpha.ovrl_c, alpha.s_c))
  out$alpha_calib = with(out, alpha_calib*alpha_nom)
  return(out)
}
### Resampling metric calculation ####
resamp_metrics = function(param.dat, resamp_param, resamp_calib, resample){

  calibrate=FALSE
  alpha_ovrl = unique(param.dat$alpha[param.dat$Subgrps==0])
  alpha_s = unique(param.dat$alpha[param.dat$Subgrps>0])
  ## Calculate calibrated alpha ##
  if (!is.null(resamp_calib)){
    alpha_c = calibrate_alpha(resamp_calib=resamp_calib, alpha_ovrl=alpha_ovrl,
                              alpha_s=alpha_s)
    alpha.ovrl_c = alpha_c$alpha_calib[alpha_c$ovrl_ind==1]
    alpha.s_c = alpha_c$alpha_calib[alpha_c$ovrl_ind==0]
    param.dat$alpha_c = with(param.dat, ifelse(Subgrps==0,alpha.ovrl_c,
                                               alpha.s_c))
    calibrate=TRUE
  }
  # Loop through estimands #
  looper = function(e){
    final_ests = param.dat[param.dat$estimand==e,]
    final_ests$est_resamp = NA
    final_ests$SE_resamp = NA
    if (resample=="Permutation"){
      final_ests$pval_perm = NA
    }
    if (resample=="Bootstrap"){
      final_ests = final_ests %>% mutate(bias.boot = NA, LCL.pct = NA, UCL.pct= NA)
      if (calibrate){
        final_ests$LCL.calib = with(final_ests, est - qt(1-alpha_c/2, df=N-1)*SE)
        final_ests$UCL.calib = with(final_ests, est + qt(1-alpha_c/2, df=N-1)*SE)
      }
    }
    for (sub in unique(final_ests$Subgrps)){
      if (sub<=0){ alpha = alpha_ovrl }
      if (sub>0 ){ alpha = alpha_s }
      hold = final_ests[final_ests$Subgrps==sub,]
      hold.R = resamp_param[resamp_param$Subgrps==sub & resamp_param$estimand==e,]
      est0 = hold$est
      est.vec = hold.R$est
      ## Permutation (est, SE, p-value) ##
      if (resample=="Permutation"){
        final_ests$est_resamp[final_ests$Subgrps==sub] = mean(hold.R$est, na.rm=TRUE)
        final_ests$SE_resamp[final_ests$Subgrps==sub] = sd( est.vec, na.rm=TRUE)
        final_ests$pval_perm[final_ests$Subgrps==sub] =
          (sum(abs(est.vec)>abs(est0), na.rm=TRUE) + 1 ) / (length(na.omit(est.vec))+1)
      }
      ## Bootstrap (smoothed est, SE, bias, pct CI) ##
      if (resample=="Bootstrap"){
        final_ests$est_resamp[final_ests$Subgrps==sub] = mean(hold.R$est, na.rm=TRUE)
        final_ests$SE_resamp[final_ests$Subgrps==sub] = sd( est.vec, na.rm=TRUE)
        final_ests$bias.boot[final_ests$Subgrps==sub] = mean(hold.R$bias, na.rm=TRUE)
        quants = as.numeric(
          quantile(est.vec, probs=c(alpha/2, (1-alpha/2)), na.rm = TRUE) )
        final_ests$LCL.pct[final_ests$Subgrps==sub] = quants[1]
        final_ests$UCL.pct[final_ests$Subgrps==sub] = quants[2]
      }
      ## Cross-validation (Cross-fitting estimate) ##
      if (resample=="CV"){
        N.tot = sum(hold.R$N)
        est.CF <- with(hold.R, weighted.mean(est, N))
        SE.CF <- with(hold.R, sqrt(  sum( (N/N.tot)^2*SE^2) ) )
        final_ests$est_resamp[final_ests$Subgrps==sub] = est.CF
        final_ests$SE_resamp[final_ests$Subgrps==sub] = SE.CF
        final_ests$LCL.CV[final_ests$Subgrps==sub] = est.CF - qnorm(1-alpha/2)*SE.CF
        final_ests$UCL.CV[final_ests$Subgrps==sub] = est.CF + qnorm(1-alpha/2)*SE.CF
      }
    }
    return(final_ests)
  }
  final_ests = lapply(unique(param.dat$estimand), looper)
  final_ests = do.call(rbind, final_ests)
  return(final_ests)
}

#### Weibull Mob Functions ###
wbreg <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...){
  survreg(y ~ 0 + x, weights = weights, dist = "weibull", ...)
}
logLik.survreg <- function(object, ...){
  structure(object$loglik[2], df = sum(object$df), class = "logLik")
}

### RMST Estimation: based on survRM2 ###
rmst_single = function(time, status, tau=NULL){
  if (is.null(tau)){
    tau = max(time)
  }
  km = survfit(Surv(time, status) ~ 1)
  keep = km$time <= tau
  time0 = sort(c(km$time[keep], tau))
  surv0 = km$surv[keep]
  nrisk0 = km$n.risk[keep]
  nevent0 = km$n.event[keep]
  time.diff <- diff(c(0, time0))
  AUC <- time.diff * c(1, surv0)
  est = sum(AUC)
  var.vec <- ifelse((nrisk0 - nevent0) == 0, 0, nevent0/(nrisk0 * 
                                                       (nrisk0 - nevent0)))
  var.vec = c(var.vec, 0)
  SE = sqrt( sum(cumsum(rev(AUC[-1]))^2 * rev(var.vec)[-1]) )
  return( list(rmst=est, rmst.se=SE))
}

## P-value converter ##
pval_convert <- function(p_value) {
  if (is.na(p_value)){ return( "" )}
  if (p_value < 0.001) return(c("p<0.001"))
  if (p_value >= 0.001) return(paste("p=",(round(p_value,3)),sep=""))
  else return("")
}
## RMST calculator: tau is RMST truncation time ##
rmst_calc <- function(time, surv, tau){
  idx = time <= tau
  id.time = sort(c(time[idx], tau))
  id.surv = surv[idx]
  time.diff <- diff(c(0, id.time))
  areas <- time.diff * c(1, id.surv)
  rmst = sum(areas)
  return(rmst)
}
