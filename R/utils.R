
.add_trt_details <- function(Y, Subgrps, trt_dat, family) {
  
  if (!("Prob(>0)" %in% names(trt_dat))) {
    trt_dat$`Prob(>0)` <- with(trt_dat, 1-pnorm(0, mean=est, sd=SE))
  }
  if (family=="survival") {
    event.tot <- sum(Y[,2])
    event.subs <- aggregate(Y[,2] ~ Subgrps, FUN="sum")
    colnames(event.subs) <- c("Subgrps", "events")
    event.dat <- rbind( data.frame(Subgrps="ovrl", events=event.tot),
                        event.subs)
    trt_dat <- left_join(trt_dat, event.dat, by="Subgrps")
  }
  return(trt_dat)
}
trt_data_cleaner <- function(trt_dat, round_est=4, round_SE=4,
                    round_CI=4) {
  trt_dat <- trt_dat[order(trt_dat$estimand, trt_dat$Subgrps),]
  trt_dat$est <- with(trt_dat, round(est,round_est) )
  trt_dat$SE <- with( trt_dat,  round(SE,round_SE) )
  trt_dat$LCL <- with( trt_dat,  round(LCL,round_CI) )
  trt_dat$UCL <- with( trt_dat,  round(UCL,round_CI) )
  return(trt_dat)
}
summary_output_submod <- function(object, round_est=4, round_SE=4,
                                  round_CI=4) {
  
  out <- NULL
  out$submod_call <- object$submod_call
  ind_pool <- ifelse(is.null(object$trt_assign), FALSE, TRUE)
  ind_resamp <- ifelse(is.null(object$trt_eff_resamp), FALSE, TRUE)
  
  # Set up Subgroups / Trt Ests #
  trt_eff_obs <- trt_data_cleaner(object$trt_eff_obs, round_est=round_est,
                                  round_SE=round_SE, round_CI=round_CI)
  if (!ind_pool) {
    if (inherits(object, "PRISM")) {
      numb.subs <- with(object, length(unique(out.train$Subgrps)))
    }
    if (inherits(object, "submod_train")) {
      numb.subs <- with(object, length(unique(Subgrps.train))) 
    }
  }
  if (ind_pool) {
    numb.subs <- paste("Before Pooling, ", length(unique(object$trt_assign$Subgrps)))
    trt_eff_pool <- trt_data_cleaner(object$trt_eff_pool)
    pool_name <- paste("Treatment Effect Estimates (", 
                       unique(object$trt_eff_pool$type),
                       " pooling)",sep="")
    trt_eff_dopt <- trt_data_cleaner(object$trt_eff_dopt)
  }
  if (ind_resamp) {
    trt_eff_resamp <- trt_data_cleaner(object$trt_eff_resamp)
  }
  out$`Number of Identified Subgroups` <- numb.subs
  
  # Variables selected #
  if (inherits(object, "PRISM")) {
    if (c("party") %in% class(object$submod.fit$mod)) {
      submod_vars <- getUsefulPredictors(object$submod.fit$mod)
      out$`Variables that Define the Subgroups` <- paste(submod_vars, collapse=", ")
    } 
  }
  if (inherits(object, "submod_train")) {
    if (c("party") %in% class(object$fit$mod)) {
      submod_vars <- getUsefulPredictors(object$fit$mod)
      out$`Variables that Define the Subgroups` <- paste(submod_vars, collapse=", ")
    }
  }
  
  if (!ind_pool) {
    out$`Treatment Effect Estimates (observed)` <- trt_eff_obs 
  }
  if (ind_pool) {
    out$`Treatment Effect Estimates (pooling)` <- trt_eff_pool
    out$`Treatment Effect Estimates (optimal trt)` <- trt_eff_dopt
    # names(out)[4] <- pool_name 
  }
  if (ind_resamp) {
    out$`Treatment Effect Estimates (bootstrap)` <- trt_eff_resamp
  }
  # if (inherits(object, "PRISM")) {
  #   trt_eff_prism <- trt_data_cleaner(object$param.dat)
  #   out$`Treatment Effect Estimates (Final)` <- trt_eff_prism
  # }
  return(out)
}
### Default Treatment Estimates for Plots (PRISM) ###
default_trt_plots <- function(obj, est.resamp=TRUE) {
  
  # Default Setup #
  param.dat <- obj$param.dat
  param.dat$est0 <- param.dat$est
  param.dat$SE0 <- param.dat$SE
  param.dat$LCL0 <- param.dat$LCL
  param.dat$UCL0 <- param.dat$UCL
  param.dat$prob.est <- param.dat$`Prob(>0)`
  # resample <- ifelse(is.null(obj$resample), "None", obj$resample)
  # bayes <- ifelse(is.null(obj$bayes.fun), FALSE, TRUE)
  # label.param <- ""
  # if (est.resamp & (resample %in% c("Bootstrap", "CV"))) {
  #   if (resample=="Bootstrap" & is.null(param.dat$LCL.calib)) {
  #     param.dat$est0 <- param.dat$est_resamp
  #     param.dat$SE0 <- param.dat$SE_resamp
  #     param.dat$LCL0 <- param.dat$LCL.pct
  #     param.dat$UCL0 <- param.dat$UCL.pct
  #     label.param <- "(Boot,Pct)"
  #   }
  #   if (resample=="Bootstrap" & !is.null(param.dat$LCL.calib)) {
  #     param.dat$SE0 <- param.dat$SE_resamp
  #     param.dat$LCL0 <- param.dat$LCL.calib
  #     param.dat$UCL0 <- param.dat$UCL.calib
  #     label.param <- "(Boot,Calib)"
  #   }
  #   if (resample=="CV"){
  #     param.dat$est0 <- param.dat$est_resamp
  #     param.dat$LCL0 <- param.dat$LCL.CV
  #     param.dat$UCL0 = param.dat$UCL.CV
  #     label.param <- "(CV)"
  #   }
  # }
  # if (obj$family=="survival") {
  #   if (obj$param=="cox") {
  #     param.dat$est0 = exp(param.dat$est0)
  #     param.dat$LCL0 = exp(param.dat$LCL0)
  #     param.dat$UCL0 = exp(param.dat$UCL0)
  #     param.dat$estimand = gsub("logHR", "HR", param.dat$estimand)
  #     param.dat$prob.est = 1-param.dat$`Prob(>0)`
  #   }
  # }
  obj$param.dat <- param.dat
  return(obj)
}
### Extract Predictors Used in Party/MOB Tree ###
getUsefulPredictors <- function(x) {
  varid <- partykit::nodeapply(x, ids = partykit::nodeids(x),
                     FUN = function(n) partykit::split_node(n)$varid)
  varid <- unique(unlist(varid))
  names(partykit::data_party(x))[varid]
}
## Extract summary statistics (linear regression: Y~A within subgroup) ##
lm_stats = function(summ.fit, Subgrps, s, alpha, noA) {
  # Sample size #
  n.s <- length(summ.fit$residuals)
  coeff <- summ.fit$coefficients
  if (dim(coeff)==1) {
    est <- coeff[1,1]
    SE <- coeff[1,2]
    LCL <- est - qt(1-alpha/2, df=n.s-1)*SE
    UCL <- est + qt(1-alpha/2, df=n.s-1)*SE
    pval <- 2*pt(-abs(est/SE), df=n.s-1)
    summ <- data.frame( Subgrps = s, N = n.s, 
                       estimand = "E(Y)", est, SE, LCL, UCL, pval)
  }
  if (dim(coeff)[1]>1) {
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
    looper <- function(s) {
      ns <- param.e$N[param.e$Subgrps==s]
      obs.mu <- param.e$est[param.e$Subgrps==s]
      obs.var <- ( param.e$SE[param.e$Subgrps==s] )^2
      ### Mean/Var ##
      if (s==0) {
        mu_B <- obs.mu
        var_B <- obs.var
        alpha <- alpha_ovrl
      }
      if (s>0) {
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
  for (e in unique(param.dat$estimand)) {
    param.B = get_posterior_ests(param.dat=param.dat, e = e)
    bayes_param = suppressWarnings( bind_rows(bayes_param, param.B) )
  }
  bayes_param$`Prob(>0)` = with(bayes_param, 1-pnorm(0, mean=est.bayes, sd=SE.bayes))
  param.dat <- suppressWarnings(
    left_join(param.dat, bayes_param, by=c("Subgrps", "estimand")))
  bayes.sim = function(mu, sd, n=100000) {
    return( rnorm(n=n, mean = mu, sd = sd)  )
  }
  return( list(param.dat=param.dat, bayes.sim=bayes.sim) )
}

## Generate bootstrap resampling indices ##
bootstrap_indices = function(n, R, strata=NULL) {
  if (is.null(strata)) {
    indices <- sample.int(n, n * R, replace = TRUE)
    dim(indices) <- c(R, n)
  }
  if (!is.null(strata)) {
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
permute_indices = function(n, R, strata=NULL) {
  if (is.null(strata)) {
    indices <- matrix(NA, R, n)
    for (r in 1:R) {
      indices[r,] <- sample(n, size=n, replace=FALSE)
    }
  }
  if (!is.null(strata)) {
    indices <- matrix(NA, R, n)
    for (s in unique(strata)) {
      sub_i <- seq_len(n)[strata == s]
      ns = length(sub_i)
      for (r in 1:R) {
        indices[r, sub_i] <- sample(sub_i, replace=FALSE)
      }
    }
  }
  return(indices)
}
## Generate CV sampling folds ###
CV_folds = function(n, R, strata=NULL) {
  if (is.null(strata)) {
    folds <- sample(rep(1:R,ceiling(n/R))[1:n])
  }
  if (!is.null(strata)) {
    folds = NULL
    for (s in unique(strata)) {
      sub_i <- seq_len(n)[strata == s]
      ns = length(sub_i)
      folds.s <- sample(rep(1:R,ceiling(ns/R))[1:ns])
      folds[sub_i] = folds.s
    }
  }
  return(folds)
}

## Coverage counter(for bootstrap calibration) ##
coverage_counter = function(param.dat, alpha.mat) {
  
  counts = NULL
  counter = function(i) {
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
calibrate_alpha = function(resamp_calib, alpha_ovrl, alpha_s) {
  
  # Average across overall/subgroups and within alphas #
  calib.dat = aggregate(ind ~ ovrl_ind*alpha, data=resamp_calib, FUN="mean")
  # Linear calibration #
  fn = function(alpha, ind) {
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
resamp_metrics = function(trt_eff, resamp_dist, resamp_calib, resample) {

  calibrate <- FALSE
  alpha_ovrl <- unique(trt_eff$alpha[trt_eff$Subgrps=="ovrl"])
  alpha_s <- unique(trt_eff$alpha[trt_eff$Subgrps!="ovrl"])
  ## Calculate calibrated alpha ##
  if (!is.null(resamp_calib)) {
    alpha_c <- calibrate_alpha(resamp_calib=resamp_calib, alpha_ovrl=alpha_ovrl,
                              alpha_s=alpha_s)
    alpha.ovrl_c <- alpha_c$alpha_calib[alpha_c$ovrl_ind==1]
    alpha.s_c <- alpha_c$alpha_calib[alpha_c$ovrl_ind==0]
    trt_eff$alpha_c <- with(trt_eff, ifelse(Subgrps==0,alpha.ovrl_c,
                                               alpha.s_c))
    calibrate <- TRUE
  }
  # Loop through estimands #
  looper = function(e) {
    
    trt_obs = trt_eff[trt_eff$estimand==e,]
    # if (calibrate) {
    #   LCL_calib = with(trt_resamp, est - qnorm(1-alpha_c/2, df=N-1)*SE)
    #   UCL_calib = with(trt_resamp, est + qnorm(1-alpha_c/2, df=N-1)*SE)
    # }
    est_out <- NULL
    for (sub in unique(trt_obs$Subgrps)) {
      if (sub=="ovrl") { alpha = alpha_ovrl }
      if (sub!="ovrl") { alpha = alpha_s }
      
      hold = trt_obs[trt_obs$Subgrps==sub,]
      hold.R = resamp_dist[resamp_dist$Subgrps==sub & resamp_dist$estimand==e,]
      est0 = hold$est
      est.vec = hold.R$est

      pval_resamp <- NA
      bias_resamp <- NA
      if (resample=="Permutation") {
        est_resamp <- mean(hold.R$est, na.rm=TRUE)
        SE_resamp <- sd( est.vec, na.rm=TRUE)
        pval_resamp <- (sum(abs(est.vec)>abs(est0), na.rm=TRUE) + 1) / 
          (length(na.omit(est.vec))+1)
      }
      if (resample=="Bootstrap") {
        est_resamp <- mean(hold.R$est, na.rm=TRUE)
        SE_resamp <- sd( est.vec, na.rm=TRUE)
        bias_resamp <- mean(hold.R$bias, na.rm=TRUE)
        quants <- as.numeric(
          quantile(est.vec, probs=c(alpha/2, (1-alpha/2)), na.rm = TRUE) )
        LCL_resamp <- quants[1]
        UCL_resamp <- quants[2]
      }
      if (resample=="CV") {
        N.tot = sum(hold.R$N)
        est.CF <- with(hold.R, weighted.mean(est, N))
        SE.CF <- with(hold.R, sqrt(  sum( (N/N.tot)^2*SE^2) ) )
        est_resamp <- est.CF
        SE_resamp <- SE.CF
        LCL_resamp <- est.CF - qnorm(1-alpha/2)*SE.CF
        UCL_resamp <- est.CF + qnorm(1-alpha/2)*SE.CF
      }
      hold <- data.frame(Subgrps = sub, estimand=e, est = est_resamp, 
                         SE=SE_resamp, LCL = LCL_resamp, UCL = UCL_resamp,
                         pval = pval_resamp, alpha = alpha)
      est_out <- rbind(est_out, hold)
    }
    return(est_out)
  }
  trt_resamp = lapply(unique(trt_eff$estimand), looper)
  trt_resamp = do.call(rbind, trt_resamp)
  trt_resamp$resample <- resample
  return(trt_resamp)
}

#### Survival Helper Functions ####

# Weibull Functions (for MOB) #
wbreg <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...) {
  survreg(y ~ 0 + x, weights = weights, dist = "weibull", ...)
}
logLik.survreg <- function(object, ...) {
  structure(object$loglik[2], df = sum(object$df), class = "logLik")
}
aft_mavg <- function(y, x, start = NULL, weights = NULL, offset=NULL,
                      ..., 
                      estfun = FALSE, object = FALSE) {
  
  # x <- x[, 2]
  dist.vec <- c("weibull", "lognormal", "loglogistic")
  obj <- list()
  fits <- NULL
  est_fun <- list()
  for (dist in dist.vec) {
    
    fit <- tryCatch(survreg(y ~ 0 + x, dist=dist, ...,
                            weights = weights), 
                    error = function(e) paste("error", e),
                    warning = function(w) paste("convergence error", w))
    
    if (is.character(fit)) {
      hold <- data.frame(dist=dist, AIC=NA, loglik=NA, df=NA,
                         est=NA, SE=NA, parm=NA)
    }
    if (is.list(fit)) {
      # tbl <- summary(fit)$table
      # est <- tbl[,1]
      # SE <- tbl[,2]
      est <- coef(fit)
      SE <- summary(fit)$table[names(est),2]
      hold <- data.frame(dist=dist, AIC=stats::AIC(fit), 
                         loglik = fit$loglik[2],
                         df = sum(fit$df),
                         est=est, SE=NA)
      hold$parm <- rownames(hold)
      rownames(hold) <- NULL
      est_fun00 <- sandwich::estfun(fit)
      est_fun <- append(est_fun, list(est_fun00))
      obj <- append(obj, list(fit))
    } 
    fits <- rbind(fits, hold)
  }
  names(obj) <- dist.vec
  fits <- fits[!is.na(fits$est),]
  # model averaging #
  out_dat <- NULL
  for (p in unique(fits$parm)) {
    ma.fit <- .model_average(fits=fits[fits$parm==p,])
    w.MA <- ma.fit$fits$w.MA
    hold <- data.frame(parm = p, est = ma.fit$est.MA)
    hold$loglik <- sum(ma.fit$fits$w.MA %*% ma.fit$fits$loglik)
    hold$AIC <- sum(ma.fit$fits$w.MA %*% ma.fit$fits$AIC)
    hold$df <- mean(ma.fit$fits$df)
    out_dat <- rbind(out_dat, hold)
  }
  obj <- append(obj, 
                list(wdat=data.frame(ma.fit$fits[,c("dist", "AIC", "w.MA")])))
  class(obj) <- "aft_mavg"
  # Obtain weighted est FUN #
  est_fun1 <- 0
  for (i in 1:length(w.MA)) {
    
    est_fun1 <- est_fun1 + w.MA[i]*est_fun[[i]]
    
  }
  coef_hold <- out_dat$est
  names(coef_hold) <- out_dat$parm
  
  out <- list(coefficients = coef_hold, 
              objfun = -mean(out_dat$loglik),
              estfun = if (estfun) est_fun1 else NULL,
              object = if (object) obj else NULL)
  return(out)
}
predict.aftmavg <- function(object, newdata) {
  
  dist.vec <- object$wdat$dist
  pred <- 0
  for (dist in dist.vec) {
    hold <- predict(object[[dist]], newdata=newdata)
    pred <- pred + object$wdat$w.MA[object$wdat$dist==dist]*hold
  }
  return(pred)
}
# RMST Estimation: based on survRM2 #
rmst_single = function(time, status, tau=NULL) {
  if (is.null(tau)) {
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
# IPC Weights (based on "Survival Ensembles 2006 Hothorn T et al") #
ipc_weights <- function(y_surv, max_w = 5) {
  
  event <- y_surv[,2]
  y_cens <- survival::Surv(y_surv[,1], 1-event)
  
  mod <- survfit(y_cens ~ 1)
  times <- y_surv[,1]
  cens_est <- rep(NA, length(times))
  ind_vec <- times %in% mod$time
  
  for (i in 1:length(cens_est)) {
    if (ind_vec[i]) {
      cens_est[i] <- mod$surv[mod$time == times[i]]
    }
    if (!ind_vec[i]) {
      before_t <- mod$time[mod$time < times[i]]
      if (length(before_t) == 0) {
        cens_est[i] <- 1
      }
      if (length(before_t) > 0) {
        cens_est[i] <- mod$surv[mod$time == max(before_t)]
      }
    }
  }
  cens_est[event == 0] <- 1
  ipc_w <- event / cens_est
  ipc_w[ipc_w > max_w] <- max_w
  
  return(ipc_w)
}

## P-value converter ##
pval_convert <- function(p_value) {
  if (is.na(p_value)) { return( "" )}
  if (p_value < 0.001) return(c("p<0.001"))
  if (p_value >= 0.001) return(paste("p=",(round(p_value,3)),sep=""))
  else return("")
}
## RMST calculator: tau is RMST truncation time ##
rmst_calc <- function(time, surv, tau) {
  idx = time <= tau
  id.time = sort(c(time[idx], tau))
  id.surv = surv[idx]
  time.diff <- diff(c(0, id.time))
  areas <- time.diff * c(1, id.surv)
  rmst = sum(areas)
  return(rmst)
}

## Probability Calculator (input desired threshold and PRISM.fit) ##
prob_calculator <- function(fit, thres=">0") {
  
  param.dat <- fit$param.dat
  thres.name <- paste("Prob(",thres, ")", sep="")
  # Stop if threshold has already been computed #
  if (thres.name %in% colnames(param.dat)) {
    return(fit)
  }
  dir <- substr(thres, 1, 1)
  numb <- substr(thres, 2, nchar(thres))
  if (suppressWarnings(is.na(as.numeric(numb)))) {
    numb <- substr(thres, 3, nchar(thres))
  }
  numb <- as.numeric(numb)
  if (fit$param=="param_cox") {
    numb <- log(numb)
    thres <- paste(dir,numb,sep=" ")
  }
  # Check resampling / bayesian #
  if (inherits(fit, "PRISM")) {
    resamp <- ifelse(is.null(fit$resample), "None", as.character(fit$resample))
    bayes <- ifelse(is.null(fit$bayes.fun), FALSE, TRUE)
  }
  if (inherits(fit, "submod_train")) {
    resamp <- "None"
    bayes <- FALSE
  }
  # Normal Approx (if no normal/bayes)
  if (resamp=="None" & !bayes) {
    prob.est <- with(param.dat, 1-pnorm(numb, mean=est, sd=SE))
    if (dir=="<") prob.est <- 1-prob.est
    param.dat$prob.est <- prob.est
  }
  # Bayes #
  if (resamp=="None" & bayes) {
    prob.est <- with(param.dat, 1-pnorm(numb, mean=est.bayes, sd=SE.bayes))
    if (dir=="<") prob.est <- 1-prob.est
    param.dat$prob.est <- prob.est
  }
  # CV #
  if (resamp=="CV") {
    prob.est <- with(param.dat, 1-pnorm(numb, mean=est_resamp, sd=SE_resamp))
    if (dir=="<") prob.est <- 1-prob.est
    param.dat$prob.est <- prob.est
  }
  # Resampling #
  if (resamp %in% c("Bootstrap", "Permutation")) {
    rdist <- fit$resamp_dist
    rdist$prob.est = eval(parse(text=paste("ifelse(rdist$est",thres, ", 1, 0)")))
    prob.dat <- aggregate(prob.est ~ Subgrps*estimand, data=rdist, FUN="mean")
    if (fit$param=="param_cox") {
      prob.dat$estimand = gsub("logHR", "HR", prob.dat$estimand)
    }
    param.dat <- param.dat[,!(colnames(param.dat) %in% c("prob.est"))]
    param.dat <- left_join(param.dat, prob.dat, by=c("Subgrps", "estimand"))
  }
  # Create column name #
  # colnames(param.dat)[which(colnames(param.dat)=="prob.est")] <- thres.name
  fit$param.dat <- param.dat
  return(fit)
}
#### Modeling Averaging Functions ####
.weights.MA = function(weights){
  w = weights
  w = w-min(w)
  norm.const = sum(exp(-w/2))
  w_out = exp(-w/2)/norm.const
  return(w_out)
}
.model_average <- function(fits, type="AIC"){
  
  fits$w.MA <- .weights.MA(fits$AIC)
  est.MA <- with(fits, as.numeric(w.MA %*% est) )
  var.MA <- with(fits, sum(w.MA*(SE^2+(est-est.MA)^2)))
  SE.MA <- sqrt(var.MA)
  
  return( list(est.MA=est.MA, SE.MA=SE.MA, var.MA=var.MA, fits=fits))
}