## T-Learner: Fit separate regression model by exposure level (A) ##
T_learner <- function(Y, A, X, Xtest, family, ple, hyper, pfit, tau=NULL, ...) {
  
  A_lvls <- unique(A)[order(unique(A))]
  Xfull <- if (is.null(Xtest)) X else rbind(X, Xtest)
  n_tr <- dim(X)[1]
  n_ts <- dim(Xfull)[1]
  ids_tr <- 1:nrow(X)
  
  fit_learner <- function(a) {
    fit <- do.call(ple, append(list(Y=Y[A==a], X=X[A==a,], 
                                    Xtest = Xfull,
                                    family=family), hyper))
    if (!is.null(fit$mu_train)) {
      mu_tr <- data.frame(ids = ids_tr[A==a], mu_a = fit$mu_train)
      if (is.null(fit$mu_test)) {
        mu_ts <- fit$pred.fun(fit$mod, Xfull[A!=a,], tau=tau)
      }
      if (!is.null(fit$mu_test)) {
        mu_ts <- fit$mu_test[A!=a]
      } 
      mu_ts <- data.frame(ids = ids_tr[A!=a], mu_a = mu_ts)
      mu_dat <- rbind(mu_tr, mu_ts)
      mu_dat <- mu_dat[order(mu_dat$ids),]
      mu_a <- mu_dat$mu_a
    }
    if (is.null(fit$mu_train)) {
      if (is.null(fit$mu_test)) {
        mu_a <- fit$pred.fun(fit$mod, Xfull, tau=tau)
      }
      if (!is.null(fit$mu_test)) {
        mu_a <- fit$mu_test
      } 
    }
    
    return(list(fit=fit, mu_a=mu_a))
  }
  res <- lapply(A_lvls, fit_learner)
  res <- do.call(cbind, res)
  fit0 <- res[1,]
  names(fit0) <- paste("fit", A_lvls, sep="_")
  fit <- list(fit0=fit0, fitp=pfit)
  
  mu_hat <- data.frame(do.call(cbind, res[2,]))
  colnames(mu_hat) <- paste("mu", A_lvls, sep="_")
  for (aa in A_lvls[-1]) {
    diff_name <- paste("diff", aa, A_lvls[1], sep="_")
    aa_name <- paste("mu", aa, sep="_")
    ref_name <- paste("mu", A_lvls[1], sep="_")
    mu_hat[diff_name] <- mu_hat[aa_name] - mu_hat[ref_name] 
  }
  pi_hat <- data.frame(fit$fitp$pred.fun(fit$fitp$mod, X))
  mu_hat <- data.frame(mu_hat, pi_hat)
  n_tr <- dim(X)[1]
  n_ts <- dim(Xfull)[1]
  mu_train <- mu_hat[1:n_tr,]
  mu_test <- if (!is.null(Xtest)) mu_hat[(n_tr+1):n_ts,] else NULL
  
  pred.fun <- function(fit, X, tau=NULL) {
    fit0 <- fit$fit0
    A_lvls <- sub(".*_", "", names(fit0))
    looper <- function(a) {
      name <- paste("fit", a, sep="_")
      mu_a <- fit0[[name]]$pred.fun(fit0[[name]]$mod, X=data.frame(X), tau=tau)
      return(mu_a)
    }
    mu_hat <- lapply(A_lvls, looper)
    mu_hat <- data.frame(do.call(cbind, mu_hat))
    colnames(mu_hat) <- paste("mu", A_lvls, sep="_")
    for (aa in A_lvls[-1]) {
      diff_name <- paste("diff", aa, A_lvls[1], sep="_")
      aa_name <- paste("mu", aa, sep="_")
      ref_name <- paste("mu", A_lvls[1], sep="_")
      mu_hat[diff_name] <- mu_hat[aa_name] - mu_hat[ref_name] 
    }
    pi_hat <- data.frame(fit$fitp$pred.fun(fit$fitp$mod, X))
    mu_hat <- data.frame(mu_hat, pi_hat)
    return(mu_hat)
  }
  res = list(mod=fit, pred.fun=pred.fun, mu_train=mu_train, mu_test=mu_test)
  class(res) = "T_learner"
  return(res)
}
## S-Learner: Fit single Regression Model (can do forced interactions, VT) ##
S_learner <- function(Y, A, X, Xtest, VT=FALSE, family, ple, hyper, pfit,
                      tau=NULL,...) {
  
  A_lvls <- unique(A)[order(unique(A))]
  
  if (!VT) {
    fit <- do.call(ple, append(list(Y=Y, X=data.frame(A=A, X),
                                    family=family), hyper))
    fit$A_lvls <- A_lvls
    fit <- list(fit0=fit, fitp=pfit)
    pred.fun <- function(fit, X, tau=NULL) {
      fit0 <- fit$fit0
      looper <- function(a) {
        mu_a <- fit0$pred.fun(fit0$mod, X=data.frame(A=a, X), tau=tau)
        return(mu_a)
      }
      mu_hat <- lapply(A_lvls, looper)
      mu_hat <- data.frame(do.call(cbind, mu_hat))
      colnames(mu_hat) <- paste("mu", A_lvls, sep="_")
      for (aa in A_lvls[-1]) {
        diff_name <- paste("diff", aa, A_lvls[1], sep="_")
        aa_name <- paste("mu", aa, sep="_")
        ref_name <- paste("mu", A_lvls[1], sep="_")
        mu_hat[diff_name] <- mu_hat[aa_name] - mu_hat[ref_name] 
      }
      pi_hat <- data.frame(fit$fitp$pred.fun(fit$fitp$mod, X))
      mu_hat <- data.frame(mu_hat, pi_hat)
      return(mu_hat)
    }
  }
  if (VT) {
    ## TO DO ##
  }
  res = list(mod=fit, pred.fun=pred.fun)
  class(res) = "S_learner"
  return(res)
}
## X-Learner ##
X_learner <- function(Y, A, X, Xtest, family, ple, hyper, pfit, tau=NULL, ...) {
  
  A_lvls <- unique(A)[order(unique(A))]
  
  # Stage 1: T-Learner (treatment-specific models) #
  stage1 <- T_learner(Y=Y, A=A, X=X, Xtest=Xtest,
                      family = family, ple=ple, hyper=hyper, 
                      pfit = pfit, tau = tau)
  
  mu_tr <- stage1$mu_train

  # Step 2: Imputed Treatment Effects #
  next_learner <- function(aa) {
    aa_name <- paste("mu", aa, sep="_")
    ref_name <- paste("mu", A_lvls[1], sep="_")
    if (is.Surv(Y)) { Y <- Y[,1]; }
    imp_aa <- (Y[A==aa]-mu_tr[ref_name])[,1]
    imp_ref <- (mu_tr[aa_name]-Y[A==A_lvls[1]])[,1]
    fit_aa <- do.call(ple, append(list(Y=imp_aa, X=X, Xtest=Xtest,
                                       family="gaussian"), hyper))
    fit_ref <- do.call(ple, append(list(Y=imp_ref, X=X, Xtest=Xtest,
                                        family="gaussian"), hyper))
    fit <- list(fit_aa=fit_aa, fit_ref=fit_ref)
    return(fit)
  }
  fitd <- lapply(A_lvls[-1], next_learner)
  fitd <- do.call(rbind, fitd)
  rownames(fitd) <- paste("fit", A_lvls[-1], A_lvls[1], sep="_")
  
  # Obtain final mu_train / mu_test #
  mu_train <- stage1$mu_train
  mu_test <- stage1$mu_test
  diff_tr <- NULL
  diff_ts <- NULL
  diff_names <- gsub("fit", "diff", rownames(fitd))
  diff_lvls <- do.call(rbind, strsplit(diff_names, split="_"))
  for (ii in 1:length(diff_names)) {
    aa <- diff_lvls[ii,2]
    ref <- diff_lvls[ii,3]
    diff_name <- paste("diff", aa, ref, sep="_")
    fitd_ii <- fitd[ii,]
    # Training #
    tau_aa <- fitd_ii$fit_aa$mu_train
    tau_ref <- fitd_ii$fit_ref$mu_train
    pi_hat <- data.frame(pfit$pred.fun(pfit$mod, X))
    pi_aa <- as.numeric(pi_hat[[paste("pi", aa, sep="_")]])
    pi_ref <- as.numeric(pi_hat[[paste("pi", ref, sep="_")]])
    tau_hat <- tau_aa*pi_ref + tau_ref*pi_aa
    mu_train[diff_name] <-  tau_hat
    if (!is.null(mu_test)) {
      tau_aa <- fitd_ii$fit_aa$mu_test
      tau_ref <- fitd_ii$fit_ref$mu_test
      pi_hat <- data.frame(pfit$pred.fun(pfit$mod, Xtest))
      pi_aa <- as.numeric(pi_hat[[paste("pi", aa, sep="_")]])
      pi_ref <- as.numeric(pi_hat[[paste("pi", ref, sep="_")]])
      tau_hat <- tau_aa*pi_ref + tau_ref*pi_aa
      mu_test[diff_name] <-  tau_hat
    }
  }
  
  fit <- list(fit0=stage1$mod$fit0, fitd=fitd, fitp=pfit)
  
  # Prediction Function #
  pred.fun <- function(fit, X, tau=NULL) {
    fit0 <- fit$fit0
    fitd <- fit$fitd
    fitp <- fit$fitp
    # Treatment-specific # 
    A_lvls <- sub(".*_", "", names(fit$fit0))
    looper1 <- function(a) {
      name <- paste("fit", a, sep="_")
      mu_a <- fit$fit0[[name]]$pred.fun(fit$fit0[[name]]$mod, 
                                        X=data.frame(X), tau=tau)
      return(mu_a)
    }
    mu_hat <- lapply(A_lvls, looper1)
    mu_hat <- data.frame(do.call(cbind, mu_hat))
    colnames(mu_hat) <- paste("mu", A_lvls, sep="_")
    # Differences #
    diff_names <- rownames(fit$fitd)
    diff_lvls <- do.call(rbind, strsplit(diff_names, split="_"))
    looper2 <- function(ii) {
      aa <- diff_lvls[ii,2]
      ref <- diff_lvls[ii,3]
      name <- paste("fit", aa, ref, sep="_")
      fitd_ii <- fit$fitd[rownames(fit$fitd)==name,]
      tau_aa <- fitd_ii$fit_aa$pred.fun(fitd_ii$fit_aa$mod, X=data.frame(X), tau=tau)
      tau_ref <- fitd_ii$fit_ref$pred.fun(fitd_ii$fit_ref$mod, X=data.frame(X), tau=tau)
      pi_hat <- data.frame(fit$fitp$pred.fun(fit$fitp$mod, X))
      pi_aa <- as.numeric(pi_hat[[paste("pi", aa, sep="_")]])
      pi_ref <- as.numeric(pi_hat[[paste("pi", ref, sep="_")]])
      tau_hat <- tau_aa*pi_ref + tau_ref*pi_aa
      tau_hat <- data.frame(tau_hat = tau_hat)
      colnames(tau_hat) <- paste("diff", aa, ref, sep="_")
      return(tau_hat)
    }
    tau_hat <- lapply(1:dim(diff_lvls)[1], looper2)
    tau_hat <- data.frame(do.call(cbind, tau_hat))
    # Combine All Estimates #
    pi_hat <- fit$fitp$pred.fun(fit$fitp$mod, X)
    mu_hat <- data.frame(mu_hat, tau_hat, pi_hat)
    return(mu_hat)
  }
  res = list(mod=fit, pred.fun=pred.fun, mu_train=mu_train, mu_test=mu_test)
  class(res) = "X_learner"
  return(res)
}
## Create Pseudo Outcomes (using tmle-update; easier for prediction) ##
.create_pseudo <- function(Y, A, mu_tr, family) {
  
  A_lvls <- unique(A)[order(unique(A))]
  pseudo_dat <- NULL
  
  # Continuous or binary #
  if (family %in% c("gaussian", "binomial")) {
    
    for (aa in A_lvls) {
      aa_mu <- paste("mu", aa, sep="_")
      aa_pi <- paste("pi", aa, sep="_")
      Y_a <- (A==aa)*(Y-mu_tr[[aa_mu]])/(mu_tr[[aa_pi]]) + mu_tr[[aa_mu]]
      pseudo_dat <- cbind(pseudo_dat, Y_a)
    }
    pseudo_dat <- data.frame(pseudo_dat)
    colnames(pseudo_dat) <- paste("mu", A_lvls, sep="_")
    
    eps_star <- NULL
    fluct_dat <- NULL
  }
  
  # Survival (TMLE clever covariate + AFT model average) #
  if (family %in% c("survival")) {
    
    # Create Clever Covariate / E(Y|A=a, X) #
    H_dat <- NULL
    mu_aa <- NA
    for (aa in A_lvls) {
      H_aa <- (A==aa) / mu_tr[[paste("pi", aa, sep="_")]]
      H_dat <- cbind(H_dat, H_aa)
      mu_aa <- ifelse(A==aa, mu_tr[[paste("mu", aa, sep="_")]], mu_aa)
    }
    mu_aa <- log(mu_aa)
    H_dat <- data.frame(H_dat)
    colnames(H_dat) <- paste("H", A_lvls, sep="_")
    Q_tr <- log(mu_tr[,paste("mu", A_lvls, sep="_")])
    
    tmle_adj <- function(Y, mu_aa, H_dat, Q_tr, dist) {
      ftmle <- survreg(Y ~ -1 + offset(mu_aa) + ., data=H_dat, dist = dist)
      eps <- as.numeric(ftmle$coefficients)
      return(list(dist=dist, AIC=AIC(ftmle), eps=eps))
    }
    dist.vec <- c("weibull", "lognormal", "loglogistic")
    out_tmle <- lapply(dist.vec, tmle_adj, Y=Y, mu_aa=mu_aa,
                       H_dat=H_dat, Q_tr=Q_tr)
    out_tmle <- do.call(rbind, out_tmle)
    
    fluct_dat <- data.frame(dist = do.call(rbind, out_tmle[,1]), 
                            AIC = do.call(rbind,  out_tmle[,2]))
    eps_dat <- data.frame(eps = do.call(rbind, out_tmle[,3]))
    colnames(eps_dat) <- paste("eps", A_lvls, sep="_")
    fluct_dat <- data.frame(fluct_dat, eps_dat)
    fluct_dat$wMA <- .weights.MA(fluct_dat$AIC)
    fluct_dat$est_loss <- weighted.mean(fluct_dat$AIC, fluct_dat$wMA)
    eps_star <- colSums(fluct_dat$wMA * eps_dat) 
    
    # Obtain tmle "fluctuated" estimates #
    pseudo_dat <- Q_tr + sweep(as.matrix(H_dat), MARGIN=2, eps_star, `*`)
    
  }
  return(list(pseudo_dat=pseudo_dat, eps_star=eps_star, fluct_dat=fluct_dat))
}
## DR Learner: Fit T-Learner ==> Regress Pseudo-outcomes in second stage ##
DR_learner <- function(Y, A, X, Xtest, family, ple, hyper, pfit, tau=NULL, ...) {
  
  A_lvls <- unique(A)[order(unique(A))]
  
  # Stage 1: T-Learner (treatment-specific models) #
  stage1 <- T_learner(Y=Y, A=A, X=X, Xtest=Xtest,
                      family = family, ple=ple, hyper=hyper, 
                      pfit = pfit, tau = tau)
  
  # Stage 2: Regress Pseudo-Outcomes #
  pseudo_res <- .create_pseudo(Y=Y, A=A, mu_tr=stage1$mu_train, 
                               family=family)
  star_dat <- pseudo_res$pseudo_dat
  
  next_learner <- function(aa) {
    aa_name <- paste("mu", aa, sep="_")
    ref_name <- paste("mu", A_lvls[1], sep="_")
    diff_dr <- star_dat[[aa_name]] - star_dat[[ref_name]]
    fit <- do.call(ple, append(list(Y=diff_dr, X=X, Xtest=Xtest,
                                    family="gaussian"), hyper))
    return(fit)
  }
  fitd <- lapply(A_lvls[-1], next_learner)
  fitd <- do.call(rbind, fitd)
  rownames(fitd) <- paste("fit", A_lvls[-1], A_lvls[1], sep="_")
  
  # Obtain final mu_train / mu_test #
  mu_train <- stage1$mu_train
  mu_test <- stage1$mu_test
  diff_tr <- NULL
  diff_ts <- NULL
  for (ii in rownames(fitd)) {
    diff_name <- gsub("fit", "diff", ii)
    mu_train[diff_name] <- fitd[ii,]$mu_train
    if (!is.null(mu_test)) {
      mu_test[diff_name] <- fitd[ii,]$mu_test
    }
  }
  
  fit <- list(fit0=stage1$mod$fit0, fitd=fitd, fitp=pfit)
  
  # Prediction Function #
  pred.fun <- function(fit, X, tau=NULL) {
    fit0 <- fit$fit0
    fitd <- fit$fitd
    fitp <- fit$fitp
    # Treatment-specific # 
    A_lvls <- sub(".*_", "", names(fit$fit0))
    looper1 <- function(a) {
      name <- paste("fit", a, sep="_")
      mu_a <- fit$fit0[[name]]$pred.fun(fit$fit0[[name]]$mod, 
                                        X=data.frame(X), tau=tau)
      colnames(mu_a) <- paste("mu", a, sep="_")
      return(mu_a)
    }
    mu_hat <- lapply(A_lvls, looper1)
    mu_hat <- do.call(cbind, mu_hat)
    # Differences #
    diff_names <- rownames(fit$fitd)
    diff_lvls <- do.call(rbind, strsplit(diff_names, split="_"))
    looper2 <- function(ii) {
      aa <- diff_lvls[ii,2]
      ref <- diff_lvls[ii,3]
      name <- paste("fit", aa, ref, sep="_")
      fitd_ii <- fit$fitd[rownames(fit$fitd)==name,]
      tau_hat <- fitd_ii$pred.fun(fitd_ii$mod, X=data.frame(X))
      colnames(tau_hat) <- paste("diff", aa, ref, sep="_")
      return(tau_hat)
    }
    tau_hat <- lapply(1:dim(diff_lvls)[1], looper2)
    tau_hat <- data.frame(do.call(cbind, tau_hat))
    # Combine All Estimates #
    pi_hat <- fit$fitp$pred.fun(fit$fitp$mod, X)
    mu_hat <- data.frame(mu_hat, tau_hat, pi_hat)
    return(mu_hat)
  }
  
  res = list(mod=fit, pred.fun=pred.fun, mu_train=mu_train, mu_test=mu_test, 
             fluct_dat = pseudo_res$fluct_dat)
  class(res) = "DR_learner"
  return(res)
}


## Propensity Learner ##
prop_learner <- function(A, X, learner, hyper, ...) {
  if (learner=="ranger" | learner=="ple_ranger") {
    mod <- ranger::ranger(A~., data=data.frame(A, X), probability = TRUE)
    mod$A_lvls <- unique(A)[order(unique(A))]
    pred.fun <- function(mod, X, ...) {
      pi_hat <- predict(mod, X)$predictions
      colnames(pi_hat) <- paste("pi", rev(mod$A_lvls), sep="_")
      return(pi_hat)
    }
  }
  if (learner=="linear" | learner=="ple_linear") {
    mod <- glm(A~., data=data.frame(A, X), family="binomial")
    mod$A_lvls <- unique(A)[order(unique(A))]
    pred.fun <- function(mod, X, ...) {
      pi_hat <- predict(mod, X, type="response")
      pi_hat <- data.frame(pi_hat, 1-pi_hat)
      colnames(pi_hat) <- paste("pi", rev(mod$A_lvls), sep="_")
      return(pi_hat)
    }
  }
  mod <- list(mod=mod, pred.fun=pred.fun)
  pred.fun <- function(mod, X, ...) {
    pi_hat <- mod$pred.fun(mod$mod, X)
    return(pi_hat)
  }
  fit <- list(mod=mod, pred.fun=pred.fun)
  return(fit)
}