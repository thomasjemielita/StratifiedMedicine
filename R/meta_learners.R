## T-Learner: Fit separate regression model by exposure level (A) ##
T_learner <- function(Y, A, X, family, ple, hyper, pfit, tau=NULL, ...) {
  A_lvls <- unique(A)[order(unique(A))]
  fit_learner <- function(a) {
    fit <- do.call(ple, append(list(Y=Y[A==a], X=X[A==a,], 
                                    family=family), hyper))
    return(list(fit=fit))
  }
  fit <- lapply(A_lvls, fit_learner)
  fit <- do.call(cbind, fit)
  names(fit) <- paste("fit", A_lvls, sep="_")
  fit <- list(fit0=fit, fitp=pfit)
  # Prediction Function #
  pred.fun <- function(fit, X, tau=NULL) {
    fit0 <- fit$fit0
    A_lvls <- sub(".*_", "", names(fit0))
    looper <- function(a) {
      name <- paste("fit", a, sep="_")
      mu_a <- fit0[[name]]$pred.fun(fit0[[name]]$mod, X=data.frame(X), tau=tau)
      colnames(mu_a) <- paste("mu", a, sep="_")
      return(mu_a)
    }
    mu_hat <- lapply(A_lvls, looper)
    mu_hat <- do.call(cbind, mu_hat)
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
  res = list(mod=fit, pred.fun=pred.fun)
  class(res) = "T_learner"
  return(res)
}
## S-Learner: Fit single Regression Model (can do forced interactions, VT) ##
S_learner <- function(Y, A, X, VT=FALSE, family, ple, hyper, pfit,
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
        colnames(mu_a) <- paste("mu", a, sep="_")
        return(mu_a)
      }
      mu_hat <- lapply(A_lvls, looper)
      mu_hat <- do.call(cbind, mu_hat)
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
X_learner <- function(Y, A, X, family, ple, hyper, pfit, tau=NULL, ...) {
  A_lvls <- unique(A)[order(unique(A))]
  # Step 1: T-Learner #
  init_learner <- function(a) {
    fit <- do.call(ple, append(list(Y=Y[A==a], X=X[A==a,], 
                                    family=family), hyper))
    mu_a <- fit$pred.fun(fit$mod, X, tau=tau)
    return(list(fit=fit, mu_a=mu_a))
  }
  fit <- lapply(A_lvls, init_learner)
  fit <- do.call(cbind, fit)
  mu_hat <- do.call(cbind, fit[2,])
  fit0 <- fit[1,]
  colnames(mu_hat) <- paste("mu", A_lvls, sep="_")
  names(fit0) <- paste("fit", A_lvls, sep="_")
  # Step 2: Imputed Treatment Effects #
  next_learner <- function(aa) {
    aa_name <- paste("mu", aa, sep="_")
    ref_name <- paste("mu", A_lvls[1], sep="_")
    if (is.Surv(Y)) { Y <- Y[,1]; }
    imp_aa <- (Y[A==aa]-mu_hat[ref_name])[,1]
    imp_ref <- (mu_hat[aa_name]-Y[A==A_lvls[1]])[,1]
    fit_aa <- do.call(ple, append(list(Y=imp_aa, X=X, 
                                       family="gaussian"), hyper))
    fit_ref <- do.call(ple, append(list(Y=imp_ref, X=X, 
                                       family="gaussian"), hyper))
    fit <- list(fit_aa=fit_aa, fit_ref=fit_ref)
    return(fit)
  }
  fitd <- lapply(A_lvls[-1], next_learner)
  fitd <- do.call(rbind, fitd)
  rownames(fitd) <- paste("fit", A_lvls[-1], A_lvls[1], sep="_")
  fit <- list(fit0=fit0, fitd=fitd, fitp=pfit)
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
      tau_aa <- fitd_ii$fit_aa$pred.fun(fitd_ii$fit_aa$mod, X=data.frame(X), tau=tau)
      tau_ref <- fitd_ii$fit_ref$pred.fun(fitd_ii$fit_ref$mod, X=data.frame(X), tau=tau)
      pi_hat <- data.frame(fit$fitp$pred.fun(fit$fitp$mod, X))
      pi_aa <- as.numeric(pi_hat[[paste("pi", aa, sep="_")]])
      pi_ref <- as.numeric(pi_hat[[paste("pi", ref, sep="_")]])
      tau_hat <- tau_aa*pi_ref + tau_ref*pi_aa
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
  res = list(mod=fit, pred.fun=pred.fun)
  class(res) = "X_learner"
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