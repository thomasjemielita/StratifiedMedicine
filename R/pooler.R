#### Pooling Function #####
# @param Y: Outcome
# @param A: Treatment
# @param X: Variables
# @param fit: Initial subgroup model fit (should include objects "Subgrps" and "trt_eff")
# @param pool: Method of pooling. Options include "trteff", "otr:logistic", and "otr:rf".
# @param otr_prob_thres: For OTR pooling, method for choosing probability threshold. 
# Options include "youden" (use Youden index; balances sensitivity/specificity). 
# Can also include numeric threshold, i.e. otr_prob_thres=0.5.
# @param wrapper: Wrapper function. Should output objects "Subgrps" and "trt_eff".
# @param delta: Clinically meaningful threshold (should be character)
# @param alpha: Alpha for CIs
# @param R_pool: Number of resamples for pooling. Used for pool="boot_trteff". 
# @param pool_owl For OTR pooling, if TRUE (default), fit weighted regressions.
# If FALSE, fit unweighted regressions. 
## Output: 
# @trt_assign: Dataframe with Subgrps (initial Subgroups), dopt (optimal treatment assignments), for all subjects
# @trt_eff: Treatment effect estimates (naive or resample based) and dopt assignments.

.pooler <- function(Y, A, X, fit, pool = "trteff",
                   wrapper=NULL, delta="> 0", R_pool=20, 
                   otr_prob_thres = "youden", pool_owl=TRUE,
                   combine = "SS",
                   alpha_ovrl=0.05, alpha_s=0.05) {
  
  # Do we already have dopt assignments? #
  if (!is.null(fit$opt_dat)) {
    
    trt_assign <- data.frame(Subgrps = fit$Subgrps, dopt = fit$dopt)
    
  }
  
  # Initial Subgroup Estimates #
  if (!is.null(fit$trt_eff)) {
    trt_eff_pool <- fit$trt_eff
    trt_eff_pool$Subgrps <- as.character(trt_eff_pool$Subgrps)
  }
  if (is.null(fit$trt_eff)) {
    stop("For pooling, fit object must include trt_eff")
  }
  # Initialize classification data (only used for OTR-based) #
  class_dat <- NULL
 
  # Treatment Effect Pooling #
  if (pool %in% c("trteff", "trteff_boot")) {
    
    if (is.null(trt_eff_pool$dopt)) {
      trt_eff_pool$dopt = eval(parse(text = paste("ifelse(trt_eff_pool$est",
                                                  delta, ", 1, 0)")))
      trt_eff_pool$dopt <- ifelse(trt_eff_pool$dopt==1, "dopt=1", "dopt=0")
    }
    trt_assign <- data.frame(Subgrps=as.character(fit$Subgrps))
    trt_assign <- dplyr::left_join(trt_assign, 
                                   trt_eff_pool[,c("Subgrps", "dopt")], 
                                   by=c("Subgrps"))
    trt_eff_pool$type <- "trteff"
    
  }
  # OTR Pooling #
  if (pool %in% c("otr:logistic", "otr:rf")) {
    
    if (is.null(fit$mu_train)) {
      stop("OTR pooling needs counterfactual estimates (fit$mu_train)")
    }
    if (!requireNamespace("pROC", quietly = TRUE)) {
      stop("Package pROC needed for OTR pooling. Please install.")
    }
    
    res_otr <- .otr_pooler(fit=fit, X=X, pool=pool, delta=delta,
                          pool_owl=pool_owl)
    trt_eff_pool <- dplyr::left_join(trt_eff_pool, 
                                     res_otr$prob_dat, by="Subgrps")
    
    if (otr_prob_thres=="youden") {
      trt_eff_pool$dopt_thres <- unique(trt_eff_pool$thres_youden)
    }
    if (otr_prob_thres=="F1") {
      trt_eff_pool$dopt_thres <- unique(trt_eff_pool$thres_f1)
    }
    if (is.numeric(otr_prob_thres)) {
      trt_eff_pool$dopt_thres <- otr_prob_thres
    }
    trt_eff_pool$dopt <- with(trt_eff_pool, ifelse(prob_dopt <= dopt_thres,
                                                   "dopt=0", "dopt=1"))
    trt_assign <- data.frame(Subgrps=as.character(fit$Subgrps))
    trt_assign <- dplyr::left_join(trt_assign, 
                                   trt_eff_pool[,c("Subgrps", "dopt")],
                                   by=c("Subgrps"))
    trt_eff_pool$type <- pool
  }
  # Last: After pooling, estimate dopt-specific estimates (weighted averages) #
  trt_eff_dopt <- NULL
  for (e in unique(trt_eff_pool$estimand)) {
    holder <- NULL
    for (d in unique(trt_assign$dopt)) {
      sub_d <- trt_eff_pool[trt_eff_pool$dopt==d & trt_eff_pool$Subgrps!="ovrl",]
      out_p <- param_combine(param.dat = sub_d, 
                             combine=combine, 
                             alpha=alpha_s)
      out_p$Subgrps <- d
      holder <- rbind(holder, out_p)
    }
    holder$alpha <- alpha_s
    hold_ovrl <- param_combine(param.dat = holder, 
                               combine=combine, alpha=alpha_ovrl)
    hold_ovrl$alpha <- alpha_ovrl
    hold_all <- rbind(hold_ovrl, holder)
    hold_all$estimand <- e
    trt_eff_dopt <- rbind(trt_eff_dopt, hold_all) 
  }
  # Re-order columns #
  hold_id <- trt_eff_dopt[,colnames(trt_eff_dopt) %in% 
                            c("Subgps", "N", "estimand")]
  hold_other <- trt_eff_dopt[,!(colnames(trt_eff_dopt) %in% 
                               c("Subgps", "N", "estimand"))]
  trt_eff_dopt <- data.frame(hold_id, hold_other)
  
  return(list(trt_assign=trt_assign, trt_eff_pool = trt_eff_pool,
              trt_eff_dopt = trt_eff_dopt))
}

## OTR Pooler Function ##
.otr_pooler <- function(fit, X, delta,
                       pool = "otr:logistic", fast = TRUE,
                       pool_owl = TRUE) {
  
  # Set up data #
  mu_hat <- fit$mu_train
  ple_name <- colnames(mu_hat)[grepl("diff", colnames(mu_hat))]
  ind_ple <- eval(parse(text=paste(paste("ifelse(mu_hat$",ple_name,sep=""),
                                   delta, ", 1, 0)")))
  # ind_ple <- factor(ind_ple, levels = c("0", "1"))
  
  # Threshold as numeric #
  delta.num <- substr(delta, 2, nchar(delta))
  if (suppressWarnings(is.na(as.numeric(delta.num)))) {
    delta.num <- substr(delta, 3, nchar(delta))
  }
  delta.num <- as.numeric(delta.num)
  w_ple <- abs(mu_hat[[ple_name]] - delta.num)
  Subgrps.mat <- data.frame(Subgrps=as.factor(fit$Subgrps))
  Subgrps.mat <- data.frame(model.matrix(~.-1, data=Subgrps.mat))
  otr_dat <- data.frame(ind_ple, Subgrps.mat)
  
  if (pool_owl==FALSE) {
    w_ple <- rep(1, length(w_ple))
  }
  
  # logistic (only subgroups) #
  if (pool=="otr:logistic") {
    mod <- suppressWarnings(glm(ind_ple ~ . -1, 
                                data = otr_dat,
                                family = "binomial", weights = w_ple))
    prob_dopt <- as.numeric(predict(mod, type="response"))
  }
  # ranger: include X again #
  if (pool=="otr:rf") {
    mod <- ranger::ranger(ind_ple ~ ., 
                          data=data.frame(otr_dat, X),
                          case.weights = w_ple)
    sub_hold <- data.frame(Subgrps = fit$Subgrps)
    hold_est <- aggregate(mod$predictions ~ sub_hold$Subgrps, FUN="mean")
    colnames(hold_est) <- c("Subgrps", "prob_dopt")
    prob_dopt <- dplyr::left_join(sub_hold[c("Subgrps")], hold_est,
                                  by="Subgrps")
    prob_dopt <- prob_dopt$prob_dopt
  }
  # Classify Patients: Metric Specific #
  class_dat <- .gen_class_metrics(outcome=ind_ple, pred=prob_dopt)
  thres_youden <- unique(class_dat$yindex)
  
  # Subgroup-Level OTR Probabilities #
  prob_dat <- aggregate(prob_dopt ~ fit$Subgrps, FUN="mean")
  colnames(prob_dat) <- c("Subgrps", "prob_dopt")
  
  # Overall #
  prob_dat0 <- data.frame(Subgrps = "ovrl", prob_dopt = mean(prob_dopt))
  
  prob_dat <- rbind(prob_dat0, prob_dat)
  prob_dat$thres_youden <- thres_youden
  
  if (fast) {
    return(list(otr_mod = NULL, class_dat=class_dat, prob_dat = prob_dat))
  }
  if (!fast) {
    return(list(otr_mod = mod, class_dat=class_dat, prob_dat = prob_dat))
  }
}

## Generate PPV, NPV, Sens, Spec, etc (across cut-points) ##
.gen_class_metrics <- function(outcome, pred) {
  
  rocobj <- suppressMessages(pROC::roc(outcome, pred))
  yindex <- suppressWarnings(pROC::coords(rocobj, "best"))$threshold
  
  delta_vec <- seq(0, 1, by=0.01)
  summ <- NULL
  for (delta in delta_vec) {
    pred_01 <- ifelse(pred>delta, 1, 0)
    pred_01 <- factor(pred_01, levels = c("0", "1"))
    tab <- table(pred_01, outcome, exclude = "no")
    acc <- tab[1,1] + tab[2,2] / length(outcome)
    spec <- tab[1,1] / sum(tab[,1])
    sens <- tab[2,2] / sum(tab[,2])
    npv <- tab[1,1] / sum(tab[1,])
    ppv <- tab[2,2] / sum(tab[2,])
    youden <- sens + spec - 1
    # confusionMatrix(as.factor(pred_01), as.factor(outcome), positive = "1")
    f1 <- 2*(ppv*sens)/(ppv+sens)
    hold <- data.frame(delta=delta, acc=acc, 
                       sens=sens, spec=spec, 
                       ppv=ppv, npv=npv, f1=f1)
    summ <- rbind(summ, hold)
  }
  summ$yindex <- yindex
  return(summ)
}