## Learners ##

### Filter (reduce covariate space) ###
## Elastic Net Filter ##
filter_glmnet = function(Y, A, X, lambda="lambda.min", family="gaussian",
                         interaction=FALSE, ...){
  
  ## Model Matrix #
  fact.vars <- sapply(X, is.factor)
  X.mat <- X
  colnames(X.mat)[fact.vars] <- paste(colnames(X.mat)[fact.vars], "_lvl_", sep="")
  X.mat <- model.matrix(~., data = X.mat)
  X.mat <- X.mat[, colnames(X.mat) != "(Intercept)"]
  W <- X.mat
  intercept <- TRUE
  
  if (interaction){
    A.mat <- model.matrix(~., data=data.frame(A))[,-1]
    X_inter = X.mat*A.mat
    colnames(X_inter) = paste(colnames(X.mat), "_trtA", sep="")
    W = cbind(X.mat, A, X_inter)
  }
  # Center and Scale #
  n <- dim(W)[1]
  W_centered <- apply(W, 2, function(x) x - mean(x))
  Ws <- apply(W_centered, 2, function(x) x / sqrt(sum(x^2) / n))
  
  # Fit Elastic Net #
  if (family=="survival") { family = "cox" }
  mod <- suppressWarnings( 
    glmnet::cv.glmnet(x = Ws, y = Y, nlambda = 100, alpha=0.5, family=family,
              intercept=intercept) )
  
  ### Extract filtered variable based on lambda ###
  VI <- coef(mod, s = lambda)[,1]
  VI <- VI[ names(VI) != "(Intercept)" ]
  # Extract variables that pass the filter ##
  filter.vars <- names(VI[VI!=0])
  filter.vars <- unique( gsub("_lvl_.*","",filter.vars) )
  if (interaction) {
    filter.vars <- unique(sub("_trtA","",filter.vars))
    filter.vars <- filter.vars[filter.vars!="A"]
  }
  # Store selected lambda in model #
  mod$sel.lambda <- lambda
  mod$sel.family <- family
  # Return model fit, filter.vars #
  return(list(mod=mod, filter.vars=filter.vars))
}
# VI Plot #
plot_vimp_glmnet <- function(mod) {
  # Extract filtered variable based on lambda #
  VI <- coef(mod, s = mod$sel.lambda)[,1]
  VI = VI[ names(VI) != "(Intercept)" ]
  
  # Importance Plot #
  plot.title <- paste("Elastic Net (", mod$sel.family, ") Importance Plot", sep="")
  VI.dat <- data.frame(covariate = names(VI),
                       est = as.numeric(VI))
  VI.dat <- VI.dat[VI.dat$est!=0,]
  VI.dat$rank <- 1
  vimp.plt <- ggplot2::ggplot(VI.dat, aes(x=reorder(.data$covariate, abs(.data$est)), 
                                          y=.data$est)) + 
    ggplot2::geom_bar(stat="identity")+
    ggplot2::ylab("Beta (standardized)")+
    ggplot2::xlab("Variable")+
    ggplot2::ggtitle(plot.title)+
    ggplot2::coord_flip()
  return(vimp.plt)
}
## Random Forest based Filtering ##
filter_ranger = function(Y, A, X, b=0.66, K=200, DF2=FALSE, FDR=FALSE, pval.thres=0.10,
                         family="gaussian", ...) {
  
  if (DF2==TRUE) { #Generate the interaction covariates #
    A_lvls <- unique(A)[order(unique(A))]
    ## Generate the interaction covariates ###
    A.mat <- model.matrix(~., data=data.frame(A))[,-1]
    X_inter = X*A.mat
    colnames(X_inter) = paste(colnames(X), A_lvls[2], sep="_")
    W = cbind(A=A.mat, X, X_inter)
    # X_inter = X*A
    # colnames(X_inter) = paste(colnames(X), "_A", sep="")
    # W = cbind(X, A, X_inter)
  }
  if (DF2==FALSE) {
    W = cbind(X)
  }
  n <- dim(W)[1]
  
  ## Calculcate observed variable importance ##
  mod0 <- ranger::ranger(Y ~ ., data = data.frame(Y,W), importance = "permutation")
  VI0 = as.numeric(mod0$variable.importance)
  
  ### Subsample function ###
  b1 = (n)^(b)
  sub_VI = function(s) {
    ## Subsample data ##
    hold = data.frame(Y, W)
    hold = hold[sample(nrow(hold), size=b1),]
    mod.s <- ranger::ranger(Y ~ ., data = hold, importance = "permutation")
    VI = c( as.numeric(mod.s$variable.importance) )
    return(VI)
  }
  res = lapply(1:K, sub_VI)
  res = do.call(rbind, res)
  VI.var = NULL
  for (j in 1:dim(res)[2]){
    VI.var = c(VI.var, b1/((n-b1)*K) * sum((res[,j] - VI0[j])^2, na.rm = TRUE ))
  }
  #### Initial Output of VI and SE(VI) for each variable ###
  out = data.frame(Variables = colnames(W), est = VI0, SE = sqrt(VI.var))
  if (DF2==FALSE) {
    out.F = out
    ### T-Statistics and one-sided p-values ###
    out.F$Tstat = with(out.F, est / SE)
    out.F$pval = with(out.F, 1-pnorm(Tstat))
    out.F$pval.fdr = p.adjust(out.F$pval, "fdr")
  }
  if (DF2==TRUE) {
    out.F = data.frame(Variables=colnames(X), est.2DF=NA, SE.2DF=NA)
    for (var in colnames(X)) {
      inter.name = paste(var, A_lvls[2], sep="_")
      # Extract VI estimates and SEs #
      est0 = out$est[out$Variables==var]; SE0 = out$SE[out$Variable==var]
      estI = out$est[out$Variables==inter.name]; SEI = out$SE[out$Variable==inter.name]
      # Calculate empirical covariance between main/interaction VIs #
      indices = match(c(var, inter.name), colnames(W))
      cov.est = b1/((length(Y)-b1)*K) *
        sum((res[,indices[1]]-est0)*(res[,indices[2]]-estI), na.rm=TRUE)
      # Calculcate "2 DF" VI Estimate and SE #
      est.F = estI+est0
      SE.F = sqrt( SE0^2 + SEI^2 + 2*cov.est )
      # Store #
      out.F$est.2DF[out.F$Variables==var] =  est.F
      out.F$SE.2DF[out.F$Variables==var] =  SE.F
    }
    # T-Statistics and One-sided Pval #
    out.F$Tstat = with(out.F, est.2DF/SE.2DF)
    out.F$pval = with(out.F, 1-pnorm(Tstat))
    out.F$pval.fdr = p.adjust(out.F$pval, "fdr")
    out.F$est <- out.F$est.2DF
    out.F$SE <- out.F$SE.2DF
  }
  
  # FDR Adjustment? #
  if (FDR) {
    out.F$filter_vars = with(out.F, ifelse(pval.fdr<pval.thres, 1, 0))
  }
  if (!FDR) {
    out.F$filter_vars = with(out.F, ifelse(pval<pval.thres, 1, 0))
  }
  out.F$LCL <- with(out.F, est-qnorm(1-pval.thres)*SE)
  out.F$UCL <- with(out.F, est+qnorm(1-pval.thres)*SE)
  filter.vars = with(out.F, as.character(Variables[filter_vars==1]))
  out.F$rank <- with(out.F, rank(pval))
  out.F$pval.thres <- pval.thres
  VI.dat <- out.F
  VI.dat$family <- family
  mod = list(VI.dat=VI.dat)
  # Return VI information and variables that we filter #
  return(list(mod=mod, filter.vars=filter.vars) )
}
# VIMP Plot #
plot_vimp_ranger <- function(mod, top_n=NULL) {
  VI.dat <- mod$VI.dat
  family <- unique(VI.dat$family)
  pval.thres <- unique(VI.dat$pval.thres)
  family <- ifelse(family=="gaussian", "regression", family)
  plot.title <- paste("Random Forest (", family,
                      ") Importance Plot", sep="")
  y.label <- paste("Variable Importance ", 
                   "(", (1-pval.thres)*100, "% CI)",sep="")
  if (!is.null(top_n)) {
    VI.dat <- VI.dat[VI.dat$rank<=top_n,]
    plot.title <- paste(plot.title, " [Top ", top_n, "]", sep="")
  }
  vimp.plt <- ggplot2::ggplot(VI.dat, aes(x=reorder(.data$Variables,-.data$pval), 
                                          y=.data$est, ymin=.data$LCL, ymax=.data$UCL)) +
    ggplot2::geom_errorbar(width=.1)+geom_point()+
    ggplot2::geom_hline(yintercept=0, color="red")+
    ggplot2::ylab(y.label)+
    ggplot2::xlab("Variable")+
    ggplot2::ggtitle(plot.title)+
    ggplot2::coord_flip()
  return(vimp.plt)
}

### PLE (regression + prediction) ###

## Random Forest (Ranger) ##
ple_ranger <- function(Y, X, family="gaussian", min.node.pct=0.10,
                       mtry = NULL, ...) {
  
  if (is.null(mtry)) {
    mtry <- floor(sqrt(dim(X)[2]))
  }
  probability = FALSE
  if (family=="binomial") {
    probability = TRUE
  }
  mod <- ranger::ranger(Y ~ ., data = data.frame(Y, X), 
                        probability = probability, 
                        min.node.size = min.node.pct*dim(X)[1], mtry = mtry)
  mod = list(mod=mod)
  pred.fun = function(mod, X, tau=NULL, ...) {
    treetype = mod[[1]]$treetype
    if (treetype!="Survival") {
      mu_hat = data.frame(mu_hat = predict(mod$mod, X)$predictions)
    }
    if (treetype=="Survival") {
      preds = predict(mod$mod, X)
      looper_rmst <- function(i, surv, time) {
        est.rmst <- rmst_calc(surv = surv[i,],
                              time = time,
                              tau=tau)
        return(est.rmst)
      }
      rmst <- lapply(1:dim(X)[1], looper_rmst, surv=preds$survival,
                     time=preds$unique.death.times)
      mu_hat <- do.call(rbind, rmst)
      mu_hat <- data.frame(mu_hat = mu_hat)
    }
    return(mu_hat)
  }
  res = list(mod=mod, pred.fun=pred.fun)
  return(res)
}
# Linear Model (OLS, GLM, or Cox) #
ple_linear <- function(Y, X, family="gaussian", ...) {
  
  if (family=="gaussian") {
    mod <- lm(Y~., data=data.frame(Y, X))
    pred.fun <- function(mod, X, ...){
      mu_hat <- data.frame(mu_hat = predict(mod, newdata=X))
      return(mu_hat)
    }
  }
  if (family=="binomial") {
    mod <- glm(Y~., data=data.frame(Y, X), family="binomial")
    pred.fun <- function(mod, X, ...){
      mu_hat <- data.frame(mu_hat = predict(mod, newdata=X, type="response"))
      return(mu_hat)
    }
  }
  if (family=="survival") {
    mod <- survival::coxph(Y~., data=data.frame(Y, X))
    pred.fun <- function(mod, X, ...){
      mu_hat <- data.frame(mu_hat = predict(mod, newdata=X))
      return(mu_hat)
    }
  }
  res = list(mod=mod, pred.fun=pred.fun)
  return(res)
}
## GLNET ##
ple_glmnet <- function(Y, X, family="gaussian", 
                        lambda="lambda.min", ...) {
  W = model.matrix(~., data = X)[,-1]
  if (family=="survival") { family = "cox"  }
  mod <- glmnet::cv.glmnet(x = W, y = Y, alpha=0.5, family=family)
  mod$sel.lambda <- lambda
  pred.fun = function(mod, X, ...) {
    lambda <- mod$sel.lambda
    X = model.matrix(~., data = X)[,-1]
    mu_hat = data.frame(mu=as.numeric(predict(mod, newx = X, s=lambda))) 
    return(mu_hat)
  }
  res = list(mod = mod, pred.fun=pred.fun)
}
## BART ##
ple_bart <- function(Y, X, family="gaussian", sparse=FALSE, ...) {
  if (!requireNamespace("BART", quietly = TRUE)) {
    stop("Package BART needed for ple_bart. Please install.")
  }
  if (family=="gaussian") {
    mod <- BART::gbart(x.train = X, y.train = Y, 
                       type = "wbart", sparse = sparse)
    pred.fun <- function(mod, X, ...){
      hold <- predict(mod, newdata=X)
      mu_hat <- data.frame(mu_hat = apply(hold, 2, mean))
      return(mu_hat)
    }
  }
  if (family=="binomial") {
    mod <- BART::gbart(x.train = X, y.train = Y,
                       type = "pbart", sparse = sparse)
    pred.fun <- function(mod, X, ...) {
      hold <- predict(mod, newdata=X)
      mu_hat <- data.frame(mu_hat = hold$prob.test.mean)
      return(mu_hat)
    }
  }
  if (family=="survival") {
    ## TO DO ##
    stop("Survival BART not yet implemented")
  }
  res = list(mod=mod, pred.fun=pred.fun)
  return(res)
}

### submod (Subgroup Identification) ###
## RPART ##
submod_rpart = function(Y, A, X, mu_train, minbucket = floor(dim(X)[1]*0.10),
                        maxdepth = 4, outcome_PLE=FALSE, family="gaussian", ...) {
  
  if (!requireNamespace("rpart", quietly = TRUE)) {
    stop("Package rpart needed for submod_rpart. Please install.")
  }
  
  ## Use PLE as outcome? #
  if (outcome_PLE==TRUE){
    Y <- mu_train$PLE
  }
  ## Fit Model ##
  mod <- rpart::rpart(Y ~ ., data = X,
                      control = 
                        rpart::rpart.control(minbucket=minbucket, maxdepth=maxdepth))
  mod <- as.party(mod)
  # Prediction Function #
  pred.fun <- function(mod, X=NULL, type="subgrp") {
    pred <- NULL
    Subgrps <- as.numeric( predict(mod, type="node", newdata = X) )
    if (type=="all"){
      ## Response Predictions ##
      pred <- data.frame(Subgrps=Subgrps,
                         mu = as.numeric( predict(mod, newdata = X, type="response")))
    }
    return(list(Subgrps=Subgrps, pred=pred))
  }
  ## Return Results ##
  res <- list(mod=mod, pred.fun=pred.fun)
  class(res) <- "submod_rpart"
  return(res)
}
## lmtree (MOB) ##
submod_lmtree = function(Y, A, X, alpha=0.10,
                         minsize = floor(dim(X)[1]*0.10),
                         maxdepth = 4, parm=NULL, ...){
  
  ## Fit Model ##
  mod <- partykit::lmtree(Y~A | ., data = X, alpha=alpha, maxdepth = maxdepth, 
                minsize=minsize, parm=parm)
  
  # Prediction Function #
  pred.fun <- function(mod, X=NULL, type="subgrp"){
    Subgrps <- NULL; pred <- NULL
    Subgrps <- as.numeric( predict(mod, type="node", newdata = X) )
    if (type=="all"){
      L.mat <- rbind( c(1,0), c(1,1) )
      pred <- data.frame(Subgrps=Subgrps, mu0 = NA, mu1 = NA)
      for (s in unique(Subgrps)){
        hold <- suppressWarnings(  summary(mod)[[as.character(s)]] )
        hold.est <- as.numeric( L.mat %*% coef(hold)  ) 
        pred$mu0[pred$Subgrps==s] <- hold.est[1]
        pred$mu1[pred$Subgrps==s] <- hold.est[2]
      }
      pred$PLE <- with(pred, mu1-mu0) 
    }
    return(list(Subgrps=Subgrps, pred=pred))
  }
  res <- list(mod=mod, pred.fun=pred.fun)
  class(res) <- "submod_lmtree"
  ## Return Results ##
  return(res)
}
## glmtree (MOB) ##
submod_glmtree = function(Y, A, X, 
                          glm.fam = binomial, link="identity", alpha=0.10,
                          minsize = floor(dim(X)[1]*0.10),
                          maxdepth = 4, parm=NULL, ...){
  
  ## Fit Model ##
  mod <- partykit::glmtree(Y~A | ., data = X, family= glm.fam(link=link), 
                 alpha=alpha, maxdepth = maxdepth, 
                 minsize=minsize, parm=parm)
  
  # Prediction Function #
  pred.fun <- function(mod, X=NULL, type="subgrp"){
    Subgrps <- NULL; pred <- NULL
    Subgrps <- as.numeric( predict(mod, type="node", newdata = X) )
    if (type=="all"){
      L.mat <- rbind( c(1,0), c(1,1) )
      pred <- data.frame(Subgrps=Subgrps, mu0 = NA, mu1 = NA)
      for (s in unique(Subgrps)){
        hold <- suppressWarnings(  summary(mod)[[as.character(s)]] )
        hold.est <- as.numeric( L.mat %*% coef(hold)  ) 
        pred$mu0[pred$Subgrps==s] <- hold.est[1]
        pred$mu1[pred$Subgrps==s] <- hold.est[2]
      }
      pred$PLE <- with(pred, mu1-mu0) 
    }
    return(list(Subgrps=Subgrps, pred=pred))
  }
  res <- list(mod=mod, pred.fun=pred.fun)
  class(res) <- "submod_glmtree"
  ## Return Results ##
  return(  res )
}
## Conditional Inference Trees ##
submod_ctree = function(Y, A, X, mu_train, alpha=0.10,
                        minbucket = floor(dim(X)[1]*0.10),
                        maxdepth = 4, outcome_PLE=FALSE, ...) {
  
  ## Use PLE as outcome? ##
  if (outcome_PLE==TRUE){
    Y <- mu_train$PLE
  }
  # Fit Model #
  mod <- partykit::ctree(Y ~ ., data = X,
               control = partykit::ctree_control(alpha=alpha, 
                                       minbucket=minbucket, maxdepth=maxdepth))
  
  # Prediction Function #
  pred.fun <- function(mod, X=NULL, type="subgrp"){
    Subgrps <- NULL
    pred <- NULL
    Subgrps <- as.numeric( predict(mod, type="node", newdata = X) )
    if (type=="all"){
      pred <- data.frame(Subgrps=Subgrps,
                         mu=as.numeric( predict(mod, newdata = X, type="response")) )
    }
    return( list(Subgrps=Subgrps, pred=pred) )
  }
  res <- list(mod=mod, pred.fun=pred.fun)
  class(res) = "submod_ctree"
  ## Return Results ##
  return(  res )
}
## OTR: I(PLE>thres) ~ X, weights = abs(PLE), with CTREE ###
submod_otr = function(Y, A, X, mu_train, alpha=0.10,
                      minbucket = floor(dim(X)[1]*0.10),
                      maxdepth = 4, thres=">0", ...){
  # Identify treatment difference name #
  ple_name <- colnames(mu_train)[grepl("diff", colnames(mu_train))]
  ple_name <- ple_name[1]
  
  ## Set up data ##
  mu_train$PLE <- mu_train[[ple_name]]
  print(summary(mu_train$PLE))
  ind_PLE <- eval(parse(text=paste("ifelse(mu_train$PLE", thres, ", 1, 0)")))
  w_PLE <- abs(mu_train$PLE)
  hold <- data.frame(ind_PLE, X)
  # Fit Model #
  mod <- suppressWarnings(partykit::ctree(ind_PLE ~ ., data = hold, weights = w_PLE,
                                 control = partykit::ctree_control(alpha=alpha,
                                                         minbucket=minbucket,
                                                         maxdepth=maxdepth)))
  # Prediction Function #
  pred.fun <- function(mod, X=NULL, type="subgrp"){
    Subgrps <- NULL; pred <- NULL;
    Subgrps = as.numeric( predict(mod, type="node", newdata = X) )
    if (type=="all"){
      pred = data.frame(Subgrps=Subgrps,
                        as.numeric( predict(mod, newdata = X, type="response")))
    }
    return( list(Subgrps=Subgrps, pred=pred) )
  }
  ## Return Results ##
  res <- list(mod=mod, pred.fun=pred.fun)
  class(res) <- "submod_otr"
  return(res)
}
## MOB: Weibull ##
submod_mob_weib = function(Y, A, X, alpha=0.10,
                          minsize = floor(dim(X)[1]*0.10),
                          maxdepth = 4, parm=NULL, ...) {
  
  ## Fit Model ##
  mod <- partykit::mob(Y ~ A | ., data = X,
             fit = wbreg, control = partykit::mob_control(alpha=alpha, minsize=minsize,
                                                maxdepth=maxdepth, parm=parm))
  # Prediction Function #
  pred.fun <- function(mod, X=NULL, type="subgrp"){
    pred <- NULL
    Subgrps <- as.numeric( predict(mod, type="node", newdata = X) )
    if (type=="all"){
      L.mat <- rbind( c(1,0), c(1,1) )
      pred <- data.frame(Subgrps=Subgrps, mu0 = NA, mu1 = NA)
      for (s in unique(Subgrps)){
        hold <- summary(mod)[[as.character(s)]]
        hold.est <- exp( as.numeric( L.mat %*% coef(hold)  ) )
        pred$mu0[pred$Subgrps==s] <- hold.est[1]
        pred$mu1[pred$Subgrps==s] <- hold.est[2]
      } 
    }
    return(list(Subgrps=Subgrps, pred=pred))
  }
  res <- list(mod=mod, pred.fun=pred.fun)
  class(res) <- "submod_weibull"
  ## Return Results ##
  return(res)
}

## Parameter Models ##
## Linear Regression ##
param_lm <- function(Y, A, alpha, ...) {
  
  n <- length(Y)
  if (is.null(A)) {
    indata <- data.frame(Y)
    lm.mod = tryCatch(lm(Y ~ 1 , data=indata), error = function(e) "param error")
    est = summary(lm.mod)$coefficients[1,1]
    SE = summary(lm.mod)$coefficients[1,2]
    LCL = est - qt(1-alpha/2, df=n-1)*SE
    UCL = est + qt(1-alpha/2, df=n-1)*SE
    pval = 2*pt(-abs(est/SE), df=n-1)
    summ = data.frame( N = n, estimand = "mu", est, SE, LCL, UCL, pval)
  }
  if (!is.null(A)) {
    A_lvls <- unique(A)[order(unique(A))]
    estimands <- c(paste("mu", A=A_lvls[1], sep="_"),
                   paste("mu", A=A_lvls[2], sep="_"))
    estimands <- c(estimands, paste(estimands[2], "-", estimands[1], sep=""))
    indata <- data.frame(Y, A)
    lm.mod = tryCatch(lm(Y ~ A , data=indata), error = function(e) "param error")
    L.mat = rbind( c(1,0), c(1,1), c(0,1) )
    est = L.mat %*% coef(lm.mod)
    SE = sqrt(  diag( L.mat %*% vcov(lm.mod) %*% t(L.mat) ) )
    LCL = est - qt(1-alpha/2, df=n-1)*SE
    UCL = est + qt(1-alpha/2, df=n-1)*SE
    pval = 2*pt(-abs(est/SE), df=n-1)
    summ = data.frame(N = n, estimand = estimands, 
                      est, SE, LCL, UCL, pval) 
  }
   return(summ)                 
}
## Doubly-Robust ##
param_dr <- function(Y, A, mu_hat, alpha, ...) {
  
  n <- length(Y)
  if (is.null(A)) {
    stop("DR parameter estimation not usable for A=NULL")
  }
  if (!is.null(A)) {
    A_lvls <- unique(A)[order(unique(A))]
    estimands <- c(paste("mu", A=A_lvls[1], sep="_"),
                   paste("mu", A=A_lvls[2], sep="_"))
    estimands <- c(estimands, paste(estimands[2], "-", estimands[1], sep=""))
    mu_A0 = estimands[1]
    mu_A1 = estimands[2]
    A_ind <- ifelse(A==A_lvls[2], 1, 0)
    pi_hat <- mu_hat[[paste("pi", A_lvls[2], sep="_")]]
    # EIF ###
    eif.0 = ( (1-A_ind)*Y + (A_ind-pi_hat)*mu_hat[,mu_A0] ) / (1-pi_hat)
    eif.1 = ( A_ind*Y - (A_ind-pi_hat)*mu_hat[,mu_A1])/ pi_hat
    eif = eif.1 - eif.0
    # Double robust estimator: Average eifs #
    est = c( mean(eif.0, na.rm=TRUE), 
             mean(eif.1, na.rm=TRUE),
             mean(eif, na.rm=TRUE) )
    # EIF for variance estimate #
    SE = sqrt( n^(-2) * c( sum( (eif.0-est[1])^2 ),
                             sum( (eif.1-est[2])^2 ),
                             sum( (eif-est[3])^2 ) ) )
    LCL = est-qt( (1-alpha/2), df=n-1 )*SE
    UCL = est+qt( (1-alpha/2), df=n-1 )*SE
    pval = 2*pt(-abs(est/SE), df=n-1)
    summ = data.frame(N = n, estimand = estimands,
                      est, SE, LCL, UCL, pval)
  }
  return(summ)                 
}
## Average patient-level estimates ##
param_ple <- function(Y, A, mu_hat, alpha, ...) {
  
  n <- length(Y)
  if (is.null(A)) {
    est = mean(mu_hat[,1])
    SE = sqrt( n^(-2)*sum((est-Y)^2) ) 
    LCL = est-qt( (1-alpha/2), df=n-1 )*SE
    UCL = est+qt( (1-alpha/2), df=n-1 )*SE
    pval = 2*pt(-abs(est/SE), df=n-1)
    summ = data.frame(N = n, estimand = "E(Y)", est, SE, LCL, UCL, pval)
  }
  if (!is.null(A)) {
    A_lvls <- unique(A)[order(unique(A))]
    estimands <- c(paste("mu", A=A_lvls[1], sep="_"),
                   paste("mu", A=A_lvls[2], sep="_"))
    estimands <- c(estimands, paste(estimands[2], "-", estimands[1], sep=""))
    mu_A0 = estimands[1]
    mu_A1 = estimands[2]
    A_ind <- ifelse(A==A_lvls[2], 1, 0)
    pi_hat <- mu_hat[[paste("pi", A_lvls[2], sep="_")]]
    # EIF ###
    eif.0 = ( (1-A_ind)*Y + (A_ind-pi_hat)*mu_hat[,mu_A0] ) / (1-pi_hat)
    eif.1 = ( A_ind*Y - (A_ind-pi_hat)*mu_hat[,mu_A1])/ pi_hat
    eif = eif.1 - eif.0
    # Double robust estimator: Average eifs #
    est = c( mean(eif.0, na.rm=TRUE), 
             mean(eif.1, na.rm=TRUE),
             mean(mu_hat[,mu_A1]-mu_hat[,mu_A0], na.rm=TRUE))
    # EIF for variance estimate #
    SE = sqrt( n^(-2) * c( sum( (eif.0-est[1])^2 ),
                           sum( (eif.1-est[2])^2 ),
                           sum( (eif-est[3])^2 ) ) )
    LCL = est-qt( (1-alpha/2), df=n-1 )*SE
    UCL = est+qt( (1-alpha/2), df=n-1 )*SE
    pval = 2*pt(-abs(est/SE), df=n-1)
    summ = data.frame(N = n, estimand = estimands,
                      est, SE, LCL, UCL, pval)
  }
  return(summ)                 
}
## Cox Regression ##
param_cox <- function(Y, A, alpha, ...) {
  
  n <- length(Y[,1])
  if (is.null(A)) {
    stop("Cox regression parameter estimation not applicable for A=NULL")
  }
  if (!is.null(A)) {
    A_lvls <- unique(A)[order(unique(A))]
    estimands <- c(paste("logHR", A=A_lvls[1], sep="_"),
                   paste("logHR", A=A_lvls[2], sep="_"))
    estimands <- paste(estimands[2], "-", estimands[1], sep="")
    indata = data.frame(Y=Y, A=A)
    cox.mod = tryCatch( survival::coxph(Y ~ A , data=indata),
                        error = function(e) "fit error",
                        warning = function(w) "convergence issues")
    if (is.character(cox.mod)) {
      summ = data.frame(N=n,est=NA, SE=NA, LCL=NA, UCL=NA, pval=NA)
    }
    if (is.list(cox.mod)) {
      est = summary(cox.mod)$coefficients[1]
      SE = summary(cox.mod)$coefficients[3]
      LCL = confint(cox.mod, level=1-alpha)[1]
      UCL = confint(cox.mod, level=1-alpha)[2]
      pval = summary(cox.mod)$coefficients[5]
      summ = data.frame(N = n, estimand = estimands, est, SE, LCL, UCL, pval)
    }
  }
  return(summ)                 
}
## RMST ##
param_rmst <- function(Y, A, alpha, ...) {
  
  if (!requireNamespace("survRM2", quietly = TRUE)) {
    stop("Package survRM2 needed for RMST parameter estimation. Please install.")
  }
  n <- length(Y[,1])
  time = Y[,1]
  status = Y[,2]
  if (is.null(A)) {
    obj = tryCatch( rmst_single(time, status, tau=NULL),
                    error = function(e) "param error" )
    if (is.character(obj)){
      est = NA; SE = NA; pval = NA; LCL = NA; UCL = NA;
    }
    if (is.list(obj)){
      est = obj$rmst
      SE = obj$rmst.se
      LCL =  est - qnorm(1-alpha/2)*SE
      UCL =  est + qnorm(1-alpha/2)*SE
      pval = 2*pt(-abs(est/SE), df=n-1)
      summ = data.frame(N = n, est, SE, LCL, UCL, pval)
    }
  }
  if (!is.null(A)) {
    A_lvls <- unique(A)[order(unique(A))]
    estimands <- c(paste("mu", A=A_lvls[1], sep="_"),
                   paste("mu", A=A_lvls[2], sep="_"))
    estimands <- paste(estimands[2], "-", estimands[1], sep="")
    arm = A
    obj = tryCatch( survRM2::rmst2(time, status, arm),
                    error = function(e) "param error" ) 
    if (is.character(obj)){
      est = NA; SE = NA; pval = NA; LCL = NA; UCL = NA;
    }
    if (is.list(obj)){
      est = obj$unadjusted.result[1,1]
      SE = sqrt( obj$RMST.arm1$rmst.var + obj$RMST.arm0$rmst.var )
      LCL =  est - qnorm(1-alpha/2)*SE
      UCL =  est + qnorm(1-alpha/2)*SE
      pval = 2*pt(-abs(est/SE), df=n-1)
      summ = data.frame(N = n, estimand=estimands, est, SE, LCL, UCL, pval)
    }
  }
  return(summ)                 
}