# PRISM (Train): Patient Response Identifier for Stratified Medicine
PRISM_train = function(Y, A, X, Xtest=NULL, family="gaussian",
                 filter="filter_glmnet", ple=NULL, submod=NULL, param=NULL,
                 alpha_ovrl=0.05, alpha_s = 0.05,
                 filter.hyper=NULL, ple.hyper=NULL, submod.hyper = NULL,
                 param.hyper = NULL, verbose=TRUE) {
  
  ### Step 1: Variable Filtering ###
  if (!(filter=="None")) {
    if (verbose) message(paste("Filtering:", filter, sep=" "))
    step1 <- tryCatch(filter_train(Y=Y, A=A, X=X, family = family, 
                                   filter=filter, hyper = filter.hyper),
                       error = function(e) e )
    filter.mod <- step1$mod
    filter.vars <- step1$filter.vars
  }
  if ( filter=="None" ){
    filter.mod <- NULL; filter.vars <- NULL;
  }
  # Drop variables depending on filter #
  if ( filter=="None" ){ X.star <- X; Xtest.star <- Xtest }
  if ( !(filter=="None") ){
    if (length(filter.vars)==0){ X.star <- X; Xtest.star <- Xtest  }
    if (length(filter.vars) >0){
      X.star <- X[, colnames(X) %in% c(filter.vars, "A", "Y"), drop=FALSE]
      Xtest.star <- Xtest[, colnames(Xtest) %in% c(filter.vars, "A", "Y"), drop=FALSE] }
  }
  ### Step 2: PLE ###
  if ( !(ple=="None") ){
    if (verbose) message( paste("PLE:", ple, sep=" " ) )
    step2 <- ple_train(Y=Y,A=A,X=X.star,Xtest=Xtest.star,family=family, ple=ple,
                      hyper = ple.hyper)
    ple.fit <- step2$fit
    mu_train <- step2$mu_train
    mu_test <- step2$mu_test
  }
  if ( ple=="None" ){
    ple.fit <- NULL
    mu_train <- NULL
    mu_test <- NULL
  }
  ### Step 3: Subgroup Identification ###
  if ( !(submod=="None") ){
    if (verbose) message( paste("Subgroup Identification:",
                                submod, sep=" "))
    step3 <- tryCatch( submod_train(Y=Y, A=A, X=X.star, Xtest=Xtest.star, 
                                   mu_train=mu_train, family = family, 
                                   submod=submod, hyper = submod.hyper),
                      error = function(e) "submod error" )
    if (is.character(step3)){
      submod.fit <- step3
      Rules <- step3
      Subgrps.train <- 1; Subgrps.test <- 1
    }
    if (is.list(step3)){
      submod.fit <- step3$fit
      Rules <- step3$Rules;
      Subgrps.train <- step3$Subgrps.train; Subgrps.test = step3$Subgrps.test
    }
  }
  if ( submod=="None" ){
    submod.fit <- NULL; Rules <- NULL;
    Subgrps.train <- NULL; Subgrps.test <- NULL;
  }
  ### Step 4: Parameter Estimation and Inference ###
  if (verbose){ message(paste("Parameter Estimation:", param, sep=" ")) }
  param.dat <- tryCatch( do.call( param, list(Y=Y, A=A, X=X.star, mu_hat = mu_train,
                                             Subgrps=Subgrps.train,
                                             alpha_ovrl=alpha_ovrl,
                                             alpha_s=alpha_s)  ),
                        error = function(e) "param error" )
  if (!is.character(param.dat)){
    if (is.null(param.dat$estimand)){
      param.dat$estimand <- "est"
    }
    if (!is.character(param.dat)){
      param.dat$alpha = ifelse(param.dat$Subgrps==0, alpha_ovrl, alpha_s)
    }
  }
  
  output <- list( mu_train = mu_train, mu_test = mu_test, filter.mod = filter.mod,
                 filter.vars = filter.vars, ple.fit = ple.fit, submod.fit = submod.fit,
                 Subgrps.train = Subgrps.train, Subgrps.test=Subgrps.test,
                 Rules = Rules, param.dat=param.dat)
  ### Return Outputs ###
  return( output )
}