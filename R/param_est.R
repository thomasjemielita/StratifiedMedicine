#' Parameter Estimation: Across Subgroups
#'
#' For each identified subgroup, obtain point-estimates and variability metrics
#' (est, SE, CI). fit separate linear regression models. Point-estimates and
#' variability metrics in the overall population are obtained by aggregating subgroup
#' specific results (adaptive weighting or sample size weighting).
#'
#' @inheritParams PRISM
#' @param mu_hat Patient-level estimates (see \code{ple_train})
#' @param Subgrps Identified subgroups. Can be pre-specified, or determined
#' adaptively (see \code{submod_train}).
#' @param alpha_ovrl Two-sided alpha level for overall population
#' @param alpha_s Two-sided alpha level at subgroup
#' @param combine Given identified subgroups and correspond point-estimates/SEs/sample sizes,
#' combine="SS" will use sample size weighting for estimates at the overall level. Not 
#' applicable for param="dr","ple". 
#' @param ... Any additional parameters, not currently passed through.
#'
#' @return Data-set with parameter estimates and corresponding
#' variability metrics, for overall and subgroups. Subgrps="ovrl" corresponds to the overall
#' population by default.
#'  \itemize{
#'   \item param.dat - Parameter estimates and variability metrics (est, SE,
#'   LCL/UCL = lower/upper confidence limits, pval = p-value).
#'   }
#' @export
#' @importFrom survival coxph
#' @importFrom stats AIC
#' @references Funk et al. Doubly Robust Estimation of Causal Effects. 
#' Am J Epidemiol 2011. 173(7): 761-767.
#' @references Andersen, P. and Gill, R. (1982). Cox’s regression model for counting 
#' processes, a large sample study. Annals of Statistics 10, 1100-1120.
#' @references Uno et al. Moving beyond the hazard ratio in quantifying the 
#' between-group difference in survival analysis. Journal of clinical Oncology 2014, 
#' 32, 2380-2385.
#' @examples
#' library(StratifiedMedicine)
#'
#' ## Continuous ##
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#'
#' ## Identify Subgroups: MOB (lmtree) ##
#' res_lmtree = submod_train(Y, A, X, submod="lmtree")
#'
#' ## Parameter-estimation ##
#' param.dat = param_est(Y, A, X, param="lm", Subgrps = res_lmtree$Subgrps.train)
#' param.dat
#'
#' @seealso \code{\link{param_combine}}
#' 
param_est <- function(Y, A, X, param, mu_hat=NULL, Subgrps, 
                     alpha_ovrl=0.05, alpha_s=0.05, combine="SS",...) {
  
  if (param %in% c("lm", "ple", "gcomp", "dr", "rmst", "cox", "aft")) {
    param <- paste("param", param, sep="_")
  }
  numb_subs <- length(unique(Subgrps))
  Subgrps <- as.character(Subgrps)
  type_mu <- ifelse(is.null(dim(mu_hat)), "vec", "df")
  # Estimation across subgroups #
  looper <- function(s, alpha, mu_hat, sname = NULL) {
    Y.s <- Y[Subgrps %in% s]
    A.s <- A[Subgrps %in% s]
    X.s <- X[Subgrps %in% s,]
    if (type_mu=="vec") {
      mu_hat.s <- mu_hat[Subgrps %in% s]
    }
    if (type_mu=="df") {
      mu_hat.s <- mu_hat[Subgrps %in% s,]
    }
    res <- do.call(param, list(Y=Y.s, A=A.s, X=X.s, mu_hat=mu_hat.s, 
                               alpha=alpha))
    if (is.null(sname)) { s_use = s}
    if (!is.null(sname)) {s_use = sname}
    res <- data.frame(Subgrps=s_use, res, alpha=alpha)
    return(res)
  }
  S_levels <- names(table(Subgrps))
  res_s <- lapply(S_levels, looper, alpha=alpha_s, mu_hat=mu_hat)
  param.dat <- do.call(rbind, res_s)
  ## If One Subgroup ##
  if (numb_subs==1) {
    param_ovrl <- param.dat
    param_ovrl$Subgrps <- "ovrl"
    param.dat = rbind(param.dat, param_ovrl)
  }
  if (numb_subs>1) {
    if (param %in% c("param_dr", "param_ple")) {
      param_ovrl <- looper(S_levels, alpha=alpha_ovrl, mu_hat=mu_hat,
                           sname = "ovrl")
    }
    if (!(param %in% c("param_dr", "param_ple"))) {
      param_ovrl <- NULL
      for (e in unique(param.dat$estimand)){
        hold <- param_combine(param.dat = param.dat[param.dat$estimand==e,],
                             alpha=alpha_ovrl, combine=combine)
        hold$estimand <- e
        hold <- hold[,c("Subgrps", "N", "estimand", "est", "SE", 
                       "LCL", "UCL", "pval")]
        param_ovrl <- rbind(param_ovrl, hold)
      }
    }
    param_ovrl$alpha <- alpha_ovrl
    param.dat <- rbind(param_ovrl, param.dat)
  }
  return(param.dat)
}
