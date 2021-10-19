#' Overall Population Estimate: Aggregating Subgroup-Specific Parameter Estimates
#'
#' Function that combines subgroup-specific estimates to obtain an overall population
#' estimate. Options including sample size weighting and max Z weighting
#'
#' @param param.dat Parameter data-set with subgroup-specific point estimates, SEs, and
#' sample sizes.
#' @param combine Method to combine subgroup-specific estimates. Default is "SS",
#' or sample size weighting. Another option is "maxZ"(see Mehrotra
#' and Marceau-West 2021).
#' @param alpha Two-sided alpha level for overall population. Default=0.05
#' @param ... Any additional parameters, not currently passed through.
#'
#' @importFrom stats cor dnorm integrate uniroot
#' @return Data-frame with point estimate, SE, and CI
#' @export

param_combine = function(param.dat, combine="SS", alpha=0.05, ...){

  # Sample Size Weighting #
  if (combine=="SS") {
    w.vec <- param.dat$N / sum(param.dat$N)
    beta0 = as.numeric( w.vec %*% param.dat$est )
    SE.beta0 = sqrt(  sum(w.vec^2*param.dat$SE^2)    )
    beta0.LCL = beta0 - qt(1-alpha/2, df=sum(param.dat$N)-1) * SE.beta0
    beta0.UCL = beta0 + qt(1-alpha/2, df=sum(param.dat$N)-1) * SE.beta0
    param.dat0 = data.frame(Subgrps="ovrl", N = sum(param.dat$N), 
                            est = beta0, SE = SE.beta0,
                            LCL = beta0.LCL, UCL = beta0.UCL)
    param.dat0$pval = with(param.dat0, 2*pt(-abs(est/SE), df=N-1) )
    param.dat0$alpha <- alpha
  }
  # Zmax #
  if (combine=="maxZ") {
    
    resz <- zmax_helper(param.dat = param.dat, alpha=alpha)
    
    param.dat0 = data.frame(Subgrps="ovrl", N = sum(param.dat$N), 
                            resz)
  }
  
  return(param.dat0)
}

# Helper Function for Combining Groups Zmax 
# See: 5-STAR, Mehrotra and Marceau-West 2021 #
zmax_helper <- function(param.dat, alternative="two-sided",
                        alpha=0.05) {
  
  est <- param.dat$est
  SE <- param.dat$SE
  N <- param.dat$N
  Subgrps <- param.dat$Subgrps
  numb_subs <- length(Subgrps)
  
  # Test Stats #
  Z_I <- sum(N*est) / sqrt(sum(N^2 * SE^2))
  
  Z_II <- sum(N*(est / SE)) / sqrt(sum(N^2))
  
  # Calculate p-values #
  if (alternative == "less") {
    p1 = pnorm(Z_I,lower.tail=TRUE)
    p2 = pnorm(Z_II,lower.tail=TRUE)
  }
  if (alternative == "greater") {
    p1 = pnorm(Z_I,lower.tail=FALSE)
    p2 = pnorm(Z_II,lower.tail=FALSE)
  } 
  if (alternative == "two-sided") {
    p1 = 2*pnorm(abs(Z_I),lower.tail=FALSE)
    p2 = 2*pnorm(abs(Z_II),lower.tail=FALSE)
  }
  
  minP = min(p1, p2)
  Zstar <- max(Z_I, Z_II)
  
  # If only one subgroup #
  if (numb_subs==1) {
    summ <- data.frame(est=est, SE=SE)
    summ$LCL <- with(summ, est - qnorm(1-alpha/2)*SE)
    summ$UCL <- with(summ, est + qnorm(1-alpha/2)*SE)
    summ$pval <- minP
    return(summ)
  }
  
  # If >1 subgroups #
  if (numb_subs>1) {
   
    maxz_pdf = function(z, rho) {
      fz = 2*dnorm(z)*pnorm((z*(rho-1))/sqrt(1-rho^2))
      return(fz)
    }
    
    rho <- sum(N^2 * SE) /
      (sqrt(sum(N^2 * SE^2))*sqrt(sum(N^2)))
    if (rho==1) {
      rho=0.999
    }
    if (Zstar==Z_I) {
      w.vec <- N 
    }
    if (Zstar==Z_II) {
      w.vec <- N / SE
    }
    est0 <- sum(w.vec * est) / sum(w.vec)
    SE0 = sqrt(sum(w.vec^2*SE^2) / sum(w.vec)^2)
    
    if (alternative=="less") {
      pstar <- integrate(function(x) maxz_pdf(x,rho),lower=-Inf,upper=Zstar)$value
      calpha = uniroot(function(x) integrate(function(z)
        maxz_pdf(z,rho),lower=-Inf,upper=x)$value - alpha, lower=-10,
        upper=0)$root
    }
    if (alternative=="greater") {
      pstar <- integrate(function(x) maxz_pdf(-x,rho),lower=Zstar,upper=Inf)$value
      calpha = uniroot(function(x) integrate(function(z)
        maxz_pdf(-z,rho),lower=x,upper=Inf)$value - alpha, lower=0,
        upper=10)$root
    }
    if (alternative=="two-sided") {
      pstar <- 2*integrate(function(x) maxz_pdf(-x,rho),
                           lower=abs(Zstar),upper=Inf)$value
      calpha = uniroot(function(x) integrate(function(z)
        maxz_pdf(-z,rho),lower=x,upper=Inf)$value - alpha/2, lower=0,
        upper=10)$root
    } 
    summ <- data.frame(est=est0, SE=SE0)
    summ$LCL <- with(summ, est - calpha*SE)
    summ$UCL <- with(summ, est + calpha*SE)
    summ$pval <- pstar
    return(summ)
  }
}
