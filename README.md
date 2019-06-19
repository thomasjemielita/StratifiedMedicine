
<!-- README.md is generated from README.Rmd. Please edit that file -->
StratifiedMedicine
==================

<!-- badges: start -->
<!-- badges: end -->
The goal of StratifiedMedicine is to develop analytic and visualization tools to aid in stratified and personalized medicine. Stratified medicine aims to find subsets or subgroups of patients with similar treatment effects, for example responders vs non-responders, while personalized medicine aims to understand treatment effects at the individual level (does a specific individual respond to the study treatment?). Development of this package is ongoing.

Currently, the main algorithm in this package is “PRISM” (Patient Response Identifiers for Stratified Medicine; Jemielita and Mehrotra 2019 in progress). Given a data-structure of (Y,A,X) (outcome, treatments, covariates), PRISM is a five step procedure:

1.  **Estimand**: Determine the question or estimand of interest. For example, *θ*<sub>0</sub> = *E*(*Y*|*A* = 1)−*E*(*Y*|*A* = 0), where A is a binary treatment variable. While this isn't an explicit step in the PRISM function, the question of interest guides how to set up PRISM.

2.  **Filter (filter)**: Reduce covariate space by removing variables unrelated to outcome/treatment.

3.  **Patient-level estimate (ple)**: Estimate counterfactual patient-level quantities, for example the individual treatment effect, *θ*(*x*)=*E*(*Y*|*X* = *x*, *A* = 1)−*E*(*Y*|*X* = *x*, *A* = 0).

4.  **Subgroup model (submod)**: Partition the data into subsets of patients (likely with similar treatment effects).

5.  **Parameter estimation and inference (param)**: For the overall population and discovered subgroups, output point estimates and variability metrics. These outputs are crucial for Go-No-Go decision making.

Installation
------------

You can install the released version of StratifiedMedicine from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("StratifiedMedicine")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("thomasjemielita/StratifiedMedicine")
```

Example: Continuous Outcome with Binary Treatment
-------------------------------------------------

Suppose the estimand or question of interest is the average treatment effect, *θ*<sub>0</sub> = *E*(*Y*|*A* = 1)−*E*(*Y*|*A* = 0). The goal is to understand whether there is any treatment heterogeneity across patients and if there are any distinct subgroups with similar responses.

``` r
library(StratifiedMedicine)
## basic example code
dat_ctns = generate_subgrp_data(family="gaussian")
Y = dat_ctns$Y
X = dat_ctns$X # 50 covariates, 46 are noise variables, X1 and X2 are truly predictive
A = dat_ctns$A # binary treatment, 1:1 randomized 

# PRISM Default: filter_glmnet, ple_ranger, submod_lmtree, param_ple #
res0 = PRISM(Y=Y, A=A, X=X)
#> Observed Data
#> Filtering: filter_glmnet
#> PLE: ple_ranger
#> Subgroup Identification: submod_lmtree
#> Parameter Estimation: param_ple
## Plot the distribution of PLEs ###
hist(res0$mu_train$PLE, main="Distribution of PLEs",
     xlab = "PLEs: E(Y|X=x, A=1)-E(Y|X=x,A=0)")
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r
## Plot of the subgroup model (lmtree) ##
plot(res0$Sub.mod)
```

<img src="man/figures/README-example-2.png" width="100%" />

``` r
## Overall/subgroup specific parameter estimates/inference
res0$param.dat
#>   Subgrps   N         est         SE          LCL       UCL       pval
#> 1       0 800 0.119163550 0.06277086 -0.004051715 0.2423788 0.05800485
#> 2       3  71 0.269953773 0.20845268 -0.145791990 0.6856995 0.19956283
#> 3       5 303 0.001823749 0.10423555 -0.203296193 0.2069437 0.98605215
#> 4       6 107 0.145690459 0.16108360 -0.173673450 0.4650544 0.36781312
#> 5       8 221 0.124661570 0.11491806 -0.101819583 0.3511427 0.27920349
#> 6       9  98 0.331351436 0.18714952 -0.040088556 0.7027914 0.07978260
## Forest plot: Overall/subgroup specific parameter estimates (CIs)
plot(res0)
```

<img src="man/figures/README-example-3.png" width="100%" />

Overall, PRISM provides information at the patient-level, the subgroup-level (if any), and the overall population. While there are defaults in place, the user can also input their own functions/model wrappers into the PRISM algorithm. For more details and more examples, we refer the reader to the vignette, [PRISM\_vignette](https://cran.r-project.org/web/packages/StratifiedMedicine/vignettes/SM_PRISM.html).
