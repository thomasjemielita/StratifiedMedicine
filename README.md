
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

Suppose the estimand or question of interest is the average treatment effect, *θ*<sub>0</sub> = *E*(*Y*|*A* = 1)−*E*(*Y*|*A* = 0). The goal is to understand whether there is any treatment heterogeneity across patients and if there are any distinct subgroups with similar responses. In this example, we simulate continuous data where roughly 30% of the patients receive no treatment-benefit for using *A* = 1 vs *A* = 0. Responders vs non-responders are defined by the continuous predictive covariates *X*<sub>1</sub> and *X*<sub>2</sub> for a total of four subgroups. Subgroup treatment effects are: *θ*<sub>1</sub> = 0 (*X*<sub>1</sub> ≤ 0, *X*<sub>2</sub> ≤ 0), *θ*<sub>2</sub> = 0.25(*X*<sub>1</sub> &gt; 0, *X*<sub>2</sub> ≤ 0), *θ*<sub>3</sub> = 0.45(*X*<sub>1</sub> ≤ 0, *X*2 &gt; 0), *θ*<sub>4</sub> = 0.65(*X*<sub>1</sub> &gt; 0, *X*<sub>2</sub> &gt; 0).

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
plot(res0, type="PLE:density") # Density plot of PLEs #
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r
plot(res0, type="PLE:waterfall") # waterfall plot of PLEs
```

<img src="man/figures/README-example-2.png" width="100%" />

``` r
## Plot of the subgroup model (lmtree) ##
plot(res0$submod.fit$mod)
```

<img src="man/figures/README-example-3.png" width="100%" />

``` r
## Overall/subgroup specific parameter estimates/inference
res0$param.dat
#>    Subgrps   N          estimand       est         SE         LCL
#> 1        0 800          E(Y|A=0) 1.6388450 0.05389889  1.53304489
#> 2        0 800          E(Y|A=1) 1.8390738 0.06026975  1.72076804
#> 3        0 800 E(Y|A=1)-E(Y|A=0) 0.2002288 0.07625826  0.05053855
#> 4        3 149          E(Y|A=0) 1.4350696 0.13346377  1.17132883
#> 5        3 149          E(Y|A=1) 1.5989985 0.14895625  1.30464272
#> 6        3 149 E(Y|A=1)-E(Y|A=0) 0.1639289 0.19134782 -0.21419779
#> 7        4 277          E(Y|A=0) 1.6262684 0.07905851  1.47063408
#> 8        4 277          E(Y|A=1) 1.7630294 0.09070726  1.58446338
#> 9        4 277 E(Y|A=1)-E(Y|A=0) 0.1367610 0.11828108 -0.09608672
#> 10       6 267          E(Y|A=0) 1.6788166 0.09052970  1.50057066
#> 11       6 267          E(Y|A=1) 1.8995220 0.09614313  1.71022365
#> 12       6 267 E(Y|A=1)-E(Y|A=0) 0.2207054 0.12863914 -0.03257507
#> 13       7 107          E(Y|A=0) 1.8554233 0.17218365  1.51405244
#> 14       7 107          E(Y|A=1) 2.2194088 0.17545332  1.87155556
#> 15       7 107 E(Y|A=1)-E(Y|A=0) 0.3639856 0.23182371 -0.09562749
#>          UCL          pval
#> 1  1.7446452 1.601252e-135
#> 2  1.9573795 3.482672e-136
#> 3  0.3499190  8.813376e-03
#> 4  1.6988103  2.746867e-20
#> 5  1.8933543  3.061527e-20
#> 6  0.5420557  3.929921e-01
#> 7  1.7819027  1.214684e-57
#> 8  1.9415953  1.303949e-53
#> 9  0.3696087  2.485840e-01
#> 10 1.8570625  7.633166e-50
#> 11 2.0888204  4.277946e-54
#> 12 0.4739859  8.738245e-02
#> 13 2.1967941  9.947318e-19
#> 14 2.5672621  6.568918e-23
#> 15 0.8235986  1.193731e-01
## Forest plot: Overall/subgroup specific parameter estimates (CIs)
plot(res0, type="forest")
```

<img src="man/figures/README-example-4.png" width="100%" />

``` r

## Heatmap of PLEs #
grid.data = expand.grid(X1 = seq(min(X$X1), max(X$X1), by=0.30),
                    X2 = seq(min(X$X2), max(X$X2), by=0.30))
plot(res0, type="heatmap", grid.data = grid.data)
#> $heatmap.est
```

<img src="man/figures/README-example-5.png" width="100%" />

    #> 
    #> $heatmap.prob

<img src="man/figures/README-example-6.png" width="100%" />

Overall, PRISM provides information at the patient-level, the subgroup-level (if any), and the overall population. While there are defaults in place, the user can also input their own functions/model wrappers into the PRISM algorithm. For more details and more examples, we refer the reader to the vignette, [PRISM\_vignette](https://CRAN.R-project.org/package=StratifiedMedicine/vignettes/SM_PRISM.html).
