
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
plot(res0) # default: tree plot 
```

![](man/figures/README-example-1.png)

``` r
plot(res0, type="PLE:waterfall") # waterfall plot of PLEs
```

![](man/figures/README-example-2.png)

``` r
## Overall/subgroup specific parameter estimates/inference
res0$param.dat
#>    Subgrps   N          estimand        est         SE         LCL
#> 1        0 800          E(Y|A=0) 1.63979818 0.04567662  1.55013782
#> 2        0 800          E(Y|A=1) 1.83953825 0.04871379  1.74391613
#> 3        0 800 E(Y|A=1)-E(Y|A=0) 0.19974008 0.06329221  0.07550143
#> 4        3 149          E(Y|A=0) 1.28297958 0.11169103  1.06226442
#> 5        3 149          E(Y|A=1) 1.31807996 0.10932336  1.10204360
#> 6        3 149 E(Y|A=1)-E(Y|A=0) 0.03510038 0.15343734 -0.26811059
#> 7        4 277          E(Y|A=0) 1.61280771 0.06705679  1.48079995
#> 8        4 277          E(Y|A=1) 1.67336273 0.07708837  1.52160684
#> 9        4 277 E(Y|A=1)-E(Y|A=0) 0.06055502 0.10098218 -0.13823814
#> 10       7  99          E(Y|A=0) 1.58968017 0.14067817  1.31050892
#> 11       7  99          E(Y|A=1) 1.92148579 0.13161892  1.66029233
#> 12       7  99 E(Y|A=1)-E(Y|A=0) 0.33180562 0.19007951 -0.04540098
#> 13       8 168          E(Y|A=0) 1.76688140 0.08927908  1.59062031
#> 14       8 168          E(Y|A=1) 2.04137779 0.09774962  1.84839357
#> 15       8 168 E(Y|A=1)-E(Y|A=0) 0.27449640 0.13104528  0.01577751
#> 16       9 107          E(Y|A=0) 2.05338724 0.13530977  1.78512246
#> 17       9 107          E(Y|A=1) 2.60314625 0.11616776  2.37283238
#> 18       9 107 E(Y|A=1)-E(Y|A=0) 0.54975901 0.17746914  0.19790919
#>          UCL          pval alpha  Prob(>0)
#> 1  1.7294585 8.032773e-169  0.05 1.0000000
#> 2  1.9351604 7.200995e-180  0.05 1.0000000
#> 3  0.3239787  1.660367e-03  0.05 0.9991998
#> 4  1.5036947  3.103578e-22  0.05 1.0000000
#> 5  1.5341163  9.481730e-24  0.05 1.0000000
#> 6  0.3383114  8.193710e-01  0.05 0.5904724
#> 7  1.7448155  1.087589e-69  0.05 1.0000000
#> 8  1.8251186  1.235131e-61  0.05 1.0000000
#> 9  0.2593482  5.492246e-01  0.05 0.7256337
#> 10 1.8688514  1.876468e-19  0.05 1.0000000
#> 11 2.1826792  2.520863e-26  0.05 1.0000000
#> 12 0.7090122  8.401201e-02  0.05 0.9595611
#> 13 1.9431425  1.190149e-45  0.05 1.0000000
#> 14 2.2343620  1.958038e-48  0.05 1.0000000
#> 15 0.5332153  3.770973e-02  0.05 0.9818998
#> 16 2.3216520  2.475816e-28  0.05 1.0000000
#> 17 2.8334601  5.223720e-42  0.05 1.0000000
#> 18 0.9016088  2.496084e-03  0.05 0.9990251
## Forest plot: Overall/subgroup specific parameter estimates (CIs)
plot(res0, type="forest")
```

![](man/figures/README-example-3.png)

``` r

## Dependence Plots (univariate and heat maps)
plot_dependence(res0, vars="X1")
#> $res.est
#> `geom_smooth()` using method = 'loess' and formula 'y ~ x'
```

![](man/figures/README-example-4.png)

``` r
plot_dependence(res0, vars="X2")
#> $res.est
#> `geom_smooth()` using method = 'loess' and formula 'y ~ x'
```

![](man/figures/README-example-5.png)

``` r
plot_dependence(res0, vars=c("X1", "X2"))
#> $heatmap.est
```

![](man/figures/README-example-6.png)

    #> 
    #> $heatmap.prob

![](man/figures/README-example-7.png)

Overall, PRISM provides information at the patient-level, the subgroup-level (if any), and the overall population. While there are defaults in place, the user can also input their own functions/model wrappers into the PRISM algorithm. For more details and more examples, we refer the reader to the following vignettes, [PRISM\_vignette](https://CRAN.R-project.org/package=StratifiedMedicine/vignettes/SM_PRISM.html), [User\_Specific\_Models](https://CRAN.R-project.org/package=StratifiedMedicine/vignettes/SM_User_Models.html).
