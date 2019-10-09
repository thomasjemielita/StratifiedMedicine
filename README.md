
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
## Plot the distribution of PLEs ###
plot(res0, type="PLE:density") # Density plot of PLEs #
```

![](man/figures/README-example-2.png)

``` r
plot(res0, type="PLE:waterfall") # waterfall plot of PLEs
```

![](man/figures/README-example-3.png)

``` r
## Overall/subgroup specific parameter estimates/inference
res0$param.dat
#>    Subgrps   N          estimand        est         SE         LCL
#> 1        0 800          E(Y|A=0) 1.63776866 0.04569877  1.54806484
#> 2        0 800          E(Y|A=1) 1.84286992 0.04866618  1.74734124
#> 3        0 800 E(Y|A=1)-E(Y|A=0) 0.20510125 0.06328619  0.08087441
#> 4        3 149          E(Y|A=0) 1.28263595 0.11176578  1.06177307
#> 5        3 149          E(Y|A=1) 1.31794039 0.10980919  1.10094399
#> 6        3 149 E(Y|A=1)-E(Y|A=0) 0.03530445 0.15439447 -0.26979794
#> 7        4 277          E(Y|A=0) 1.60937815 0.06712763  1.47723094
#> 8        4 277          E(Y|A=1) 1.67799950 0.07693142  1.52655259
#> 9        4 277 E(Y|A=1)-E(Y|A=0) 0.06862135 0.10088192 -0.12997442
#> 10       7  99          E(Y|A=0) 1.59590822 0.14128296  1.31553678
#> 11       7  99          E(Y|A=1) 1.92893256 0.13048087  1.66999752
#> 12       7  99 E(Y|A=1)-E(Y|A=0) 0.33302434 0.18994596 -0.04391723
#> 13       8 168          E(Y|A=0) 1.76751942 0.08920057  1.59141332
#> 14       8 168          E(Y|A=1) 2.04252025 0.09776124  1.84951308
#> 15       8 168 E(Y|A=1)-E(Y|A=0) 0.27500084 0.13063657  0.01708886
#> 16       9 107          E(Y|A=0) 2.04080610 0.13504592  1.77306443
#> 17       9 107          E(Y|A=1) 2.60756289 0.11565560  2.37826441
#> 18       9 107 E(Y|A=1)-E(Y|A=0) 0.56675679 0.17669715  0.21643749
#>          UCL          pval alpha
#> 1  1.7274725 1.879424e-168  0.05
#> 2  1.9383986 1.723772e-180  0.05
#> 3  0.3293281  1.241074e-03  0.05
#> 4  1.5034988  3.314673e-22  0.05
#> 5  1.5349368  1.324710e-23  0.05
#> 6  0.3404068  8.194458e-01  0.05
#> 7  1.7415254  1.972697e-69  0.05
#> 8  1.8294464  5.334423e-62  0.05
#> 9  0.2672171  4.969388e-01  0.05
#> 10 1.8762797  1.916696e-19  0.05
#> 11 2.1878676  1.077698e-26  0.05
#> 12 0.7099659  8.268452e-02  0.05
#> 13 1.9436255  1.028800e-45  0.05
#> 14 2.2355274  1.856343e-48  0.05
#> 15 0.5329128  3.678032e-02  0.05
#> 16 2.3085478  3.359840e-28  0.05
#> 17 2.8368614  3.054161e-42  0.05
#> 18 0.9170761  1.770848e-03  0.05
## Forest plot: Overall/subgroup specific parameter estimates (CIs)
plot(res0, type="forest")
```

![](man/figures/README-example-4.png)

``` r

## Heatmap of PLEs #
grid.data = expand.grid(X1 = seq(min(X$X1), max(X$X1), by=0.30),
                    X2 = seq(min(X$X2), max(X$X2), by=0.30))
plot(res0, type="heatmap", grid.data = grid.data)
#> $heatmap.est
```

![](man/figures/README-example-5.png)

    #> 
    #> $heatmap.prob

![](man/figures/README-example-6.png)

Overall, PRISM provides information at the patient-level, the subgroup-level (if any), and the overall population. While there are defaults in place, the user can also input their own functions/model wrappers into the PRISM algorithm. For more details and more examples, we refer the reader to the following vignettes, [PRISM\_vignette](https://CRAN.R-project.org/package=StratifiedMedicine/vignettes/SM_PRISM.html), [User\_Specific\_Models](https://CRAN.R-project.org/package=StratifiedMedicine/vignettes/SM_User_Models.html).
