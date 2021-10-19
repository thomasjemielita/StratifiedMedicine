# StratifiedMedicine 1.0.4
* Fixed bugs for plot() for type="resample" and tree.plots="density", such that density plots match up with right tree nodes.
* Added pool="trteff","trteff_boot" (pools subgroups based on naive or bootstrap resampling based treatment effect estimates)
* Added submod="ctree_cate", CATE~ctree(X), and submod="rpart_cate", CATE~rpart(X).

# StratifiedMedicine 1.0.3
* Fixed bugs for plot() (when family="survival"), such that Kaplan-meier plots match up with right tree nodes
* Updated documentation and plot labels

# StratifiedMedicine 1.0.2

* Fixed bugs relating to binary outcome data (family="binomial")
* Fixed bugs relating to resampling estimates with OTR pooling for non-default delta vlues.
(ex: PRISM(Y, A, X, resample="Bootstrap", pool="otr:logistic", delta=">0.10"))

# StratifiedMedicine 1.0.1

* Fixed mapping issue for boxplots for function plot(PRISM.fit,type="tree", tree.plots="outcome")
* Added propensity-based estimation for functions ple_train and PRISM

# StratifiedMedicine 1.0.0

* Added a `NEWS.md` file to track changes to the package.
* Version 1.0.0: Added new features and  improved functionality. See vignette for more details.
Changes included improved plotting features, treatment difference estimates through base-learners and meta-learners, partial dependence plots, variable importance plots, and enhanced functionality on individual functions. 
