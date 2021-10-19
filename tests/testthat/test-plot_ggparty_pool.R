test_that("Test whether plot_ggparty works with pooling (ctns)", {
  
  skip_on_cran()
  
  library(ggparty)
  library(partykit)
  library(rpart)
  
  dat_ctns = generate_subgrp_data(family="gaussian")
  Y = dat_ctns$Y
  X = dat_ctns$X
  A = dat_ctns$A
  family <- "gaussian"
  
  # Fit #
  res0 <- PRISM(Y, A, X, param="lm", pool = "trteff")
  Subgrps <- res0$out.train$Subgrps
  
  p0_est <- res0$param.dat
  p0_est <- p0_est[p0_est$Subgrps!="ovrl",]
  p0_est$label <- with(p0_est, paste( sprintf("%.2f", round(est,2)),
                                      " [",
                                      sprintf("%.2f", round(LCL,2)), ",",
                                      sprintf("%.2f", round(UCL,2)), "]", sep=""))
  p0_est$prob.est <- with(p0_est, sprintf("%.2f", round(`Prob(>0)`,2)))
  sub_map <- unique(res0$trt_assign)
  colnames(sub_map) <- c("Subgrps0", "Subgrps")
  p0_est <- dplyr::left_join(p0_est, sub_map, by="Subgrps")
  
  p0 <- plot(res0, type="tree")
  
  # Make sure the node estimates match up #
  pdat0_est <- data.frame(Subgrps = p0$data$id,
                          estimand = p0$data$estimand,
                          prob.est = p0$data$prob.est, 
                          label = p0$data$label)
  pdat0_est$Subgrps0 <- as.character(pdat0_est$Subgrps)
  
  p0_check <- dplyr::left_join(p0_est, pdat0_est, by=c("Subgrps0"))
  
  p0_label_check <- all.equal(p0_check$label.x, p0_check$label.y)
  p0_prob_check <- all.equal(p0_check$prob.est.x, p0_check$prob.est.y)
  p0_estimand_check <- ifelse(unique(na.omit(pdat0_est$estimand))=="E(Y|A=1)-E(Y|A=0)", 
                              TRUE, FALSE)
  
  # Check #
  expect_equal(p0_prob_check, TRUE)
  expect_equal(p0_label_check, TRUE)
  expect_equal(p0_estimand_check, TRUE)
  
})

test_that("Test whether plot_ggparty works with pooling (binomial)", {
  
  skip_on_cran()
  
  library(ggparty)
  library(partykit)
  library(rpart)
  
  dat_bin = generate_subgrp_data(family="binomial")
  Y = dat_bin$Y
  X = dat_bin$X
  A = dat_bin$A
  family <- "binomial"
  
  # Fit #
  res0 <- PRISM(Y, A, X, param="lm", pool = "trteff")
  Subgrps <- res0$out.train$Subgrps
  
  p0_est <- res0$param.dat
  p0_est <- p0_est[p0_est$Subgrps!="ovrl",]
  p0_est$label <- with(p0_est, paste( sprintf("%.2f", round(est,2)),
                                      " [",
                                      sprintf("%.2f", round(LCL,2)), ",",
                                      sprintf("%.2f", round(UCL,2)), "]", sep=""))
  p0_est$prob.est <- with(p0_est, sprintf("%.2f", round(`Prob(>0)`,2)))
  sub_map <- unique(res0$trt_assign)
  colnames(sub_map) <- c("Subgrps0", "Subgrps")
  p0_est <- dplyr::left_join(p0_est, sub_map, by="Subgrps")
  
  p0 <- plot(res0, type="tree")
  
  # Make sure the node estimates match up #
  pdat0_est <- data.frame(Subgrps = p0$data$id,
                          estimand = p0$data$estimand,
                          prob.est = p0$data$prob.est, 
                          label = p0$data$label)
  pdat0_est$Subgrps0 <- as.character(pdat0_est$Subgrps)
  
  p0_check <- dplyr::left_join(p0_est, pdat0_est, by=c("Subgrps0"))
  
  p0_label_check <- all.equal(p0_check$label.x, p0_check$label.y)
  p0_prob_check <- all.equal(p0_check$prob.est.x, p0_check$prob.est.y)
  p0_estimand_check <- ifelse(unique(na.omit(pdat0_est$estimand))=="P(Y=1|A=1)-P(Y=1|A=0)", 
                              TRUE, FALSE)
  
  # Check #
  expect_equal(p0_prob_check, TRUE)
  expect_equal(p0_label_check, TRUE)
  expect_equal(p0_estimand_check, TRUE)
  
})

test_that("Test whether plot_ggparty works with pooling (survival)", {
  
  skip_on_cran()
  
  library(ggparty)
  library(partykit)
  library(rpart)
  library(survival)
  require(TH.data); require(coin)
  data("GBSG2", package = "TH.data")
  surv.dat = GBSG2
  # Design Matrices ###
  Y = with(surv.dat, Surv(time, cens))
  X = surv.dat[,!(colnames(surv.dat) %in% c("time", "cens")) ]
  set.seed(513)
  A = rbinom(n = dim(X)[1], size=1, prob=0.5)
  
  family <- "survival"
  
  # Fit #
  res0 <- PRISM(Y, A, X, family = family, param="cox", ple="None",
                pool = "trteff")
  Subgrps <- res0$out.train$Subgrps
  
  p0_est <- res0$param.dat
  p0_est <- p0_est[p0_est$Subgrps!="ovrl",]
  p0_est$label <- with(p0_est, paste( sprintf("%.2f", round(exp(est),2)),
                                      " [",
                                      sprintf("%.2f", round(exp(LCL),2)), ",",
                                      sprintf("%.2f", round(exp(UCL),2)), "]", sep=""))
  p0_est$prob.est <- with(p0_est, sprintf("%.2f", round(1-`Prob(>0)`,2)))
  sub_map <- unique(res0$trt_assign)
  colnames(sub_map) <- c("Subgrps0", "Subgrps")
  p0_est <- dplyr::left_join(p0_est, sub_map, by="Subgrps")
  
  p0 <- plot(res0, type="tree")
  
  # Make sure the node estimates match up #
  pdat0_est <- data.frame(Subgrps = p0$data$id,
                          estimand = p0$data$estimand,
                          prob.est = p0$data$prob.est, 
                          label = p0$data$label)
  pdat0_est$Subgrps0 <- as.character(pdat0_est$Subgrps)
  
  p0_check <- dplyr::left_join(p0_est, pdat0_est, by=c("Subgrps0"))
  
  p0_label_check <- all.equal(p0_check$label.x, p0_check$label.y)
  p0_prob_check <- all.equal(p0_check$prob.est.x, p0_check$prob.est.y)
  p0_estimand_check <- ifelse(unique(na.omit(pdat0_est$estimand))=="HR_1-HR_0", 
                              TRUE, FALSE)
  
  # Check #
  expect_equal(p0_prob_check, TRUE)
  expect_equal(p0_label_check, TRUE)
  expect_equal(p0_estimand_check, TRUE)
  
})