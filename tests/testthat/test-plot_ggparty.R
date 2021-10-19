test_that("Test whether plot_ggparty works (ctns)", {
  
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
  res0 <- PRISM(Y, A, X, param="lm")
  Subgrps <- res0$out.train$Subgrps
  
  p0_est <- res0$param.dat
  p0_est <- p0_est[p0_est$Subgrps!="ovrl",]
  p0_est$label <- with(p0_est, paste( sprintf("%.2f", round(est,2)),
                                            " [",
                                            sprintf("%.2f", round(LCL,2)), ",",
                                            sprintf("%.2f", round(UCL,2)), "]", sep=""))
  p0_est$prob.est <- with(p0_est, sprintf("%.2f", round(`Prob(>0)`,2)))
  
  p0 <- plot(res0, type="tree")
  
  # Make sure the node estimates match up #
  pdat0_est <- data.frame(Subgrps = p0$data$id,
                          estimand = p0$data$estimand,
                          prob.est = p0$data$prob.est, 
                          label = p0$data$label)
  pdat0_est$Subgrps <- as.character(pdat0_est$Subgrps)
  
  p0_check <- dplyr::left_join(p0_est, pdat0_est, by=c("Subgrps"))
  
  p0_label_check <- all.equal(p0_check$label.x, p0_check$label.y)
  p0_prob_check <- all.equal(p0_check$prob.est.x, p0_check$prob.est.y)
  p0_estimand_check <- ifelse(unique(na.omit(pdat0_est$estimand))=="E(Y|A=1)-E(Y|A=0)", 
                              TRUE, FALSE)
  
  # Make sure the outcome plots are correct #
  p0_outdat <- NULL
  
  for (s in unique(p0_check$Subgrps)) {
    hold <- p0$data$info_list[[as.numeric(s)]]$plot.dat
    hold <- data.frame(Subgrps = s, est_Y = hold$est,
                       Y = Y[Subgrps==s])
    p0_outdat <- rbind(hold, p0_outdat)
  }
  p0_outplot_check <- all.equal(p0_outdat$est_Y, p0_outdat$Y)
    
  
  ### Run with resampling: Check Density ###
  
  # Fit #
  res1 <- PRISM(Y, A, X, param="lm", ple="None",
                resample = "Bootstrap", filter.resamp = "None",
                R = 50)
  Subgrps <- res1$out.train$Subgrps
  
  p1_est <- res1$param.dat
  p1_est <- p1_est[p1_est$Subgrps!="ovrl",]
  p1_est$label <- with(p1_est, paste( sprintf("%.2f", round(est_resamp,2)),
                                      " [",
                                      sprintf("%.2f", round(LCL.pct,2)), ",",
                                      sprintf("%.2f", round(UCL.pct,2)), "]", sep=""))
  p1_est$prob.est <- with(p1_est, sprintf("%.2f", round(`Prob(>0)`,2)))
  
  p1 <- plot(res1, type="tree", tree.plots = "density")
  
  # Make sure the node estimates match up #
  pdat1_est <- data.frame(Subgrps = p1$data$id,
                          estimand = p1$data$estimand,
                          prob.est = p1$data$prob.est, 
                          label = p1$data$label)
  pdat1_est$Subgrps <- as.character(pdat1_est$Subgrps)
  
  p1_check <- dplyr::left_join(p1_est, pdat1_est, by=c("Subgrps"))
  
  p1_label_check <- all.equal(p1_check$label.x, p1_check$label.y)
  p1_prob_check <- all.equal(p1_check$prob.est.x, p1_check$prob.est.y)
  p1_estimand_check <- ifelse(unique(na.omit(pdat1_est$estimand))=="E(Y|A=1)-E(Y|A=0)", 
                              TRUE, FALSE)
  
  # Make sure the density plots are correct #
  p1_densdat <- NULL
  
  for (s in unique(p1_check$Subgrps)) {
    
    plot.s <- p1$data$info_list[[as.numeric(s)]]$dat.dens
    hold.s = res1$resamp.dist[res1$resamp.dist$Subgrps==s,]
    dat.s = with(density(hold.s$est, na.rm=T), data.frame(x, y))
    hold <- data.frame(plot.s, dat.s)
    p1_densdat <- rbind(hold, p1_densdat)
  }
  p1_densx_check <- all.equal(p1_densdat$x, p1_densdat$x.1)
  p1_densy_check <- all.equal(p1_densdat$y, p1_densdat$y.1)
  
  # Check #
  expect_equal(p0_prob_check, TRUE)
  expect_equal(p0_label_check, TRUE)
  expect_equal(p0_outplot_check, TRUE)
  expect_equal(p0_estimand_check, TRUE)
  
  expect_equal(p1_prob_check, TRUE)
  expect_equal(p1_label_check, TRUE)
  expect_equal(p1_densx_check, TRUE)
  expect_equal(p1_densy_check, TRUE)
  expect_equal(p1_estimand_check, TRUE)
  
})

test_that("Test whether plot_ggparty works (binomial)", {
  
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
  res0 <- PRISM(Y, A, X, param="lm")
  Subgrps <- res0$out.train$Subgrps
  
  p0_est <- res0$param.dat
  p0_est <- p0_est[p0_est$Subgrps!="ovrl",]
  p0_est$label <- with(p0_est, paste( sprintf("%.2f", round(est,2)),
                                      " [",
                                      sprintf("%.2f", round(LCL,2)), ",",
                                      sprintf("%.2f", round(UCL,2)), "]", sep=""))
  p0_est$prob.est <- with(p0_est, sprintf("%.2f", round(`Prob(>0)`,2)))
  
  p0 <- plot(res0, type="tree")
  
  # Make sure the node estimates match up #
  pdat0_est <- data.frame(Subgrps = p0$data$id,
                          estimand = p0$data$estimand,
                          prob.est = p0$data$prob.est, 
                          label = p0$data$label)
  pdat0_est$Subgrps <- as.character(pdat0_est$Subgrps)
  
  p0_check <- dplyr::left_join(p0_est, pdat0_est, by=c("Subgrps"))
  
  p0_label_check <- all.equal(p0_check$label.x, p0_check$label.y)
  p0_prob_check <- all.equal(p0_check$prob.est.x, p0_check$prob.est.y)
  p0_estimand_check <- ifelse(unique(na.omit(pdat0_est$estimand))=="P(Y=1|A=1)-P(Y=1|A=0)", 
                              TRUE, FALSE)
  
  # Make sure the outcome plots are correct #
  p0_outdat <- NULL
  
  for (s in unique(p0_check$Subgrps)) {
    hold <- p0$data$info_list[[as.numeric(s)]]$plot.dat
    hold <- data.frame(Subgrps = s, est_Y = hold$est,
                       Y = Y[Subgrps==s])
    p0_outdat <- rbind(hold, p0_outdat)
  }
  p0_outplot_check <- all.equal(p0_outdat$est_Y, p0_outdat$Y)
  
  
  ### Run with resampling: Check Density ###
  
  # Fit #
  res1 <- PRISM(Y, A, X, param="lm", ple="None",
                resample = "Bootstrap", filter.resamp = "None",
                R = 50)
  Subgrps <- res1$out.train$Subgrps
  
  p1_est <- res1$param.dat
  p1_est <- p1_est[p1_est$Subgrps!="ovrl",]
  p1_est$label <- with(p1_est, paste( sprintf("%.2f", round(est_resamp,2)),
                                      " [",
                                      sprintf("%.2f", round(LCL.pct,2)), ",",
                                      sprintf("%.2f", round(UCL.pct,2)), "]", sep=""))
  p1_est$prob.est <- with(p1_est, sprintf("%.2f", round(`Prob(>0)`,2)))
  
  p1 <- plot(res1, type="tree", tree.plots = "density")
  
  # Make sure the node estimates match up #
  pdat1_est <- data.frame(Subgrps = p1$data$id,
                          estimand = p1$data$estimand,
                          prob.est = p1$data$prob.est, 
                          label = p1$data$label)
  pdat1_est$Subgrps <- as.character(pdat1_est$Subgrps)
  
  p1_check <- dplyr::left_join(p1_est, pdat1_est, by=c("Subgrps"))
  
  p1_label_check <- all.equal(p1_check$label.x, p1_check$label.y)
  p1_prob_check <- all.equal(p1_check$prob.est.x, p1_check$prob.est.y)
  p1_estimand_check <- ifelse(unique(na.omit(pdat1_est$estimand))=="P(Y=1|A=1)-P(Y=1|A=0)", 
                              TRUE, FALSE)
  
  # Make sure the density plots are correct #
  p1_densdat <- NULL
  
  for (s in unique(p1_check$Subgrps)) {
    
    plot.s <- p1$data$info_list[[as.numeric(s)]]$dat.dens
    hold.s = res1$resamp.dist[res1$resamp.dist$Subgrps==s,]
    dat.s = with(density(hold.s$est, na.rm=T), data.frame(x, y))
    hold <- data.frame(plot.s, dat.s)
    p1_densdat <- rbind(hold, p1_densdat)
  }
  p1_densx_check <- all.equal(p1_densdat$x, p1_densdat$x.1)
  p1_densy_check <- all.equal(p1_densdat$y, p1_densdat$y.1)
  
  # Check #
  expect_equal(p0_prob_check, TRUE)
  expect_equal(p0_label_check, TRUE)
  expect_equal(p0_outplot_check, TRUE)
  expect_equal(p0_estimand_check, TRUE)
  
  expect_equal(p1_prob_check, TRUE)
  expect_equal(p1_label_check, TRUE)
  expect_equal(p1_densx_check, TRUE)
  expect_equal(p1_densy_check, TRUE)
  expect_equal(p1_estimand_check, TRUE)
  
})

test_that("Test whether plot_ggparty works (survival)", {
  
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
  res0 <- PRISM(Y, A, X, family = family, param="cox", ple="None")
  Subgrps <- res0$out.train$Subgrps
  
  p0_est <- res0$param.dat
  p0_est <- p0_est[p0_est$Subgrps!="ovrl",]
  p0_est$label <- with(p0_est, paste( sprintf("%.2f", round(exp(est),2)),
                                      " [",
                                      sprintf("%.2f", round(exp(LCL),2)), ",",
                                      sprintf("%.2f", round(exp(UCL),2)), "]", sep=""))
  p0_est$prob.est <- with(p0_est, sprintf("%.2f", round(1-`Prob(>0)`,2)))
  
  p0 <- plot(res0, type="tree")
  
  # Make sure the node estimates match up #
  pdat0_est <- data.frame(Subgrps = p0$data$id,
                          estimand = p0$data$estimand,
                          prob.est = p0$data$prob.est, 
                          label = p0$data$label)
  pdat0_est$Subgrps <- as.character(pdat0_est$Subgrps)
  
  p0_check <- dplyr::left_join(p0_est, pdat0_est, by=c("Subgrps"))
  
  p0_label_check <- all.equal(p0_check$label.x, p0_check$label.y)
  p0_prob_check <- all.equal(p0_check$prob.est.x, p0_check$prob.est.y)
  p0_estimand_check <- ifelse(unique(na.omit(pdat0_est$estimand))=="HR_1-HR_0", 
                              TRUE, FALSE)
  
  # Make sure the outcome plots are correct #
  p0_outdat <- NULL
  
  for (s in unique(p0_check$Subgrps)) {
    hold <- p0$data$info_list[[as.numeric(s)]]$plot.dat
    hold.dat <- data.frame(Y, A)[Subgrps==as.numeric(s),]
    pred.s <- NULL
    for (a in c(0,1) ) {
      km_a <- survfit(Y ~ 1, data=hold.dat[hold.dat$A==a,])
      surv_a <- data.frame(A=a, time=c(0,km_a$time), 
                           surv=c(1, km_a$surv))
      pred.s <- rbind(pred.s, surv_a)
    }
    pred.s$A <- factor(pred.s$A, levels = c("0", "1"))
    pred.s <- data.frame(id=s, pred.s)
    hold <- data.frame(hold, pred.s)
    p0_outdat <- rbind(hold, p0_outdat)
  }
  p0_outplot_check <- all.equal(p0_outdat$surv, p0_outdat$surv.1)
  
  
  ### Run with resampling: Check Density ###
  
  # Fit #
  res1 <- PRISM(Y, A, X, param="cox", ple="None",
                resample = "Bootstrap", filter.resamp = "None",
                R = 50)
  Subgrps <- res1$out.train$Subgrps
  
  p1_est <- res1$param.dat
  p1_est <- p1_est[p1_est$Subgrps!="ovrl",]
  p1_est$label <- with(p1_est, paste( sprintf("%.2f", round(exp(est_resamp),2)),
                                      " [",
                                      sprintf("%.2f", round(exp(LCL.pct),2)), ",",
                                      sprintf("%.2f", round(exp(UCL.pct),2)), "]", sep=""))
  p1_est$prob.est <- with(p1_est, sprintf("%.2f", round(1-`Prob(>0)`,2)))
  
  p1 <- plot(res1, type="tree", tree.plots = "density")
  
  # Make sure the node estimates match up #
  pdat1_est <- data.frame(Subgrps = p1$data$id,
                          estimand = p1$data$estimand,
                          prob.est = p1$data$prob.est, 
                          label = p1$data$label)
  pdat1_est$Subgrps <- as.character(pdat1_est$Subgrps)
  
  p1_check <- dplyr::left_join(p1_est, pdat1_est, by=c("Subgrps"))
  
  p1_label_check <- all.equal(p1_check$label.x, p1_check$label.y)
  p1_prob_check <- all.equal(p1_check$prob.est.x, p1_check$prob.est.y)
  p1_estimand_check <- ifelse(unique(na.omit(pdat1_est$estimand))=="HR_1-HR_0", 
                              TRUE, FALSE)
  
  # Make sure the density plots are correct #
  p1_densdat <- NULL
  
  for (s in unique(p1_check$Subgrps)) {
    
    plot.s <- p1$data$info_list[[as.numeric(s)]]$dat.dens
    hold.s = res1$resamp.dist[res1$resamp.dist$Subgrps==s,]
    dat.s = with(density(exp(hold.s$est), na.rm=T), data.frame(x, y))
    hold <- data.frame(plot.s, dat.s)
    p1_densdat <- rbind(hold, p1_densdat)
  }
  p1_densx_check <- all.equal(p1_densdat$x, p1_densdat$x.1)
  p1_densy_check <- all.equal(p1_densdat$y, p1_densdat$y.1)
  
  # Check #
  expect_equal(p0_prob_check, TRUE)
  expect_equal(p0_label_check, TRUE)
  expect_equal(p0_outplot_check, TRUE)
  expect_equal(p0_estimand_check, TRUE)
  
  expect_equal(p1_prob_check, TRUE)
  expect_equal(p1_label_check, TRUE)
  expect_equal(p1_densx_check, TRUE)
  expect_equal(p1_densy_check, TRUE)
  expect_equal(p1_estimand_check, TRUE)
  
})