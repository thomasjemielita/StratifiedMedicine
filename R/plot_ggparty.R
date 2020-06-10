## ggparty: Incorporate parameter estimates/CIs/etc directly onto tree plot ##
plot_ggparty = function(object, plots, prob.thres, width_out, nudge_out,
                     width_dens, nudge_dens) {
  # Extract #
  n <- dim(object$out.train)[1]
  family <- object$family
  prob.label <- paste("Prob(",prob.thres, ")", sep="")
  if (family=="survival"){
    n.events <- object$param.dat$events[object$param.dat$Subgrps==0] 
  }
  alpha_s <- object$alpha_s
  resample <- ifelse(is.null(object$resample), "None", object$resample)
  ct <- object$submod.fit$mod
  # Parameter Estimates #
  param.dat <- object$param.dat
  # param.dat$est0 <- param.dat$est
  # param.dat$SE0 <- param.dat$SE
  # param.dat$LCL0 <- param.dat$LCL
  # param.dat$UCL0 <- param.dat$UCL
  param.dat$prob.est <- param.dat$`Prob(>0)`
  Subgrps <- as.numeric(unique(param.dat$Subgrps[param.dat$Subgrps!="ovrl"]))
  
  param.dat$label <- with(param.dat, paste( sprintf("%.2f", round(est0,2)),
                                            " [",
                                            sprintf("%.2f", round(LCL0,2)), ",",
                                            sprintf("%.2f", round(UCL0,2)), "]", sep=""))
  param.dat$prob.est <- with(param.dat, sprintf("%.2f", round(prob.est,2)))
  param.dat$id <- param.dat$Subgrps
  param.ovrl <- param.dat[param.dat$Subgrps=="ovrl",]
  param.subs <- param.dat[param.dat$Subgrps!="ovrl",]
  estimand <- unique(param.dat$estimand)
  
  # A levels (if not NULL) #
  if (!is.null(object$out.train$A)) {
    A_lvls <- with(object$out.train, unique(A)[order(unique(A))])
  }
  ## Prep Data for Outcome Plots ##
  if (family!="survival") {
    if (!is.null(object$out.train$A)) {
      mu_A0 <- paste("mu", A_lvls[1], sep="_")
      mu_A1 <- paste("mu", A_lvls[2], sep="_")
      E_diff <- paste(mu_A1, "-", mu_A0, sep="")
      param.subs <- param.subs[param.subs$estimand==E_diff,]
      param.subs$estimand <- as.character(param.subs$estimand)
      estimand <- E_diff
      if (object$ple!="None") {
        mu_hat <- object$mu_train
        plot.dat <- rbind( data.frame(estimand=mu_A0, est=mu_hat[,mu_A0], 
                                      Subgrps = object$out.train$Subgrps),
                           data.frame(estimand=mu_A1, est=mu_hat[,mu_A1],
                                      Subgrps = object$out.train$Subgrps))
        plot.dat$id <- as.character(plot.dat$Subgrps)
      }
      if (object$ple=="None") {
        plot.dat <- object$out.train[,colnames(object$out.train) %in% 
                                       c("Y", "A", "Subgrps")]
        colnames(plot.dat)[1] <- c("est")
        plot.dat$estimand <- with(plot.dat, paste("mu",A,sep="_"))
        plot.dat$estimand <- factor(plot.dat$estimand, levels=c(mu_A0, mu_A1))
        plot.dat$id <- as.character(plot.dat$Subgrps)
      }
    }
    # If no Treatment #
    if (is.null(object$out.train$A)) {
      plot.dat <- object$out.train[,colnames(object$out.train) %in% 
                                     c("Y", "Subgrps")]
      colnames(plot.dat)[1] <- c("est")
      plot.dat$estimand <- "E(Y)"
      plot.dat$id <- plot.dat$Subgrps
      plot.dat$id <- as.character(plot.dat$Subgrps)
    } 
  }
  if (family=="survival") {
    # Create KM survival prediction data-set #
    out.train <- object$out.train
    pred.surv <- NULL
    for (s in unique(out.train$Subgrps)) {
      hold.dat <- out.train[out.train$Subgrps==s,]
      if (is.null(out.train$A)){
        km_mod <- survfit(Y ~ 1, data=hold.dat)
        pred.s <- data.frame(A = "",
                             time=c(0,km_mod$time), surv=c(1, km_mod$surv))
      }
      if (!is.null(out.train$A)) {
        pred.s <- NULL
        for (a in 1:length(A_lvls)) {
          km_a <- survfit(Y ~ 1, data=hold.dat[hold.dat$A==A_lvls[a],])
          surv_a <- data.frame(A=A_lvls[a], time=c(0,km_a$time), 
                               surv=c(1, km_a$surv))
          pred.s <- rbind(pred.s, surv_a)
        }
        pred.s$A <- factor(pred.s$A, levels = A_lvls)
      }
      pred.s <- data.frame(id=s, pred.s)
      pred.surv <- suppressWarnings( bind_rows(pred.surv, pred.s) )
    } 
    plot.dat <- pred.surv
  }
  
  ## Create density data ##
  post.prob <- object$resamp.dist
  gen_post <- FALSE
  if (is.null(object$resample)) {
    gen_post <- TRUE
  }
  if (!is.null(object$resample)) {
    if (object$resample=="CV") gen_post <- TRUE
  }
  if (!gen_post & object$param=="cox") {
    post.prob$estimand <- gsub("logHR", "HR", post.prob$estimand)
  }
  if (gen_post) {
    post.prob <- NULL
    if (object$param=="cox") {
      param.subs$est0 <- log(param.subs$est0)
    }
    for (s in unique(Subgrps)) {
      dat.s <- rnorm(10000, mean = param.subs$est0[param.subs$Subgrps==s],
                     sd = param.subs$SE0[param.subs$Subgrps==s])
      dat.s <- data.frame(Subgrps=s, estimand=estimand, est=dat.s)
      post.prob <- rbind(post.prob, dat.s)
    }
  }
  post.prob <- post.prob[post.prob$Subgrps!="ovrl",]
  post.prob <- post.prob[post.prob$estimand==estimand,]
  if (family=="survival" & object$param=="cox") {
    post.prob$est <- exp(post.prob$est)
  }
  dat.dens <- NULL
  for (s in unique(post.prob$Subgrps)){
    hold.s = post.prob %>% filter(Subgrps==s)
    dat.s = with(density(hold.s$est, na.rm=T), data.frame(x, y))
    dat.s = data.frame(id = s, dat.s)
    dat.dens = rbind(dat.dens, dat.s)
  }
  max.dens <- max(dat.dens$y)*1.05
  dat.dens$id <- as.character(dat.dens$id)
  
  # If Applicable: Map "Pooled Subgroups" back to original tree #
  if (!is.null(object$pool.dat)) {
    pool.dat <- object$pool.dat
    pool.dat <- pool.dat[,c("Subgrps", "pred_opt")]
    colnames(pool.dat) <- c("Subgrps0", "Subgrps")
    # Merge with param #
    param.dat <- left_join(param.dat, pool.dat, by="Subgrps")
    param.dat$Subgrps <- with(param.dat, ifelse(Subgrps=="ovrl", "ovrl",
                                                as.character(Subgrps0)))
    # Merge with plot datas #
    colnames(pool.dat) <- c("Subgrps0", "id")
    plot.dat <- left_join(plot.dat, pool.dat, by="id")
    plot.dat$trt_assign <- plot.dat$id
    plot.dat$id <- plot.dat$Subgrps0 
    dat.dens <- left_join(dat.dens, pool.dat, by="id")
    dat.dens$trt_assign <- dat.dens$id
    dat.dens$id <- dat.dens$Subgrps0
    Subgrps <- as.numeric(unique(dat.dens$Subgrps0))
    param.subs <- left_join(param.subs, pool.dat, by="id")
    param.subs$trt_assign <- param.subs$id
    param.subs$Subgrps <- param.subs$Subgrps0
  }
  # Add estimates into tree #
  ct_node <- as.list(ct$node)
  for (s in Subgrps) {
    ct_node[[s]]$info$label <- param.subs$label[param.subs$Subgrps==s] 
    ct_node[[s]]$info$N <- param.subs$N[param.subs$Subgrps==s] 
    ct_node[[s]]$info$estimand <- param.subs$estimand[param.subs$Subgrps==s]
    ct_node[[s]]$info$prob.est <- param.subs$prob.est[param.subs$Subgrps==s]
    if (family=="survival") {
      ct_node[[s]]$info$events <- param.subs$events[param.subs$Subgrps==s] 
    }
  }
  # Round Edge Labels #
  for (ii in 1:length(ct_node)) {
    brk = ct_node[[ii]]$split$breaks
    if (!is.null(brk)){
      ct_node[[ii]]$split$breaks <- ifelse(is.numeric(brk), round(brk,2), brk)
    }
  }
  ct$node <- as.partynode(ct_node)
  
  
  ## Plot setup ##
  add_vars_list <- list(N = "$node$info$N", est = "$node$info$label", 
                        estimand = "$node$info$estimand", 
                        prob.est = "$node$info$prob.est",
                        label = "$node$info$label")
  if (family=="survival"){
    add_vars_list[["events"]] = "$node$info$events"
  }
  
  plt <- ggparty::ggparty(ct, horizontal = TRUE, add_vars = add_vars_list) +
    ggparty::geom_edge() +
    ggparty::geom_edge_label() +
    ggparty::geom_node_label(line_list = list(ggplot2::aes(label = splitvar),
                                     ggplot2::aes(label = pval_convert(p.value))),
                    line_gpar = list(list(size = 8, col = "black", fontface = "bold"),
                                     list(size = 8)),ids = "inner")
  # Set up prob cutoff params #
  dir <- substr(prob.thres, 1, 1)
  thres <- substr(prob.thres, 2, nchar(prob.thres))
  if (suppressWarnings(is.na(as.numeric(thres)))) {
    thres <- substr(prob.thres, 3, nchar(prob.thres))
  }
  thres <- as.numeric(thres)
  fill_L <- ifelse(dir=="<", "green", "red")
  fill_R <- ifelse(dir==">", "green", "red")
  ## Continuous/Binary Plots ##
  if (family!="survival") {
    # Add Node Info #
    plt <- plt + 
      ggparty::geom_node_label(
      line_list = list(
        ggplot2::aes(label = paste("Node ", id) ),
        ggplot2::aes(label = paste("N = ", N, " (", round(N/n*100,1), "%)", sep="")),
        ggplot2::aes(label = paste(estimand, " [", (1-alpha_s)*100,"% CI]", sep="")),
        ggplot2::aes(label = label),
        ggplot2::aes(label = paste(prob.label, "=", prob.est))
      ),
      line_gpar = list(list(size = 8, fontface="bold"),
                       list(size = 8),
                       list(size = 8, fontface="bold"),
                       list(size = 8),
                       list(size = 8)),
      ids = "terminal")
    # Outcome Plot #
    out_plot <- ggparty::geom_node_plot(
      gglist = list(ggplot2::geom_boxplot(data=plot.dat,
                                          ggplot2::aes(x=estimand, y=est, fill=estimand)),
                                      ggplot2::scale_x_discrete(expand=c(0.4, 0.4)), 
                                      ggplot2::xlab(""),
                                      ggplot2::ylab("Estimate"),
                                      ggplot2::theme_bw(), 
                                      ggplot2::theme(axis.text.y = element_blank()),
                                      ggplot2::coord_flip()),
                               nudge_x = nudge_out, width = width_out,
                               shared_axis_labels = TRUE,
                               legend_separator = TRUE, scales = "fixed")
    # Density Plot #
    dens_plot <- ggparty::geom_node_plot(
      gglist = list( ggplot2::geom_area(data=dat.dens, 
                     ggplot2::aes(x = ifelse(x < thres , x, thres),y=y), fill = fill_L),
                     ggplot2::geom_area(data=dat.dens, aes(x = ifelse(x > thres , x, thres),
                                                  y=y), fill = fill_R),
                     ggplot2::geom_line(data=dat.dens, aes(x=x,y=y)),
                     ggplot2::ylim(0,max.dens),
                     ggplot2::theme_bw(),
                     ggplot2::ggtitle(estimand),
                     ggplot2::theme(plot.title = element_text(size = 8, face = "bold")),
                     ggplot2::xlab("Density"), ylab("") ), shared_axis_labels = TRUE,
      width=width_dens, nudge_x = nudge_dens)
    # Outcome Only #
    if (plots=="outcome") {
      plt.out <- plt + out_plot 
    }
    if (plots=="density") {
      plt.out <- plt + dens_plot
    }
    if (plots=="both") {
      plt.out <- plt + out_plot + dens_plot
    }
  }
  ### Survival Plot ###
  if (family=="survival") {
    # Add Node Info #
    plt <- plt +
      ggparty::geom_node_label(
        line_list = list(
          ggplot2::aes(label = paste("Node ", id) ),
          ggplot2::aes(label = paste("N = ", N, " (", round(N/n*100,1), "%)", sep="")),
          ggplot2::aes(label = paste("Events = ", events,
                            " (", round(events/n.events*100,1), "%)", sep="")),
          ggplot2::aes(label = paste(estimand, " [", (1-alpha_s)*100,"% CI]", sep="")),
          ggplot2::aes(label = label),
          ggplot2::aes(label = paste(prob.label, "=", prob.est))
        ),
        line_gpar = list(list(size = 8, fontface="bold"),
                         list(size = 8),
                         list(size = 8),
                         list(size = 8, fontface="bold"),
                         list(size = 8), 
                         list(size = 8)),
        ids = "terminal")
    if (is.null(out.train$A)) {
      gg_line <- ggplot2::geom_line(data=plot.dat,aes(x=time, y=surv), size=1.5)
    }
    if (!is.null(out.train$A)) {
      gg_line <- ggplot2::geom_line(data=plot.dat,aes(x=time, y=surv, col=A), size=1.5)
      plot.dat$A <- factor(plot.dat$A, levels = A_lvls)
    }
    out_plot <- ggparty::geom_node_plot(gglist = list(gg_line,
                                             ggplot2::xlab("Time"),
                                             ggplot2::ylab( "Survival Probability" ),
                                             ggplot2::theme_bw(),
                                             ggplot2::theme(legend.title = element_blank())),
                               shared_axis_labels = TRUE,
                               legend_separator = TRUE, scales = "fixed", 
                               width=width_out, nudge_x=nudge_out)
    # Density #
    dens_plot <- ggparty::geom_node_plot(
      gglist = list( ggplot2::geom_area(data=dat.dens, 
                     ggplot2::aes(x = ifelse(x < thres , x, thres),y=y), fill = fill_L),
                     ggplot2::geom_area(data=dat.dens, aes(x = ifelse(x > thres , x, thres),
                                                  y=y), fill = fill_R),
                     ggplot2::geom_line(data=dat.dens, aes(x=x,y=y)),
                     ggplot2::ylim(0,max.dens), 
                     ggplot2::theme_bw(), ggplot2::xlab("Density"), 
                     ggplot2::ylab(""),
                     ggtitle(estimand),
                     ggplot2::theme(plot.title = element_text(size = 8, face = "bold"))),
      shared_axis_labels = TRUE, width = width_dens, nudge_x = nudge_dens)
    if (plots=="outcome") {
      plt.out <- plt + out_plot
    }
    if (plots=="density") {
      plt.out <- plt + dens_plot
    }
    if (plots=="both") {
      plt.out <- plt + out_plot + dens_plot
    }
  }
  plt.out
  return(plt.out)
}