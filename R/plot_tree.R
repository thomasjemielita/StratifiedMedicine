## PRISM submod plot (with resampling) ##
plot_tree = function(object, plots, prob.thres, width_out, nudge_out,
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
  Subgrps <- unique(param.dat$Subgrps[param.dat$Subgrps>0])
  
  param.dat$label <- with(param.dat, paste( sprintf("%.2f", round(est0,2)),
                                            " [",
                                            sprintf("%.2f", round(LCL0,2)), ",",
                                            sprintf("%.2f", round(UCL0,2)), "]", sep=""))
  param.dat$prob.est <- with(param.dat, sprintf("%.2f", round(prob.est,2)))
  param.dat$id <- param.dat$Subgrps
  param.ovrl <- param.dat[param.dat$Subgrps==0,]
  param.subs <- param.dat[param.dat$Subgrps>0,]
  estimand <- unique(param.dat$estimand)
  
  # A levels (if not NULL) #
  if (!is.null(object$out.train$A)) {
    A_lvls <- with(object$out.train, unique(A)[order(unique(A))])
  }
  
  if (family!="survival" & !is.null(object$out.train$A)) {
    E_A0 <- paste("E(Y|A=", A_lvls[1], ")", sep="")
    E_A1 <- paste("E(Y|A=", A_lvls[2], ")", sep="")
    E_diff <- paste(E_A1, "-", E_A0, sep="")
    mu_A0 <- paste("mu", A_lvls[1], sep="_")
    mu_A1 <- paste("mu", A_lvls[2], sep="_")
    param.subs <- param.dat[param.dat$estimand==E_diff,]
    param.subs$estimand <- as.character(param.subs$estimand)
    estimand <- E_diff
    if (object$ple!="None") {
      mu_hat <- object$mu_train
      plot.dat <- rbind( data.frame(estimand=E_A0, est=mu_hat[,mu_A0], 
                                    Subgrps = object$out.train$Subgrps),
                         data.frame(estimand=E_A1, est=mu_hat[,mu_A1],
                                    Subgrps = object$out.train$Subgrps))
      plot.dat$id <- plot.dat$Subgrps
    }
    if (object$ple=="None") {
      plot.dat <- object$out.train[,colnames(object$out.train) %in% 
                                     c("Y", "A", "Subgrps")]
      colnames(plot.dat)[1] <- c("est")
      plot.dat$estimand <- with(plot.dat, paste("E(Y|A=",A,")",sep=""))
      plot.dat$estimand <- factor(plot.dat$estimand, levels=c(E_A0, E_A1))
      plot.dat$id <- plot.dat$Subgrps
    }
  }
  # If no Treatment #
  if (is.null(object$out.train$A)) {
    plot.dat <- object$out.train[,colnames(object$out.train) %in% 
                                   c("Y", "Subgrps")]
    colnames(plot.dat)[1] <- c("est")
    plot.dat$estimand <- "E(Y)"
    plot.dat$id <- plot.dat$Subgrps
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
  
  ## Create density data ##
  post.prob <- object$resamp.dist
  gen_post <- FALSE
  if (is.null(object$resample)) {
    gen_post <- TRUE
  }
  if (!is.null(object$resample)) {
    if (object$resample=="CV") gen_post <- TRUE
  }
  if (!gen_post & object$param=="param_cox") {
    post.prob$estimand <- gsub("logHR", "HR", post.prob$estimand)
  }
  if (gen_post) {
    post.prob <- NULL
    if (object$param=="param_cox") {
      param.subs$est0 <- log(param.subs$est0)
    }
    for (s in unique(Subgrps)) {
      dat.s <- rnorm(10000, mean = param.subs$est0[param.subs$Subgrps==s],
                     sd = param.subs$SE0[param.subs$Subgrps==s])
      dat.s <- data.frame(Subgrps=s, estimand=estimand, est=dat.s)
      post.prob <- rbind(post.prob, dat.s)
    }
  }
  post.prob <- post.prob[post.prob$Subgrps>0,]
  post.prob <- post.prob[post.prob$estimand==estimand,]
  if (family=="survival" & object$param=="param_cox") {
    post.prob$est <- exp(post.prob$est)
  }
  dat.dens <- NULL
  for (s in unique(post.prob$Subgrps)){
    hold.s = post.prob %>% filter(Subgrps==s)
    dat.s = with(density(hold.s$est), data.frame(x, y))
    dat.s = data.frame(id = s, dat.s)
    dat.dens = rbind(dat.dens, dat.s)
  }
  max.dens <- max(dat.dens$y)*1.05
  
  ## Plot setup ##
  add_vars_list <- list(N = "$node$info$N", est = "$node$info$label", 
                        estimand = "$node$info$estimand", 
                        prob.est = "$node$info$prob.est",
                        label = "$node$info$label")
  if (family=="survival"){
    add_vars_list[["events"]] = "$node$info$events"
  }
  
  plt <- ggparty(ct, horizontal = TRUE, add_vars = add_vars_list) +
    geom_edge() +
    geom_edge_label() +
    geom_node_label(line_list = list(aes(label = splitvar),
                                     aes(label = pval_convert(p.value))),
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
    plt <- plt + geom_node_label(
      line_list = list(
        aes(label = paste("Node ", id) ),
        aes(label = paste("N = ", N, " (", round(N/n*100,1), "%)", sep="")),
        aes(label = paste(estimand, " [", (1-alpha_s)*100,"% CI]", sep="")),
        aes(label = label),
        aes(label = paste(prob.label, "=", prob.est))
      ),
      line_gpar = list(list(size = 8, fontface="bold"),
                       list(size = 8),
                       list(size = 8, fontface="bold"),
                       list(size = 8),
                       list(size = 8)),
      ids = "terminal")
    # Outcome Plot #
    out_plot <- geom_node_plot(gglist =
                                 list(geom_boxplot(data=plot.dat,
                                                   aes(x=estimand, y=est, fill=estimand)),
                                      scale_x_discrete(expand=c(0.4, 0.4)), 
                                      xlab(""),
                                      ylab("Estimate"),
                                      theme_bw(), theme( axis.text.y = element_blank()),
                                      coord_flip()),
                               nudge_x = nudge_out, width = width_out,
                               shared_axis_labels = TRUE,
                               legend_separator = TRUE, scales = "fixed")
    # Density Plot #
    dens_plot <- geom_node_plot(
      gglist = list( geom_area(data=dat.dens, 
                               aes(x = ifelse(x < thres , x, thres),y=y), fill = fill_L),
                     geom_area(data=dat.dens, aes(x = ifelse(x > thres , x, thres),
                                                  y=y), fill = fill_R),
                     geom_line(data=dat.dens, aes(x=x,y=y)),
                     ylim(0,max.dens),
                     theme_bw(),
                     ggtitle(estimand),
                     theme(plot.title = element_text(size = 8, face = "bold")),
                     xlab("Density"), ylab("") ), shared_axis_labels = TRUE,
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
      geom_node_label(
        line_list = list(
          aes(label = paste("Node ", id) ),
          aes(label = paste("N = ", N, " (", round(N/n*100,1), "%)", sep="")),
          aes(label = paste("Events = ", events,
                            " (", round(events/n.events*100,1), "%)", sep="")),
          aes(label = paste(estimand, " [", (1-alpha_s)*100,"% CI]", sep="")),
          aes(label = label),
          aes(label = paste(prob.label, "=", prob.est))
        ),
        line_gpar = list(list(size = 8, fontface="bold"),
                         list(size = 8),
                         list(size = 8),
                         list(size = 8, fontface="bold"),
                         list(size = 8), 
                         list(size = 8)),
        ids = "terminal")
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
    if (is.null(out.train$A)) {
      gg_line <- geom_line(data=pred.surv,aes(x=time, y=surv), size=1.5)
    }
    if (!is.null(out.train$A)) {
      gg_line <- geom_line(data=pred.surv,aes(x=time, y=surv, col=A), size=1.5)
      pred.surv$A <- factor(pred.surv$A, levels = A_lvls)
    }
    out_plot <- geom_node_plot(gglist = list(gg_line,
                                             xlab("Time"),
                                             ylab( "Survival Probability" ),
                                             theme_bw(),
                                             theme(legend.title = element_blank())),
                               shared_axis_labels = TRUE,
                               legend_separator = TRUE, scales = "fixed", 
                               width=width_out, nudge_x=nudge_out)
    # Density #
    dens_plot <- geom_node_plot(
      gglist = list( geom_area(data=dat.dens, 
                               aes(x = ifelse(x < thres , x, thres),y=y), fill = fill_L),
                     geom_area(data=dat.dens, aes(x = ifelse(x > thres , x, thres),
                                                  y=y), fill = fill_R),
                     geom_line(data=dat.dens, aes(x=x,y=y)),
                     ylim(0,max.dens), theme_bw(), xlab("Density"), ylab(""),
                     ggtitle(estimand),
                     theme(plot.title = element_text(size = 8, face = "bold"))),
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