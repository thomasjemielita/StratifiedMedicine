### Subgroup Stability ###
# @fit: submod_train or PRISM output
subgrp_stability <- function(fit) {
  
  # fit <- res0
  if (inherits(fit, "PRISM")) {
    Subgrps0 <- fit$out.train$Subgrps
  }
  if (inherits(fit, "submod_train")) {
    Subgrps0 <- fit$Subgrps.train
  }
  uniq_subs <- unique(Subgrps0) 
  resamp_subgrps <- fit$resamp_subgrps
  id <- resamp_subgrps$id
  r_subs <- resamp_subgrps[,!(colnames(resamp_subgrps) %in% "id")]
  
  prob_dat <- NULL
  for (ss in uniq_subs) {
    hold <- r_subs
    hold[hold == ss] <- "1"
    hold[hold != "1"] <- "0"
    hold <- as.data.frame(sapply(hold, as.numeric))
    prob_hold <- rowMeans(hold, na.rm = TRUE, dims = 1)
    prob_hold <- data.frame(id=id, sub_0 = Subgrps0,
                            sub_b = ss, 
                            est = prob_hold)
    prob_dat <- rbind(prob_dat, prob_hold)
  }
  
  # Create summary Statistics #
  summ_dat <- NULL
  for (ss in uniq_subs) {
    
    id_s <- id[Subgrps0==ss]
    prob_s <- prob_dat[prob_dat$id %in% id_s,]
    hold <- NULL
    for (jj in uniq_subs) {
      hold_jj <- na.omit(prob_s$est[prob_s$sub_b==jj])
      mu_jj <- mean(hold_jj)
      min_jj <- min(hold_jj)
      max_jj <- max(hold_jj)
      lab_jj = paste( sprintf("%.2f", round(mu_jj,2)),
                      " [",
                      sprintf("%.2f", round(min_jj,2)), ",",
                      sprintf("%.2f", round(max_jj,2)), "]", sep="")
      hold <- cbind(hold, lab_jj)
    }
    hold <- data.frame(org_subgrp = ss, hold)
    summ_dat <- rbind(summ_dat, hold)
  }
  colnames(summ_dat)[-c(1)] <- paste("Prob", uniq_subs, sep="_")
  return(list(prob_dat=prob_dat, prob_summ = summ_dat))
}

plot_stability <- function(fit, type="density", display="instability") {
  
  sub_res <- subgrp_stability(fit)
  prob_dat <- sub_res$prob_dat
  plot_dat <- prob_dat
  
  if (display=="instability") {
    plot_dat <- plot_dat[plot_dat$sub_0!=plot_dat$sub_b,]
    xlabel <- "Prob(Instability)"
    title <- "Instability Probability Estimates by Original Subgroups"
  }
  if (display=="stability") {
    plot_dat <- plot_dat[plot_dat$sub_0==plot_dat$sub_b,]
    xlabel <- "Prob(Stability)"
    title <- "Stability Probability Estimates by Original Subgroups"
  }
  
  mu_dat <- aggregate(est ~ sub_b*sub_0, data=plot_dat, FUN="mean")
  res <- NULL
  if (type=="density") {
    caption_use <- "Note: Dashed lines correspond to average probability of switching to bootstrap subgroup"
    res = ggplot2::ggplot(plot_dat, ggplot2::aes_string(x="est", fill="sub_b")) +
      ggplot2::geom_density(alpha=0.30) +
      ggplot2::geom_vline(data=mu_dat, ggplot2::aes_string(xintercept="est", color="sub_b"),
                          linetype="dashed", size=1) + 
      ggplot2::facet_wrap(~sub_0) + 
      ggplot2::ggtitle(paste("Density Plot:", title, sep=" ")) +
      ggplot2::ylab("") + 
      ggplot2::xlab(xlabel)+
      ggplot2::theme(plot.title=ggplot2::element_text(size=16,face="bold"),
                     axis.text.y=ggplot2::element_blank(),
                     axis.ticks.y=ggplot2::element_blank(),
                     axis.text.x=ggplot2::element_text(face="bold"),
                     axis.title=ggplot2::element_text(size=12,face="bold"))+
      ggplot2::labs(fill = "Bootstrap Subgroup", color = "Bootstrap Subgroup",
                    caption = caption_use) + 
      ggplot2::theme_bw()
    
  }

  
  return(res)
} 

