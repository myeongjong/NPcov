###########################################################################
###
###    Author : MYEONGJONG KANG
###    E-mail : kmj.stat@gmail.com
###
###    Description : Visualization - regression with white noise errors
###
###########################################################################

library(NPcov)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)

###########################################################################
###
###########################################################################

rm(list = ls())
load(file = "out/simout_regression_with_white_noise_errors.RData")
rm(ptm1, ptm2)

###########################################################################
###
###########################################################################

idx     <- 1
x       <- result_wave_complex$history_data[[idx]]$x
y       <- result_wave_complex$history_data[[idx]]$y
y_err   <- result_wave_complex$history_data[[idx]]$y_err
fit_lc  <- result_wave_complex$history_fit[[idx]]$fit_lc
fit_ll  <- result_wave_complex$history_fit[[idx]]$fit_ll
fit_us  <- result_wave_complex$history_fit[[idx]]$fit_us

visout_cwave <- data.frame(type = c(rep("Local constant", length(as.numeric(t(fit_lc$eval)))),
                                    rep("Local linear", length(as.numeric(t(fit_ll$eval)))),
                                    rep("Our isotropic", length(fit_us$eval)),
                                    rep("Our monotone", length(fit_us$eval))),
                           x = c(as.numeric(t(fit_lc$eval)),
                                 as.numeric(t(fit_ll$eval)),
                                 fit_us$eval,
                                 fit_us$eval),
                           y = c(fit_lc$mean, fit_ll$mean, fit_us$yhat, rep(NA, length(fit_us$eval))))

visout_cwave$type <- factor(visout_cwave$type, levels = c("Our isotropic", "Our monotone", "Local constant", "Local linear"))

p1 <- ggplot() +
  geom_point(data = data.frame(x = x, y = y_err), aes(x = x, y = y), color = "gray75") +
  geom_line(data = data.frame(x = x, y = y), aes(x = x, y = y), color = "black", lwd = 1.25) +
  geom_line(data = visout_cwave, aes(x = x, y = y, col = type, lty = type), lwd = 1.2, alpha = 0.75) +
  scale_color_manual(name = NULL, labels = c("Our isotropic" = "Our isotropic", "Our monotone" = "Our monotone", "Local constant" = "Local constant", "Local linear" = "Local linear"), values = c("Our isotropic" = "#E41A1C", "Our monotone" = "#4DAF4A", "Local constant" = "#377EB8", "Local linear" = "#e0e002")) +
  scale_linetype_manual(name = NULL, labels = c("Our isotropic" = "Our isotropic", "Our monotone" = "Our monotone", "Local constant" = "Local constant", "Local linear" = "Local linear"), values = c("Our isotropic" = "twodash", "Our monotone" = "twodash", "Local constant" = "dashed", "Local linear" = "dotted")) +
  theme_bw() + xlab(NULL) + ylab(NULL) +
  theme(legend.direction = 'horizontal', legend.key.width = unit(2, "cm"))

rm(idx, x, y, y_err, fit_lc, fit_ll, fit_us)

###########################################################################
###
###########################################################################

idx     <- 6
x       <- result_spherical_complex$history_data[[idx]]$x
y       <- result_spherical_complex$history_data[[idx]]$y
y_err   <- result_spherical_complex$history_data[[idx]]$y_err
fit_lc  <- result_spherical_complex$history_fit[[idx]]$fit_lc
fit_ll  <- result_spherical_complex$history_fit[[idx]]$fit_ll
fit_us  <- result_spherical_complex$history_fit[[idx]]$fit_us

visout_csphe <- data.frame(type = c(rep("Local constant", length(as.numeric(t(fit_lc$eval)))),
                                    rep("Local linear", length(as.numeric(t(fit_ll$eval)))),
                                    rep("Our isotropic", length(fit_us$eval)),
                                    rep("Our monotone", length(fit_us$eval))),
                           x = c(as.numeric(t(fit_lc$eval)),
                                 as.numeric(t(fit_ll$eval)),
                                 fit_us$eval,
                                 fit_us$eval),
                           y = c(fit_lc$mean, fit_ll$mean, rep(NA, length(fit_us$eval)), fit_us$yhat))

visout_csphe$type <- factor(visout_csphe$type, levels = c("Our isotropic", "Our monotone", "Local constant", "Local linear"))

p2 <- ggplot() +
  geom_point(data = data.frame(x = x, y = y_err), aes(x = x, y = y), color = "gray75") +
  geom_line(data = data.frame(x = x, y = y), aes(x = x, y = y), color = "black", lwd = 1.25, alpha = 0.75) +
  geom_line(data = visout_csphe, aes(x = x, y = y, col = type, lty = type), lwd = 1.2, alpha = 0.75) +
  scale_color_manual(name = NULL, labels = c("Our isotropic" = "Our isotropic", "Our monotone" = "Our monotone", "Local constant" = "Local constant", "Local linear" = "Local linear"), values = c("Our isotropic" = "#E41A1C", "Our monotone" = "#4DAF4A", "Local constant" = "#377EB8", "Local linear" = "#e0e002")) +
  scale_linetype_manual(name = NULL, labels = c("Our isotropic" = "Our isotropic", "Our monotone" = "Our monotone", "Local constant" = "Local constant", "Local linear" = "Local linear"), values = c("Our isotropic" = "twodash", "Our monotone" = "twodash", "Local constant" = "dashed", "Local linear" = "dotted")) +
  theme_bw() + xlab(NULL) + ylab(NULL) +
  theme(legend.direction = 'horizontal', legend.key.width = unit(2, "cm"))

tmp       <- ggplot_gtable(ggplot_build(p1))
leg       <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
mylegend  <- tmp$grobs[[leg]]

### Merge the two plots ###

p.merge    <- grid.arrange(mylegend, arrangeGrob(p1 + theme(legend.position = "none"),
                                                 p2 + theme(legend.position = "none"),
                                                 ncol = 2), nrow = 2, heights = c(1, 10))

ggplot2::ggsave("img/simulation_regression.pdf", p.merge, width = 7.5, height = 3.5)

