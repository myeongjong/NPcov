###########################################################################
###
###    Author : MYEONGJONG KANG
###    E-mail : kmj.stat@gmail.com
###
###    Description : Simulation - covariance function estimation
###
###########################################################################

library(NPcov)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(gridExtra)
library(grid)

###########################################################################
###
###########################################################################

rm(list = ls())
load(file = "rdata/simout_covariance_function_estimation.RData")

###########################################################################
###
###########################################################################

visdat_swave  <- data.frame(x = simgrf_swave$coords[, 1], y = simgrf_swave$coords[, 2], z = simgrf_swave$data)

p1            <- ggplot() +
  geom_point(data = visdat_swave, aes(x = x, y = y, col = z), size = 2) +
  viridis::scale_color_viridis(direction = 1, option = "turbo") +
  theme_bw() + xlab("spatial domain dimension 1") + ylab("spatial domain dimension 2") +
  theme(legend.title = element_blank(), legend.position = 'top', legend.direction = 'horizontal', legend.key.width = unit(1.25, "cm"))

visdat_swave  <- data.frame(type = c(rep("Our isotropic", length(fit_us_swave$eval)),
                                     rep("Shapiro-Botha", length(sbout_swave$svm$fit$u))),
                            x = c(fit_us_swave$eval, sbout_swave$svm$fit$u),
                            y = c(fit_us_swave$yhat * vhat_swave, sbout_swave$svm$sill - sbout_swave$svm$fit$fitted.sv))

p2            <- ggplot() +
  geom_point(data = data.frame(x = input_swave$x, y = input_swave$y_err), aes(x = x, y = y), color = "gray75") +
  geom_line(data = data.frame(x = input_swave$x, y = input_swave$y), aes(x = x, y = y), color = "black", lwd = 1.25, alpha = 0.75) +
  geom_line(data = visdat_swave, aes(x = x, y = y, col = type, lty = type), lwd = 1.2, alpha = 1) +
  scale_color_manual(name = NULL, labels = c("Our isotropic" = "Our isotropic", "Shapiro-Botha" = "Shapiro-Botha"), values = c("Our isotropic" = "#E41A1C", "Shapiro-Botha" = "#377EB8")) +
  scale_linetype_manual(name = NULL, labels = c("Our isotropic" = "Our isotropic", "Shapiro-Botha" = "Shapiro-Botha"), values = c("Our isotropic" = "twodash", "Shapiro-Botha" = "dotted")) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  coord_cartesian(ylim = c(-0.25, 1.05) * sbout_swave$svm$sill) +
  theme_bw() + xlab("distance") + ylab("covariance") +
  theme(legend.position = "top", legend.direction = 'horizontal', legend.key.width = unit(2, "cm"))

###########################################################################
###
###########################################################################

visdat_sexpo  <- data.frame(x = simgrf_sexpo$coords[, 1], y = simgrf_sexpo$coords[, 2], z = simgrf_sexpo$data)

p3            <- ggplot() +
  geom_point(data = visdat_sexpo, aes(x = x, y = y, col = z), size = 2) +
  viridis::scale_color_viridis(direction = 1, option = "turbo") +
  theme_bw() + xlab("spatial domain dimension 1") + ylab("spatial domain dimension 2") +
  theme(legend.title = element_blank(), legend.position = 'top', legend.direction = 'horizontal', legend.key.width = unit(1.25, "cm"))

visdat_sexpo  <- data.frame(type = c(rep("Our monotone", length(fit_us_sexpo$eval)),
                                     rep("Shapiro-Botha", length(sbout_sexpo$svm$fit$u))),
                            x = c(fit_us_sexpo$eval, sbout_sexpo$svm$fit$u),
                            y = c(fit_us_sexpo$yhat * vhat_sexpo, sbout_sexpo$svm$sill - sbout_sexpo$svm$fit$fitted.sv))

p4            <- ggplot() +
  geom_point(data = data.frame(x = input_sexpo$x, y = input_sexpo$y_err), aes(x = x, y = y), color = "gray75") +
  geom_line(data = data.frame(x = input_sexpo$x, y = input_sexpo$y), aes(x = x, y = y), color = "black", lwd = 1.25, alpha = 0.75) +
  geom_line(data = visdat_sexpo, aes(x = x, y = y, col = type, lty = type), lwd = 1.2, alpha = 1) +
  scale_color_manual(name = NULL, labels = c("Our monotone" = "Our monotone", "Shapiro-Botha" = "Shapiro-Botha"), values = c("Our monotone" = "#4DAF4A", "Shapiro-Botha" = "#377EB8")) +
  scale_linetype_manual(name = NULL, labels = c("Our monotone" = "Our monotone", "Shapiro-Botha" = "Shapiro-Botha"), values = c("Our monotone" = "twodash", "Shapiro-Botha" = "dotted")) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  coord_cartesian(ylim = c(-0.8, 1.2) * vhat_sexpo) +
  theme_bw() + xlab("distance") + ylab("covariance") +
  theme(legend.position = "top", legend.direction = 'horizontal', legend.key.width = unit(2, "cm"))

###########################################################################
###
###########################################################################

p.merge.1     <- grid.arrange(p1, p2, nrow = 1, widths = c(2.75, 4.75))
ggplot2::ggsave("plot/simulation_covariance_wave.pdf", p.merge.1, device = cairo_pdf, width = 7.5, height = 3.5)

p.merge.2     <- grid.arrange(p3, p4, nrow = 1, widths = c(2.75, 4.75))
ggplot2::ggsave("plot/simulation_covariance_expo.pdf", p.merge.2, device = cairo_pdf, width = 7.5, height = 3.5)
