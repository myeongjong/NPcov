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
load(file = "out/simout_covariance_function_estimation.RData")

###########################################################################
###
###########################################################################

visdat_swave  <- data.frame(x = simgrf_swave$coords[, 1], y = simgrf_swave$coords[, 2], z = simgrf_swave$data)

p1            <- ggplot() +
  geom_point(data = visdat_swave, aes(x = x, y = y, col = z), size = 2) +
  viridis::scale_color_viridis(direction = 1, option = "turbo") +
  theme_bw() + xlab(NULL) + ylab(NULL) +
  theme(legend.title = element_blank(), legend.position = 'top', legend.direction = 'horizontal', legend.key.width = unit(1.25, "cm"))

p2            <- ggplot() +
  geom_point(data = data.frame(x = input_swave$x, y = input_swave$y_err), aes(x = x, y = y), color = "gray75") +
  geom_line(data = data.frame(x = input_swave$x, y = input_swave$y), aes(x = x, y = y), color = "black", lwd = 1.25, alpha = 0.75) +
  geom_line(data = data.frame(x = fit_us_swave$eval, y = fit_us_swave$yhat * output_swave$vhat), aes(x = x, y = y), color = "#E41A1C", lwd = 1.2) +
  theme_bw() + xlab(NULL) + ylab(NULL) +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  coord_cartesian(ylim = c(-0.75, 1.25) * output_swave$vhat)

###########################################################################
###
###########################################################################

visdat_sexpo  <- data.frame(x = simgrf_sexpo$coords[, 1], y = simgrf_sexpo$coords[, 2], z = simgrf_sexpo$data)

p3            <- ggplot() +
  geom_point(data = visdat_sexpo, aes(x = x, y = y, col = z), size = 2) +
  viridis::scale_color_viridis(direction = 1, option = "turbo") +
  theme_bw() + xlab(NULL) + ylab(NULL) +
  theme(legend.title = element_blank(), legend.position = 'top', legend.direction = 'horizontal', legend.key.width = unit(1.25, "cm"))

p4            <- ggplot() +
  geom_point(data = data.frame(x = input_sexpo$x, y = input_sexpo$y_err), aes(x = x, y = y), color = "gray75") +
  geom_line(data = data.frame(x = input_sexpo$x, y = input_sexpo$y), aes(x = x, y = y), color = "black", lwd = 1.25, alpha = 0.75) +
  geom_line(data = data.frame(x = fit_us_sexpo$eval, y = fit_us_sexpo$yhat * output_sexpo$vhat), aes(x = x, y = y), color = "#E41A1C", lwd = 1.2) +
  theme_bw() + xlab(NULL) + ylab(NULL) +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  coord_cartesian(ylim = c(-0.75, 1.25) * output_sexpo$vhat)

###########################################################################
###
###########################################################################

p.merge.1     <- grid.arrange(p1, p2, nrow = 1, widths = c(2.75, 4.75))
ggplot2::ggsave("img/simulation_covariance_wave.pdf", p.merge.1, width = 7.5, height = 3.5)

p.merge.2     <- grid.arrange(p3, p4, nrow = 1, widths = c(2.75, 4.75))
ggplot2::ggsave("img/simulation_covariance_expo.pdf", p.merge.2, width = 7.5, height = 3.5)
