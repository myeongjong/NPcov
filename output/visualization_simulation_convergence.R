###########################################################################
###
###    Author : MYEONGJONG KANG
###    E-mail : kmj.stat@gmail.com
###
###    Description : Visualization
###
###########################################################################

library(NPcov)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)

###########################################################################
### simulation 1: simple wave
###########################################################################

rm(list = ls())
set.seed(10182022)

x           <- sort(runif(200, min = 0.1, max = 10))
y           <- ftn_wave(x, 0, 1, 1)
y           <- 1 - y
y_err       <- y + rnorm(length(y), 0, 0.2)

load(file = "rdata/simout_gaussian_wave_simple.RData")

#####

# plot(x, y_err, xlab = "x", ylab = "y")
# lines(x, y, col = "black")
# lines(bestout$eval, bestout$yhat, col = "red")

#####

visout1     <- data.frame(x = bestout$eval, y = bestout$yhat)
p1          <- ggplot() +
  geom_point(data = data.frame(x = x, y = y_err), aes(x = x, y = y), color = "gray75") +
  geom_line(data = data.frame(x = x, y = y), aes(x = x, y = y), color = "black", lwd = 1.25, alpha = 0.75) +
  geom_line(data = visout1, aes(x = x, y = y), color = "#E41A1C", lty = "twodash", lwd = 1.2) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 10.25)) +
  theme_bw() + xlab("x") + ylab("y") +
  theme(legend.direction = 'horizontal', legend.key.width = unit(2, "cm"))

visout2     <- rbind(data.frame(type = "min", x = seq(bestout$iteration + 1), y = bestout$history_obj_min),
                     data.frame(type = "maxsel", x = seq(bestout$iteration + 1), y = bestout$history_obj_max_sel),
                     data.frame(type = "mean", x = seq(bestout$iteration + 1), y = bestout$history_obj_mean),
                     data.frame(type = "maxall", x = seq(bestout$iteration + 1), y = bestout$history_obj_max))
p2          <- ggplot() +
  geom_line(data = visout2, aes(x = x, y = y, lty = type), color = "black", lwd = 1) +
  scale_linetype_manual(values = c("min" = "twodash", "maxsel" = "solid", "mean" = "dashed", "maxall" = "dotted")) +
  theme_bw() + xlab("iteration") + ylab("objective value") + coord_cartesian(ylim = range(visout2$y[visout2$type != "maxall"])) +
  theme(legend.position = "none")

visout3     <- data.frame(x = seq(bestout$iteration + 1), y = bestout$history_KL)
p3          <- ggplot() +
  geom_line(data = visout3, aes(x = x, y = y), color = "black", lwd = 1) +
  theme_bw() + xlab("iteration") + ylab("KL divergence") +
  theme(legend.position = "none")

p.merge     <- grid.arrange(p1, p2, p3, nrow = 1)

ggplot2::ggsave("plot/simulation_convergence_wave_simple.pdf", p.merge, device = cairo_pdf, width = 7.5, height = 3.5)

###########################################################################
### simulation 2: complex wave
###########################################################################

rm(list = ls())
set.seed(10182022)

x           <- sort(runif(200, min = 0.1, max = 10))
y           <- ftn_wave(x, 0, 1, 0.5)
y           <- 1 - y
y_err       <- y + rnorm(length(y), 0, 0.2)

load(file = "rdata/simout_gaussian_wave_complex.RData")

#####

# plot(x, y_err, xlab = "x", ylab = "y")
# lines(x, y, col = "black")
# lines(bestout$eval, bestout$yhat, col = "red")

#####

visout1     <- data.frame(x = bestout$eval, y = bestout$yhat)
p1          <- ggplot() +
  geom_point(data = data.frame(x = x, y = y_err), aes(x = x, y = y), color = "gray75") +
  geom_line(data = data.frame(x = x, y = y), aes(x = x, y = y), color = "black", lwd = 1.25, alpha = 0.75) +
  geom_line(data = visout1, aes(x = x, y = y), color = "#E41A1C", lty = "twodash", lwd = 1.2) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 10.25)) +
  theme_bw() + xlab("x") + ylab("y") +
  theme(legend.direction = 'horizontal', legend.key.width = unit(2, "cm"))

visout2     <- rbind(data.frame(type = "min", x = seq(bestout$iteration + 1), y = bestout$history_obj_min),
                     data.frame(type = "maxsel", x = seq(bestout$iteration + 1), y = bestout$history_obj_max_sel),
                     data.frame(type = "mean", x = seq(bestout$iteration + 1), y = bestout$history_obj_mean),
                     data.frame(type = "maxall", x = seq(bestout$iteration + 1), y = bestout$history_obj_max))
p2          <- ggplot() +
  geom_line(data = visout2, aes(x = x, y = y, lty = type), color = "black", lwd = 1) +
  scale_linetype_manual(values = c("min" = "twodash", "maxsel" = "solid", "mean" = "dashed", "maxall" = "dotted")) +
  theme_bw() + xlab("iteration") + ylab("objective value") + coord_cartesian(ylim = range(visout2$y[visout2$type != "maxall"])) +
  theme(legend.position = "none")

visout3     <- data.frame(x = seq(bestout$iteration + 1), y = bestout$history_KL)
p3          <- ggplot() +
  geom_line(data = visout3, aes(x = x, y = y), color = "black", lwd = 1) +
  theme_bw() + xlab("iteration") + ylab("KL divergence") +
  theme(legend.position = "none")

p.merge     <- grid.arrange(p1, p2, p3, nrow = 1)

ggplot2::ggsave("plot/simulation_convergence_wave_complex.pdf", p.merge, device = cairo_pdf, width = 7.5, height = 3.5)

###########################################################################
### simulation 3: simple spherical case
###########################################################################

rm(list = ls())
set.seed(10182022)

x           <- sort(runif(200, min = 0.1, max = 10))
y           <- ftn_spherical(x, 0, 1, 1)
y           <- 1 - y
y_err       <- y + rnorm(length(y), 0, 0.2)

load(file = "rdata/simout_gaussian_sphe_simple.RData")

#####

# plot(x, y_err, xlab = "x", ylab = "y")
# lines(x, y, col = "black")
# lines(bestout$eval, bestout$yhat, col = "red")

#####

visout1     <- data.frame(x = bestout$eval, y = bestout$yhat)
p1          <- ggplot() +
  geom_point(data = data.frame(x = x, y = y_err), aes(x = x, y = y), color = "gray75") +
  geom_line(data = data.frame(x = x, y = y), aes(x = x, y = y), color = "black", lwd = 1.25, alpha = 0.75) +
  geom_line(data = visout1, aes(x = x, y = y), color = "#4DAF4A", lty = "twodash", lwd = 1.2) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 10.25)) +
  theme_bw() + xlab("x") + ylab("y") +
  theme(legend.direction = 'horizontal', legend.key.width = unit(2, "cm"))

visout2     <- rbind(data.frame(type = "min", x = seq(bestout$iteration + 1), y = bestout$history_obj_min),
                     data.frame(type = "maxsel", x = seq(bestout$iteration + 1), y = bestout$history_obj_max_sel),
                     data.frame(type = "mean", x = seq(bestout$iteration + 1), y = bestout$history_obj_mean),
                     data.frame(type = "maxall", x = seq(bestout$iteration + 1), y = bestout$history_obj_max))
p2          <- ggplot() +
  geom_line(data = visout2, aes(x = x, y = y, lty = type), color = "black", lwd = 1) +
  scale_linetype_manual(values = c("min" = "twodash", "maxsel" = "solid", "mean" = "dashed", "maxall" = "dotted")) +
  theme_bw() + xlab("iteration") + ylab("objective value") + coord_cartesian(ylim = range(visout2$y[visout2$type != "maxall"])) +
  theme(legend.position = "none")

visout3     <- data.frame(x = seq(bestout$iteration + 1), y = bestout$history_KL)
p3          <- ggplot() +
  geom_line(data = visout3, aes(x = x, y = y), color = "black", lwd = 1) +
  theme_bw() + xlab("iteration") + ylab("KL divergence") +
  theme(legend.position = "none")

p.merge     <- grid.arrange(p1, p2, p3, nrow = 1)

ggplot2::ggsave("plot/simulation_convergence_sphe_simple.pdf", p.merge, device = cairo_pdf, width = 7.5, height = 3.5)

###########################################################################
### simulation 4: complex spherical case
###########################################################################

rm(list = ls())
set.seed(10182022)

x           <- sort(runif(200, min = 0.1, max = 10))
y           <- ftn_spherical(x, 0, 0.8, 2)
y           <- 1 - y
y_err       <- y + rnorm(length(y), 0, 0.2)

load(file = "rdata/simout_gaussian_sphe_complex.RData")

#####

# plot(x, y_err, xlab = "x", ylab = "y")
# lines(x, y, col = "black")
# lines(bestout$eval, bestout$yhat, col = "red")

#####

visout1     <- data.frame(x = bestout$eval, y = bestout$yhat)
p1          <- ggplot() +
  geom_point(data = data.frame(x = x, y = y_err), aes(x = x, y = y), color = "gray75") +
  geom_line(data = data.frame(x = x, y = y), aes(x = x, y = y), color = "black", lwd = 1.25, alpha = 0.75) +
  geom_line(data = visout1, aes(x = x, y = y), color = "#4DAF4A", lty = "twodash", lwd = 1.2) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 10.25)) +
  theme_bw() + xlab("x") + ylab("y") +
  theme(legend.direction = 'horizontal', legend.key.width = unit(2, "cm"))

visout2     <- rbind(data.frame(type = "min", x = seq(bestout$iteration + 1), y = bestout$history_obj_min),
                     data.frame(type = "maxsel", x = seq(bestout$iteration + 1), y = bestout$history_obj_max_sel),
                     data.frame(type = "mean", x = seq(bestout$iteration + 1), y = bestout$history_obj_mean),
                     data.frame(type = "maxall", x = seq(bestout$iteration + 1), y = bestout$history_obj_max))
p2          <- ggplot() +
  geom_line(data = visout2, aes(x = x, y = y, lty = type), color = "black", lwd = 1) +
  scale_linetype_manual(values = c("min" = "twodash", "maxsel" = "solid", "mean" = "dashed", "maxall" = "dotted")) +
  theme_bw() + xlab("iteration") + ylab("objective value") + coord_cartesian(ylim = range(visout2$y[visout2$type != "maxall"]) + c(0, 5)) +
  theme(legend.position = "none")

visout3     <- data.frame(x = seq(bestout$iteration + 1), y = bestout$history_KL)
p3          <- ggplot() +
  geom_line(data = visout3, aes(x = x, y = y), color = "black", lwd = 1) +
  theme_bw() + xlab("iteration") + ylab("KL divergence") +
  theme(legend.position = "none")

p.merge     <- grid.arrange(p1, p2, p3, nrow = 1)

ggplot2::ggsave("plot/simulation_convergence_sphe_complex.pdf", p.merge, device = cairo_pdf, width = 7.5, height = 3.5)

###########################################################################
### simulation 5: simple exponential case
###########################################################################

rm(list = ls())
set.seed(10182022)

x           <- sort(runif(200, min = 0.1, max = 10))
y           <- ftn_exponential(x, 0, 1, 1)
y           <- 1 - y
y_err       <- y + rnorm(length(y), 0, 0.2)

load(file = "rdata/simout_gaussian_expo_simple.RData")

#####

# plot(x, y_err, xlab = "x", ylab = "y")
# lines(x, y, col = "black")
# lines(bestout$eval, bestout$yhat, col = "red")

#####

visout1     <- data.frame(x = bestout$eval, y = bestout$yhat)
p1          <- ggplot() +
  geom_point(data = data.frame(x = x, y = y_err), aes(x = x, y = y), color = "gray75") +
  geom_line(data = data.frame(x = x, y = y), aes(x = x, y = y), color = "black", lwd = 1.25, alpha = 0.75) +
  geom_line(data = visout1, aes(x = x, y = y), color = "#4DAF4A", lty = "twodash", lwd = 1.2) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 10.25)) +
  theme_bw() + xlab("x") + ylab("y") +
  theme(legend.direction = 'horizontal', legend.key.width = unit(2, "cm"))

visout2     <- rbind(data.frame(type = "min", x = seq(bestout$iteration + 1), y = bestout$history_obj_min),
                     data.frame(type = "maxsel", x = seq(bestout$iteration + 1), y = bestout$history_obj_max_sel),
                     data.frame(type = "mean", x = seq(bestout$iteration + 1), y = bestout$history_obj_mean),
                     data.frame(type = "maxall", x = seq(bestout$iteration + 1), y = bestout$history_obj_max))
p2          <- ggplot() +
  geom_line(data = visout2, aes(x = x, y = y, lty = type), color = "black", lwd = 1) +
  scale_linetype_manual(values = c("min" = "twodash", "maxsel" = "solid", "mean" = "dashed", "maxall" = "dotted")) +
  theme_bw() + xlab("iteration") + ylab("objective value") + coord_cartesian(ylim = range(visout2$y[visout2$type != "maxall"])) +
  theme(legend.position = "none")

visout3     <- data.frame(x = seq(bestout$iteration + 1), y = bestout$history_KL)
p3          <- ggplot() +
  geom_line(data = visout3, aes(x = x, y = y), color = "black", lwd = 1) +
  theme_bw() + xlab("iteration") + ylab("KL divergence") +
  theme(legend.position = "none")

p.merge     <- grid.arrange(p1, p2, p3, nrow = 1)

ggplot2::ggsave("plot/simulation_convergence_expo_simple.pdf", p.merge, device = cairo_pdf, width = 7.5, height = 3.5)

###########################################################################
### simulation 6: complex exponential case
###########################################################################

rm(list = ls())
set.seed(10182022)

x           <- sort(runif(200, min = 0.1, max = 10))
y           <- ftn_exponential(x, 0, 1, 4)
y           <- 1 - y
y_err       <- y + rnorm(length(y), 0, 0.2)

load(file = "rdata/simout_gaussian_expo_complex.RData")

#####

# plot(x, y_err, xlab = "x", ylab = "y")
# lines(x, y, col = "black")
# lines(bestout$eval, bestout$yhat, col = "red")

#####

visout1     <- data.frame(x = bestout$eval, y = bestout$yhat)
p1          <- ggplot() +
  geom_point(data = data.frame(x = x, y = y_err), aes(x = x, y = y), color = "gray75") +
  geom_line(data = data.frame(x = x, y = y), aes(x = x, y = y), color = "black", lwd = 1.25, alpha = 0.75) +
  geom_line(data = visout1, aes(x = x, y = y), color = "#4DAF4A", lty = "twodash", lwd = 1.2) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 10.25)) +
  theme_bw() + xlab("x") + ylab("y") +
  theme(legend.direction = 'horizontal', legend.key.width = unit(2, "cm"))

visout2     <- rbind(data.frame(type = "min", x = seq(bestout$iteration + 1), y = bestout$history_obj_min),
                     data.frame(type = "maxsel", x = seq(bestout$iteration + 1), y = bestout$history_obj_max_sel),
                     data.frame(type = "mean", x = seq(bestout$iteration + 1), y = bestout$history_obj_mean),
                     data.frame(type = "maxall", x = seq(bestout$iteration + 1), y = bestout$history_obj_max))
p2          <- ggplot() +
  geom_line(data = visout2, aes(x = x, y = y, lty = type), color = "black", lwd = 1) +
  scale_linetype_manual(values = c("min" = "twodash", "maxsel" = "solid", "mean" = "dashed", "maxall" = "dotted")) +
  theme_bw() + xlab("iteration") + ylab("objective value") + coord_cartesian(ylim = range(visout2$y[visout2$type != "maxall"])) +
  theme(legend.position = "none")

visout3     <- data.frame(x = seq(bestout$iteration + 1), y = bestout$history_KL)
p3          <- ggplot() +
  geom_line(data = visout3, aes(x = x, y = y), color = "black", lwd = 1) +
  theme_bw() + xlab("iteration") + ylab("KL divergence") +
  theme(legend.position = "none")

p.merge     <- grid.arrange(p1, p2, p3, nrow = 1)

ggplot2::ggsave("plot/simulation_convergence_expo_complex.pdf", p.merge, device = cairo_pdf, width = 7.5, height = 3.5)
