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
###
###########################################################################

rm(list = ls())
set.seed(10182022)

x           <- sort(runif(200, min = 0.1, max = 10))
y           <- ftn_wave(x, 0, 1, 0.5)
y           <- 1 - y
y_err       <- y + rnorm(length(y), 0, 0.2)

load(file = "out/simout_gaussian_wave_complex.RData")
cvout_gaussian_cwave <- cvout
bestout_gaussian_cwave <- bestout

load(file = "out/simout_epan_wave_complex.RData")
cvout_epan_cwave <- cvout
bestout_epan_cwave <- bestout

load(file = "out/simout_uniform_wave_complex.RData")
cvout_uniform_cwave <- cvout
bestout_uniform_cwave <- bestout

rm(cvout, bestout, ptm_cv, ptm_best)

###

plot(x, y_err, xlab = "x", ylab = "y")
lines(x, y, col = "black")
lines(bestout_uniform_cwave$eval, bestout_uniform_cwave$yhat, col = "green")
lines(bestout_epan_cwave$eval, bestout_epan_cwave$yhat, col = "blue")
lines(bestout_gaussian_cwave$eval, bestout_gaussian_cwave$yhat, col = "red")

###

visout_cwave <- data.frame(type = c(rep("Uniform", length(bestout_uniform_cwave$eval)),
                                    rep("Epanechnikov", length(bestout_epan_cwave$eval)),
                                    rep("Gaussian", length(bestout_gaussian_cwave$eval))),
                           x = c(bestout_uniform_cwave$eval,
                                 bestout_epan_cwave$eval,
                                 bestout_gaussian_cwave$eval),
                           y = c(bestout_uniform_cwave$yhat,
                                 bestout_epan_cwave$yhat,
                                 bestout_gaussian_cwave$yhat))

visout_cwave$type <- factor(visout_cwave$type, levels = c("Gaussian", "Epanechnikov", "Uniform"))

p1 <- ggplot() +
  geom_point(data = data.frame(x = x, y = y_err), aes(x = x, y = y), color = "gray75") +
  geom_line(data = data.frame(x = x, y = y), aes(x = x, y = y), color = "black", lwd = 1.25) +
  geom_line(data = visout_cwave, aes(x = x, y = y, col = type, lty = type), lwd = 1.2, alpha = 0.75) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 10.25)) +
  scale_color_manual(name = NULL, labels = c("Gaussian" = "Gaussian", "Epanechnikov" = "Epanechnikov", "Uniform" = "Uniform"), values = c("Gaussian" = "#E41A1C", "Epanechnikov" = "#4DAF4A", "Uniform" = "#377EB8")) +
  scale_linetype_manual(name = NULL, labels = c("Gaussian" = "Gaussian", "Epanechnikov" = "Epanechnikov", "Uniform" = "Uniform"), values = c("Gaussian" = "twodash", "Epanechnikov" = "dashed", "Uniform" = "dotted")) +
  theme_bw() + xlab("x") + ylab("y") +
  theme(legend.direction = 'horizontal', legend.key.width = unit(2, "cm"))

###########################################################################
###
###########################################################################

set.seed(10182022)

x           <- sort(runif(200, min = 0.1, max = 10))
y           <- ftn_spherical(x, 0, 0.8, 2)
y           <- 1 - y
y_err       <- y + rnorm(length(y), 0, 0.2)

load(file = "out/simout_gaussian_sphe_complex.RData")
cvout_gaussian_csphe <- cvout
bestout_gaussian_csphe <- bestout

load(file = "out/simout_epan_sphe_complex.RData")
cvout_epan_csphe <- cvout
bestout_epan_csphe <- bestout

load(file = "out/simout_uniform_sphe_complex.RData")
cvout_uniform_csphe <- cvout
bestout_uniform_csphe <- bestout

rm(cvout, bestout, ptm_cv, ptm_best)

###

plot(x, y_err, xlab = "x", ylab = "y")
lines(x, y, col = "black")
lines(bestout_uniform_csphe$eval, bestout_uniform_csphe$yhat, col = "green")
lines(bestout_epan_csphe$eval, bestout_epan_csphe$yhat, col = "blue")
lines(bestout_gaussian_csphe$eval, bestout_gaussian_csphe$yhat, col = "red")

###

visout_csphe <- data.frame(type = c(rep("Uniform", length(bestout_uniform_csphe$eval)),
                                    rep("Epanechnikov", length(bestout_epan_csphe$eval)),
                                    rep("Gaussian", length(bestout_gaussian_csphe$eval))),
                           x = c(bestout_uniform_csphe$eval,
                                 bestout_epan_csphe$eval,
                                 bestout_gaussian_csphe$eval),
                           y = c(bestout_uniform_csphe$yhat,
                                 bestout_epan_csphe$yhat,
                                 bestout_gaussian_csphe$yhat))

visout_csphe$type <- factor(visout_csphe$type, levels = c("Gaussian", "Epanechnikov", "Uniform"))

p2 <- ggplot() +
  geom_point(data = data.frame(x = x, y = y_err), aes(x = x, y = y), color = "gray75") +
  geom_line(data = data.frame(x = x, y = y), aes(x = x, y = y), color = "black", lwd = 1.25) +
  geom_line(data = visout_csphe, aes(x = x, y = y, col = type, lty = type), lwd = 1.2, alpha = 0.75) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 10.25)) +
  scale_color_manual(name = NULL, labels = c("Gaussian" = "Gaussian", "Epanechnikov" = "Epanechnikov", "Uniform" = "Uniform"), values = c("Gaussian" = "#E41A1C", "Epanechnikov" = "#377EB8", "Uniform" = "#4DAF4A")) +
  scale_linetype_manual(name = NULL, labels = c("Gaussian" = "Gaussian", "Epanechnikov" = "Epanechnikov", "Uniform" = "Uniform"), values = c("Gaussian" = "twodash", "Epanechnikov" = "dashed", "Uniform" = "dotted")) +
  theme_bw() + xlab("x") + ylab("y") +
  theme(legend.direction = 'horizontal', legend.key.width = unit(2, "cm"))

###########################################################################
###
###########################################################################

tmp       <- ggplot_gtable(ggplot_build(p1))
leg       <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
mylegend  <- tmp$grobs[[leg]]

### Merge the two plots ###

p.merge    <- grid.arrange(mylegend, arrangeGrob(p1 + theme(legend.position = "none"),
                                                 p2 + theme(legend.position = "none"),
                                                 ncol = 2), nrow = 2, heights = c(1, 10))

ggplot2::ggsave("img/simulation_kernels.pdf", p.merge, width = 7.5, height = 3.5)
