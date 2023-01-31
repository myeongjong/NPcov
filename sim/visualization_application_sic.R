###########################################################################
###
###    Author : MYEONGJONG KANG
###    E-mail : kmj.stat@gmail.com
###
###    Description : Visualization
###
###########################################################################

library(geoR)
library(sf)
library(ggplot2)
library(gridExtra)
library(grid)

###########################################################################
###
###########################################################################

data(SIC)

p1          <- ggplot() +
  geom_sf(data = st_linestring(sic.borders), lwd = 2, col = "grey20") +
  geom_sf(data = st_as_sf(data.frame(x = sic.100$coords[, "V2"],
                                     y = sic.100$coords[, "V3"],
                                     z = sic.100$data),
                          coords = c("x", "y")),
          pch = 21, aes(size = z), fill = alpha("#377EB8", 0.7), col = "grey20") +
  theme_bw() + xlab(NULL) + ylab(NULL) +
  guides(size = guide_legend(title="(1/10th of a mm)", title.position = "right")) +
  theme(legend.text = element_text(size = 8), legend.title = element_text(size = 8), legend.direction = 'horizontal', legend.position = "top", legend.spacing.x = unit(0.1, 'cm'))

# p2          <- ggplot(data = data.frame(x = sic.100$data), ) +
#   geom_histogram(aes(x), bins = 12) + geom_freqpoly(aes(x), bins = 12)
#
# p.merge.1     <- grid.arrange(p1, p2, nrow = 1)

# ggplot2::ggsave("img/application_sic100.pdf", p1, width = 7.5, height = 3.5)

###########################################################################
###
###########################################################################

load(file = "out/appout_sic100_gaussian.RData")

visdat <- data.frame(type = c(rep("Our isotropic", length(bestout$eval)),
                              rep("Our monotone", length(bestout_mon$eval))),
                     x = c(bestout$eval / output_sic100$xmax * output_sic100$xmax0,
                           bestout_mon$eval / output_sic100$xmax * output_sic100$xmax0),
                     y = c(bestout$yhat * output_sic100$vhat,
                           bestout_mon$yhat * output_sic100$vhat))
p3 <- ggplot() +
  geom_point(data = data.frame(x = output_sic100$r / output_sic100$xmax * output_sic100$xmax0, y = output_sic100$chat * output_sic100$vhat), aes(x = x, y = y), color = "gray75") +
  geom_line(data = visdat, aes(x = x, y = y, col = type, lty = type), lwd = 1.2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_color_manual(name = NULL, labels = c("Our isotropic" = "Our isotropic", "Our monotone" = "Our monotone"), values = c("Our isotropic" = "#E41A1C", "Our monotone" = "#4DAF4A")) +
  scale_linetype_manual(name = NULL, labels = c("Our isotropic" = "Our isotropic", "Our monotone" = "Our monotone"), values = c("Our isotropic" = "solid", "Our monotone" = "solid")) +
  theme_bw() + xlab(NULL) + ylab(NULL) +
  coord_cartesian(ylim = c(-0.75, 1.25) * output_sic100$vhat) +
  theme(legend.direction = 'horizontal', legend.key.width = unit(1.5, "cm"), legend.position = "top")

p.merge       <- grid.arrange(p1, p3, nrow = 1)
ggplot2::ggsave("img/application_sic100.pdf", p.merge, width = 7.5, height = 3.5)
