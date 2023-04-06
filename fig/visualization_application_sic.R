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
  theme_bw() + xlab("coordinate 1") + ylab("coordinate 2") +
  guides(size = guide_legend(title="(1/10th of a mm)", title.position = "right", nrow = 2, byrow = TRUE)) +
  theme(legend.text = element_text(size = 8), legend.title = element_text(size = 8), legend.direction = 'horizontal', legend.position = "top", legend.spacing.x = unit(0.1, 'cm'))

###########################################################################
###
###########################################################################

load(file = "out/appout_sic100_gaussian.RData")

visdat <- data.frame(type = c(rep("Our isotropic", length(bestout_iso$eval)),
                              rep("Our monotone", length(bestout_mon$eval)),
                              rep("Shapiro-Botha", length(sbout$svm$fit$u))),
                     x = c(bestout_iso$eval / output_sic100$xmax * output_sic100$xmax0,
                           bestout_mon$eval / output_sic100$xmax * output_sic100$xmax0,
                           sbout$svm$fit$u),
                     y = c(bestout_iso$yhat * vhat_iso,
                           bestout_mon$yhat * vhat_mon,
                           sbout$svm$sill - sbout$svm$fit$fitted.sv))

formatting <- function(l) {

  l <- format(l, scientific = TRUE) # turn in to character string in scientific notation
  # l <- gsub("^(.*)e", "'\\1'e", l) # quote the part before the exponent to keep all the digits
  # l <- gsub("e", "%*%10^", l) # turn the 'e+' into plotmath format
  # parse(text = l) # return this as an expression
}

p3 <- ggplot() +
  geom_point(data = data.frame(x = output_sic100$r / output_sic100$xmax * output_sic100$xmax0, y = output_sic100$chat), aes(x = x, y = y), color = "gray75") +
  geom_line(data = visdat, aes(x = x, y = y, col = type, lty = type), lwd = 1.2, alpha = 0.9) +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_color_manual(name = NULL, labels = c("Our isotropic" = "Our isotropic", "Our monotone" = "Our monotone", "Shapiro-Botha" = "Shapiro-Botha"), values = c("Our isotropic" = "#E41A1C", "Our monotone" = "#4DAF4A", "Shapiro-Botha" = "#377EB8")) +
  scale_linetype_manual(name = NULL, labels = c("Our isotropic" = "Our isotropic", "Our monotone" = "Our monotone", "Shapiro-Botha" = "Shapiro-Botha"), values = c("Our isotropic" = "twodash", "Our monotone" = "twodash", "Shapiro-Botha" = "dotted")) +
  # scale_y_continuous(labels = formatting) +
  theme_bw() + xlab("distance") + ylab("covariance") +
  coord_cartesian(ylim = c(-0.25, 1.05) * output_sic100$vhat) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.direction = 'horizontal', legend.key.width = unit(1.5, "cm"), legend.position = "top", legend.justification = "right")

p.merge       <- grid.arrange(p1, p3, nrow = 1)
ggplot2::ggsave("img/application_sic100.pdf", p.merge, width = 7.5, height = 3.5)
