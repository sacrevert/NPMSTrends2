library(ggplot2)

outDF <-data.frame(mu_se)
outDF$Year <- 2015:2019
outDF$low <- outDF$Mean - outDF$Naive.SE
outDF$high <- outDF$Mean + outDF$Naive.SE
outDF$species <- "Achillea millefolium"
outDF$species <- "Holcus lanatus"
outDF$species <- "Pastinaca sativa"
outDF$species <- "Sonchus arvensis"
outDF$species <- "Gymnadenia conopsea"
#outDF$species <- "Trifolium campestre"

ggplot(outDF, aes_string(x = "Year", y = "Mean")) +
  theme_bw() +
  geom_ribbon(data = outDF,
              aes_string(group = 1, ymin = "low", ymax = "high"),
              alpha = 0.2) +
  geom_line(size = 1, col = "black") +
  #geom_point(size = 4, aes(col = rhat_threshold)) +
  #scale_color_manual(name = 'Rhat', values = c('Bad (>1.1)' = 'red','Good (<1.1)' = 'blue')) +
  ylab("Mean ZI cover") +
  xlab("Time period") +
  scale_y_continuous(limits = c(0, 0.20)) +
  ggtitle(outDF$species) + 
  theme(plot.title = element_text(lineheight = .8, face = "bold"),
        legend.position = 'bottom')