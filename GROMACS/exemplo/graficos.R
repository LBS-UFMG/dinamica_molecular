library(ggplot2)

rmsd = read.table(file = "rmsd.xvg", header = TRUE, sep = "", dec = ".")
ggplot(rmsd, aes(x = x, y = y)) + 
  geom_line(size = 0.75) + 
  xlab("Tempo (ns)") +
  ylab("RMSD (nm)") +
  labs(title ="myo.d_50ns.A1") +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14, colour = "black"),
        panel.border = element_rect(size = 2, fill = NA),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        axis.ticks = element_line(colour = "black", size = 1),
        axis.ticks.length.y = unit(.25, "cm"),
        axis.ticks.length.x = unit(.25, "cm"),
        plot.margin = margin(0.5,1,0.5,0.5, "cm"),
        plot.title = element_text(size = rel(2), hjust = 0.5)) +
  coord_cartesian(ylim = c(0, 0.7))
ggsave("rmsd.png", dpi = 600, height = 5, width = 8)

rmsf = read.table(file = "rmsf_residue.xvg", header = TRUE, sep = "", dec = ".")
ggplot(rmsf, aes(x = x, y = y)) + 
  geom_line(size = 0.75) +
  xlab("Res?duo") +
  ylab("RMSF (nm)") +
  labs(title ="myo.d_50ns.A1") +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14, colour = "black"),
        panel.border = element_rect(size = 2, fill = NA),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        axis.ticks = element_line(colour = "black", size = 1),
        axis.ticks.length.y = unit(.25, "cm"),
        axis.ticks.length.x = unit(.25, "cm"),
        plot.margin = margin(0.5,1,0.5,0.5, "cm"),
        plot.title = element_text(size = rel(2), hjust = 0.5)) +
  coord_cartesian(ylim = c(0, 0.7))
ggsave("rmsf.png", dpi = 600, height = 5, width = 8)
