library(ggplot2)


dat <- read.table("tssd.dat", sep = "|",
                  col.names = c("pos", "d", "d.smooth", "side"))
dat$sep <- ifelse(dat$side == "combined", "combined", "separate")


p <- ggplot(dat) +
  geom_point(aes(pos, d, col=side)) +
  geom_line(aes(pos, d.smooth, col=side)) +
  facet_wrap(~sep, ncol = 1, scale = "free_y") +
  scale_color_brewer(palette = "Set1") +
  theme_bw(16) + theme(legend.position = "top")
print(p)
