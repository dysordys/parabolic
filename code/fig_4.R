library(tidyverse)
library(patchwork)


dat <-
  read_rds("../data/evol_data.rds") |>
  summarize(numSpecies = sum(simplex + 2*duplex > 3e-5), .by = c(time, a, r)) |>
  select(time, resource = r, numTypes = numSpecies, assocRate = a) |>
  mutate(across(time | resource | assocRate, \(x) round(x, 5))) |>
  mutate(time = time / 1e6)

tmax <- max(dat$time)


p1 <- dat |>
  slice(1:400 * 100) |> # Reduce resolution (same quality, smaller file size)
  ggplot(aes(x = time, y = assocRate)) +
  geom_line(color = viridis::plasma(1)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(limits = c(0, NA), breaks = c(0, 50, 100)) +
  labs(x = NULL, y = "Association\nrate") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())

p2 <- dat |>
  slice(1:400 * 100) |> # Reduce resolution (same quality, smaller file size)
  ggplot(aes(x = time, y = 0, color = assocRate, fill = assocRate)) +
  geom_tile() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_gradient(low = "firebrick3", high = "skyblue1") +
  scale_fill_gradient(low = "firebrick3", high = "skyblue1") +
  labs(x = NULL, y = NULL) +
  guides(color = "none", fill = "none") +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  annotate(geom="text", x=0.175*tmax, y=0, label="WARM", color="black", size=3) +
  annotate(geom="text", x=0.500*tmax, y=0, label="COOL", color="black", size=3) +
  annotate(geom="text", x=0.825*tmax, y=0, label="WARM", color="black", size=3)

p3 <- dat |>
  ggplot(aes(x = time, y = numTypes)) +
  geom_line(color = viridis::plasma(1)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  labs(x = NULL, y = "Number of species") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())

p4 <- (dat |>
         slice(1:400 * 100) |> # Reduce resolution (same quality, smaller file size)
         ggplot(aes(x = time, y = 0, color = resource, fill = resource)) +
         geom_tile() +
         scale_x_continuous(expand = c(0, 0)) +
         scale_y_continuous(expand = c(0, 0)) +
         scale_color_viridis_c(begin = 0.7, direction = -1) +
         scale_fill_viridis_c(begin = 0.7, direction = -1) +
         labs(x = NULL, y = NULL) +
         guides(color = "none", fill = "none") +
         theme_minimal() +
         theme(axis.text = element_blank(), axis.ticks = element_blank())) |>
  (\(plt) reduce(1:20, \(p, i) {
    x <- seq(0.025*tmax, 0.975*tmax, l = 20)[i]
    lab <- ifelse(i %% 2 == 0, "S", "A")
    p + annotate(geom = "text", x = x, y = 0, label = lab, color = "black", size = 3)
  }, .init = plt))()


p5 <- dat |>
  slice(1:400 * 100) |> # Reduce resolution (same quality, smaller file size)
  ggplot(aes(x = time, y = resource)) +
  geom_line(color = viridis::plasma(1)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(0, 25, 50)) +
  labs(x = expression(paste("Time (", phantom() %*% 10^6, ")")),
       y = "Resource\nconcentration") +
  theme_bw() +
  theme(panel.grid = element_blank())


p <- (p1 / p2 / p3 / p4 / p5 +
        plot_layout(heights = c(1, 0.25, 4, 0.25, 1), tag_level = "new"))

show(p)
#ggsave("../figures/Fig4.pdf", device = cairo_pdf, width = 6, height = 7)
