library(tidyverse)
library(patchwork)


dat <- read_tsv("../data/fig4_data.tsv", col_types = "ddid-",
                col_names = c("time", "resource", "numTypes", "assocRate")) |>
  mutate(time = time / 1e6)


p1 <- dat |>
  ggplot(aes(x = time, y = assocRate)) +
  geom_line(color = viridis::plasma(1)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Association\nrate") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())

p2 <- dat |>
  ggplot(aes(x = time, y = 0, color = assocRate, fill = assocRate)) +
  geom_tile() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_gradient2(midpoint = 50) +
  scale_fill_gradient2(midpoint = 50) +
  labs(x = NULL, y = NULL) +
  guides(color = "none", fill = "none") +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  annotate(geom = "text", x = 0.25, y = 0, label = "warm", color = "black", size = 4) +
  annotate(geom = "text", x = 1.00, y = 0, label = "cool", color = "white", size = 4) +
  annotate(geom = "text", x = 1.75, y = 0, label = "warm", color = "black", size = 4)

p3 <- dat |>
  ggplot(aes(x = time, y = numTypes)) +
  geom_line(color = viridis::plasma(1)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "No. of coexisting species") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())

p4 <- (dat |>
         ggplot(aes(x = time, y = 0, color = resource, fill = resource)) +
         geom_tile() +
         scale_x_continuous(expand = c(0, 0)) +
         scale_y_continuous(expand = c(0, 0)) +
         scale_color_viridis_c(begin = 0.6, direction = -1) +
         scale_fill_viridis_c(begin = 0.6, direction = -1) +
         labs(x = NULL, y = NULL) +
         guides(color = "none", fill = "none") +
         theme_minimal() +
         theme(axis.text = element_blank(), axis.ticks = element_blank())) |>
  (\(plt) reduce(1:20, \(p, i) {
    x <- seq(0.05, 1.95, l = 20)[i]
    lab <- ifelse(i %% 2 == 0, "S", "A")
    p + annotate(geom = "text", x = x, y = 0, label = lab, color = "black", size = 4)
  }, .init = plt))()


p5 <- dat |>
  ggplot(aes(x = time, y = resource)) +
  geom_line(color = viridis::plasma(1)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Resource\nconcentration") +
  theme_bw() +
  theme(panel.grid = element_blank())


p <- (p1 / p2 / p3 / p4 / p5 +
        plot_layout(heights = c(1, 0.25, 4, 0.25, 1), tag_level = "new"))

show(p)
#ggsave("../figures/Fig4.pdf", device = cairo_pdf, width = 6, height = 7)
