library(tidyverse)
library(patchwork)



dat <- read_tsv("../data/B3/alldata_tidy_B3.tsv", col_types = "cldicdddi") |>
  filter(conc > 0) |>
  mutate(
    species = as_factor(species),
    type = fct_relevel(type, "S-species", "E-species"),
    param = fct_rev(as_factor(ifelse(
      param == "D",
      "'Dilution rate,'~italic('D')",
      "'Resource concentration in the inflow,'~italic(rho)"
    ))),
    expType = ifelse(expType, "With E-species", "Only S-species")
  )


p1 <- dat |>
  ggplot(aes(x = paramValue, y = conc, color = growthRate,
             group = species, linetype = type)) +
  geom_line() +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_viridis_c(option = "C", end = 0.9) +
  labs(x = NULL, y = "Relative concentration",
       linetype = NULL, color = expression(paste("Growth rate, ", lambda[i]))) +
  facet_grid(expType ~ param, labeller = labeller(.cols = label_parsed),
             scales = "free_x") +
  coord_cartesian(ylim = c(1e-4, NA)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_blank())


p2 <- dat |>
  mutate(resourceConc = resourceConc * 2) |> # Because secondary y-axis will be rescaled
  pivot_longer(resourceConc | numTypes) |>
  mutate(name = fct_relevel(name, "numTypes", "resourceConc")) |>
  ggplot(aes(x = paramValue, y = value, color = name, linetype = expType)) +
  geom_line(linewidth = 0.6, alpha = 0.6) +
  facet_wrap(~ param, labeller = label_parsed,
             strip.position = "bottom", scales = "free_x") +
  scale_x_log10() +
  scale_y_continuous(
    name = "Number of species",
    breaks = (0:5) * 2,
    limits = c(0, 10),
    sec.axis = dup_axis(name = "Equilibrium resource concentration",
                        transform = ~ . / 2)
  ) +
  scale_color_manual(
    values = c("numTypes" = "#46A4E9", "resourceConc" = "darkgreen"),
    guide = "none"
  ) +
  labs(x = NULL, color = NULL, linetype = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.placement = "outside",
        axis.title.y.left = element_text(color = "#46A4E9"),
        axis.text.y.left = element_text(color = "#46A4E9"),
        axis.ticks.y.left = element_line(color = "#46A4E9"),
        axis.title.y.right = element_text(color = "darkgreen"),
        axis.text.y.right = element_text(color = "darkgreen"),
        axis.ticks.y.right = element_line(color = "darkgreen"))



p <- (p1 / p2 +
        plot_layout(heights = c(2, 1), tag_level = "new") +
        plot_annotation(tag_levels = list(c("A", "B"), "1")))

show(p)
#ggsave("../figures/Fig_B3.pdf", device = cairo_pdf, width = 8, height = 7)
