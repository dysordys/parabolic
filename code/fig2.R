library(tidyverse)
library(patchwork)



dat <- read_tsv("../data/alldata_tidy.tsv", col_types = "cldicdddi") |>
  filter(conc > 0) |>
  mutate(
    species = as_factor(species),
    type = fct_relevel(type, "S-species", "E-species"),
    param = as_factor(ifelse(param == "m",
                             "'Varying target replicator concentration'~italic('m')",
                             "'Varying resource concentration'~italic('r')")),
    expType = ifelse(expType, "With E-species", "Only S-species")
  )


p1 <- dat |>
  ggplot(aes(x = paramValue, y = conc, color = type,
             group = species, alpha = growthRate)) +
  geom_line() +
  facet_grid(expType ~ param, labeller = labeller(.cols = label_parsed)) +
  scale_x_log10(labels = scales::label_log()) +
  scale_y_log10(labels = scales::label_log()) +
  scale_color_manual(values = c("S-species" = "steelblue",
                                "E-species" = "goldenrod")) +
  labs(x = "Parameter value", y = "Relative concentration",
       color = NULL, alpha = "Growth rate") +
  guides(alpha = guide_legend(reverse = TRUE)) +
  coord_cartesian(ylim = c(1e-4, NA)) +
  theme_bw() +
  theme(panel.grid = element_blank())


p2 <- dat |>
  pivot_longer(excessProd | numTypes) |>
  mutate(name = ifelse(name == "numTypes", "Number of coexisting species",
                       "Relative excess production")) |>
  ggplot(aes(x = paramValue, y = value, color = expType, linetype = name)) +
  geom_line(linewidth = 0.6, alpha = 0.6) +
  facet_wrap(~ param, labeller = label_parsed) +
  scale_x_log10(labels = scales::label_log()) +
  scale_y_continuous(name = expression(paste(
    "No. of coexisting species /\nRelative excess production, ", hat(phi)
  )), breaks = (0:5) * 2, limits = c(0, 10)) +
  scale_color_manual(values = c("Only S-species" = "steelblue",
                                "With E-species" = "goldenrod")) +
  labs(x = "Parameter value", color = NULL, linetype = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank())



(p1 / p2 +
    plot_layout(heights = c(2, 1), tag_level = "new") +
    plot_annotation(tag_levels = list(c("a)", "b)"), "1"))) |>
  ggsave(filename = "../figures/new_fig2.pdf",
         device = cairo_pdf, width = 8, height = 8)

