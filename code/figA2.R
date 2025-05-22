library(tidyverse)


read_tsv("../data/alldata_tidy.tsv", col_types = "cldicdddi") |>
  filter(conc > 0) |>
  mutate(
    type = fct_relevel(type, "S-species", "E-species"),
    param = as_factor(ifelse(param == "m",
                             "'Target replicator concentration,'~italic('m')",
                             "'Resource concentration,'~italic('r')")),
    expType = ifelse(expType, "With E-species", "Only S-species")
  ) |>
  ggplot(aes(x = paramValue, y = conc, color = growthRate,
             group = species, linetype = type)) +
  geom_line() +
  scale_x_log10(labels = scales::label_log()) +
  scale_color_viridis_c(option = "C", end = 0.9) +
  labs(x = NULL, y = "Relative concentration",
       linetype = NULL, color = expression(paste("Growth rate, ", lambda[i]))) +
  facet_grid(expType ~ param, labeller = labeller(.cols = label_parsed), switch = "x") +
  coord_cartesian(ylim = c(1e-4, NA)) +
  theme_bw() +
  theme(panel.grid = element_blank(), strip.placement = "outside")

#ggsave("../figures/Fig_A2.pdf", device = cairo_pdf, width = 8, height = 5)
