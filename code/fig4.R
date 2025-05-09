library(tidyverse)


read_tsv("../data/fig4_data.tsv", col_types = "ddid-",
         col_names = c("time", "resource", "numTypes", "assocRate")) |>
  mutate(time = time / 1e6) |>
  pivot_longer(!time) |>
  mutate(name = case_match(name,
                           "resource"  ~ "Resource concentration",
                           "numTypes"  ~ "No. of coexisting species",
                           "assocRate" ~ "Association rate")) |>
  ggplot(aes(x = time, y = value)) +
  geom_line(color = "steelblue") +
  labs(x = expression(paste("time (", phantom()%*%10^6, ")")), y = NULL) +
  facet_wrap(~ name, ncol = 1, strip.position = "left") +
  theme_bw() +
  theme(panel.grid = element_blank(), strip.text = element_text(size = 10.5),
        strip.placement = "outside", strip.background = element_blank())

#ggsave("../figures/new_fig4.pdf", device = cairo_pdf, width = 6, height = 7)
