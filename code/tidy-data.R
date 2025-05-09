library(tidyverse)


tidyTabs <- function(file, param = "r", expType = FALSE) {
  read_tsv(file, col_types = "d-didddddddddd-",
           col_names = c("paramValue", "excessProd", "numTypes", 1:10)) |>
    pivot_longer(matches("[[:digit:]]"), names_to = "species", values_to = "conc") |>
    mutate(type = ifelse(species == "5" & expType, "E-species", "S-species"))
}


tibble(file = Sys.glob("../data/dat_*")) |>
  mutate(param = str_sub(file, -1L)) |>
  mutate(expType = str_detect(file, "exp")) |>
  mutate(data = pmap(list(file, param, expType), tidyTabs)) |>
  unnest(data) |>
  select(!file) |>
  left_join(read_tsv("../data/growth_rates.tsv", col_types = "cd"),
            by = join_by(species)) |>
  mutate(species = as_factor(as.integer(species))) |>
  relocate(param, expType, paramValue, species, type, conc, growthRate, excessProd) |>
  arrange(param, expType, paramValue, species) |>
  write_tsv("../data/alldata_tidy.tsv")
