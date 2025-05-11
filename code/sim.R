library(tidyverse)
library(deSolve)
library(withr)


parabDyn <- function(t, state, p) {
  x <- state[1:p$S]
  y <- state[(p$S + 1):(2*p$S)]
  phi <- sum(p$c*x) * p$r / p$m
  dxdt <- 2*p$b*y - 2*p$a*x^2 - p$r*p$c*x - phi*x
  dydt <- p$a*x^2 - p$b*y + p$r*p$c*x - phi*y
  list(c(dxdt, dydt))
}


integrateEqs <- function(ic, params, tseq = seq(0, 1e6, l = 1001),
                         func = parabDyn, method = "bdf", ...) {
  sol <- ode(ic, tseq, func, params, method, ...)
  as_tibble(as.data.frame(sol)) |>
    pivot_longer(!time, names_to = "species", values_to = "conc") |>
    mutate(species = as.integer(species)) |>
    filter(time == max(time)) |>
    select(!time) |>
    mutate(type = ifelse(species <= params$S, "simplex", "duplex")) |>
    mutate(species = ((as.integer(species) - 1L) %% params$S + 1L)) |>
    summarize(conc = conc[type == "simplex"] + 2*conc[type == "duplex"], .by = species)
}


paramList <- function(r, m,
                      na = c(9,  9,  9,  9,  11, 11, 11, 13, 15, 20),
                      nb = c(20, 4,  8,  11, 15, 19, 20, 4,  2,  8),
                      nc = c(17, 20, 13, 20, 5,  3,  17, 8,  17, 15)) {
  if ((length(na)!=length(nb)) || (length(nb)!=length(nc))) stop("Input lengths differ")
  list(a = 50 + 2.5*na, b = 5 + 0.25*nb, c = 1 + 0.2*nc, r = r, m = m, S = length(na))
}



sol <- crossing(r = 1, m = exp(seq(log(9.668819e-05), log(1e2), l = 271))) |>
  mutate(ic = list(with_seed(259751L, runif(20, min = 0, max = 1)))) |>
  mutate(params = map2(r, m, paramList)) |>
  rowid_to_column("id") |>
  mutate(sol = pmap(list(ic, params, id), \(ic, params, id) {
    cat(str_c(id, ".\n"))
    integrateEqs(ic, params)
  } )) |>
  unnest(sol)

sol |>
  filter(conc > 0) |>
  mutate(rel_conc = conc / sum(conc), .by = c(id, r, m)) |>
  ggplot(aes(x = m, y = rel_conc, group = species)) +
  geom_line(color = "steelblue") +
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(ylim = c(1e-4, NA)) +
  theme_bw()
