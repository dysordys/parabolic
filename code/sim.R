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
    mutate(species = ((as.integer(species) - 1L) %% params$S + 1L))
}


paramList <- function(r, m,
                      ai = c(72.5, 72.5, 72.5, 72.5, 77.5, 77.5, 77.5, 82.5, 87.5, 100),
                      bi = c(5.0,  6.0,  7.0,  7.75, 8.75, 9.75, 10.0, 6.0,  5.5, 7.0),
                      ci = c(3.4,  5.0,  3.6,  5.0,  2.0,  1.6,  4.4,  2.8,  4.4, 4.0)) {
  if ((length(ai)!=length(bi)) || (length(bi)!=length(ci))) stop("Input lengths differ")
  list(a = ai, b = bi, c = ci, r = r, m = m, S = length(ai))
}



sol <- crossing(r = 1, m = exp(seq(log(9.668819e-05), log(1e2), l = 271))) |>
  mutate(ic = list(with_seed(259751L, runif(20, min = 0, max = 1)))) |>
  mutate(params = map2(r, m, paramList)) |>
  mutate(sol = map2(ic, params, integrateEqs, .progress = TRUE)) |>
  unnest(sol)

sol |>
  summarize(conc = conc[type == "simplex"] + 2*conc[type == "duplex"],
            .by = c(r, m, species)) |>
  filter(conc > 0) |>
  mutate(rel_conc = conc / sum(conc), .by = c(r, m)) |>
  ggplot(aes(x = m, y = rel_conc, group = species)) +
  geom_line(color = "steelblue") +
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(ylim = c(1e-4, NA)) +
  theme_bw()
