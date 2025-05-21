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


paramList <- function(r, m, ai, bi, ci) {
  if ((length(ai)!=length(bi)) || (length(bi)!=length(ci))) stop("Input lengths differ")
  list(a = ai, b = bi, c = ci, r = r, m = m, S = length(ai))
}


standardParams <- function(r, m, expType = FALSE) {
  ai <- c(77.5, 72.5, 72.5, 77.5,  72.5, 82.5, 77.5,  72.5,  100,   87.5)
  bi <- c(10.0,  5.0,  6.0,  9.75,  7.0,  6.0,  8.75,  7.75,   7.0,  5.5)
  ci <- c( 4.4,  3.4,  5.0,  1.6,   3.6,  2.8,  2.0,   5.0,    4.0,  4.4)
  if (expType) ai[5] <- 0
  paramList(r, m, ai, bi, ci)
}


standardSeq <- function(from = log(9.668819e-05), to = log(1e2), l = 271) {
  exp(seq(from, to, length.out = l))
}



sol <-
  crossing(param = "m", r = 1, m = standardSeq(), expType = c(FALSE, TRUE)) |>
  bind_rows(crossing(param = "r", r = standardSeq(), m = 2, expType = c(FALSE, TRUE))) |>
  mutate(ic = list(rep(1 / 20, times = 20))) |>
  mutate(params = pmap(list(r, m, expType), standardParams)) |>
  mutate(sol = map2(ic, params, integrateEqs, .progress = TRUE)) |>
  unnest(sol)

sol |>
  summarize(conc = conc[type == "simplex"] + 2 * conc[type == "duplex"],
            .by = c(param, r, m, species, expType)) |>
  mutate(type = ifelse(species == "5" & expType, "E-species", "S-species")) |>
  mutate(paramValue = ifelse(param == "m", m, r)) |>
  select(!r & !m) |>
  mutate(species = as.character(species)) |>
  left_join(read_tsv("../data/growth_rates.tsv", col_types = "cd"),
            by = join_by(species)) |>
  mutate(numTypes = sum(conc > 3e-5), .by = c(param, expType, paramValue)) |>
  relocate(param, expType, paramValue, species, type, conc, numTypes) |>
  write_tsv("../data/alldata_reproduced.tsv")
