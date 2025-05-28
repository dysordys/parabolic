library(tidyverse)


parabDyn <- function(t, state, p) {
  x <- state[1:p$S]
  y <- state[(p$S + 1):(2*p$S)]
  phi <- sum(p$c*x) * p$r / p$m
  dxdt <- 2*p$b*y - 2*p$a*x^2 - p$r*p$c*x - phi*x
  dydt <- p$a*x^2 - p$b*y + p$r*p$c*x - phi*y
  list(c(dxdt, dydt))
}


integrateEqs <- function(params, tseq = seq(0, 1e6, l = 1001),
                         func = parabDyn, method = "bdf", ...) {
  sol <- deSolve::ode(params$initCond, tseq, func, params, method, ...)
  as_tibble(as.data.frame(sol)) |>
    pivot_longer(!time, names_to = "species", values_to = "conc") |>
    mutate(species = as.integer(species)) |>
    filter(time == max(time)) |>
    select(!time) |>
    mutate(type = ifelse(species <= params$S, "simplex", "duplex")) |>
    mutate(species = ((as.integer(species) - 1L) %% params$S + 1L))
}


paramList <- function(r, m, ai, bi, ci, initCond = NULL) {
  if ((length(ai)!=length(bi)) || (length(bi)!=length(ci))) stop("Input lengths differ")
  S <- length(ai)
  init <- if (length(initCond) != 2*S) rep(1, times = 2*S) else initCond
  list(a = ai, b = bi, c = ci, r = r, m = m, S = S, initCond = init)
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


addMutant <- function(params, a, b, c) {
  S <- params$S
  # New mutant with simplex concentration x = 1e-5 and duplex concentration y = 0:
  init <- c(params$initCond[1:S], 1e-5, params$initCond[(S + 1):(2*S)], 0)
  paramList(r = params$r, m = params$m, ai = c(params$a, a), bi = c(params$b, b),
            ci = c(params$c, c), initCond = init)
}


evoDyn <- function(params, iter, tseq = seq(0, 1e6, l = 1001)) {
  tibble(n = 0:iter, params = list(params))
}



sol <-
  crossing(param = "m", r = 1, m = standardSeq(), expType = c(FALSE, TRUE)) |>
  bind_rows(crossing(param = "r", r = standardSeq(), m = 2, expType = c(FALSE, TRUE))) |>
  mutate(params = pmap(list(r, m, expType), standardParams)) |>
  mutate(sol = map(params, integrateEqs, .progress = TRUE)) |>
  unnest(sol)

sol |>
  mutate(excessProd = sum(params[[1]]$c * conc[type == "simplex"]) *
           params[[1]]$r / params[[1]]$m,
         .by = c(param, r, m, expType)) |>
  summarize(conc = conc[type == "simplex"] + 2 * conc[type == "duplex"],
            .by = c(param, r, m, species, expType, excessProd)) |>
  mutate(type = ifelse(species == "5" & expType, "E-species", "S-species")) |>
  mutate(paramValue = ifelse(param == "m", m, r)) |>
  select(!r & !m) |>
  left_join(read_tsv("../data/growth_rates.tsv", col_types = "id"),
            by = join_by(species)) |>
  mutate(numTypes = sum(conc > 3e-5), .by = c(param, expType, paramValue)) |>
  relocate(param, expType, paramValue, species, type, conc, growthRate, excessProd) |>
  write_tsv("../data/alldata.tsv")
