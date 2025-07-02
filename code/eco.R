library(tidyverse)


parabDyn <- function(t, state, p) {
  S <- length(p$species)
  x <- state[1:S]
  y <- state[(S + 1):(2*S)]
  phi <- sum(p$c*x) * p$r / p$m
  dxdt <- 2*p$b*y - 2*p$a*x^2 - p$r*p$c*x - phi*x
  dydt <- p$a*x^2 - p$b*y + p$r*p$c*x - phi*y
  list(c(dxdt, dydt))
}


integrateEqs <- function(params, tseq = seq(0, 1e6, l = 101),
                         func = parabDyn, method = "bdf", ...) {
  numSpecies <- nrow(params)
  initCond <- c(params$simplex, params$duplex)
  sol <- deSolve::ode(initCond, tseq, func, as.list(params), method, ...)
  as_tibble(as.data.frame(sol)) |>
    pivot_longer(!time, names_to = "species", values_to = "conc") |>
    filter(time == max(time)) |>
    select(!time) |>
    mutate(species = as.integer(species)) |>
    mutate(type = ifelse(species <= numSpecies, "simplex", "duplex")) |>
    mutate(species = ((as.integer(species) - 1L) %% numSpecies + 1L)) |>
    pivot_wider(names_from = type, values_from = conc) |>
    left_join(select(params, !simplex & !duplex), by = join_by(species))
}


paramTab <- function(r, m, a, b, c, initCond = NULL) {
  if ((length(a) != length(b)) || (length(b) != length(c))) stop("Input lengths differ")
  S <- length(a)
  init <- if (length(initCond) != 2*S) rep(1, times = 2*S) else initCond
  tibble(species = 1:S, simplex = init[1:S], duplex = init[(S+1):(2*S)],
         a = a, b = b, c = c, r = r, m = m)
}


standardParams <- function(r, m, expType = FALSE) {
  a <- c(77.5, 72.5, 72.5, 77.5,  72.5, 82.5, 77.5,  72.5,  100,   87.5)
  b <- c(10.0,  5.0,  6.0,  9.75,  7.0,  6.0,  8.75,  7.75,   7.0,  5.5)
  c <- c( 4.4,  3.4,  5.0,  1.6,   3.6,  2.8,  2.0,   5.0,    4.0,  4.4)
  if (expType) a[5] <- 0
  paramTab(r, m, a, b, c)
}


standardSeq <- function(from = log(9.668819e-05), to = log(1e2), l = 271) {
  exp(seq(from, to, length.out = l))
}



eco <-
  crossing(param = "m", r = 1, m = standardSeq(), expType = c(FALSE, TRUE)) |>
  bind_rows(crossing(param = "r", r = standardSeq(), m = 2, expType = c(FALSE, TRUE))) |>
  mutate(params = pmap(list(r, m, expType), standardParams)) |>
  mutate(sol = map(params, integrateEqs, .progress = TRUE)) |>
  select(!r & !m) |>
  unnest(sol)

eco |>
  mutate(growthRate = (sqrt(b^2 + c^2*r[1]^2 + 6*b*c*r[1]) - (b + c*r[1])) / 2) |>
  mutate(excessProd = sum(c * simplex) * r[1] / m[1],
         .by = c(param, r, m, expType)) |>
  summarize(conc = simplex + 2*duplex,
            .by = c(param, r, m, species, expType, excessProd, growthRate)) |>
  mutate(type = ifelse(species == "5" & expType, "E-species", "S-species")) |>
  mutate(paramValue = ifelse(param == "m", m, r)) |>
  select(!r & !m) |>
  mutate(numTypes = sum(conc > 3e-5), .by = c(param, expType, paramValue)) |>
  relocate(param, expType, paramValue, species, type, conc, growthRate, excessProd) |>
  write_rds("../data/eco_data.rds", compress = "xz")
