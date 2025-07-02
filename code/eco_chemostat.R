library(tidyverse)


parabDynChemo <- function(t, state, p) {
  S <- length(p$species)
  x <- state[1:S]
  y <- state[(S + 1):(2*S)]
  r <- state[2*S + 1]
  dxdt <- 2*p$b*y - 2*p$a*x^2 - p$c*r*x - p$D*x
  dydt <- p$a*x^2 - p$b*y + p$c*r*x - p$D*y
  drdt <- (p$rho[1] - r)*p$D[1] - r*sum(p$c*x)
  list(c(dxdt, dydt, drdt))
}


integrateEqs <- function(params, tseq = seq(0, 1e6, l = 101),
                         func = parabDynChemo, method = "bdf", ...) {
  numSpecies <- nrow(params)
  initCond <- c(params$simplex, params$duplex, params$resourceConc[1])
  sol <- deSolve::ode(initCond, tseq, func, as.list(params), method, ...)
  as_tibble(as.data.frame(sol)) |>
    pivot_longer(1:(2L*numSpecies) + 1L, names_to = "species", values_to = "conc") |>
    filter(time == max(time)) |>
    select(!time) |>
    mutate(species = as.integer(species)) |>
    mutate(type = ifelse(species <= numSpecies, "simplex", "duplex")) |>
    mutate(species = ((as.integer(species) - 1L) %% numSpecies + 1L)) |>
    rename(resourceConc = matches("\\d+")) |>
    pivot_wider(names_from = type, values_from = conc) |>
    left_join(select(params, !simplex & !duplex & !resourceConc), by = join_by(species))
}


paramTab <- function(D, rho, a, b, c, initCond = NULL) {
  if ((length(a) != length(b)) || (length(b) != length(c))) stop("Input lengths differ")
  S <- length(a)
  init <- if (length(initCond) != 2L*S + 1L) rep(1, times = 2*S + 1) else initCond
  tibble(species = 1:S, simplex = init[1:S], duplex = init[(S+1):(2*S)],
         resourceConc = init[2*S + 1], a = a, b = b, c = c, D = D, rho = rho)
}


standardParams <- function(D, rho, expType = FALSE) {
  a <- c(77.5, 72.5, 72.5, 77.5,  72.5, 82.5, 77.5,  72.5,  100,   87.5)
  b <- c(10.0,  5.0,  6.0,  9.75,  7.0,  6.0,  8.75,  7.75,   7.0,  5.5)
  c <- c( 4.4,  3.4,  5.0,  1.6,   3.6,  2.8,  2.0,   5.0,    4.0,  4.4)
  if (expType) a[5] <- 0
  paramTab(D, rho, a, b, c)
}



eco <-
  bind_rows(
    crossing(param = "rho", D = 1, rho = seq(0.26, 100, by = 0.01),
             expType = c(FALSE, TRUE)),
    crossing(param = "D", D = exp(seq(-3.229, 1.352, by = 0.01613)),
             rho = 2, expType = c(FALSE, TRUE))
  ) |>
  mutate(params = pmap(list(D, rho, expType), standardParams)) |>
  mutate(sol = map(params, integrateEqs, .progress = TRUE)) |>
  select(!D & !rho) |>
  unnest(sol)

eco |>
  mutate(growthRate = (sqrt(b^2 + c^2*resourceConc[1]^2 + 6*b*c*resourceConc[1]) -
                         (b + c*resourceConc[1])) / 2) |>
  mutate(conc = simplex + 2*duplex) |>
  mutate(type = ifelse(species == "5" & expType, "E-species", "S-species")) |>
  mutate(paramValue = ifelse(param == "D", D, rho)) |>
  mutate(numTypes = sum(conc > 3e-5),
         .by = c(param, expType, paramValue, resourceConc)) |>
  select(param, expType, paramValue, species, type, conc, growthRate,
         resourceConc, numTypes) |>
  write_tsv("../data/eco_data_chemostat.tsv")
