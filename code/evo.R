library(tidyverse)


assocFun = function(t) 100 - 90 * abs(t / 1e6 - 1)


resourceFun = function(t) 25 * (sin(2*pi*t / 2e5) + 1)


parabDyn <- function(t, state, p) {
  S <- length(p$species)
  a <- assocFun(t)
  r <- resourceFun(t)
  m <- p$m[1]
  x <- state[1:S]
  y <- state[(S + 1):(2*S)]
  phi <- sum(p$c*x) * r/m
  dxdt <- 2*p$b*y - 2*a*x^2 - r*p$c*x - phi*x
  dydt <- a*x^2 - p$b*y + r*p$c*x - phi*y
  list(c(dxdt, dydt))
}


integrateEqs <- function(params, tstart, tseq = seq(0.1, 100, by = 0.1),
                         func = parabDyn) {
  numSpecies <- nrow(params)
  initCond <- c(params$simplex, params$duplex)
  sol <- deSolve::ode(initCond, tstart + tseq, func, as.list(params))
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


removeExtinct <- function(params, threshold = 3e-5) {
  params |>
    filter(simplex + 2*duplex >= threshold) |>
    mutate(species = row_number())
}


addMutant <- function(params) {
  crossing(b = 5 + 0.25 * 0:20, c = 1 + 0.2 * 0:20) |>
    anti_join(select(params, b, c), by = join_by(b, c)) |>
    slice_sample(n = 1) |>
    mutate(species = max(params$species) + 1L, simplex = 1e-5, duplex = 0,
           a = params$a[1], r = params$r[1], m = params$m[1]) |>
    add_row(.data = params)
}


evoDyn <- function(params, timeline) {
  tibble(time = timeline, params = list(params)) |>
    mutate(params = accumulate2(params, time[-1], \(acc, p, t) {
      if (t %% 1000 == 0) cat(str_c("time: ", t, "\n"))
      acc |>
        addMutant() |>
        integrateEqs(t) |>
        removeExtinct()
    } ))
}



# Evolutionary simulations:
evo <- paramTab(r = 25, m = 2, a = 10, b = 7.5, c = 3) |>
  evoDyn(timeline = seq(0, 2e6, by = 100)) |>
  unnest(params) |>
  mutate(a = assocFun(time), r = resourceFun(time))

write_rds(evo, "../data/evo_data.rds", compress = "xz")
