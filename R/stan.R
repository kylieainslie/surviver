#' Build Stan input for the time-varying VE model
#'
#' Aggregates person-days into cells of age group x vaccination status x day
#' x day since vaccination, and attaches the force-of-infection offset.
#'
#' Days since vaccination are binned by day: bin 1 is unvaccinated, bin
#' `t + 1` is day `t` after vaccination. Days beyond `max_days_since_vax`
#' share the last bin.
#'
#' @param person_days Output of [expand_person_days()].
#' @param foi Output of [force_of_infection()].
#' @param max_days_since_vax Last day since vaccination with its own bin.
#'   Defaults to the longest observed.
#'
#' @return A list with `data` (the list passed to Stan) and `cells` (the
#'   aggregated data frame, for inspection).
#' @export
make_stan_data <- function(person_days, foi,
                           max_days_since_vax = NULL) {
  if (is.null(max_days_since_vax)) {
    max_days_since_vax <- max(c(1, person_days$days_since_vax), na.rm = TRUE)
  }
  cells <- person_days |>
    dplyr::mutate(
      tvax_bin = dplyr::if_else(
        .data$vax_status == 1L,
        as.integer(pmin(.data$days_since_vax, max_days_since_vax)) + 1L,
        1L
      )
    ) |>
    dplyr::group_by(.data$age_group, .data$vax_status, .data$day,
                    .data$tvax_bin) |>
    dplyr::summarise(Y = sum(.data$event), n_persdays = dplyr::n(),
                     .groups = "drop") |>
    dplyr::left_join(foi, by = c("day", "age_group"))

  if (anyNA(cells$Ibar)) {
    stop("`foi` does not cover every day and age group in `person_days`.",
         call. = FALSE)
  }

  data <- list(
    N            = nrow(cells),
    K            = nlevels(cells$age_group),
    T_vax        = as.integer(max_days_since_vax) + 1L,
    Y            = as.integer(cells$Y),
    log_Ibar     = log(cells$Ibar),
    log_persdays = log(cells$n_persdays),
    vax          = as.integer(cells$vax_status),
    tvax_bin     = as.integer(cells$tvax_bin)
  )
  list(data = data, cells = cells)
}

#' Path to a bundled Stan model
#'
#' @param model `"noncentered"` (default) or `"centered"` parameterisation of
#'   the random walk on log hazard ratios.
#' @return File path to the `.stan` file.
#' @export
ve_stan_file <- function(model = c("noncentered", "centered")) {
  model <- match.arg(model)
  system.file("stan", paste0("ve_", model, ".stan"), package = "surviver",
              mustWork = TRUE)
}

#' Fit the time-varying VE model
#'
#' Fits a Poisson approximation to a Cox model in which the log hazard ratio
#' for vaccination, \eqn{\gamma(t)}, follows a random walk over days since
#' vaccination:
#' \deqn{\log \mu = \log\beta + \log \bar I + \log(\text{person-days}) +
#'   \gamma(t)\,\text{vax}.}
#' VE is \eqn{1 - \exp(\gamma(t))}. Requires
#' [cmdstanr](https://mc-stan.org/cmdstanr/) and CmdStan. The compiled model
#' is cached in `tools::R_user_dir("surviver", "cache")`.
#'
#' @param stan_data Output of [make_stan_data()], or its `data` element.
#' @param model Stan model; see [ve_stan_file()].
#' @param chains,iter_warmup,iter_sampling,adapt_delta,max_treedepth,seed
#'   Passed to `$sample()`.
#' @param ... Further arguments to `$sample()`.
#'
#' @return A `CmdStanMCMC` fit.
#' @export
fit_ve <- function(stan_data,
                   model         = c("noncentered", "centered"),
                   chains        = 4,
                   iter_warmup   = 1000,
                   iter_sampling = 1000,
                   adapt_delta   = 0.95,
                   max_treedepth = 10,
                   seed          = NULL,
                   ...) {
  rlang::check_installed("cmdstanr", reason = "to fit the VE model.")
  if (!is.null(stan_data$data)) stan_data <- stan_data$data
  model <- match.arg(model)

  exe_dir <- tools::R_user_dir("surviver", "cache")
  dir.create(exe_dir, recursive = TRUE, showWarnings = FALSE)
  mod <- cmdstanr::cmdstan_model(ve_stan_file(model), dir = exe_dir)

  mod$sample(
    data            = stan_data,
    chains          = chains,
    parallel_chains = chains,
    iter_warmup     = iter_warmup,
    iter_sampling   = iter_sampling,
    adapt_delta     = adapt_delta,
    max_treedepth   = max_treedepth,
    seed            = seed,
    ...
  )
}

#' Summarise VE by day since vaccination
#'
#' @param fit A fit from [fit_ve()].
#' @param probs Lower and upper quantiles for the credible interval.
#'
#' @return A tibble with `days_since_vax`, `median`, `lower`, `upper` (VE
#'   scale), plus `rhat` and `ess_bulk`.
#' @export
extract_ve <- function(fit, probs = c(0.05, 0.95)) {
  s <- fit$summary(
    "VE",
    median = stats::median,
    ~ stats::quantile(.x, probs = probs, names = FALSE),
    rhat     = posterior::rhat,
    ess_bulk = posterior::ess_bulk
  )
  names(s)[3:4] <- c("lower", "upper")
  bin <- as.integer(sub("VE\\[(\\d+)\\]", "\\1", s$variable))
  dplyr::tibble(
    days_since_vax = bin - 1L,
    median         = s$median,
    lower          = s$lower,
    upper          = s$upper,
    rhat           = s$rhat,
    ess_bulk       = s$ess_bulk
  ) |>
    dplyr::filter(.data$days_since_vax > 0)
}

#' True VE curve used in the HERMES simulations
#'
#' Linear rise from 0 to `vaccine_effectiveness` over `ve_ramp_days`, then
#' exponential decay with half-life `ve_half_life` towards `ve_floor`. Matches
#' `compute_time_varying_ve()` in the
#' [HERMES waning-VE extension](https://github.com/kylieainslie/ibm_hermes/tree/waning-vaccine-effectiveness).
#'
#' @param t Days since vaccination.
#' @param params A list with `vaccine_effectiveness`, `ve_ramp_days`,
#'   `ve_half_life`, `ve_floor`, e.g. `ve_sim_waning$params`.
#' @return Numeric vector of VE values.
#' @export
#' @examples
#' true_ve(c(7, 14, 59), ve_sim_waning$params)
true_ve <- function(t, params) {
  ramp <- t <= params$ve_ramp_days
  ve <- numeric(length(t))
  ve[ramp] <- params$vaccine_effectiveness * t[ramp] / params$ve_ramp_days
  ve[!ramp] <- params$ve_floor +
    (params$vaccine_effectiveness - params$ve_floor) *
    exp(-log(2) / params$ve_half_life * (t[!ramp] - params$ve_ramp_days))
  ve
}

#' Simulate event counts from the VE model
#'
#' Replaces the outcome in Stan input with Poisson counts generated under
#' known parameters, keeping the real covariate structure. Used to check
#' that [fit_ve()] recovers known values.
#'
#' @param stan_data Output of [make_stan_data()], or its `data` element.
#' @param log_beta True log transmission scaling.
#' @param gamma True log hazard ratio per bin (length `T_vax`; element 1 is
#'   the unvaccinated bin and is ignored).
#' @param seed Optional random seed, applied locally.
#' @return `stan_data` (same form as supplied) with `Y` replaced.
#' @export
simulate_counts <- function(stan_data, log_beta, gamma, seed = NULL) {
  wrapped <- !is.null(stan_data$data)
  d <- if (wrapped) stan_data$data else stan_data
  if (length(gamma) != d$T_vax) {
    stop("`gamma` must have length T_vax (", d$T_vax, ").", call. = FALSE)
  }
  mu <- exp(log_beta + d$log_Ibar + d$log_persdays + gamma[d$tvax_bin] * d$vax)
  draw <- function() as.integer(stats::rpois(d$N, mu))
  d$Y <- if (is.null(seed)) draw() else withr::with_seed(seed, draw())
  if (wrapped) {
    stan_data$data <- d
    stan_data$cells$Y <- d$Y
    stan_data
  } else {
    d
  }
}

#' Plot estimated VE by day since vaccination
#'
#' @param ve Output of [extract_ve()]. Several fits can be compared by
#'   binding their outputs with a `model` (or other) column named in `colour`.
#' @param truth Optional HERMES parameter list (e.g. `ve_sim_waning$params`);
#'   the true VE(t) is drawn as a dashed line.
#' @param colour Optional column of `ve` to colour by.
#' @param max_days Upper limit of the x axis.
#' @return A ggplot object.
#' @export
plot_ve <- function(ve, truth = NULL, colour = NULL, max_days = NULL) {
  rlang::check_installed("ggplot2", reason = "to plot VE.")
  mapping <- if (is.null(colour)) {
    ggplot2::aes(.data$days_since_vax, .data$median)
  } else {
    ggplot2::aes(.data$days_since_vax, .data$median,
                 colour = .data[[colour]], fill = .data[[colour]])
  }
  ribbon_args <- list(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                      alpha = 0.2, colour = NA)
  line_args <- list(linewidth = 0.9)
  if (is.null(colour)) {
    ribbon_args$fill <- "steelblue"
    line_args$colour <- "steelblue"
  }
  p <- ggplot2::ggplot(ve, mapping) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
    do.call(ggplot2::geom_ribbon, ribbon_args) +
    do.call(ggplot2::geom_line, line_args)
  if (!is.null(truth)) {
    days <- seq(0, max(ve$days_since_vax))
    p <- p + ggplot2::geom_line(
      data = data.frame(days_since_vax = days, median = true_ve(days, truth)),
      ggplot2::aes(.data$days_since_vax, .data$median),
      inherit.aes = FALSE, linetype = "dashed", colour = "firebrick",
      linewidth = 0.9
    )
  }
  p <- p +
    ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
    ggplot2::labs(x = "Days since vaccination", y = "Vaccine effectiveness",
                  colour = NULL, fill = NULL) +
    ggplot2::theme_minimal()
  if (!is.null(max_days)) p <- p + ggplot2::coord_cartesian(xlim = c(0, max_days))
  p
}

#' Observed and model-expected case counts
#'
#' Sums observed cases and the posterior distribution of expected cases
#' (\eqn{\mu} from the fitted model) over groups of Stan cells, e.g. by
#' calendar day and vaccination status. Comparing the two shows how well
#' the model fits the case data.
#'
#' @param fit A fit from [fit_ve()].
#' @param stan_data The output of [make_stan_data()] used for the fit (with
#'   its `cells`).
#' @param by Columns of `stan_data$cells` to group by.
#' @param probs Lower and upper quantiles for the interval.
#'
#' @return A tibble with the `by` columns, `observed`, and the posterior
#'   `mean`, `lower` and `upper` of the expected count.
#' @export
expected_cases <- function(fit, stan_data, by = c("day", "vax_status"),
                           probs = c(0.05, 0.95)) {
  rlang::check_installed("posterior", reason = "to summarise draws.")
  cells <- stan_data$cells
  if (is.null(cells)) {
    stop("`stan_data` must be the full output of make_stan_data().", call. = FALSE)
  }
  log_mu <- posterior::as_draws_matrix(fit$draws("log_mu"))
  if (ncol(log_mu) != nrow(cells)) {
    stop("`fit` and `stan_data` do not match (different numbers of cells).",
         call. = FALSE)
  }
  group <- interaction(cells[by], drop = TRUE, lex.order = TRUE)
  # indicator matrix: cells x groups, so draws %*% ind sums mu within groups
  ind <- indicator_matrix(group)
  mu_by_group <- exp(unclass(log_mu)) %*% ind

  keys <- cells[!duplicated(group), by, drop = FALSE]
  keys <- keys[order(group[!duplicated(group)]), , drop = FALSE]
  out <- dplyr::as_tibble(keys)
  out$observed <- as.vector(tapply(cells$Y, group, sum))
  out$mean  <- colMeans(mu_by_group)
  out$lower <- apply(mu_by_group, 2, stats::quantile, probs = probs[1], names = FALSE)
  out$upper <- apply(mu_by_group, 2, stats::quantile, probs = probs[2], names = FALSE)
  out
}

indicator_matrix <- function(group) {
  m <- matrix(0, nrow = length(group), ncol = nlevels(group))
  m[cbind(seq_along(group), as.integer(group))] <- 1
  m
}
