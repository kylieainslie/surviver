#' Default age groups
#'
#' Age-group breaks and labels used throughout surviver unless overridden.
#' Groups are left-closed: `[0, 5)`, `[5, 18)`, `[18, 50)`, `[50, 65)`,
#' `[65, Inf)`.
#'
#' @export
default_age_breaks <- c(0, 5, 18, 50, 65, Inf)

#' @rdname default_age_breaks
#' @export
default_age_labels <- c("0-4", "5-17", "18-49", "50-64", "65+")

#' Build a person-level table from VE input data
#'
#' Joins demographic, vaccination and reported-case data into one row per
#' person.
#'
#' @param x Either a directory containing `demographic_data.rds`,
#'   `vaccination_data.rds` and `incidence_data_reduced.rds` (the HERMES
#'   export format, optionally with `scenario_params.rds`), or a list with
#'   elements `demographic`, `vaccination` and `incidence` such as
#'   [ve_sim_waning].
#' @param study_end Last day of follow-up. Defaults to `x$study_end` (or the
#'   value in `scenario_params.rds`); if neither exists, the last whole
#'   reporting day, with a message.
#'
#' @return A data frame with one row per person: the demographic columns,
#'   `time_of_vaccination` and `time_of_reporting` (`NA` if none), and logical
#'   `vaccinated` and `infected`. The study end is stored in the
#'   `"study_end"` attribute.
#' @export
#' @examples
#' pop <- ve_population(ve_sim_waning)
#' attr(pop, "study_end")
ve_population <- function(x, study_end = NULL) {
  if (is.character(x) && length(x) == 1) {
    read <- function(f) {
      d <- readRDS(file.path(x, f))
      d[setdiff(names(d), "time")]
    }
    demo <- read("demographic_data.rds")
    vax  <- read("vaccination_data.rds")
    inc  <- read("incidence_data_reduced.rds")
    params_file <- file.path(x, "scenario_params.rds")
    if (is.null(study_end) && file.exists(params_file)) {
      study_end <- readRDS(params_file)$study_end
    }
  } else if (is.list(x) && all(c("demographic", "vaccination", "incidence") %in% names(x))) {
    demo <- x$demographic
    vax  <- x$vaccination
    inc  <- x$incidence
    if (is.null(study_end)) study_end <- x$study_end
  } else {
    stop("`x` must be a data directory or a list with `demographic`, ",
         "`vaccination` and `incidence` elements.", call. = FALSE)
  }

  vax <- vax[c("ind_id", "time_of_vaccination")]
  inc <- inc[!is.na(inc$time_of_reporting), c("ind_id", "time_of_reporting")]
  for (d in list(demo = demo, vax = vax, inc = inc)) {
    if (anyDuplicated(d$ind_id)) {
      stop("Duplicated `ind_id` in input data.", call. = FALSE)
    }
  }

  pop <- demo |>
    dplyr::left_join(vax, by = "ind_id") |>
    dplyr::left_join(inc, by = "ind_id") |>
    dplyr::mutate(
      vaccinated = !is.na(.data$time_of_vaccination),
      infected   = !is.na(.data$time_of_reporting)
    )

  if (is.null(study_end)) {
    study_end <- floor(max(inc$time_of_reporting))
    message("`study_end` not supplied; using last reporting day (", study_end, ").")
  }
  attr(pop, "study_end") <- study_end
  pop
}

#' Add an age-group column
#'
#' @param data Data frame with an `age` column.
#' @param age_breaks,age_labels Passed to [cut()] (left-closed intervals).
#' @return `data` with a factor column `age_group`.
#' @export
add_age_group <- function(data,
                          age_breaks = default_age_breaks,
                          age_labels = default_age_labels) {
  data$age_group <- cut(data$age, breaks = age_breaks, labels = age_labels,
                        right = FALSE)
  data
}

#' Randomly sample a cohort
#'
#' @param pop Person-level data from [ve_population()].
#' @param n Cohort size. If `n >= nrow(pop)` the whole population is returned.
#' @param seed Optional random seed, applied locally.
#' @return A subset of `pop` (attributes preserved).
#' @export
sample_cohort <- function(pop, n, seed = NULL) {
  if (n >= nrow(pop)) return(pop)
  draw <- function() sample.int(nrow(pop), n)
  idx <- if (is.null(seed)) draw() else withr::with_seed(seed, draw())
  out <- pop[idx, , drop = FALSE]
  attr(out, "study_end") <- attr(pop, "study_end")
  out
}

#' Daily cases and population at risk
#'
#' Counts reported cases and the population at risk on each study day by age
#' group and vaccination status.
#'
#' A person is at risk on day `t` if they have not been reported before day
#' `t`. They count as vaccinated on day `t` if vaccinated before day `t`, and
#' a case counts as vaccinated if vaccination preceded reporting.
#'
#' @param pop Person-level data from [ve_population()].
#' @param study_end Last day of follow-up.
#' @inheritParams add_age_group
#'
#' @return A tibble with one row per `day` x `age_group` x `vax_status`:
#'   `n_at_risk`, `n_cases` and `incidence_rate` (`n_cases / n_at_risk`).
#' @export
daily_incidence <- function(pop,
                            study_end  = attr(pop, "study_end"),
                            age_breaks = default_age_breaks,
                            age_labels = default_age_labels) {
  check_study_end(study_end)
  pop <- add_age_group(pop, age_breaks, age_labels)

  at_risk <- lapply(seq_len(study_end), function(t) {
    risk <- is.na(pop$time_of_reporting) | pop$time_of_reporting >= t
    vax  <- !is.na(pop$time_of_vaccination) & pop$time_of_vaccination < t
    tab  <- table(
      age_group  = pop$age_group[risk],
      vax_status = factor(as.integer(vax[risk]), levels = 0:1)
    )
    out <- as.data.frame(tab, responseName = "n_at_risk",
                         stringsAsFactors = FALSE)
    out$day <- t
    out
  })
  at_risk <- do.call(rbind, at_risk)
  at_risk$vax_status <- as.integer(at_risk$vax_status)

  cases <- pop |>
    dplyr::filter(.data$infected, .data$time_of_reporting <= study_end) |>
    dplyr::mutate(
      day        = as.integer(.data$time_of_reporting),
      vax_status = as.integer(.data$vaccinated &
                                .data$time_of_vaccination < .data$time_of_reporting),
      age_group  = as.character(.data$age_group)
    ) |>
    dplyr::count(.data$day, .data$age_group, .data$vax_status, name = "n_cases")

  at_risk |>
    dplyr::left_join(cases, by = c("day", "age_group", "vax_status")) |>
    dplyr::mutate(
      n_cases        = dplyr::coalesce(.data$n_cases, 0L),
      age_group      = factor(.data$age_group, levels = age_labels),
      incidence_rate = dplyr::if_else(.data$n_at_risk > 0,
                                      .data$n_cases / .data$n_at_risk,
                                      NA_real_)
    ) |>
    dplyr::select("day", "age_group", "vax_status", "n_at_risk", "n_cases",
                  "incidence_rate") |>
    dplyr::arrange(.data$day, .data$age_group, .data$vax_status) |>
    dplyr::as_tibble()
}

#' Daily number infectious
#'
#' Approximates the number of people infectious on each day from reported
#' cases, assuming each case was infectious for `duration_infection` days
#' starting `reporting_delay` days before they were reported.
#'
#' @param incidence Output of [daily_incidence()].
#' @param duration_infection Days infectious.
#' @param reporting_delay Days from start of infectiousness to reporting.
#'
#' @return A tibble with one row per `day` x `age_group` x `vax_status`
#'   (vaccination status at reporting) and column `n_infectious`, from day 0
#'   (or earlier, if a case was infectious before day 0) to the last study
#'   day.
#' @export
daily_prevalence <- function(incidence,
                             duration_infection = 7,
                             reporting_delay    = 2) {
  cases <- incidence[incidence$n_cases > 0, ]
  offsets <- seq_len(duration_infection) - 1L - reporting_delay
  study_end <- max(incidence$day)

  spread <- cases[rep(seq_len(nrow(cases)), each = length(offsets)), ]
  spread$day <- spread$day + rep(offsets, times = nrow(cases))

  grid <- tidyr::expand_grid(
    day        = seq(min(c(0L, spread$day)), study_end),
    age_group  = factor(levels(incidence$age_group),
                        levels = levels(incidence$age_group)),
    vax_status = 0:1
  )

  spread |>
    dplyr::filter(.data$day <= study_end) |>
    dplyr::group_by(.data$day, .data$age_group, .data$vax_status) |>
    dplyr::summarise(n_infectious = sum(.data$n_cases), .groups = "drop") |>
    dplyr::right_join(grid, by = c("day", "age_group", "vax_status")) |>
    dplyr::mutate(n_infectious = dplyr::coalesce(.data$n_infectious, 0L)) |>
    dplyr::arrange(.data$day, .data$age_group, .data$vax_status)
}

#' Counting-process data with time-varying vaccination
#'
#' Splits each person's follow-up at vaccination, giving `(tstart, tstop]`
#' intervals with `vax_status` 0 before and 1 after vaccination. Follow-up
#' ends at reporting (`event = 1`) or `study_end` (censored); cases reported
#' after `study_end` are censored at `study_end`.
#'
#' People reported before their vaccination date are already handled: their
#' follow-up ends at the report, so they never count as vaccinated. People
#' infected before vaccination but reported after it would otherwise count
#' as vaccinated cases. `censor_at_vaccination` lets these be censored at
#' their vaccination date instead (follow-up ends unvaccinated, no event).
#' Which people to censor is the analyst's choice; see
#' [reported_soon_after_vaccination()] for a rule based only on observed
#' dates.
#'
#' @param cohort Person-level data from [ve_population()] or
#'   [sample_cohort()].
#' @param study_end Last day of follow-up.
#' @param censor_at_vaccination Optional vector of `ind_id`s to censor at
#'   their vaccination date. Ignored for people who are unvaccinated or whose
#'   follow-up ends before vaccination.
#' @inheritParams add_age_group
#'
#' @return A data frame from [survival::tmerge()] with `ind_id`, `age`,
#'   `age_group`, `sex`, `time_of_vaccination`, `tstart`, `tstop`, `event`,
#'   `vax_status`.
#' @export
counting_process <- function(cohort,
                             study_end  = attr(cohort, "study_end"),
                             censor_at_vaccination = NULL,
                             age_breaks = default_age_breaks,
                             age_labels = default_age_labels) {
  check_study_end(study_end)

  base <- cohort |>
    add_age_group(age_breaks, age_labels) |>
    dplyr::mutate(
      event = as.integer(.data$infected & .data$time_of_reporting <= study_end),
      tstop = dplyr::if_else(.data$event == 1L,
                             .data$time_of_reporting, as.numeric(study_end))
    )

  if (!is.null(censor_at_vaccination)) {
    cens <- base$ind_id %in% censor_at_vaccination & base$vaccinated &
      base$time_of_vaccination < base$tstop
    base$tstop[cens] <- base$time_of_vaccination[cens]
    base$event[cens] <- 0L
    # follow-up now ends at vaccination, so no vaccinated interval is created
    base$vaccinated[cens] <- FALSE
  }

  if (any(base$tstop <= 0)) {
    warning(sum(base$tstop <= 0), " people reported on or before day 0 ",
            "were dropped.", call. = FALSE)
    base <- base[base$tstop > 0, ]
  }

  keep <- base[c("ind_id", "age", "age_group", "sex", "time_of_vaccination",
                 "tstop", "event")]
  cp <- survival::tmerge(keep, keep, id = ind_id,
                         event = event(tstop, event))

  vaccinated <- base[base$vaccinated, c("ind_id", "time_of_vaccination")]
  cp <- survival::tmerge(cp, vaccinated, id = ind_id,
                         vax_status = tdc(time_of_vaccination))
  cp$vax_status <- as.integer(dplyr::coalesce(cp$vax_status, 0))
  cp
}

#' Expand counting-process data to person-days
#'
#' Turns each `(tstart, tstop]` interval into one row per day, with
#' `event = 1` only on the day of reporting. Optionally attaches the daily
#' incidence rate for the person's age group and vaccination status.
#'
#' @param cp Output of [counting_process()].
#' @param incidence Optional output of [daily_incidence()].
#'
#' @return A tibble with one row per person-day: `ind_id`, `age`,
#'   `age_group`, `sex`, `day`, `tstart`, `tstop` (= `day`), `event`,
#'   `vax_status`, `time_of_vaccination`, `days_since_vax` (`NA` when
#'   unvaccinated) and, if `incidence` is given, `incidence_rate`.
#' @export
expand_person_days <- function(cp, incidence = NULL) {
  if (any(cp$tstart != round(cp$tstart)) || any(cp$tstop != round(cp$tstop))) {
    stop("Interval endpoints must be whole days.", call. = FALSE)
  }
  n_days <- as.integer(cp$tstop - cp$tstart)
  idx    <- rep(seq_len(nrow(cp)), n_days)

  pd <- cp[idx, c("ind_id", "age", "age_group", "sex", "vax_status",
                  "time_of_vaccination")]
  pd$day    <- as.integer(cp$tstart[idx]) + sequence(n_days)
  pd$tstart <- pd$day - 1L
  pd$tstop  <- pd$day
  pd$event  <- as.integer(cp$event[idx] == 1L & pd$day == cp$tstop[idx])
  pd$days_since_vax <- ifelse(pd$vax_status == 1L,
                              pd$day - pd$time_of_vaccination, NA_real_)

  pd <- dplyr::as_tibble(pd)
  if (!is.null(incidence)) {
    pd <- dplyr::left_join(
      pd,
      incidence[c("day", "age_group", "vax_status", "incidence_rate")],
      by = c("day", "age_group", "vax_status")
    )
  }
  dplyr::relocate(pd, "ind_id", "age", "age_group", "sex", "day", "tstart",
                  "tstop", "event", "vax_status")
}

check_study_end <- function(study_end) {
  if (is.null(study_end) || length(study_end) != 1 || !is.finite(study_end) ||
      study_end < 1) {
    stop("`study_end` must be a single positive number. Supply it, or build ",
         "the data with ve_population() so it is stored as an attribute.",
         call. = FALSE)
  }
  invisible(study_end)
}

#' People reported soon after vaccination
#'
#' A rule, based only on observed dates, for identifying vaccinated people
#' who may have been infected before vaccination: those reported within
#' `window` days after their vaccination date. Infection precedes reporting
#' by the incubation period plus the reporting delay, so a case reported
#' shortly after vaccination may have been infected before it. Pass the
#' result to `counting_process(censor_at_vaccination = )`.
#'
#' Censoring also removes some genuine early vaccine failures, so VE is not
#' estimable for the first `window` days after vaccination. Longer windows
#' catch more prior infections but discard more genuine cases.
#'
#' @param pop Person-level data from [ve_population()].
#' @param window Days after vaccination within which a report triggers
#'   censoring (reports on days `1` to `window` after vaccination).
#' @return Character vector of `ind_id`s.
#' @export
#' @examples
#' pop <- ve_population(ve_sim_waning)
#' ids <- reported_soon_after_vaccination(pop, window = 7)
#' length(ids)
#' cp <- counting_process(pop, censor_at_vaccination = ids)
reported_soon_after_vaccination <- function(pop, window) {
  gap <- pop$time_of_reporting - pop$time_of_vaccination
  as.character(pop$ind_id[pop$vaccinated & pop$infected & !is.na(gap) &
                            gap > 0 & gap <= window])
}
