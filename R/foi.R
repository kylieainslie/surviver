#' Population size by age group
#'
#' @param pop Person-level data from [ve_population()].
#' @inheritParams add_age_group
#' @return A named integer vector of group sizes.
#' @export
age_group_sizes <- function(pop,
                            age_breaks = default_age_breaks,
                            age_labels = default_age_labels) {
  pop <- add_age_group(pop, age_breaks, age_labels)
  sizes <- table(pop$age_group)
  stats::setNames(as.integer(sizes), names(sizes))
}

#' Daily force-of-infection offset
#'
#' Computes the infectious pressure \eqn{\bar I(t, a)} on each study day for
#' each age group, used as an offset in the VE model.
#'
#' With age structure, for recipient age group \eqn{a} on day \eqn{t}:
#' \deqn{\bar I(t, a) = \sum_{a'} \frac{C_{a'a}}{N_a} I^{\mathrm{eff}}_{a'}(t - \ell)}
#' where \eqn{I^{\mathrm{eff}}_{a'} = I^u_{a'} + (1 - \theta) I^v_{a'}} is the
#' number infectious in group \eqn{a'} (unvaccinated plus vaccinated, the
#' latter down-weighted by \eqn{\theta}) and \eqn{\ell} is `lag`.
#'
#' Without age structure every group gets
#' \eqn{\bar I(t) = \sum_{a'} I^{\mathrm{eff}}_{a'}(t - \ell) / N}.
#'
#' @param prevalence Output of [daily_prevalence()].
#' @param pop_size Named vector of population size by age group, e.g. from
#'   [age_group_sizes()].
#' @param contact_matrix Square matrix of contacts with rows and columns in
#'   the order of `names(pop_size)`; `contact_matrix[a', a]` is the mean daily
#'   contacts a person in `a'` has with group `a`. Ignored when
#'   `age_structured = FALSE`.
#' @param theta Relative reduction in infectiousness of vaccinated cases
#'   (0 = as infectious as unvaccinated, 1 = not infectious). The HERMES
#'   simulations assume vaccination does not affect infectiousness, i.e.
#'   `theta = 0`.
#' @param age_structured Use the contact matrix (`TRUE`) or overall
#'   prevalence (`FALSE`).
#' @param lag Days between prevalence and the infections it causes.
#' @param study_end Last day to compute; defaults to the last prevalence day.
#' @param floor Lower bound applied so the log offset is finite.
#'
#' @return A tibble with `day` (1 to `study_end`), `age_group` and `Ibar`.
#' @export
force_of_infection <- function(prevalence,
                               pop_size,
                               contact_matrix = surviver::polymod_nl,
                               theta          = 0,
                               age_structured = TRUE,
                               lag            = 1L,
                               study_end      = max(prevalence$day),
                               floor          = 1e-10) {
  groups <- names(pop_size)
  if (is.null(groups) || !setequal(groups, levels(prevalence$age_group))) {
    stop("`pop_size` must be named by the age groups in `prevalence`.",
         call. = FALSE)
  }
  days <- seq_len(study_end)

  # effective number infectious by (day, age group): rows = days 0..study_end
  prev_days <- seq(0, study_end)
  i_eff <- matrix(0, nrow = length(prev_days), ncol = length(groups),
                  dimnames = list(NULL, groups))
  keep <- prevalence$day %in% prev_days
  w    <- ifelse(prevalence$vax_status[keep] == 1L, 1 - theta, 1)
  agg  <- stats::aggregate(
    list(n = prevalence$n_infectious[keep] * w),
    by = list(day = prevalence$day[keep],
              age_group = as.character(prevalence$age_group[keep])),
    FUN = sum
  )
  i_eff[cbind(match(agg$day, prev_days), match(agg$age_group, groups))] <- agg$n

  # lagged prevalence for each study day (zero before the data start)
  lag_rows <- match(days - lag, prev_days)
  i_lag <- i_eff[lag_rows, , drop = FALSE]
  i_lag[is.na(i_lag)] <- 0

  if (age_structured) {
    if (!all(groups %in% rownames(contact_matrix)) ||
        !all(groups %in% colnames(contact_matrix))) {
      stop("`contact_matrix` dimnames must match the age groups.", call. = FALSE)
    }
    cm <- contact_matrix[groups, groups]
    ibar <- sweep(i_lag %*% cm, 2, pop_size[groups], "/")
  } else {
    total <- rowSums(i_lag) / sum(pop_size)
    ibar  <- matrix(total, nrow = length(days), ncol = length(groups))
  }

  dplyr::tibble(
    day       = rep(days, times = length(groups)),
    age_group = factor(rep(groups, each = length(days)), levels = groups),
    Ibar      = pmax(as.vector(ibar), floor)
  ) |>
    dplyr::arrange(.data$day, .data$age_group)
}
