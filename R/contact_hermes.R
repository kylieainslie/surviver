#' Contact matrix implied by a HERMES population
#'
#' Computes the age-by-age contact matrix that HERMES effectively uses for
#' transmission, from the household, school-class and workplace structure of
#' a simulated population.
#'
#' Each day, HERMES gives every infectious person \eqn{i} and susceptible
#' person \eqn{j} a contact weight
#' \deqn{w_{ij} = c_\text{community} + c_\text{household} 1[\text{same household}] +
#'   c_\text{school} 1[\text{same class}] + c_\text{workplace} 1[\text{same workplace}],}
#' with \eqn{c_\text{community} = 1 - \exp(-k/N)} for `k` community contacts
#' per day, and infection probability `transmission_prob` \eqn{\times w_{ij}}.
#' This function returns
#' \eqn{C_{a'a} = \frac{1}{N_{a'}} \sum_{i \in a'} \sum_{j \in a, j \ne i} w_{ij}},
#' the mean daily weighted contacts a person in group \eqn{a'} has with group
#' \eqn{a}, in the same orientation as [polymod_nl]. It is reciprocal:
#' \eqn{N_{a'} C_{a'a} = N_a C_{aa'}}.
#'
#' Averaging over each age group ignores clustering: a household member of a
#' case is at much higher risk than the age-group average suggests.
#'
#' @param demographic Data frame with `age`, `hh_id`, `classroom_id` and
#'   `workplace_id` (`NA` when not in school or work), e.g.
#'   `ve_sim_waning$demographic`.
#' @param params HERMES parameter list with `contact_prob_household`,
#'   `contact_prob_school`, `contact_prob_workplace` and
#'   `num_contacts_community_day`, e.g. `ve_sim_waning$params`.
#' @inheritParams add_age_group
#'
#' @return A square numeric matrix with dimnames `age_group` (rows) and
#'   `contact_age_group` (columns).
#' @export
#' @examples
#' cm <- contact_matrix_hermes(ve_sim_waning$demographic, ve_sim_waning$params)
#' round(cm, 2)
contact_matrix_hermes <- function(demographic, params,
                                  age_breaks = default_age_breaks,
                                  age_labels = default_age_labels) {
  d <- add_age_group(demographic, age_breaks, age_labels)
  k <- length(age_labels)
  n_age <- tabulate(as.integer(d$age_group), nbins = k)

  # within-setting pairs: sum over groups of (n_a' x n_a), minus self-pairs
  pairs_within <- function(group_id) {
    keep <- !is.na(group_id)
    if (!any(keep)) return(matrix(0, k, k))
    x <- table(factor(group_id[keep]), d$age_group[keep])
    x <- matrix(as.numeric(x), nrow = nrow(x))
    m <- crossprod(x)
    m - diag(tabulate(as.integer(d$age_group[keep]), nbins = k))
  }

  n <- nrow(d)
  c_comm <- 1 - exp(-params$num_contacts_community_day / n)
  community <- c_comm * (outer(n_age, n_age) - diag(n_age))

  total <- community +
    params$contact_prob_household * pairs_within(d$hh_id) +
    params$contact_prob_school    * pairs_within(d$classroom_id) +
    params$contact_prob_workplace * pairs_within(d$workplace_id)

  cm <- total / n_age   # divide each row a' by N_a'
  dimnames(cm) <- list(age_group = age_labels, contact_age_group = age_labels)
  cm
}
