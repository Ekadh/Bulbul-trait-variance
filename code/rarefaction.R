#libs
library(tidyverse)

bulbul <- read.csv("../data/Bulbul_AVONET_data_integrated_trimmed.csv")

bulbul_trimmed <- bulbul %>%
  group_by(AviList) %>%
  filter(n() > 50) %>%
  ungroup() %>%
  select(AviList, Wing.Length, Beak.Length_Nares, Beak.Depth, Tarsus.Length)

## Rarefaction analysis to determine how many specimens are needed to get a reliable estimate of CV for each trait

# The CV calculated from all available (non-missing) measurements for a species
# is treated as its reference, or "true", CV.  At every sample size, 1,000
# random nested samples are drawn.  A species is considered to have converged
# at the first sample size for which the *median* rarefied CV is within +/- 5%
# of its reference CV for that sample size and the next two sample sizes.  The
# three consecutive sizes prevent a single chance agreement from determining
# the required sample size.

set.seed(20260804)

traits <- c(
  "Wing.Length",
  "Beak.Length_Nares",
  "Beak.Depth",
  "Tarsus.Length"
)
n_reps <- 1000L
n_consecutive <- 3L
relative_tolerance <- 0.05

# Calculate CVs at every sample size from random permutations.  Taking the
# cumulative sums makes each replicate a nested rarefaction curve and avoids
# repeatedly drawing a fresh sample at every n.
rarefy_one_species <- function(values, n_reps) {
  n_total <- length(values)
  sample_sizes <- 2:n_total

  cv_matrix <- vapply(seq_len(n_reps), function(i) {
    sampled_values <- sample(values, size = n_total, replace = FALSE)
    sums <- cumsum(sampled_values)
    sum_squares <- cumsum(sampled_values^2)

    # Sample SD / sample mean for each prefix of the random permutation.
    prefix_n <- sample_sizes
    prefix_variance <-
      (sum_squares[prefix_n] - sums[prefix_n]^2 / prefix_n) /
      (prefix_n - 1)
    # Floating-point rounding can make a true zero variance very slightly
    # negative; variance itself cannot be negative.
    prefix_sd <- sqrt(pmax(prefix_variance, 0))
    prefix_sd / (sums[prefix_n] / prefix_n)
  }, numeric(length(sample_sizes)))

  if (is.null(dim(cv_matrix))) {
    cv_matrix <- matrix(cv_matrix, ncol = 1)
  }

  tibble(
    n_specimens = sample_sizes,
    median_rarefied_cv = apply(cv_matrix, 1, median),
    mean_rarefied_cv = rowMeans(cv_matrix),
    lower_95_ci = apply(cv_matrix, 1, quantile, probs = 0.025),
    upper_95_ci = apply(cv_matrix, 1, quantile, probs = 0.975)
  )
}

find_convergence_n <- function(rarefaction_curve, true_cv,
                               relative_tolerance, n_consecutive) {
  within_tolerance <- abs(
    rarefaction_curve$median_rarefied_cv / true_cv - 1
  ) <= relative_tolerance

  # A run must stay within the tolerance for n, n + 1, and n + 2.
  run_starts <- which(vapply(
    seq_len(length(within_tolerance) - n_consecutive + 1),
    function(i) all(within_tolerance[i:(i + n_consecutive - 1)]),
    logical(1)
  ))

  if (length(run_starts) == 0) NA_integer_ else
    rarefaction_curve$n_specimens[run_starts[1]]
}

rarefaction_curves <- list()
species_requirements <- list()
counter <- 1L

for (trait in traits) {
  for (species in unique(bulbul_trimmed$AviList)) {
    values <- bulbul_trimmed %>%
      filter(AviList == species, !is.na(.data[[trait]])) %>%
      pull(all_of(trait))

    # All retained species have at least 12 usable observations per trait,
    # but this guard makes the code safe if the input data are updated.
    if (length(values) < n_consecutive + 1 || mean(values) == 0) next

    true_cv <- sd(values) / mean(values)
    curve <- rarefy_one_species(values, n_reps) %>%
      mutate(
        species = species,
        trait = trait,
        true_cv = true_cv,
        relative_error = median_rarefied_cv / true_cv - 1
      )

    required_n <- find_convergence_n(
      curve,
      true_cv = true_cv,
      relative_tolerance = relative_tolerance,
      n_consecutive = n_consecutive
    )

    rarefaction_curves[[counter]] <- curve
    species_requirements[[counter]] <- tibble(
      species = species,
      trait = trait,
      n_available = length(values),
      true_cv = true_cv,
      specimens_required = required_n
    )
    counter <- counter + 1L
  }
}

rarefaction_curves <- bind_rows(rarefaction_curves)
species_requirements <- bind_rows(species_requirements)

# These are the four requested values: mean requirement across the 52 species
# for each trait.  n_species confirms whether any species could not converge.
trait_specimen_requirements <- species_requirements %>%
  group_by(trait) %>%
  summarise(
    mean_specimens_required = mean(specimens_required, na.rm = TRUE),
    median_specimens_required = median(specimens_required, na.rm = TRUE),
    n_species = n(),
    n_converged = sum(!is.na(specimens_required)),
    .groups = "drop"
  )

print(trait_specimen_requirements)

# Save the final four values and the species-level/curve-level results needed
# for reporting or plotting the sensitivity analysis.
write_csv(trait_specimen_requirements,
          "../results/Rarefaction test/rarefaction_trait_specimen_requirements.csv")
write_csv(species_requirements,
          "../results/Rarefaction test/rarefaction_species_specimen_requirements.csv")
write_csv(rarefaction_curves,
          "../results/Rarefaction test/rarefaction_curves.csv")

##Now plot specimens_required against total sample size for each species

species_requirements <- read_csv("../results/Rarefaction test/rarefaction_species_specimen_requirements.csv")

#Calculate true CV for each species and trait
true_cv_df <- bulbul_trimmed %>%
  group_by(AviList) %>%
  summarise(
    true_cv_wing = sd(Wing.Length, na.rm = TRUE) / mean(Wing.Length, na.rm = TRUE),
    true_cv_beak_length = sd(Beak.Length_Nares, na.rm = TRUE) / mean(Beak.Length_Nares, na.rm = TRUE),
    true_cv_beak_depth = sd(Beak.Depth, na.rm = TRUE) / mean(Beak.Depth, na.rm = TRUE),
    true_cv_tarsus_length = sd(Tarsus.Length, na.rm = TRUE) / mean(Tarsus.Length, na.rm = TRUE)
  )

## Make a plot of specimens_required vs CV for each species and trait

# Merge species_requirements with true_cv_df
species_requirements_merged <- species_requirements %>%
  left_join(
    true_cv_df %>% rename(AviList = AviList),
    by = c("species" = "AviList")
  )

# Create a plot of specimens_required vs CV for each trait
plot_data <- species_requirements_merged %>%
  pivot_longer(
    cols = starts_with("true_cv_"),
    names_to = "cv_type",
    values_to = "cv_value"
  ) %>%
  filter(
    (trait == "Wing.Length" & cv_type == "true_cv_wing") |
    (trait == "Beak.Length_Nares" & cv_type == "true_cv_beak_length") |
    (trait == "Beak.Depth" & cv_type == "true_cv_beak_depth") |
    (trait == "Tarsus.Length" & cv_type == "true_cv_tarsus_length")
  )

ggplot(plot_data, aes(x = cv_value, y = log(specimens_required))) +
  geom_point() +
  facet_wrap(~trait) +
  labs(x = "True CV", y = "Log Specimens Required", title = "Specimens Required vs CV by Trait") +
  theme_minimal()

## Run an LM for each trait to see if there is a relationship between CV and specimens_required

lm_beak_length <- lm(log(specimens_required) ~ true_cv_beak_length, data = species_requirements_merged)
summary(lm_beak_length)

lm_beak_depth <- lm(log(specimens_required) ~ true_cv_beak_depth, data = species_requirements_merged)
summary(lm_beak_depth)

lm_wing_length <- lm(log(specimens_required) ~ true_cv_wing, data = species_requirements_merged)
summary(lm_wing_length)

lm_tarsus_length <- lm(log(specimens_required) ~ true_cv_tarsus_length, data = species_requirements_merged)
summary(lm_tarsus_length)
