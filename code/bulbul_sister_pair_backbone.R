# Sister-pair range-matched variance analysis backbone
#
# Edit the CONFIG block below for a new dataset, then run from the project root:
#   Rscript code/bulbul_sister_pair_backbone.R
#
# Intended final workflow:
# 1. Add Latitude and Longitude columns to the specimen table.
# 2. Fill in the sister_pairs table.
# 3. Use extent_rule = "restricted_convex_hull".
# 4. Run repeated rarefaction and paired models.

rm(list = ls())
set.seed(123454)

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
})

#===========================================================
# CONFIG: edit these values when swapping datasets
#===========================================================

data_path <- "data/Bulbul_AVONET_data_integrated.csv"
output_dir <- "results/Sister_pair_analysis"

species_col <- "AviList"
sex_col <- "Sex"
source_col <- "Source"
specimen_id_col <- "Specimen.number"
locality_col <- "Locality"
country_col <- "Country_WRI"

# These columns are placeholders in the current integrated dataset.
# Add them later using decimal degrees, WGS84:
#   Latitude:  north positive, south negative
#   Longitude: east positive, west negative
lat_col <- "Latitude"
lon_col <- "Longitude"

trait_cols <- c(
  "Beak.Length_Culmen",
  "Beak.Length_Nares",
  "Beak.Width",
  "Beak.Depth",
  "Tarsus.Length",
  "Wing.Length",
  "Kipps.Distance",
  "Tail.Length"
)

# Main final option: "restricted_convex_hull".
# Temporary pilot option before coordinates exist: "restricted_countries".
extent_rule <- "restricted_countries"

# Sample sizes to stress-test. Each pair uses the largest requested n that both
# species can support inside the matched extent; lower n values are also run.
target_sample_sizes <- c(20, 30, 40)
n_repeats <- 500

# Optional nuisance balancing. The sampler preserves these variables exactly
# only where both species have enough specimens in each stratum. If that is too
# strict for a pair, it falls back to unstratified rarefaction and records it.
balance_vars <- c("Sex", "Source")

# If a restricted species has too few unique coordinate points for a convex
# hull, use a point/line buffer around its occupied records.
extent_buffer_km <- 50

# Fill these with your chosen sister pairs. The status labels are generated
# from restricted_species and widespread_species.
sister_pairs <- tribble(
  ~pair_id, ~restricted_species, ~widespread_species,
  # "pair_01", "Restricted species name", "Widespread sister name",
  # "pair_02", "Restricted species name", "Widespread sister name",
  # "pair_03", "Restricted species name", "Widespread sister name",
  # "pair_04", "Restricted species name", "Widespread sister name",
  # "pair_05", "Restricted species name", "Widespread sister name",
  # "pair_06", "Restricted species name", "Widespread sister name",
  # "pair_07", "Restricted species name", "Widespread sister name",
  # "pair_08", "Restricted species name", "Widespread sister name",
  # "pair_09", "Restricted species name", "Widespread sister name",
  # "pair_10", "Restricted species name", "Widespread sister name"
)

#===========================================================
# Helpers
#===========================================================

make_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}

clean_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

required_columns <- function() {
  unique(c(
    species_col, sex_col, source_col, specimen_id_col, locality_col,
    country_col, trait_cols
  ))
}

ensure_coordinate_placeholders <- function(dat) {
  if (!lat_col %in% names(dat)) {
    dat[[lat_col]] <- NA_real_
  }
  if (!lon_col %in% names(dat)) {
    dat[[lon_col]] <- NA_real_
  }
  dat[[lat_col]] <- clean_numeric(dat[[lat_col]])
  dat[[lon_col]] <- clean_numeric(dat[[lon_col]])
  dat
}

validate_inputs <- function(dat, pairs) {
  missing_cols <- setdiff(required_columns(), names(dat))
  if (length(missing_cols) > 0) {
    stop(
      "The dataset is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (nrow(pairs) == 0) {
    species_counts <- dat |>
      count(.data[[species_col]], sort = TRUE, name = "n_specimens")

    write_csv(
      species_counts,
      file.path(output_dir, "species_counts_for_pair_selection.csv")
    )

    stop(
      "Fill in the sister_pairs table near the top of this script before ",
      "running the analysis. I wrote species counts to ",
      file.path(output_dir, "species_counts_for_pair_selection.csv"),
      call. = FALSE
    )
  }

  pair_species <- unique(c(pairs$restricted_species, pairs$widespread_species))
  missing_species <- setdiff(pair_species, unique(dat[[species_col]]))
  if (length(missing_species) > 0) {
    stop(
      "These sister-pair species are not present in ", species_col, ": ",
      paste(missing_species, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

has_complete_coordinates <- function(dat) {
  all(c(lat_col, lon_col) %in% names(dat)) &&
    any(!is.na(dat[[lat_col]]) & !is.na(dat[[lon_col]]))
}

make_specimen_sf <- function(dat) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop(
      "Package 'sf' is required for extent_rule = 'restricted_convex_hull'. ",
      "Install or repair sf/GDAL before running coordinate-based extents.",
      call. = FALSE
    )
  }

  dat |>
    filter(!is.na(.data[[lat_col]]), !is.na(.data[[lon_col]])) |>
    sf::st_as_sf(
      coords = c(lon_col, lat_col),
      crs = 4326,
      remove = FALSE
    )
}

make_restricted_extent <- function(restricted_sf) {
  if (nrow(restricted_sf) == 0) {
    stop("No georeferenced restricted-species specimens for this pair.")
  }

  restricted_projected <- restricted_sf |>
    sf::st_transform(6933)

  unique_points <- restricted_projected |>
    sf::st_coordinates() |>
    as_tibble() |>
    distinct(X, Y)

  if (nrow(unique_points) >= 3) {
    restricted_projected |>
      sf::st_union() |>
      sf::st_convex_hull() |>
      sf::st_as_sf() |>
      sf::st_set_crs(6933) |>
      sf::st_transform(4326)
  } else {
    restricted_projected |>
      sf::st_union() |>
      sf::st_buffer(dist = extent_buffer_km * 1000) |>
      sf::st_as_sf() |>
      sf::st_set_crs(6933) |>
      sf::st_transform(4326)
  }
}

filter_inside_extent <- function(pair_dat, restricted_species) {
  specimen_sf <- make_specimen_sf(pair_dat)
  restricted_sf <- specimen_sf |>
    filter(.data[[species_col]] == restricted_species)
  extent <- make_restricted_extent(restricted_sf)
  inside <- lengths(sf::st_intersects(specimen_sf, extent)) > 0
  sf::st_drop_geometry(specimen_sf[inside, ])
}

filter_restricted_countries <- function(pair_dat, restricted_species) {
  restricted_countries <- pair_dat |>
    filter(.data[[species_col]] == restricted_species) |>
    pull(.data[[country_col]]) |>
    discard(~ is.na(.x) || .x == "") |>
    unique()

  pair_dat |>
    filter(.data[[country_col]] %in% restricted_countries)
}

prepare_pair_extent <- function(dat, pair_row) {
  pair_dat <- dat |>
    filter(.data[[species_col]] %in% c(
      pair_row$restricted_species,
      pair_row$widespread_species
    )) |>
    mutate(
      pair_id = pair_row$pair_id,
      range_status = if_else(
        .data[[species_col]] == pair_row$restricted_species,
        "restricted",
        "widespread"
      )
    )

  if (extent_rule == "restricted_convex_hull") {
    if (!has_complete_coordinates(pair_dat)) {
      stop(
        "extent_rule is 'restricted_convex_hull', but no usable Latitude/",
        "Longitude values are present for pair ", pair_row$pair_id, ".",
        call. = FALSE
      )
    }
    filter_inside_extent(pair_dat, pair_row$restricted_species)
  } else if (extent_rule == "restricted_countries") {
    filter_restricted_countries(pair_dat, pair_row$restricted_species)
  } else {
    stop("Unknown extent_rule: ", extent_rule, call. = FALSE)
  }
}

usable_balance_vars <- function(dat) {
  balance_vars[balance_vars %in% names(dat)]
}

sample_unstratified_by_status <- function(dat, status_col, n_target) {
  counts <- dat |>
    count(.data[[status_col]], name = "n_available")

  if (!all(c("restricted", "widespread") %in% counts[[status_col]]) ||
      min(counts$n_available) < n_target) {
    return(NULL)
  }

  dat |>
    group_by(.data[[status_col]]) |>
    slice_sample(n = n_target) |>
    ungroup()
}

sample_stratified <- function(dat, status_col, n_target, vars) {
  if (length(vars) == 0) {
    return(sample_unstratified_by_status(dat, status_col, n_target))
  }

  strata_counts <- dat |>
    filter(!if_any(all_of(vars), is.na)) |>
    count(.data[[status_col]], across(all_of(vars)), name = "n")

  if (nrow(strata_counts) == 0) {
    return(sample_unstratified_by_status(dat, status_col, n_target))
  }

  strata_wide <- strata_counts |>
    pivot_wider(
      names_from = all_of(status_col),
      values_from = n,
      values_fill = 0
    )

  if (!all(c("restricted", "widespread") %in% names(strata_wide))) {
    return(sample_unstratified_by_status(dat, status_col, n_target))
  }

  strata_plan <- strata_wide |>
    mutate(n_available = pmin(restricted, widespread)) |>
    filter(n_available > 0)

  if (sum(strata_plan$n_available) < n_target) {
    return(sample_unstratified_by_status(dat, status_col, n_target))
  }

  strata_plan <- strata_plan |>
    mutate(
      ideal = n_target * n_available / sum(n_available),
      n_draw = floor(ideal),
      remainder = ideal - n_draw
    )

  leftover <- n_target - sum(strata_plan$n_draw)
  if (leftover > 0) {
    add_to <- order(strata_plan$remainder, decreasing = TRUE)[seq_len(leftover)]
    strata_plan$n_draw[add_to] <- strata_plan$n_draw[add_to] + 1
  }

  strata_plan <- strata_plan |>
    filter(n_draw > 0) |>
    select(all_of(vars), n_draw)

  dat |>
    inner_join(strata_plan, by = vars) |>
    group_by(across(all_of(c(status_col, vars)))) |>
    group_modify(~ slice_sample(.x, n = unique(.x$n_draw))) |>
    ungroup() |>
    select(-n_draw)
}

draw_rarefied_pair <- function(pair_dat, trait, n_target) {
  trait_dat <- pair_dat |>
    filter(!is.na(.data[[trait]])) |>
    mutate(
      range_status = factor(
        range_status,
        levels = c("restricted", "widespread")
      )
    )

  counts <- trait_dat |>
    count(range_status, name = "n_available")

  if (!all(c("restricted", "widespread") %in% counts$range_status)) {
    return(NULL)
  }

  if (min(counts$n_available) < n_target) {
    return(NULL)
  }

  vars <- usable_balance_vars(trait_dat)
  sampled <- trait_dat |>
    sample_stratified("range_status", n_target, vars)

  if (is.null(sampled)) {
    return(NULL)
  }

  if (nrow(sampled) != n_target * 2) {
    sampled <- trait_dat |>
      group_by(range_status) |>
      slice_sample(n = n_target) |>
      ungroup()

    sampled$balanced_sampling <- FALSE
  } else {
    sampled$balanced_sampling <- length(vars) > 0
  }

  sampled
}

summarise_variance <- function(sampled_dat, trait, replicate_id, n_target) {
  sampled_dat |>
    group_by(pair_id, range_status, .data[[species_col]]) |>
    summarise(
      trait = trait,
      replicate_id = replicate_id,
      target_n = n_target,
      n_used = sum(!is.na(.data[[trait]])),
      trait_mean = mean(.data[[trait]], na.rm = TRUE),
      trait_sd = sd(.data[[trait]], na.rm = TRUE),
      trait_var = var(.data[[trait]], na.rm = TRUE),
      trait_cv = trait_sd / trait_mean,
      log_sd = log(trait_sd),
      log_cv = log(trait_cv),
      balanced_sampling = all(balanced_sampling),
      .groups = "drop"
    ) |>
    rename(species = all_of(species_col))
}

run_rarefaction <- function(extent_data) {
  results <- vector("list", length = 0)
  skipped <- vector("list", length = 0)
  result_i <- 1
  skipped_i <- 1

  for (pair_id_i in unique(extent_data$pair_id)) {
    pair_dat <- extent_data |>
      filter(pair_id == pair_id_i)

    for (trait in trait_cols) {
      for (n_target in target_sample_sizes) {
        for (rep_i in seq_len(n_repeats)) {
          sampled <- draw_rarefied_pair(pair_dat, trait, n_target)

          if (is.null(sampled)) {
            skipped[[skipped_i]] <- tibble(
              pair_id = pair_id_i,
              trait = trait,
              target_n = n_target,
              reason = "insufficient specimens after extent/trait filtering"
            )
            skipped_i <- skipped_i + 1
            break
          }

          results[[result_i]] <- summarise_variance(
            sampled,
            trait = trait,
            replicate_id = rep_i,
            n_target = n_target
          )
          result_i <- result_i + 1
        }
      }
    }
  }

  list(
    metrics = bind_rows(results),
    skipped = bind_rows(skipped) |> distinct()
  )
}

make_pair_contrasts <- function(metrics) {
  metrics |>
    select(
      pair_id, trait, replicate_id, target_n, range_status,
      species, n_used, trait_sd, trait_cv, log_sd, log_cv
    ) |>
    pivot_wider(
      names_from = range_status,
      values_from = c(species, n_used, trait_sd, trait_cv, log_sd, log_cv)
    ) |>
    mutate(
      delta_log_sd = log_sd_widespread - log_sd_restricted,
      delta_log_cv = log_cv_widespread - log_cv_restricted,
      ratio_sd = trait_sd_widespread / trait_sd_restricted,
      ratio_cv = trait_cv_widespread / trait_cv_restricted
    )
}

fit_trait_models <- function(metrics) {
  metrics |>
    filter(is.finite(log_cv), n_used > 1) |>
    group_by(trait, target_n) |>
    group_modify(~ {
      dat <- .x

      if (n_distinct(dat$pair_id) < 2 ||
          n_distinct(dat$range_status) < 2) {
        return(tibble(
          term = "range_statuswidespread",
          estimate = NA_real_,
          std_error = NA_real_,
          model = "not enough pairs"
        ))
      }

      fit <- lmer(
        log_cv ~ range_status + log(n_used) + (1 | pair_id) +
          (1 | pair_id:replicate_id),
        data = dat,
        REML = FALSE
      )

      coef_table <- coef(summary(fit))

      tibble(
        term = rownames(coef_table),
        estimate = coef_table[, "Estimate"],
        std_error = coef_table[, "Std. Error"],
        model = "lmer_log_cv"
      )
    }) |>
    ungroup()
}

plot_contrasts <- function(contrasts) {
  p <- contrasts |>
    ggplot(aes(x = factor(target_n), y = delta_log_cv)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_boxplot(outlier.alpha = 0.2) +
    facet_grid(trait ~ pair_id, scales = "free_y") +
    labs(
      x = "Rarefied specimens per species",
      y = "log(CV widespread subset) - log(CV restricted species)"
    ) +
    theme_bw(base_size = 10)

  ggsave(
    file.path(output_dir, "paired_delta_log_cv_by_pair.pdf"),
    p,
    width = 12,
    height = 10
  )
}

#===========================================================
# Main workflow
#===========================================================

make_dir(output_dir)

bulbuls <- read.csv(data_path, stringsAsFactors = FALSE, check.names = FALSE) |>
  ensure_coordinate_placeholders()

for (trait in trait_cols[trait_cols %in% names(bulbuls)]) {
  bulbuls[[trait]] <- clean_numeric(bulbuls[[trait]])
}

validate_inputs(bulbuls, sister_pairs)

if (extent_rule == "restricted_countries") {
  message(
    "Using country-level matched extents because coordinates are not yet ",
    "available. Switch extent_rule to 'restricted_convex_hull' after adding ",
    lat_col, " and ", lon_col, "."
  )
}

extent_data <- sister_pairs |>
  split(seq_len(nrow(sister_pairs))) |>
  map_dfr(~ prepare_pair_extent(bulbuls, .x))

extent_summary <- extent_data |>
  pivot_longer(
    cols = all_of(trait_cols),
    names_to = "trait",
    values_to = "trait_value"
  ) |>
  group_by(pair_id, range_status, .data[[species_col]], trait) |>
  summarise(
    n_in_extent = n(),
    n_with_trait = sum(!is.na(trait_value)),
    n_georeferenced = sum(!is.na(.data[[lat_col]]) & !is.na(.data[[lon_col]])),
    .groups = "drop"
  ) |>
  rename(species = all_of(species_col))

write_csv(
  extent_summary,
  file.path(output_dir, "extent_specimen_counts_by_pair.csv")
)

rarefaction <- run_rarefaction(extent_data)

if (nrow(rarefaction$metrics) == 0) {
  write_csv(
    rarefaction$skipped,
    file.path(output_dir, "skipped_rarefactions.csv")
  )

  stop(
    "No rarefaction draws were possible. Check ",
    file.path(output_dir, "extent_specimen_counts_by_pair.csv"),
    " and lower target_sample_sizes if needed.",
    call. = FALSE
  )
}

write_csv(
  rarefaction$metrics,
  file.path(output_dir, "rarefied_trait_variance_metrics.csv")
)

write_csv(
  rarefaction$skipped,
  file.path(output_dir, "skipped_rarefactions.csv")
)

pair_contrasts <- make_pair_contrasts(rarefaction$metrics)

write_csv(
  pair_contrasts,
  file.path(output_dir, "paired_variance_contrasts.csv")
)

contrast_summary <- pair_contrasts |>
  group_by(pair_id, trait, target_n) |>
  summarise(
    n_repeats = n(),
    mean_delta_log_cv = mean(delta_log_cv, na.rm = TRUE),
    median_delta_log_cv = median(delta_log_cv, na.rm = TRUE),
    q025_delta_log_cv = quantile(delta_log_cv, 0.025, na.rm = TRUE),
    q975_delta_log_cv = quantile(delta_log_cv, 0.975, na.rm = TRUE),
    mean_ratio_cv = mean(ratio_cv, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  contrast_summary,
  file.path(output_dir, "paired_variance_contrast_summary.csv")
)

model_summary <- fit_trait_models(rarefaction$metrics)

write_csv(
  model_summary,
  file.path(output_dir, "paired_lmer_log_cv_models.csv")
)

plot_contrasts(pair_contrasts)

message("Done. Outputs written to: ", output_dir)
