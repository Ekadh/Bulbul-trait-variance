rm(list = ls())
set.seed(123454)

#===========================================================
# LIBRARIES
#===========================================================

library(tidyverse)
library(ape)
library(cmdstanr)
library(Matrix)
library(scales)
library(patchwork)

#===========================================================
# LOAD DATA
#===========================================================

bulbul <- read.csv(
  "../data/Bulbul_AVONET_data_integrated_trimmed.csv"
)

bulbul_density <- read.csv("../data/Density/Bulbul_density_data.csv")
bulbul_ESI <- read.csv("../data/BIRDBASE/Bulbul_Birdbase_data.csv")
bulbul_range <- read.csv("../data/Bulbul_rangesize.csv")

# Add Density and ESI to the data frame, matching by species
bulbul <- bulbul %>%
  left_join(bulbul_density %>% select(Species, Predicted.Density.n.km2), by = c("AviList" = "Species")) %>%
  left_join(bulbul_ESI %>% select(AviList.v1.2025, ESI), by = c("AviList" = "AviList.v1.2025")) %>%
  left_join(bulbul_range %>% select(AviList, Range.Size), by = c("AviList" = "AviList"))

#===========================================================
# TRAITS TO ANALYSE
#===========================================================

traits <- tibble(

  trait = c(

    "Wing.Length",
    "Beak.Length_Culmen",
    "Tail.Length",
    "Tarsus.Length",
    "Beak.Depth",
    "Beak.Width",
    "Beak.Length_Nares"

  ),

  short_name = c(

    "Wing",
    "Culmen",
    "Tail",
    "Tarsus",
    "Beak Depth",
    "Beak Width",
    "Beak Length (Nares)"

  )

)

#===========================================================
# BASIC FILTERING
#===========================================================

bulbul <- bulbul %>%
  filter(
    !is.na(Data.type)
  )

bulbul$species <- factor(
  bulbul$AviList
)

bulbul$source <- factor(
  bulbul$Data.type
)


#===========================================================
# PHYLOGENY
#===========================================================

phylogeny <- read.tree(
  "../data/Phylogeny/mcc_dated_clements_McTavish_2025.nex"
)

phylogeny <- phylogeny[[6]]

phylogeny$tip.label <-
  gsub(
    "_",
    " ",
    phylogeny$tip.label
  )

#===========================================================
# COMPILE STAN MODEL
#===========================================================

mod <- cmdstan_model(
  "latent_variance_rangesize.stan"
)

run_trait_model <- function(
  trait,
  short_name
){

  cat("\n")
  cat("---------------------------------------\n")
  cat("Running:", short_name, "\n")
  cat("---------------------------------------\n")

  dat <-
    bulbul %>%
    filter(
      !is.na(.data[[trait]])
    ) %>%
    mutate(
      trait_value = .data[[trait]]
    )

  trait_mean <- mean(dat$trait_value)

  species_levels <-
    levels(
      factor(dat$species)
    )

  phy_trim <-
    drop.tip(
      phylogeny,
      setdiff(
        phylogeny$tip.label,
        species_levels
      )
    )

  A <-
    ape::vcv(
      phy_trim,
      corr = TRUE
    )

  A <-
    A[
      species_levels,
      species_levels
    ]

  L_A <-
    t(
      chol(A)
    )

  species_df <-

    tibble(
      species = species_levels
    ) %>%

    left_join(

      dat %>%

        group_by(species) %>%

        summarise(

          log_range =
            first(log(Range.Size)),

          density =
            first(Predicted.Density.n.km2),

          ESI =
            first(ESI),

          .groups = "drop"

        ),

      by = "species"

    )

  dat$species <-

    factor(
      dat$species,
      levels = species_levels
    )

  dat$source <-

    factor(
      dat$source
    )

  has_range <-
    !is.na(
      species_df$log_range
    )

  has_density <-
    !is.na(
      species_df$density
    )

  has_ESI <-
    !is.na(
      species_df$ESI
    )

  log_range <-
    species_df$log_range

  density <-
    species_df$density

  ESI <-
    species_df$ESI

  log_range[
    !has_range
  ] <- 0

  density[
    !has_density
  ] <- 0

  ESI[
    !has_ESI
  ] <- 0

  stan_data <- list(

    N =
      nrow(dat),

    S =
      nlevels(dat$species),

    wing =
      dat$trait_value,

    species =
      as.integer(dat$species),

    n_source =
      nlevels(dat$source),

    source =
      as.integer(dat$source),

    log_range =
      log_range,

    density =
      density,

    ESI =
      ESI,
    
    trait_mean =
      trait_mean,

    has_range =
      as.integer(has_range),

    has_density =
      as.integer(has_density),

    has_ESI =
      as.integer(has_ESI),

    L_A =
      L_A

  )

  fit <-

    mod$sample(

      data = stan_data,

      chains = 4,

      parallel_chains = 4,

      iter_warmup = 2000,

      iter_sampling = 4000,

      adapt_delta = 0.99,

      refresh = 100,

      max_treedepth = 15

    )

  CV_summary <-

    fit$summary(
      "CV_species"
    )

  CV_summary$species <-
    species_levels

  raw_cv <-

    dat %>%

    group_by(species) %>%

    summarise(

      CV_raw =
        sd(trait_value) /
        mean(trait_value),

      .groups = "drop"

    )

  cv_compare <-

    left_join(

      CV_summary %>%

        select(
          species,
          CV_stan = mean
        ),

      raw_cv,

      by = "species"

    )

  plot_df <-

    left_join(

      CV_summary %>%

        select(
          species,
          CV = mean
        ),

      species_df,

      by = "species"

    ) %>%

    mutate(

      Range_1000km2 =
        exp(log_range) /
        1000

    )

  list(

    trait = trait,

    short_name = short_name,

    fit = fit,

    species_df = species_df,

    CV_summary = CV_summary,

    cv_compare = cv_compare,

    plot_df = plot_df,

    raw_data = dat

  )

}

plot_trait_results <- function(result){

  fit <- result$fit

  trait <- result$trait

  short_name <- result$short_name

  species_df <- result$species_df

  CV_summary <- result$CV_summary

  cv_compare <- result$cv_compare

  plot_df <- result$plot_df

  ###########################################################
  ## Raw CV vs latent CV
  ###########################################################

  p_cv <- ggplot(

    cv_compare,

    aes(
      x = CV_raw,
      y = CV_stan
    )

  ) +

    geom_point(size = 2) +

    geom_smooth(
      method = "lm",
      se = TRUE
    ) +

    labs(

      title = short_name,

      x = "Raw CV",

      y = "Estimated latent CV"

    ) +

    theme_classic()

  ###########################################################
  ## CV vs Range size
  ###########################################################

  p_range <- ggplot(

    plot_df,

    aes(

      x = Range_1000km2,

      y = CV

    )

  ) +

    geom_point(size = 2) +

    geom_smooth(
      method = "lm",
      se = TRUE
    ) +

    scale_x_log10(
      labels = label_comma()
    ) +

    labs(

      title = short_name,

      x = expression(
        paste(
          "Range size (",
          1000,
          " km"^2,
          ")"
        )
      ),

      y = "Estimated latent CV"

    ) +

    theme_classic()

  ###########################################################
  ## Posterior beta
  ###########################################################

  posterior <-

    fit$draws(
      "beta_CV_range"
    )

  p_beta <- ggplot(

    data.frame(
      beta = as.vector(posterior)
    ),

    aes(
      x = beta
    )

  ) +

    geom_density(
      fill = "lightblue"
    ) +

    geom_vline(
      xintercept = 0,
      linetype = 2
    ) +

    labs(

      title = short_name,

      x = expression(beta[CV]),

      y = "Density"

    ) +

    theme_classic()

  ###########################################################

  list(

    cv_plot = p_cv,

    range_plot = p_range,

    beta_plot = p_beta

  )

}

############################################################
## RUN ALL TRAITS
############################################################

results <-

  vector(

    "list",

    nrow(traits)

  )

names(results) <-

  traits$short_name

for(i in seq_len(nrow(traits))){

  res <-

    run_trait_model(

      trait = traits$trait[i],

      short_name = traits$short_name[i]

    )

  results[[i]] <- res

  ##########################################################
  ## Print summaries
  ##########################################################

  cat("\n")

  print(

    res$fit$summary(

      variables = c(

        "beta_CV_range",

        "beta_CV_density",

        "beta_CV_ESI"

      )

    )

  )

  ##########################################################
  ## Save model
  ##########################################################

  saveRDS(

    res$fit,

    file =

      paste0(

        "../results/stan/",

        res$short_name,

        "_fit_Main.rds"

      )

  )
}

####POST-HOC CHECKS
wing <- results[[1]]
wing$fit$summary(
  variables = c(
    "beta_CV_range",
    "beta_CV_density",
    "beta_CV_ESI"
  )
)

#Descending
wing$CV_summary |>
  arrange(desc(mean)) |>
  select(species, mean, q5, q95) |>
  head(20)

#Ascending
wing$CV_summary |>
  arrange(mean) |>
  select(species, mean, q5, q95) |>
  head(20)

cor.test(
  wing$cv_compare$CV_raw,
  wing$cv_compare$CV_stan
)

culmen <- results[[2]]
culmen$fit$summary(
  variables = c(
    "beta_CV_range",
    "beta_CV_density",
    "beta_CV_ESI"
  )
)

#Descending
culmen$CV_summary |>
  arrange(desc(mean)) |>
  select(species, mean, q5, q95) |>
  head(20)

#Ascending
culmen$CV_summary |>
  arrange(mean) |>
  select(species, mean, q5, q95) |>
  head(20)

cor.test(
  culmen$cv_compare$CV_raw,
  culmen$cv_compare$CV_stan
)

tail <- results[[3]]
tail$fit$summary(
  variables = c(
    "beta_CV_range",
    "beta_CV_density",
    "beta_CV_ESI"
  )
)

#Descending
tail$CV_summary |>
  arrange(desc(mean)) |>
  select(species, mean, q5, q95) |>
  head(20)

#Ascending
tail$CV_summary |>
  arrange(mean) |>
  select(species, mean, q5, q95) |>
  head(20)

cor.test(
  tail$cv_compare$CV_raw,
  tail$cv_compare$CV_stan
)

tarsus <- results[[4]]
tarsus$fit$summary(
  variables = c(
    "beta_CV_range",
    "beta_CV_density",
    "beta_CV_ESI"
  )
)

#Descending
tarsus$CV_summary |>
  arrange(desc(mean)) |>
  select(species, mean, q5, q95) |>
  head(20)

#Ascending
tarsus$CV_summary |>
  arrange(mean) |>
  select(species, mean, q5, q95) |>
  head(20)

cor.test(
  tarsus$cv_compare$CV_raw,
  tarsus$cv_compare$CV_stan
)

beak_depth <- results[[5]]
beak_depth$fit$summary(
  variables = c(
    "beta_CV_range",
    "beta_CV_density",
    "beta_CV_ESI"
  )
)

beak_width <- results[[6]]
beak_width$fit$summary(
  variables = c(
    "beta_CV_range",
    "beta_CV_density",
    "beta_CV_ESI"
  )
)

beak_nares <- results[[7]]
beak_nares$fit$summary(
  variables = c(
    "beta_CV_range",
    "beta_CV_density",
    "beta_CV_ESI"
  )
)


avonet_means <- bulbul |>
  group_by(AviList) |>
  summarise(
    Tail.Length = mean(Tail.Length, na.rm = TRUE),
    .groups = "drop"
  )

avonet_points <- merge(
  bulbul,
  avonet_means,
  by="AviList",
  suffixes=c("_specimen","_mean")
  )

avonet_long <- avonet_points %>%
  transmute(
    AviList,
    Trait = "Tail.Length",
    specimen = Tail.Length_specimen,
    mean = Tail.Length_mean
  )

ggplot() +

  geom_point(
    data = avonet_long,
    aes(mean, specimen, colour = AviList),
    alpha = 0.4
  ) +

  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +

  theme_classic() +
  theme(legend.position = "none") +

  labs(
    x = "AVONET species mean",
    y = "Specimen measurement"
  )


####Test analysis: Check distributions of latent CVs across traits
summary(results[[7]]$CV_summary$mean)
sd(results[[4]]$CV_summary$mean)
IQR(results[[4]]$CV_summary$mean)

####To check if the most variable species are the same across traits, we can check correlations of the latent CVs across traits
cor(results$Tarsus$CV_summary$mean,
    results$Wing$CV_summary$mean)

cor(results$Tarsus$CV_summary$mean,
    results$Tail$CV_summary$mean)

cor(results$Culmen$CV_summary$mean,
    results$Tarsus$CV_summary$mean)

#I need to now load the rds files for all the traits and then make a faceted plot, one facet for each response variable (range size, density and ESI) and then plot the posterior distributions of the beta coefficients for each trait, as well as a dot and credible interval just below each distribution to show the mean and 95% credible interval for each trait. I will also add a dashed vertical line at 0 to show which traits have a significant effect on the latent CVs (the ones that don't include 0 in their credible interval, will be opaque color, while non-significant will be slightly transparent).

# ============================================================
# Posterior distributions of CV effects across morphological traits
# ============================================================

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(posterior)

# ------------------------------------------------------------
# 1. Find all CmdStanR fit files
# ------------------------------------------------------------

rds_files <- list.files(
  path = "../results/stan/",
  pattern = "_fit\\.rds$",
  full.names = TRUE
)

# Sanity check
rds_files


# ------------------------------------------------------------
# 2. Function to extract the three beta coefficients
#    from each CmdStanMCMC fit
# ------------------------------------------------------------

read_beta_draws <- function(f) {

  # Read fit
  fit <- readRDS(f)

  # Check that this is a CmdStanR MCMC fit
  if (!inherits(fit, "CmdStanMCMC")) {
    warning(
      "Skipping ", basename(f),
      ": object is not a CmdStanMCMC fit."
    )
    return(NULL)
  }

  # Get clean trait name from filename
  # e.g. "Beak Depth_fit.rds" -> "Beak Depth"
  trait <- tools::file_path_sans_ext(basename(f)) %>%
    sub("_fit$", "", .)

  # Parameters of interest
  beta_vars <- c(
    "beta_CV_range",
    "beta_CV_density",
    "beta_CV_ESI"
  )

  # Extract posterior draws directly from CmdStanR object
  draws <- fit$draws(
    variables = beta_vars,
    format = "draws_df"
  ) %>%
    as.data.frame()

  # Convert to long format
  draws_long <- draws %>%
    select(
      all_of(beta_vars),
      any_of(c(".chain", ".iteration", ".draw"))
    ) %>%
    pivot_longer(
      cols = all_of(beta_vars),
      names_to = "parameter",
      values_to = "beta"
    ) %>%
    mutate(
      trait = trait
    )

  return(draws_long)
}


# ------------------------------------------------------------
# 3. Extract draws from all seven trait models
# ------------------------------------------------------------

all_draws <- map_dfr(
  rds_files,
  read_beta_draws
)


# ------------------------------------------------------------
# 4. Clean response names
# ------------------------------------------------------------

all_draws <- all_draws %>%
  mutate(
    response = recode(
      parameter,
      beta_CV_range   = "Range size",
      beta_CV_density = "Population density",
      beta_CV_ESI     = "Ecological specialisation"
    )
  )


# ------------------------------------------------------------
# 5. Set desired trait and facet ordering
# ------------------------------------------------------------

trait_order <- c(
  "Culmen",
  "Beak Length (Nares)",
  "Beak Width",
  "Beak Depth",
  "Tarsus",
  "Wing",
  "Tail"
)

response_order <- c(
  "Range size",
  "Population density",
  "Ecological specialisation"
)

all_draws <- all_draws %>%
  mutate(
    trait = factor(
      trait,
      levels = trait_order
    ),
    response = factor(
      response,
      levels = response_order
    )
  )


# ------------------------------------------------------------
# 6. Calculate posterior mean and 95% credible interval
# ------------------------------------------------------------

posterior_summary <- all_draws %>%
  group_by(trait, response) %>%
  summarise(
    mean_beta = mean(beta),
    lower_95 = quantile(beta, 0.05),
    upper_95 = quantile(beta, 0.95),
    .groups = "drop"
  ) %>%

  # 95% CrI excludes zero if:
  # both bounds > 0 OR both bounds < 0
  mutate(
    excludes_zero = (
      lower_95 > 0 |
      upper_95 < 0
    ),
    effect_status = if_else(
      excludes_zero,
      "95% CrI excludes 0",
      "95% CrI overlaps 0"
    )
  )


# Inspect numerical results
posterior_summary


# ------------------------------------------------------------
# 7. Calculate density curves explicitly
#
# This gives us full control over:
# - vertical placement of each trait
# - density height
# - space below each density for mean + CrI
# ------------------------------------------------------------

density_data <- all_draws %>%
  filter(
    !is.na(beta),
    !is.na(trait),
    !is.na(response)
  ) %>%
  group_by(trait, response) %>%
  group_modify(~ {

    dens <- density(.x$beta, n = 512)

    tibble(
      beta = dens$x,
      density = dens$y
    )

  }) %>%
  ungroup()


# ------------------------------------------------------------
# 8. Create numeric y positions for traits
# ------------------------------------------------------------

trait_positions <- tibble(
  trait = factor(
    trait_order,
    levels = trait_order
  ),

  # Reverse so first trait appears at top
  y_base = rev(seq_along(trait_order))
)


# Add positions to density curves
density_data <- density_data %>%
  left_join(
    trait_positions,
    by = "trait"
  )


# Add positions to posterior summaries
posterior_summary <- posterior_summary %>%
  left_join(
    trait_positions,
    by = "trait"
  )


# ------------------------------------------------------------
# 9. Add significance / credible-interval status to densities
# ------------------------------------------------------------

density_data <- density_data %>%
  left_join(
    posterior_summary %>%
      select(
        trait,
        response,
        excludes_zero,
        effect_status
      ),
    by = c("trait", "response")
  )


# ------------------------------------------------------------
# 10. Scale each density independently
#
# Each posterior gets the same maximum visual height,
# regardless of how narrow or broad the posterior is.
# ------------------------------------------------------------

density_data <- density_data %>%
  group_by(trait, response) %>%
  mutate(
    density_scaled =
      density / max(density, na.rm = TRUE) * 0.55,

    y_density =
      y_base + density_scaled
  ) %>%
  ungroup()


# ------------------------------------------------------------
# 11. Position mean + 95% CrI just below each density
# ------------------------------------------------------------

posterior_summary <- posterior_summary %>%
  mutate(
    y_interval = y_base - 0.12
  )


# ------------------------------------------------------------
# 12. Plot
# ------------------------------------------------------------

p <- ggplot() +

  # Zero-effect reference line
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6,
    colour = "grey35"
  ) +

  # Posterior density distributions
  geom_ribbon(
    data = density_data,
    aes(
      x = beta,
      ymin = y_base,
      ymax = y_density,
      group = trait,
      fill = response,
      alpha = effect_status
    ),
    colour = NA
  ) +

  # Density outlines
  geom_line(
    data = density_data,
    aes(
      x = beta,
      y = y_density,
      group = trait,
      colour = response,
      alpha = effect_status
    ),
    linewidth = 0.7
  ) +

  # 95% credible intervals below each density
  geom_segment(
    data = posterior_summary,
    aes(
      x = lower_95,
      xend = upper_95,
      y = y_interval,
      yend = y_interval,
      colour = response,
      alpha = effect_status
    ),
    linewidth = 1
  ) +

  # Posterior mean
  geom_point(
    data = posterior_summary,
    aes(
      x = mean_beta,
      y = y_interval,
      colour = response,
      alpha = effect_status
    ),
    size = 2.5
  ) +

  # One separate facet per response variable
  facet_wrap(
    ~ response,
    nrow = 1,
    scales = "free_x"
  ) +

  # Trait labels
  scale_y_continuous(
    breaks = trait_positions$y_base,
    labels = as.character(trait_positions$trait),
    expand = expansion(
      mult = c(0.04, 0.08)
    )
  ) +

  # Opaque when 95% CrI excludes zero;
  # transparent when it overlaps zero
  scale_alpha_manual(
    values = c(
      "95% CrI excludes 0" = 1,
      "95% CrI overlaps 0" = 0.3
    )
  ) +

  labs(
    x = expression(beta["CV"] ~ "posterior estimate"),
    y = NULL,
    fill = "Response variable",
    colour = "Response variable",
    alpha = NULL
  ) +

  scale_fill_manual(
  values = c(
    "Range size" = "#0072B2",
    "Population density" = "#D55E00",
    "Ecological specialisation" = "#009E73"
  )
) +

scale_colour_manual(
  values = c(
    "Range size" = "#0072B2",
    "Population density" = "#D55E00",
    "Ecological specialisation" = "#009E73"
  )
) +

  # Remove alpha legend (95% CrI) and adjust legend text sizes
  guides(alpha = "none") +

  # Gives each facet a clean rectangular panel
  theme_bw(base_size = 12) +

  theme(
    # Explicit rectangular border around every facet
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.7
    ),

    # Remove internal grid lines so posterior shapes stay clean
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),

    # Facet headers as separate strips
    strip.background = element_rect(
      fill = "grey92",
      colour = "black",
      linewidth = 0.7
    ),

    strip.text = element_text(
      face = "bold",
      size = 18
    ),

    axis.text.x = element_text(
      size = 16
    ),

    axis.text.y = element_text(
      size = 16
    ),

    axis.title.x = element_text(
      size = 20,
      margin = margin(t = 8)
    ),

    axis.title.y = element_text(
      size = 20
    ),

    # Visible gap between rectangular facets
    panel.spacing = unit(
      1.2,
      "lines"
    ),

    legend.position = "none"
  )


p

ggsave(
  filename = "../results/stan/CV_beta_posteriors_all_traits.pdf",
  plot = p,
  width = 14,
  height = 8)
