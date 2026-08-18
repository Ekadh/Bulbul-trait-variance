rm(list = ls())
set.seed(123453)

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

bulbul_ESI$ESI <- 2-bulbul_ESI$ESI

# Add Density and ESI to the data frame, matching by species
bulbul <- bulbul %>%
  left_join(bulbul_density %>% dplyr::select(Species, Predicted.Density.n.km2), by = c("AviList" = "Species")) %>%
  left_join(bulbul_ESI %>% dplyr::select(AviList.v1.2025, ESI), by = c("AviList" = "AviList.v1.2025")) %>%
  left_join(bulbul_range %>% dplyr::select(AviList, Range.Size), by = c("AviList" = "AviList"))

#===========================================================
# TRAITS TO ANALYSE
#===========================================================

traits <- tibble(

  trait = c(

    "Wing.Length",
    #"Beak.Length_Culmen",
    #"Tail.Length",
    "Tarsus.Length",
    "Beak.Depth",
    #"Beak.Width",
    "Beak.Length_Nares"

  ),

  short_name = c(

    "Wing",
    #"Culmen",
    #"Tail",
    "Tarsus",
    "Beak Depth",
    #"Beak Width",
    "Beak Length"

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

tarsus <- results[[2]]
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

beak_depth <- results[[3]]
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

beak_nares <- results[[4]]
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

# ============================================================
# Publication figure: posterior CV slopes for four focal traits
# ============================================================
# Each trait is a figure row, with Range size and ESI shown as adjacent
# slope/posterior pairs.  Slope ribbons are 95% posterior credible intervals
# for the expected response, calculated from paired intercept and slope draws.

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)
library(ggbeeswarm)

trait_order <- c(
  "Beak Depth", "Beak Length (Nares)", "Tarsus", "Wing"
)

species_responses <- bulbul |>
  group_by(AviList) |>
  summarise(
    Range.Size = first(Range.Size[!is.na(Range.Size)], default = NA_real_),
    ESI = first(ESI[!is.na(ESI)], default = NA_real_),
    .groups = "drop"
  ) |>
  rename(species = AviList)

extract_trait_figure_data <- function(result, trait) {
  cv_data <- result$CV_summary |>
    transmute(
      species = as.character(species),
      latent_cv = mean,
      log_latent_cv = log(mean)
    ) |>
    filter(is.finite(latent_cv), latent_cv > 0, is.finite(log_latent_cv))

  if (nrow(cv_data) < 2) {
    stop("Insufficient positive latent-CV estimates for ", trait, ".")
  }

  # Retained to ensure the fitted line is restricted to species represented
  # in the observed data, although no species points are plotted.
  observed_data <- cv_data |>
    left_join(species_responses, by = "species")

  model_draws <- result$fit$draws(
    variables = c(
      "alpha_range", "beta_CV_range", "alpha_ESI", "beta_CV_ESI"
    ),
    format = "draws_df"
  ) |>
    as.data.frame()

  beta_draws <- bind_rows(
    transmute(model_draws, alpha = alpha_range, beta = beta_CV_range,
              response = "Range size"),
    # Reflect the entire ESI prediction through zero.  Negating both paired
    # parameters preserves every posterior fitted interval exactly, while the
    # slope and its posterior draws are shown in the niche-breadth direction.
    transmute(model_draws, alpha = -alpha_ESI, beta = -beta_CV_ESI,
              response = "ESI")
  ) |>
    mutate(trait = trait)

  prediction_grid <- tibble(
    trait = trait,
    log_latent_cv = seq(
      min(observed_data$log_latent_cv),
      max(observed_data$log_latent_cv),
      length.out = 200
    )
  )

  list(beta_draws = beta_draws, prediction_grid = prediction_grid)
}

focal_results <- results[trait_order]

if (length(focal_results) != length(trait_order) || any(map_lgl(focal_results, is.null))) {
  stop(
    "The results object must contain: ",
    paste(trait_order, collapse = ", "), "."
  )
}

trait_figure_data <- map2(
  focal_results,
  trait_order,
  extract_trait_figure_data
)

all_draws <- map_dfr(trait_figure_data, "beta_draws") |>
  mutate(
    trait = factor(trait, levels = trait_order),
    response = factor(response, levels = c("Range size", "ESI"))
  )

prediction_grid <- map_dfr(trait_figure_data, "prediction_grid") |>
  mutate(trait = factor(trait, levels = trait_order))

# Draws are deliberately thinned only for plotting; all draws are retained for
# summaries and credible intervals.
plot_draws <- all_draws |>
  group_by(trait, response) |>
  group_modify(~ slice_sample(.x, n = min(500, nrow(.x)))) |>
  ungroup()

posterior_summary <- all_draws |>
  group_by(trait, response) |>
  summarise(
    mean_beta = mean(beta),
    lower_90 = quantile(beta, 0.05),
    upper_90 = quantile(beta, 0.95),
    lower_95 = quantile(beta, 0.025),
    upper_95 = quantile(beta, 0.975),
    .groups = "drop"
  ) |>
  mutate(
    significance = case_when(
      lower_95 > 0 | upper_95 < 0 ~ "**",
      lower_90 > 0 | upper_90 < 0 ~ "*",
      TRUE ~ "NS"
    ),
    significant = significance != "NS"
  )

# CV_summary stores CV on its natural positive scale.  The Stan likelihood
# uses log(CV), so the x grid and fitted relationship use log latent CV.
# The ribbon is the fitted-model 95% credible interval.
slope_data <- all_draws |>
  group_by(trait, response) |>
  group_modify(~ slice_sample(.x, n = min(2000, nrow(.x)))) |>
  ungroup() |>
  inner_join(prediction_grid, by = "trait") |>
  mutate(expected_response = alpha + beta * log_latent_cv) |>
  group_by(trait, response, log_latent_cv) |>
  summarise(
    fitted_mean = mean(expected_response),
    fitted_lower = quantile(expected_response, 0.025),
    fitted_upper = quantile(expected_response, 0.975),
    .groups = "drop"
  )

response_colours <- c(
  "Range size" = "#3B6FB6",
  "ESI" = "#159A78"
)

publication_theme <- theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.65),
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold", size = 15),
    axis.text = element_text(colour = "black", size = 12),
    plot.title = element_text(face = "bold", size = 18, hjust = 0),
    plot.margin = margin(7, 7, 7, 7)
  )

make_response_pair <- function(
  focal_trait, focal_response, show_x_title
) {
  trait_summary <- posterior_summary |>
    filter(trait == focal_trait, response == focal_response) |>
    mutate(
      annotation_y = max(c(abs(lower_95), abs(upper_95))) * 1.28,
      annotation_line_y = annotation_y * 0.91
    )

  trait_slopes <- slope_data |>
    filter(trait == focal_trait, response == focal_response)

  trait_draws <- plot_draws |>
    filter(trait == focal_trait, response == focal_response) |>
    left_join(
      trait_summary |>
        dplyr::select(response, mean_beta, lower_95, upper_95, significance,
                      significant, annotation_y, annotation_line_y),
      by = "response"
    )

  slope_plot <- ggplot(trait_slopes, aes(log_latent_cv, fitted_mean, colour = response, fill = response)) +
    geom_ribbon(aes(ymin = fitted_lower, ymax = fitted_upper), alpha = 0.20, colour = NA) +
    geom_line(linewidth = 1) +
    scale_colour_manual(values = response_colours) +
    scale_fill_manual(values = response_colours) +
    labs(
      title = if (focal_response == "Range size") focal_trait else NULL,
      x = if (show_x_title) "log(Latent variability)" else NULL,
      y = NULL
    ) +
    publication_theme +
    theme(legend.position = "none")

posterior_plot <-

ggplot(

  trait_draws,

  aes(

    x = 0,

    y = beta,

    colour = response,

    fill = response

  )

) +

############################################################
## Zero line
############################################################

geom_hline(

  yintercept = 0,

  linetype = "dashed",

  colour = "black",

  linewidth = 0.6

) +

############################################################
## Posterior draws
############################################################

geom_quasirandom(

  width = 0.06,

  alpha = 0.18,

  size = 1.5,

  varwidth = TRUE,

  stroke = 0,

  show.legend = FALSE

) +

############################################################
## 95% credible interval
############################################################

geom_linerange(

  aes(

    x = 0,

    ymin = lower_95,

    ymax = upper_95,

    colour = response

  ),

  data = trait_summary,

  inherit.aes = FALSE,

  linewidth = 1.4,

  colour = "black"

) +

############################################################
## Mean
############################################################

geom_point(

  aes(

    x = 0,

    y = mean_beta,

    fill = response

  ),

  data = trait_summary,

  inherit.aes = FALSE,

  shape = 21,

  colour = "black",

  stroke = 0.35,

  size = 3.2

) +

############################################################
## Colours
############################################################

scale_colour_manual(

  values = response_colours

) +

scale_fill_manual(

  values = response_colours

) +

############################################################
## Axes
############################################################

coord_cartesian(

  xlim = c(-0.18, 0.18),

  clip = "off"

) +

labs(

  x = NULL,

  y = expression(beta[CV] ~ "posterior draws")

) +

publication_theme +

theme(

  legend.position = "none",

  axis.text.x = element_blank(),

  axis.ticks.x = element_blank()

)

  (slope_plot | posterior_plot) +
    plot_layout(widths = c(3.5, 1))
}

make_trait_panel <- function(focal_trait) {
  range_pair <- make_response_pair(
    focal_trait, "Range size",
    show_x_title = focal_trait == tail(trait_order, 1)
  )
  esi_pair <- make_response_pair(
    focal_trait, "ESI",
    show_x_title = focal_trait == tail(trait_order, 1)
  )

  wrap_plots(esi_pair, range_pair, ncol = 2)
}

trait_panels <- wrap_plots(
  map(trait_order, make_trait_panel),
  ncol = 1
)

response_headers <- wrap_plots(
  ggplot() + theme_void() + labs(title = "Ecological specialisation") + publication_theme,
  ggplot() + theme_void() + labs(title = "Range size") + publication_theme,
  ncol = 2
)

main_panels <- response_headers / trait_panels +
  plot_layout(heights = c(0.035, 1))

global_y_label <- wrap_elements(
  full = grid::textGrob(
    "log(Standardised response variable)", rot = 90,
    gp = grid::gpar(fontsize = 16, fontface = "bold")
  )
)

final_plot <- wrap_plots(
  global_y_label,
  main_panels,
  ncol = 2,
  widths = c(0.035, 1)
)

final_plot

ggsave(
  filename = "../results/stan/CV_beta_slope_panels_Main.pdf",
  plot = final_plot,
  width = 12,
  height = 10,
  units = "in",
  dpi = 300
)

##Save the results object for later use
saveRDS(
  results,
  file = "../results/stan/CV_trait_results_Main.rds"
)

#Correlation between bulbul$ESI and bulbul$Range.Size (IMPORTANT)
cor.test(
  bulbul$ESI,
  bulbul$Range.Size,
  method = "pearson",
  use = "complete.obs"
)

# Random check to see which row of Iole charlottae has a $Lon values higher than 108,
bulbul %>%
  filter(
    AviList == "Rubigula squamata",
    Lat > 15
  ) %>%
  dplyr::select(
    AviList,
    Specimen.number,
    Lat,
    Lon,
    Data.type
  )

  #################Estimating phylogenetic signal as variance explained by phylogeny (lambda) for each trait using brms
  phylogenetic_variance <- purrr::map_dfr(results, function(res){

  draws <- res$fit$draws(
    variables = "sigma_logCV",
    format = "draws_df"
  )

  tibble(
    Trait = res$short_name,
    Phylogenetic_SD = mean(draws$sigma_logCV),
    SD_lower = quantile(draws$sigma_logCV, 0.025),
    SD_upper = quantile(draws$sigma_logCV, 0.975),

    Phylogenetic_variance = mean(draws$sigma_logCV^2),
    Variance_lower = quantile(draws$sigma_logCV^2, 0.025),
    Variance_upper = quantile(draws$sigma_logCV^2, 0.975)
  )

})

phylogenetic_variance
write.csv(
  phylogenetic_variance,
  "../results/stan/phylogenetic_variance_estimates_Main.csv",
  row.names = FALSE
)


################Calculating the proportion of posterior draws of beta_CV_range and beta_CV_ESI that are lying on one side of zero (i.e. the proportion of posterior draws that are positive or negative) for each trait
proportion_positive <- purrr::map_dfr(results, function(res){

  draws <- res$fit$draws(
    variables = c("beta_CV_range", "beta_CV_ESI"),
    format = "draws_df"
  )

  tibble(
    Trait = res$short_name,
    Proportion_positive_beta_CV_range = mean(draws$beta_CV_range > 0),
    Proportion_positive_beta_CV_ESI = mean(draws$beta_CV_ESI > 0)
  )

})

proportion_positive


####NOT NEEDED ANYMORE
proportion_negative <- purrr::map_dfr(results, function(res){

  draws <- res$fit$draws(
    variables = c("beta_CV_range", "beta_CV_ESI"),
    format = "draws_df"
  )

  tibble(
    Trait = res$short_name,
    Proportion_negative_beta_CV_range = mean(draws$beta_CV_range < 0),
    Proportion_negative_beta_CV_ESI = mean(draws$beta_CV_ESI < 0)
  )

})

proportion_negative







################Need to also check whether sampling depth (as number of specimens per species per trait) is correlated to the latent CV estimates means in CV_summary

############################################################
## Sampling depth vs latent CV
############################################################

library(dplyr)
library(purrr)
library(ggplot2)

############################################################
## Number of specimens per species for each trait
############################################################

specimen_counts <-

  bulbul %>%

  group_by(AviList) %>%

  summarise(

    Wing = sum(!is.na(Wing.Length)),

    Tarsus = sum(!is.na(Tarsus.Length)),

    Beak_Depth = sum(!is.na(Beak.Depth)),

    Beak_Nares = sum(!is.na(Beak.Length_Nares)),

    .groups = "drop"

  )

############################################################
## Extract latent CV estimates
############################################################

bulbul_latent_cv <-

  purrr::map_dfr(results, function(res){

    res$CV_summary %>%

      transmute(

        Trait = res$short_name,

        Species = as.character(species),

        latent_cv = mean

      )

  })

############################################################
## Make trait names match specimen count names
############################################################

bulbul_latent_cv <-

  bulbul_latent_cv %>%

  mutate(

    Trait = recode(

      Trait,

      "Wing.Length" = "Wing",

      "Tarsus.Length" = "Tarsus",

      "Beak.Depth" = "Beak_Depth",

      "Beak.Length_Nares" = "Beak_Nares",

      ## If your short names are already these,
      ## these lines simply leave them unchanged.

      "Wing" = "Wing",

      "Tarsus" = "Tarsus",

      "Beak_Depth" = "Beak_Depth",

      "Beak_Nares" = "Beak_Nares"

    )

  )

############################################################
## Merge specimen counts
############################################################

latent_cv_with_counts <-

  bulbul_latent_cv %>%

  left_join(

    specimen_counts,

    by = c("Species" = "AviList")

  )

############################################################
## Select the correct specimen count for each trait
############################################################

plot_data <-

  latent_cv_with_counts %>%

  mutate(

    Specimen_Count = case_when(

      Trait == "Wing" ~

        Wing,

      Trait == "Tarsus" ~

        Tarsus,

      Trait == "Beak Depth" ~

        Beak_Depth,

      Trait == "Beak Length" ~

        Beak_Nares,

      TRUE ~

        NA_real_

    )

  ) %>%

  filter(

    !is.na(Specimen_Count),

    !is.na(latent_cv)

  )

############################################################
## Correlations
############################################################

correlation_results <-

  plot_data %>%

  group_by(Trait) %>%

  summarise(

    n_species = n(),

    correlation = cor(

      latent_cv,

      log(Specimen_Count),

      use = "complete.obs"

    ),

    p_value = cor.test(

      latent_cv,

      log(Specimen_Count)

    )$p.value,

    .groups = "drop"

  )

print(correlation_results)

############################################################
## Plot
############################################################

plot_specimen_cv <-

ggplot(

  plot_data,

  aes(

    x = log(Specimen_Count),

    y = latent_cv

  )

) +

geom_point(

  size = 2.5,

  alpha = 0.7

) +

geom_smooth(

  method = "lm",

  se = TRUE,

  colour = "#2C7BB6",

  linewidth = 0.9

) +

facet_wrap(

  ~Trait,

  scales = "free_x"

) +

theme_classic(base_size = 15) +

labs(

  x = "log(Number of specimens)",

  y = "Estimated latent CV",

)

print(plot_specimen_cv)

############################################################
## Save outputs
############################################################

write.csv(

  correlation_results,

  "../results/Latent_CV_sampling_depth_correlations.csv",

  row.names = FALSE

)

ggsave(

  "../results/Latent_CV_sampling_depth.pdf",

  plot_specimen_cv,

  width = 10,

  height = 8,

  dpi = 600

)
