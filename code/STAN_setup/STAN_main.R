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
  "../data/Bulbul_plotting.csv"
)

bulbul_density <- read.csv("../data/Density/Bulbul_density_data_Callaghan.csv")
bulbul_ESI <- read.csv("../data/BIRDBASE/Bulbul_Birdbase_data.csv")

# Add Density and ESI to the data frame, matching by species
bulbul <- bulbul %>%
  left_join(bulbul_density %>% select(Scientific.name, Density), by = c("AviList" = "Scientific.name")) %>%
  left_join(bulbul_ESI %>% select(AviList.v1.2025, ESI), by = c("AviList" = "AviList.v1.2025"))

#===========================================================
# TRAITS TO ANALYSE
#===========================================================

traits <- tibble(

  trait = c(

    "Wing.Length",
    "Beak.Length_Culmen",
    "Tail.Length",
    "Tarsus.Length"
    # "Beak.Depth",
    # "Beak.Width",
    # "Beak.Length_Nares"

  ),

  short_name = c(

    "Wing",
    "Culmen",
    "Tail",
    "Tarsus"

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
            first(log(Density)),

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

      iter_warmup = 500,

      iter_sampling = 500,

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

        "_fit.rds"

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
summary(results[[1]]$CV_summary$mean)
sd(results[[4]]$CV_summary$mean)
IQR(results[[4]]$CV_summary$mean)

####To check if the most variable species are the same across traits, we can check correlations of the latent CVs across traits
cor(results$Tarsus$CV_summary$mean,
    results$Wing$CV_summary$mean)

cor(results$Tarsus$CV_summary$mean,
    results$Tail$CV_summary$mean)

cor(results$Culmen$CV_summary$mean,
    results$Tarsus$CV_summary$mean)

