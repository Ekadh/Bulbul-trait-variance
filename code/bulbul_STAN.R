rm(list = ls())
set.seed(123454)
# libs
library(tidyverse)
library(ape)
library(cmdstanr)
library(Matrix)


bulbul <- read.csv("../data/Bulbul_plotting.csv")

# Use raw wing length (mm)
bulbul$Wing_raw <- bulbul$Wing.Length

bulbul <- bulbul %>%
  filter(
    !is.na(Wing_raw),
    !is.na(Sex),
    !is.na(Data.source),
    !is.na(Range.Size)
  )

bulbul$species  <- factor(bulbul$AviList)
bulbul$sex      <- factor(bulbul$Sex)
bulbul$source   <- factor(bulbul$Data.source)

#PHYLOGENY

phylogeny <- read.tree(
  "../data/Phylogeny/mcc_dated_clements_McTavish_2025.nex"
)

phylogeny <- phylogeny[[6]]

phylogeny$tip.label <-
  gsub("_", " ", phylogeny$tip.label)

species_levels <- levels(bulbul$species)

phylogeny <- drop.tip(
  phylogeny,
  setdiff(
    phylogeny$tip.label,
    species_levels
  )
)

A <- ape::vcv(phylogeny, corr = TRUE)

A <- A[
  species_levels,
  species_levels
]

#===========================================================
# CHOLESKY DECOMPOSITION
#===========================================================

#It seems that Stan requires the Cholesky transformation of the phylogenetic covariance matrix, but I'm not sure if it has any unwanted consequences

L_A <- t(
  chol(A)
)

#Range size data frame

species_df <- data.frame(
  species = levels(bulbul$species)
) %>%
  left_join(
    bulbul %>%
      group_by(species) %>%
      summarise(
        log_range = first(log(Range.Size)),
        .groups = "drop"
      ),
    by = "species"
  )

head(species_df)


species_id  <- as.integer(bulbul$species)
sex_id      <- as.integer(bulbul$sex)
source_id   <- as.integer(bulbul$source)

#All data needed for Stan 

stan_data_A <- list(

  N = nrow(bulbul),
  S = nlevels(bulbul$species),

  wing = bulbul$Wing_raw,

  log_range = species_df$log_range,

  species = species_id,

  n_sex = nlevels(bulbul$sex),
  sex = sex_id,

  n_source = nlevels(bulbul$source),
  source = source_id,

  L_A = L_A
)
#Final checks
length(stan_data_A$log_range)
sum(is.na(stan_data_A$log_range))

#Compiling Stan

mod <- cmdstan_model(
  "latent_variance_rangesize.stan"
)

fit <- mod$sample(

  data = stan_data_A,

  chains = 4,
  parallel_chains = 4,

  iter_warmup = 500,
  #Just trialing a shorter sampling phase to test the model, but will increase for final run, ESS and Rhat seem ok already
  iter_sampling = 500,

  adapt_delta = 0.95,

  refresh = 100
)

#Post-hoc checks

fit$summary(
  variables = c(
    "beta_CV",
    "sigma_logCV",
    "sigma_R"
  )
)

fit$summary("mean_CV")

CV_summary <- fit$summary("logCV")

CV_summary$species <- species_levels

raw_cv <- bulbul |>
  group_by(species) |>
  summarise(
    CV = sd(Wing_raw)/mean(Wing_raw)
  )

compare <- left_join(
  CV_summary[,c("species","mean")],
  raw_cv,
  by="species"
)

cor.test(
  exp(compare$mean),
  compare$CV
)

stan_cv <- data.frame(
  species = species_levels,
  CV_stan = exp(CV_summary$mean)
)

raw_cv <- bulbul |>
  group_by(species) |>
  summarise(
    CV_raw = sd(Wing_raw) /
      mean(Wing_raw)
  )

cv_compare <- left_join(
  stan_cv,
  raw_cv,
  by = "species"
)

ggplot(
  cv_compare,
  aes(
    x = CV_raw,
    y = CV_stan
  )
) +
  geom_point() +
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  theme_classic()


#Looking at mean CVs
CV_summary <- fit$summary("CV_species")

CV_summary$species <- species_levels

CV_summary |>
  arrange(mean)

summary(CV_summary$mean)

CV_summary |>
  arrange(mean) |>
  select(species, mean, q5, q95) |>
  head(20)

##PLOT
plot_df <- left_join(
  CV_summary |> select(species, CV = mean),
  species_df |> select(species, log_range),
  by = "species"
)

plot_df$Range_1000km2 <- exp(plot_df$log_range) / 1000

p <-ggplot(
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
    x = expression(paste("Range size (", 1000, " km"^2, ")")),
    y = "Estimated species CV"
  ) +
  theme_classic()
  
p

ggsave(
  "../results/stan/Bulbul_CulmenCV_range_size_model.pdf",
  width = 8,
  height = 6
)

q <-ggplot(
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
    x = expression(paste("Range size (", 1000, " km"^2, ")")),
    y = "Estimated species CV"
  ) +
  theme_classic()
  
q

#Combine p and q into one plot
library(patchwork)
combined_plot <- p + q
ggsave(
  "../results/stan/Bulbul_CulmenWingCV_range_size_model_combined.pdf",
  plot = combined_plot,
  width = 16,
  height = 10
)

#Save fitted model

saveRDS(
  fit,
  "../results/wing_latent_variance_fit.rds"
)

#Posterior distribution of beta_CV
posterior <- fit$draws("beta_CV")
ggplot(
  data.frame(beta_CV = as.vector(posterior)),
  aes(x = beta_CV)
) +
  geom_density(fill = "lightblue") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "Posterior distribution of beta_CV",
    x = "beta_CV",
    y = "Density"
  ) +
  theme_classic()

