#Start fresh
rm(list=ls())
set.seed(123454)
#libs
library(tidyverse)
library(openxlsx)
library(ape)
library(brms)
library(lme4)

#load bulbul data including range size etc.
bulbul <- read.csv("../data/Bulbul_AVONET_data_integrated_trimmed.csv")
str(bulbul)

#only keep species with > 10 specimens (rows) in the dataset to ensure more reliable CV estimates
bulbul <- bulbul %>%
  group_by(AviList) %>%
  filter(n() > 10) %>%
  ungroup()




#METHOD 1: A potential method for using CV as a predictor such that range_size ~ CV
#STEP 1:Regressing specimen-level trait data against confounding effects like Sex, Data.Source and Measurer and taking the residuals as the 'corrected' trait values.
wing <- lmer(
  Wing.Length ~ Sex + Measurer + Data.source + (1|AviList),
  data = bulbul,
  na.action = na.exclude
)

#STEP 2: Use the residuals from the model to calculate CV for each species using the adjusted SD and mean values.
resid_df <- data.frame(
  species = bulbul$AviList,
  resid = residuals(wing)
)

species_sd <- resid_df |>
  group_by(species) |>
  summarise(
    resid_sd = sd(resid, na.rm = TRUE),
    n = n()
  )

bulbul$fitted <- fitted(wing)

species_mean <- bulbul |>
  group_by(AviList) |>
  summarise(
    fitted_mean = mean(fitted, na.rm = TRUE)
  )

species_stats <- left_join(species_sd, species_mean, by = c("species" = "AviList"))

species_stats$CV <- species_stats$resid_sd / species_stats$fitted_mean
species_stats$logCV <- log(species_stats$CV)

#STEP 3: Regressing range size against CV while accounting for phylogeny
#load phylogeny
phylogeny <- read.tree("../data/Phylogeny/mcc_dated_clements_McTavish_2025.nex")

phylogeny <- phylogeny[[6]]
phylogeny$tip.label <- gsub("_", " ", phylogeny$tip.label)
phylogeny <- drop.tip(phylogeny, setdiff(phylogeny$tip.label, unique(bulbul$AviList)))
A <- ape::vcv(phylogeny, corr = TRUE)

#merge range size data with CV data
species_range <- bulbul %>%
  select(AviList, Range.Size) %>%
  distinct(AviList, .keep_all = TRUE)

range_cv <- left_join(species_stats, species_range, by = c("species" = "AviList"))

#Add standard error of CV based on sample sizes
range_cv$CV_se = range_cv$CV / sqrt(2*range_cv$n)
range_cv$logCV_se = 1 / sqrt(2 * range_cv$n)

#Now scaling SE and logCV for the model
range_cv$logCV_mean <- mean(range_cv$logCV)
range_cv$logCV_sd <- sd(range_cv$logCV)

range_cv$logCV_z <- (
  range_cv$logCV - range_cv$logCV_mean
) / range_cv$logCV_sd

range_cv$logCV_se_z <- (
  range_cv$logCV_se / range_cv$logCV_sd
)

#STEP 4: Fit a brms model with range size as the response variable, CV as the predictor and a phylogenetic random effect
model <- brm(
  log(Range.Size) ~ me(logCV_z, logCV_se_z) + (1 | gr(species, cov = A)),
  data = range_cv,
  data2 = list(A = A),
  family = gaussian(),
  chains = 4,
  iter = 4000,
  warmup = 500,
  control = list(adapt_delta = 0.99)
)

summary(model)

#PLOT THE MODEL
range_cv$predicted <- predict(model, newdata = range_cv, re_formula = NA)[, "Estimate"]
range_cv$lower <- predict(model, newdata = range_cv, re_formula = NA)[, "Q2.5"]
range_cv$upper <- predict(model, newdata = range_cv, re_formula = NA)[, "Q97.5"]    

ggplot(range_cv, aes(x = logCV_z, y = log(Range.Size))) +
  geom_point() +
  geom_errorbar(aes(ymin = log(Range.Size) - logCV_se_z, ymax = log(Range.Size) + logCV_se_z), width = 0.1) +
  geom_line(aes(y = predicted), color = "blue") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  labs(x = "Standardized Log(Coefficient of Variation)", y = "Log(Range Size)") +
  theme_minimal()

ggsave("../results/brms/Bulbul_CV_range_size_model_BeakLength.pdf", width = 10, height = 6)




#METHOD 2: Regressing trait variance as SD against range size such that SD ~ range_size (the reverse of method 1 but the easiest method)
#STEP 1: Log-transform trait values from the start
bulbul$logWing <- log(bulbul$Wing.Length)

#STEP 2: Calculate species means for logged traits to use in the SD ~ range_size model to approximate CV
species_means <- bulbul %>%
  group_by(AviList) %>%
  summarise(
    species_mean = mean(logWing, na.rm = TRUE)
  )

bulbul <- left_join(
  bulbul,
  species_means,
  by = "AviList"
)

bulbul$RangeSize_z <- scale(log(bulbul$Range.Size))[,1]
bulbul$species_mean_z <- scale(bulbul$species_mean)[,1]
bulbul$ESI_z <- scale(bulbul$ESI)[,1]

model2 <- brm(
  bf(logWing ~ Sex + Data.type + (1|gr(AviList, cov = A)),
    sigma ~ species_mean_z + RangeSize_z + ESI_z + (1|gr(AviList, cov = A))
  ),
  data = bulbul,
  data2 = list(A = A),
  family = gaussian(),
  chains = 4,
  iter = 2000,
  warmup = 500,

  control = list(adapt_delta = 0.95, max_treedepth = 15)
)

summary(model2)



##METHOD 3: Bivariate model to test the covariance between wing-length variance and range size while accounting for phylogeny

#STEP 1: Create scaled species-level predictors

bulbul$logWing <- log(bulbul$Wing.Length)

species_means <- bulbul %>%
  group_by(AviList) %>%
  summarise(
    species_mean = mean(logWing, na.rm = TRUE),
    n = n()
  )

bulbul <- left_join(
  bulbul,
  species_means,
  by = "AviList"
)

#Scale predictors
bulbul$RangeSize_z <- scale(log(bulbul$Range.Size))[,1]
bulbul$species_mean_z <- scale(bulbul$species_mean)[,1]
bulbul$log_n_z <- scale(log(bulbul$n))[,1]



#STEP 2: Trait model
#This estimates:
#- mean wing length
#- species-specific residual SDs (sigma)
#- while accounting for phylogeny and nuisance effects

trait_bf <- bf(

  logWing ~
    Sex +
    Data.source +
    (1|Measurer) +
    (1|p|gr(AviList, cov = A)),

  sigma ~
    species_mean_z +
    log_n_z +
    (1|p|gr(AviList, cov = A))

)



#STEP 3: Species-level range-size model
#This estimates species-level range size while sharing correlated
#phylogenetic random effects with the sigma structure above.

range_bf <- bf(

  RangeSize_z ~
    1 +
    (1|p|gr(AviList, cov = A))

)



#STEP 4: Fit the bivariate model
#The shared "p" term allows correlated species effects across:
#- range size
#- latent trait variability

model3 <- brm(

  trait_bf +
    range_bf +
    set_rescor(FALSE),

  data = bulbul,
  data2 = list(A = A),

  family = gaussian(),

  chains = 4,
  iter = 4000,
  warmup = 1000,

  control = list(
    adapt_delta = 0.95,
    max_treedepth = 15
  )

)

summary(model3)
