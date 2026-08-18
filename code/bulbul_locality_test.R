############################################################
## Pairwise distance-decay analysis
## Part 1
## Configuration + master geographic pairwise table
############################################################

rm(list = ls())

set.seed(12345)

############################################################
## Packages
############################################################

suppressPackageStartupMessages({

library(tidyverse)

library(geosphere)

library(broom)

library(MASS)

})

############################################################
## Paths
############################################################

data_path <-
  "../data/Bulbul_AVONET_data_integrated_trimmed.csv"

output_dir <-
  "../results/Distance_decay_2spec"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

############################################################
## Settings
############################################################

maximum_specimens <- 100

minimum_specimens <- 2

traits <- c(

"Beak.Depth",

"Beak.Length_Nares",

"Tarsus.Length",

"Wing.Length"

)

############################################################
## Read data
############################################################

bulbuls <-

read.csv(

data_path,

check.names = FALSE,

stringsAsFactors = FALSE

)

############################################################
## Keep only specimens with coordinates
############################################################

bulbuls <-

bulbuls %>%

filter(

!is.na(Lat),

!is.na(Lon)

)

############################################################
## Create globally unique specimen ID
############################################################

bulbuls <-

bulbuls %>%

mutate(

specimen_id = paste(

AviList,

Source,

Specimen.number,

sep="_"

)

)

############################################################
## Build master pairwise table
############################################################

pairwise_species <- function(dat){

  ##########################################################
  ## Skip poorly sampled species
  ##########################################################

  if(nrow(dat) < minimum_specimens)
    return(NULL)

  ##########################################################
  ## Random subsampling
  ##########################################################

  if(nrow(dat) > maximum_specimens){

    dat <-

      dat %>%

      slice_sample(
        n = maximum_specimens
      )

  }

  ##########################################################
  ## All pairwise combinations
  ##########################################################

  pairs <- combn(
    nrow(dat),
    2
  )

  ##########################################################
  ## Allocate output list
  ##########################################################

  out <-

    vector(

      "list",

      ncol(pairs)

    )

  ##########################################################
  ## Loop
  ##########################################################

  for(i in seq_len(ncol(pairs))){

    id1 <- pairs[1,i]

    id2 <- pairs[2,i]

    ########################################################
    ## Geographic distance
    ########################################################

    distance_km <-

      distHaversine(

        c(

          dat$Lon[id1],

          dat$Lat[id1]

        ),

        c(

          dat$Lon[id2],

          dat$Lat[id2]

        )

      ) / 1000

    ########################################################
    ## Source pair
    ########################################################

    SourcePair <-

      case_when(

      dat$Data.type[id1]== 2 &
      dat$Data.type[id2]== 2

      ~ "Field-Field",

      dat$Data.type[id1]== 1 &
      dat$Data.type[id2]== 1

      ~ "Museum-Museum",

      TRUE

      ~ "Museum-Field"

      )

    ########################################################
    ## Store
    ########################################################

    out[[i]] <-

      tibble(

        Species = dat$AviList[id1],

        specimen_1 = dat$specimen_id[id1],

        specimen_2 = dat$specimen_id[id2],

        source_1 = dat$Data.type[id1],

        source_2 = dat$Data.type[id2],

        SourcePair = SourcePair,

        Lat1 = dat$Lat[id1],

        Lon1 = dat$Lon[id1],

        Lat2 = dat$Lat[id2],

        Lon2 = dat$Lon[id2],

        distance_km = distance_km,

        log_distance = log(distance_km + 1)

      )

  }

  bind_rows(out)

}

############################################################
## Build pairwise geographic table
############################################################

pairwise_list <- list()

counter <- 1

for(sp in sort(unique(bulbuls$AviList))){

cat(sp,"\n")

dat <-

bulbuls %>%

filter(

AviList==sp

)

tmp <-

pairwise_species(dat)

if(is.null(tmp))
next

pairwise_list[[counter]] <- tmp

counter <- counter + 1

}

pairwise_master <-

bind_rows(pairwise_list)

write.csv(

pairwise_master,

file.path(

output_dir,

"Pairwise_master_geographic_table.csv"

),

row.names=FALSE

)

############################################################
## Sanity checks
############################################################

summary(pairwise_master$distance_km)

length(unique(pairwise_master$Species))

nrow(pairwise_master)

table(pairwise_master$SourcePair)

############################################################
## PART 2
## Build trait-specific pairwise datasets
############################################################

trait_pairwise_list <- list()

############################################################
## Loop over traits
############################################################

for(current_trait in traits){

  cat("\n")
  cat("=============================================\n")
  cat("Processing:", current_trait, "\n")
  cat("=============================================\n")

  ##########################################################
  ## Build lookup vector
  ##########################################################

  lookup <-

    bulbuls %>%

    dplyr::select(
      specimen_id,
      all_of(current_trait)
    ) %>%

    filter(!is.na(.data[[current_trait]]))

  lookup_vector <- lookup[[current_trait]]
  names(lookup_vector) <- lookup$specimen_id

  ##########################################################
  ## Copy master table
  ##########################################################

  trait_data <- pairwise_master

  ##########################################################
  ## Retrieve measurements
  ##########################################################

  trait_data$value_1 <-

    lookup_vector[
      trait_data$specimen_1
    ]

  trait_data$value_2 <-

    lookup_vector[
      trait_data$specimen_2
    ]

  ##########################################################
  ## Calculate pairwise differences
  ##########################################################

  trait_data <-

    trait_data %>%

    mutate(

      Trait = current_trait,

      trait_difference =
        abs(value_1 - value_2),

      proportional_difference =
        abs(value_1 - value_2) /
        ((value_1 + value_2)/2)

    ) %>%

    filter(

      !is.na(trait_difference),

      is.finite(proportional_difference)

    )

  ##########################################################
  ## Save
  ##########################################################

  trait_pairwise_list[[current_trait]] <- trait_data

}

############################################################
## Combine all traits
############################################################

pairwise_traits <-

  bind_rows(
    trait_pairwise_list
  )

############################################################
## Save pairwise trait table
############################################################

write.csv(

  pairwise_traits,

  file.path(

    output_dir,

    "Pairwise_trait_dataset.csv"

  ),

  row.names = FALSE

)

############################################################
## Geographic distance distribution
############################################################

p1 <-

ggplot(

pairwise_master,

aes(log_distance)

)+

geom_histogram(

bins = 50,

fill = "grey70",

colour = "black"

)+

theme_bw(base_size = 12)+

labs(

x = "log(distance + 1 km)",

y = "Pairwise comparisons"

)

ggsave(

file.path(

output_dir,

"Geographic_distance_distribution.pdf"

),

plot = p1,

width = 7,

height = 5

)

############################################################
## Trait difference distributions
############################################################

p2 <-

ggplot(

pairwise_traits,

aes(trait_difference)

)+

geom_histogram(

bins = 50,

fill = "grey70",

colour = "black"

)+

facet_wrap(

~Trait,

scales = "free",

labeller = labeller(
  Trait = c(
    "Beak.Depth" = "Beak Depth",
    "Beak.Length_Nares" = "Beak Length (Nares)",
    "Tarsus.Length" = "Tarsus Length",
    "Wing.Length" = "Wing Length"
  )
)

)+

theme_bw(base_size = 12)

ggsave(

file.path(

output_dir,

"Trait_difference_distributions.pdf"

),

plot = p2,

width = 12,

height = 8
)

############################################################
## Proportional differences
############################################################

p3 <-

ggplot(

pairwise_traits,

aes(proportional_difference)

)+

geom_histogram(

bins = 50,

fill = "grey70",

colour = "black"

)+

facet_wrap(

~Trait,

scales = "free"

)+

theme_bw(base_size = 12)

ggsave(

file.path(

output_dir,

"Proportional_difference_distributions.pdf"

),

plot = p3,

width = 12,

height = 8

)

############################################################
## Sanity checks
############################################################

pairwise_traits %>%

group_by(Trait) %>%

summarise(

n_pairs = n(),

mean_difference =
  mean(trait_difference),

median_difference =
  median(trait_difference),

max_difference =
  max(trait_difference)

) %>%

print(n = Inf)

############################################################
## PART 3
## Species-specific regressions
############################################################

############################################################
## Containers
############################################################

species_results <- list()

model_list <- list()

counter <- 1

############################################################
## Loop over traits
############################################################

for(current_trait in traits){

  cat("\n")
  cat("========================================\n")
  cat("Regression:", current_trait, "\n")
  cat("========================================\n")

  trait_dat <-

    pairwise_traits %>%

    filter(
      Trait == current_trait
    )

  ##########################################################
  ## Loop over species
  ##########################################################

  for(sp in unique(trait_dat$Species)){

    dat <-

      trait_dat %>%

      filter(
        Species == sp
      )

    ########################################################
    ## Skip species with insufficient variation
    ########################################################

    if(nrow(dat) < 2)
      next

    ## No geographic variation

    if(max(dat$distance_km, na.rm = TRUE) == 0){

      cat("Skipping", sp, "- all specimens at same locality\n")

      next

    }

    ## No variation in predictor

    if(sd(dat$log_distance, na.rm = TRUE) == 0){

      cat("Skipping", sp, "- no variation in geographic distance\n")

      next

    }

    ## No variation in response

    if(sd(dat$trait_difference, na.rm = TRUE) == 0){

      cat("Skipping", sp, "- no variation in trait\n")

      next

    }

    ########################################################
    ## Linear model
    ########################################################

    lm_fit <-

      lm(

        trait_difference ~

          log_distance,

        data = dat

      )

    ########################################################
    ## Robust regression
    ########################################################

    rlm_fit <-

      tryCatch(

        MASS::rlm(

          trait_difference ~

            log_distance,

          data = dat

        ),

        error = function(e) NULL

      )

    ########################################################
    ## Store models
    ########################################################

    model_list[[paste(
      current_trait,
      sp,
      sep = "_"
    )]] <- lm_fit

    ########################################################
    ## Tidy output
    ########################################################

    lm_tidy <-

      broom::tidy(lm_fit)

    if(is.null(rlm_fit)){

      slope_rlm <- NA_real_

    }else{

      rlm_tidy <-

        broom::tidy(rlm_fit)

      slope_rlm <-

        rlm_tidy %>%

        filter(term == "log_distance") %>%

        pull(estimate)

    }

    ########################################################
    ## Extract LM slope
    ########################################################

    lm_row <-

      lm_tidy %>%

      filter(
        term == "log_distance"
      )

    ########################################################
    ## Store summary
    ########################################################

    species_results[[counter]] <-

      tibble(

        Species = sp,

        Trait = current_trait,

        n_pairs = nrow(dat),

        mean_distance =
          mean(dat$distance_km),

        max_distance =
          max(dat$distance_km),

        slope_lm =
          lm_row$estimate,

        SE =
          lm_row$std.error,

        t =
          lm_row$statistic,

        p =
          lm_row$p.value,

        R2 =
          summary(lm_fit)$r.squared,

        Adj_R2 =
          summary(lm_fit)$adj.r.squared,

        slope_rlm =
          slope_rlm

      )

    counter <- counter + 1

  }

}

############################################################
## Combine regression summaries
############################################################

species_slopes <-

  bind_rows(
    species_results
  )

############################################################
## Save regression summaries
############################################################

write.csv(

  species_slopes,

  file.path(

    output_dir,

    "Species_distance_decay_results.csv"

  ),

  row.names = FALSE

)

############################################################
## Histogram of slopes
############################################################

p4 <-

ggplot(

species_slopes,

aes(slope_lm)

)+

geom_density(

fill = "gray70",

colour = "black"

)+

facet_wrap(

~Trait,

scales = "free",

labeller = labeller(
  Trait = c(
    "Beak.Depth" = "Beak Depth",
    "Beak.Length_Nares" = "Beak Length (Nares)",
    "Tarsus.Length" = "Tarsus Length",
    "Wing.Length" = "Wing Length"
  )

)

)+

labs(x = "Distance-decay slope", y = "Density") +

theme_bw(base_size = 12) +

theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      strip.text = element_text(size = 20),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 15)) 

ggsave(

file.path(

output_dir,

"Slope_distribution.pdf"

),

plot = p4,

width = 12,

height = 8

)

p5 <-

ggplot(

species_slopes,

aes(R2)

)+

geom_density(

fill = "gray70",

colour = "black"

)+

facet_wrap(

~Trait,

scales = "free",

labeller = labeller(
  Trait = c(
    "Beak.Depth" = "Beak Depth",
    "Beak.Length_Nares" = "Beak Length (Nares)",
    "Tarsus.Length" = "Tarsus Length",
    "Wing.Length" = "Wing Length"
  )
)

)+

labs(x = "Distance-decay R² values", y = "Density") +

facet_wrap(

~Trait,

scales = "free"

)+

theme_bw(base_size = 12)

ggsave(

file.path(

output_dir,

"R2_distribution.pdf"

),

plot = p5,

width = 12,

height = 8

)

############################################################
## Number of significant species
############################################################

species_slopes %>%

group_by(Trait) %>%

summarise(

  n_species = n(),

  significant = sum(p < 0.05),

  positive = sum(slope_lm > 0),

  negative = sum(slope_lm < 0),

  mean_slope = mean(slope_lm)

) %>%

print(n = Inf)

############################################################
## PART 4
## Merge with range size
############################################################

############################################################
## Read species-level dataset
############################################################

## Replace with your own file

species_data <- read.csv(

  "../data/Bulbul_rangesize.csv",

  check.names = FALSE,

  stringsAsFactors = FALSE

)

############################################################
## Merge
############################################################

species_slopes <-

species_slopes %>%

left_join(

species_data,

by = c("Species" = "AviList")

)

############################################################
## Check merge
############################################################

sum(is.na(species_slopes$Range.Size))

############################################################
## Range size regressions
############################################################

range_models <- list()

counter <- 1

for(current_trait in traits){

  cat(current_trait,"\n")

  dat <-

    species_slopes %>%

    filter(

      Trait == current_trait,

      !is.na(Range.Size)

    )

  if(nrow(dat) < 10)
    next

  fit <-

    lm(

      slope_lm ~

        log10(Range.Size),

      data = dat

    )

  tidy_fit <-

    broom::tidy(fit)

  slope_row <-

    tidy_fit %>%

    filter(

      term == "log10(Range.Size)"

    )

  range_models[[counter]] <-

    tibble(

      Trait = current_trait,

      Estimate = slope_row$estimate,

      SE = slope_row$std.error,

      t = slope_row$statistic,

      P = slope_row$p.value,

      R2 = summary(fit)$r.squared,

      Adj_R2 = summary(fit)$adj.r.squared

    )

  counter <- counter + 1

}

range_model_summary <-

bind_rows(range_models)

############################################################
## Save results
############################################################

write.csv(

range_model_summary,

file.path(

output_dir,

"Range_size_regressions.csv"

),

row.names = FALSE

)

############################################################
## Scatterplots
############################################################

## Combined faceted plot for all traits vs range size

dat_all <-
  species_slopes %>%
  filter(
    Trait %in% traits,
    !is.na(Range.Size)
  )

p_all <-
  ggplot(dat_all, aes(log10(Range.Size), slope_lm)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, colour = "red") +
  facet_wrap(~ Trait, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(x = "log10(range size)", y = "Distance-decay slope", title = "Trait regressions vs Range size")

ggsave(
  file.path(output_dir, "AllTraits_RangeSize.pdf"),
  plot = p_all,
  width = 10,
  height = 8
)

############################################################
## Overall summary
############################################################

species_slopes %>%

group_by(Trait) %>%

summarise(

  Mean_slope = mean(slope_lm),

  SD_slope = sd(slope_lm),

  Median_slope = median(slope_lm),

  Mean_R2 = mean(R2),

  Significant = sum(p < 0.05),

  Total_species = n(),
  
  mean_p = mean(p)

) %>%

write.csv(

file.path(

output_dir,

"Overall_summary.csv"

),

row.names = FALSE

)

############################################################
## Finished
############################################################

cat("\n")

cat("=======================================\n")

cat("Distance-decay analysis complete.\n")

cat("Outputs written to:\n")

cat(output_dir,"\n")

cat("=======================================\n")










############################################################
## Species-specific geographic maps
############################################################

## Aggregate specimen counts by unique coordinates per species

species_maps_data <- bulbuls %>%
  group_by(AviList, Lat, Lon) %>%
  summarise(
    n_specimens = n(),
    .groups = "drop"
  )

## Create map for each species

species_list <- sort(unique(species_maps_data$AviList))

pdf(
  "../data/Species_geographic_maps.pdf",
  width = 10,
  height = 8
)

# load country outlines once
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

for (sp in species_list) {

  cat("Mapping:", sp, "\n")

  sp_data <- species_maps_data %>%
    filter(AviList == sp)

  # remove rows with missing coordinates to avoid st_as_sf error
  sp_data <- sp_data %>%
    filter(!is.na(Lat) & !is.na(Lon))

  if (nrow(sp_data) == 0) {
    warning(sprintf("No valid coordinates for species '%s'; skipping map.", sp))
    next
  }

  ## Scale circle size by number of specimens
  max_specimens <- max(sp_data$n_specimens)

  ## convert points to sf so coord_sf can correctly constrain plotting extent
  sp_sf <- sf::st_as_sf(sp_data, coords = c("Lon", "Lat"), crs = 4326)

  p_map <- ggplot() +
    # country outlines
    geom_sf(data = world, inherit.aes = FALSE, fill = NA, colour = "gray40") +
    geom_sf(
      data = sp_sf,
      aes(size = n_specimens),
      alpha = 0.6,
      colour = "darkblue",
      inherit.aes = FALSE
    ) +

    scale_size_continuous(
      name = "Specimen count",
      range = c(2, 15),
      # force integer labels in legend
      labels = function(x) as.character(round(x, 0))
    ) +

    theme_bw(base_size = 12) +

    labs(
      title = sp,
      x = "Longitude",
      y = "Latitude"
    ) +

    ggspatial::annotation_scale(
      location = "bl",
      width_hint = 0.25
    ) +

    # use a different north arrow style
    ggspatial::annotation_north_arrow(
      location = "tr",
      which_north = "true",
      style = ggspatial::north_arrow_fancy_orienteering()
    ) +

    # compute bbox from sf object and pad by 1 degree
    coord_sf(xlim = c(sf::st_bbox(sp_sf)$xmin - 1, sf::st_bbox(sp_sf)$xmax + 1),
             ylim = c(sf::st_bbox(sp_sf)$ymin - 1, sf::st_bbox(sp_sf)$ymax + 1),
             expand = FALSE)

  print(p_map)

}

dev.off()


