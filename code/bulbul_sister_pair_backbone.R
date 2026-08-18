############################################################
## Local population variance analysis
## Bulbul sister species
############################################################

rm(list = ls())
set.seed(12345)

suppressPackageStartupMessages({

  library(tidyverse)

  library(sf)
  library(dbscan)

  library(lme4)
  library(lmerTest)

  library(broom.mixed)

  library(ggplot2)

})

############################################################
## CONFIGURATION
############################################################

## Input/output

data_path <- "../data/Bulbul_AVONET_data_integrated_trimmed.csv"

output_dir <- "../results/Sister_pair_analysis_5_10km"

if(!dir.exists(output_dir)){
  dir.create(output_dir,
             recursive = TRUE)
}

############################################################
## Dataset columns
############################################################

species_col <- "AviList"

lat_col <- "Lat"
lon_col <- "Lon"

specimen_col <- "Specimen.number"

locality_col <- "Locality"


############################################################
## Coordinate system
############################################################

## Equal-area projection in metres

projection_crs <- 6933

############################################################
## Sister pairs
############################################################

sister_pairs <- tribble(

  ~pair_id, ~restricted_species, ~widespread_species,

  "pair_01", "Pycnonotus taivanus", "Pycnonotus sinensis",

  "pair_02", "Chlorocichla falkensteini", "Chlorocichla flaviventris",

  "pair_03", "Arizelocichla tephrolaema", "Arizelocichla masukuensis",

  "pair_04", "Criniger olivaceus", "Criniger calurus",

  "pair_05", "Phyllastrephus viridiceps", "Phyllastrephus albigularis",

  "pair_06", "Phyllastrephus fischeri", "Phyllastrephus cabanisi",

  "pair_07", "Alophoixus frater", "Alophoixus tephrogenys",

  "pair_08", "Hypsipetes mindorensis", "Hypsipetes philippinus",

  "pair_09", "Hypsipetes olivaceus", "Hypsipetes madagascariensis",

  "pair_10", "Hypsipetes siquijorensis", "Hypsipetes philippinus",

  "pair_11", "Atimastillas flavicollis", "Atimastillas flavigula",

  "pair_12", "Phyllastrephus albigula", "Phyllastrephus debilis",

  "pair_13", "Hypsipetes thompsoni", "Hypsipetes leucocephalus",

  "pair_14", "Microtarsus fuscoflavescens", "Microtarsus melanocephalos",

  "pair_15", "Phyllastrephus poliocephalus", "Phyllastrephus albigularis",

  "pair_16", "Rubigula gularis", "Rubigula flaviventris",

  "pair_17", "Rubigula montis", "Rubigula dispar",

  "pair_18", "Pycnonotus davisoni", "Pycnonotus finlaysoni",

  "pair_19", "Pycnonotus capensis", "Pycnonotus barbatus",

  "pair_20", "Bleda notatus", "Bleda canicapillus",

  "pair_21", "Rubigula squamata", "Rubigula erythropthalmos",

  "pair_22", "Pycnonotus snouckaerti", "Pycnonotus bimaculatus",

  "pair_23", "Pycnonotus leucogenys", "Pycnonotus leucotis",

  "pair_24", "Pycnonotus nigricans", "Pycnonotus barbatus",

  "pair_25", "Alophoixus flaveolus", "Alophoixus pallidus",

  "pair_26", "Hemixos cinereus", "Hemixos flavala"

)

## Mark pair rows containing a species that occurs in more than one pair.
## This is retained in all outputs to document non-independent pairs.
duplicate_species <- sister_pairs %>%
  pivot_longer(
    cols = c(restricted_species, widespread_species),
    values_to = "AviList"
  ) %>%
  count(AviList, name = "n_pairs") %>%
  filter(n_pairs > 1) %>%
  pull(AviList)

sister_pairs <- sister_pairs %>%
  mutate(
    Dupe = as.integer(
      restricted_species %in% duplicate_species |
        widespread_species %in% duplicate_species
    )
  )

write.csv(
  sister_pairs,
  file.path(output_dir, "Sister_pairs_with_Dupe.csv"),
  row.names = FALSE
)

############################################################
## ANALYSIS SETTINGS
############################################################

## Radius (km) used by DBSCAN: choose between 10, 25 and 50 km

cluster_radius_km <- 10

## Minimum number of specimens required for a cluster
## to contribute to the analysis: choose between 5 and 10 specimens depending on radius size

minimum_cluster_size <- 5

## Minimum number of neighbours for DBSCAN

dbscan_minPts <- 2

############################################################
## MODEL SETTINGS
############################################################

response_transformation <- "logCV"

############################################################
## TRAITS
############################################################

trait_cols <- c(

  "Beak.Length_Nares",
  "Beak.Depth",
  "Tarsus.Length",
  "Wing.Length"

)

############################################################
## OUTPUTS
############################################################

save_cluster_maps <- TRUE

save_cluster_summary <- TRUE

save_model_results <- TRUE

save_diagnostic_plots <- TRUE

############################################################
## Read data
############################################################

bulbuls <- read.csv(
  data_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

############################################################
## Keep only focal species
############################################################

species_keep <- unique(c(
  sister_pairs$restricted_species,
  sister_pairs$widespread_species
))

############################################################
## Create lookup table
############################################################

pair_lookup <-

  sister_pairs %>%

  pivot_longer(
    cols = c(restricted_species,
             widespread_species),
    names_to = "status",
    values_to = "AviList"
  ) %>%

  mutate(

    range_status = ifelse(
      status == "restricted_species",
      "Restricted",
      "Widespread"
    )

  ) %>%

  select(
    AviList,
    pair_id,
    range_status,
    Dupe
  )

############################################################
## Prepare specimen table
############################################################

specimens <-

  bulbuls %>%

  filter(
    AviList %in% species_keep,
    !is.na(Lat),
    !is.na(Lon)
  ) %>%

  left_join(
    pair_lookup,
    by="AviList",
    relationship = "many-to-many"
  )

############################################################
## Convert to sf
############################################################

specimens_sf <-

  st_as_sf(

    specimens,

    coords=c("Lon","Lat"),

    crs=4326

  ) %>%

  st_transform(projection_crs)

############################################################
## Run DBSCAN separately for each species within each pair.
## This prevents a reused species from having its specimen records duplicated
## before cluster sizes and CVs are calculated.
############################################################

cluster_list <- list()

for(pair in unique(specimens_sf$pair_id)){

  for(sp in unique(specimens_sf$AviList[specimens_sf$pair_id == pair])){

  cat(pair, ": ", sp, "\n", sep = "")

  dat <-

    specimens_sf %>%

    filter(pair_id == pair, AviList == sp)

  coords <-

    st_coordinates(dat)

  db <-

    dbscan(

      coords,

      eps = cluster_radius_km*1000,

      minPts = dbscan_minPts

    )

  dat$cluster <-

    db$cluster

  cluster_list[[paste(pair, sp, sep = "__")]] <- dat

  }

}

cluster_data <-

  bind_rows(cluster_list)

#To account for the cluster == 0 for singleton specimens
cluster_data <-

cluster_data %>%

mutate(

cluster = ifelse(
cluster==0,
NA,
cluster
),

#For debugging, create a unique cluster ID
cluster_id = paste(
    AviList,
    cluster,
    sep = "_"
)


)

cluster_sizes <-

cluster_data %>%

st_drop_geometry() %>%

filter(!is.na(cluster)) %>%

count(

pair_id,

AviList,

range_status,

Dupe,

cluster,

name="n_specimens"

)

write.csv(

cluster_sizes,

file.path(

output_dir,

"Cluster_sizes.csv"

),

row.names=FALSE

)

#Remove clusters below threshold

cluster_data <-

cluster_data %>%

left_join(

cluster_sizes,

by=c(

"pair_id",

"AviList",

"range_status",

"Dupe",

"cluster"

)

) %>%

filter(

is.na(cluster) |

n_specimens >= minimum_cluster_size

)

## Check which species have no clusters remaining after filtering

cluster_species_after_filter <- cluster_data %>%
  filter(!is.na(cluster)) %>%
  distinct(pair_id, AviList, range_status, Dupe)

cluster_species_absent_after_filter <- cluster_sizes %>%
  distinct(pair_id, AviList, range_status, Dupe) %>%
  anti_join(
    cluster_species_after_filter,
    by = c("pair_id", "AviList", "range_status", "Dupe")
  )

write.csv(
  cluster_species_absent_after_filter,
  file.path(output_dir, "Species_absent_after_cluster_size_filter.csv"),
  row.names = FALSE
)


############################################################
## Build convex hulls
############################################################

cluster_hulls <-

cluster_data %>%

filter(!is.na(cluster)) %>%

group_by(

pair_id,

AviList,

range_status,

Dupe,

cluster,

cluster_id

) %>%

summarise(

n_specimens = first(n_specimens),

geometry = st_combine(geometry) %>%
           st_convex_hull(),

.groups="drop"

)

cluster_hulls$area_km2 <-

as.numeric(

st_area(cluster_hulls)

)/1e6

#Save the hulls
st_write(

cluster_hulls,

file.path(

output_dir,

"Cluster_hulls.gpkg"

),

delete_dsn=TRUE,

quiet=TRUE
)

#Write the hull summaries
write.csv(

st_drop_geometry(cluster_hulls),

file.path(

output_dir,

"Cluster_hull_summary.csv"

),

row.names=FALSE
)

############################################################
## SECTION 3: Calculate summary statistics for each cluster
############################################################

############################################################
## Helper functions
############################################################

safe_cv <- function(x){

  x <- x[!is.na(x)]

  if(length(x) < 2)
    return(NA_real_)

  m <- mean(x)

  s <- sd(x)

  if(is.na(m) || m <= 0)
    return(NA_real_)

  s/m

}

safe_log_cv <- function(x){

  cv <- safe_cv(x)

  if(is.na(cv) || cv <= 0)
    return(NA_real_)

  log(cv)

}

############################################################
## Calculate cluster statistics
############################################################

cluster_summary <- list()

for(trait in trait_cols){

  cat("Processing:", trait, "\n")

  tmp <-

    cluster_data %>%

    st_drop_geometry() %>%

    filter(!is.na(cluster_id))

  ## Convert trait to numeric

  tmp[[trait]] <- suppressWarnings(
    as.numeric(tmp[[trait]])
  )

  stats <-

    tmp %>%

    group_by(

      pair_id,

      AviList,

      range_status,

      Dupe,

      cluster,

      cluster_id

    ) %>%

    summarise(

      trait = trait,

      n = sum(!is.na(.data[[trait]])),

      mean = mean(.data[[trait]], na.rm = TRUE),

      sd = sd(.data[[trait]], na.rm = TRUE),

      cv = safe_cv(.data[[trait]]),

      log_cv = safe_log_cv(.data[[trait]]),

      .groups = "drop"

    )

  cluster_summary[[trait]] <- stats

}

cluster_summary <-

  bind_rows(cluster_summary)

#Add hull information to this
cluster_summary <-

cluster_summary %>%

left_join(

st_drop_geometry(

cluster_hulls %>%

select(

cluster_id,

area_km2

)

),

by="cluster_id"

)

#Remove unusable clusters
cluster_summary <-

cluster_summary %>%

filter(

!is.na(log_cv),

n >= minimum_cluster_size

)

write.csv(

cluster_summary,

file.path(

output_dir,

"Cluster_summary.csv"

),

row.names = FALSE
)

############################################################
## Diagnostic
############################################################

p_n <-

ggplot(

cluster_summary,

aes(

log(n),

log_cv,

colour = range_status

)

)+

geom_point(

size = 2,

alpha = 0.8

)+

geom_smooth(

method = "lm",

se = TRUE

)+

facet_wrap(

~trait,

scales = "free_y"

)+

scale_colour_manual(

values = c(

Restricted = "#2C7BB6",

Widespread = "#D95F02"

)

)+

theme_bw(base_size = 12)+

labs(

x = "log(cluster size)",

y = "log(CV)",

colour = ""

)

ggsave(

file.path(

output_dir,

"logCV_vs_sample_size.pdf"

),

plot = p_n,

width = 12,

height = 8
)

#Removing 1 obvious beak width outlier
cluster_summary %>%
  filter(
    trait == "Beak.Depth"
  ) %>%
  arrange(desc(log_cv)) %>%
  select(
    pair_id,
    AviList,
    cluster_id,
    n,
    log_cv,
  )

clusters_to_remove <- c(
  "Hypsipetes philippinus_26"
)

cluster_summary <-

  cluster_summary %>%

  filter(
    !cluster_id %in% clusters_to_remove
  )


############################################################
## SECTION 4
## Final mixed-effects model
############################################################

model_results <- list()

for(current_trait in unique(cluster_summary$trait)){

  dat <- cluster_summary %>%
    filter(trait == current_trait) %>%
    mutate(range_status = factor(
      range_status,
      levels = c("Restricted", "Widespread")
    ))

  if(nrow(dat) < 10) {
    message("Skipping ", current_trait, ": fewer than 10 clusters.")
    next
  }

  # AviList is crossed with pair_id: repeated species retain one shared
  # random intercept across all pairs in which they occur.
  model <- tryCatch(
    lmer(
      log_cv ~ log(n) + range_status + (1 | pair_id) + (1 | AviList),
      data = dat,
      REML = TRUE
    ),
    error = function(e) {
      message("Model failed for ", current_trait, ": ", e$message)
      NULL
    }
  )

  if(is.null(model)) next

  cat("\n", strrep("-", 60), "\n", current_trait, "\n", sep = "")
  print(summary(model))

  model_results[[current_trait]] <- broom.mixed::tidy(
    model,
    effects = "fixed"
  ) %>% mutate(trait = current_trait, AIC = AIC(model))

}

############################################################
## Combine tables
############################################################

model_results <- bind_rows(model_results)

write.csv(

  model_results,

  file.path(

    output_dir,

    "cluster_model_results.csv"

  ),

  row.names = FALSE

)

############################################################
## FINAL PLOT by range status
############################################################

library(ggforce)
library(dplyr)

############################################################
## Four-panel figure by range status
############################################################

trait_labels <- c(
  "Beak.Depth" = "Beak depth",
  "Beak.Length_Nares" = "Beak length",
  "Tarsus.Length" = "Tarsus length",
  "Wing.Length" = "Wing length"
)

############################################################
## Calculate box-plot and median summaries
############################################################

summary_df <-

  cluster_summary %>%

  group_by(trait, range_status) %>%

  summarise(

    n = sum(!is.na(log_cv)),

    median = median(log_cv, na.rm = TRUE),

    .groups = "drop"

  )

############################################################
## Prepare plotting data
############################################################

trait_order <- names(trait_labels)

summary_df <-

  summary_df %>%

  mutate(

    trait = factor(trait, levels = trait_order),

    range_status = factor(
      range_status,
      levels = c("Restricted", "Widespread")
    ),

    trait_label = factor(
      unname(trait_labels[as.character(trait)]),
      levels = unname(trait_labels)
    ),

    x_position = as.integer(range_status)

  )

plot_df <-

  cluster_summary %>%

  mutate(

    trait = factor(trait, levels = trait_order),

    range_status = factor(
      range_status,
      levels = c("Restricted", "Widespread")
    ),

    trait_label = factor(
      unname(trait_labels[as.character(trait)]),
      levels = unname(trait_labels)
    ),

    x_position = as.integer(range_status)

  )

############################################################
## Trait-level model significance labels
############################################################

trait_y_ranges <- plot_df %>%
  group_by(trait, trait_label) %>%
  summarise(
    y_max = max(log_cv, na.rm = TRUE),
    y_min = min(log_cv, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    y_span = pmax(y_max - y_min, 0.1),
    bracket_y = y_max + y_span * 0.12,
    bracket_bottom = y_max + y_span * 0.07,
    label_y = y_max + y_span * 0.17
  )

significance_df <- tibble(trait = factor(trait_order, levels = trait_order)) %>%
  left_join(
    model_results %>%
      filter(term == "range_statusWidespread") %>%
      transmute(trait = factor(trait, levels = trait_order), p_value = p.value),
    by = "trait"
  ) %>%
  left_join(trait_y_ranges, by = "trait") %>%
  mutate(
    significance = case_when(
      is.na(p_value) ~ "NS",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "NS"
    ),
    x_start = 1,
    x_end = 2,
    x_centre = 1.5
  )

############################################################
## Plot
############################################################

p1 <-

ggplot(

  plot_df,

  aes(

    x = x_position,

    y = log_cv,

    fill = range_status,

    colour = range_status,

    # Points are calculated separately for every trait x range-status group.
    group = interaction(trait, range_status)

  )

) +


############################################################
## Sina points
############################################################

geom_sina(

  maxwidth = 0.72,

  method = "density",

  scale = "width",

  alpha = 0.34,

  size = 2,

  stroke = 1,

  show.legend = FALSE

) +

############################################################
## Broad internal box plots and a clear median bar
############################################################

geom_boxplot(
  width = 0.64,
  alpha = 0.30,
  outlier.shape = NA,
  linewidth = 0.7,
  colour = "grey15",
  show.legend = FALSE
) +

geom_segment(
  data = summary_df,
  aes(
    x = x_position - 0.27,
    xend = x_position + 0.27,
    y = median,
    yend = median
  ),
  inherit.aes = FALSE,
  colour = "grey10",
  linewidth = 1.25
) +

# A transparent square layer supplies a clean square key for the top legend.
geom_point(
  data = summary_df,
  aes(x = x_position, y = median, fill = range_status),
  inherit.aes = FALSE,
  shape = 22,
  size = 5,
  alpha = 0,
  colour = "grey15",
  show.legend = TRUE
) +

############################################################
## Trait labels and significance brackets
############################################################

geom_segment(

  data = significance_df,

  aes(x = x_start, xend = x_end, y = bracket_y, yend = bracket_y),

  inherit.aes = FALSE,

  linewidth = 0.7,

  colour = "black"

) +

geom_segment(

  data = significance_df,

  aes(x = x_start, xend = x_start, y = bracket_bottom, yend = bracket_y),

  inherit.aes = FALSE,

  linewidth = 0.7,

  colour = "black"

) +

geom_segment(

  data = significance_df,

  aes(x = x_end, xend = x_end, y = bracket_bottom, yend = bracket_y),

  inherit.aes = FALSE,

  linewidth = 0.7,

  colour = "black"

) +

geom_text(

  data = significance_df,

  aes(x = x_centre, y = label_y, label = significance),

  inherit.aes = FALSE,

  size = 5,

  fontface = "bold"

) +

scale_x_continuous(
  breaks = NULL,
  limits = c(0.4, 3.45),
  expand = expansion(mult = c(0, 0))

) +

facet_wrap(~trait_label, ncol = 2, scales = "free_y") +

############################################################
## Colours
############################################################

scale_fill_manual(

  values = c(

    Restricted = "#2C7BB6",

    Widespread = "#D95F02"

  )

) +

scale_colour_manual(

  values = c(

    Restricted = "#2C7BB6",

    Widespread = "#D95F02"

  ),
  guide = "none"

) +

guides(
  fill = guide_legend(
    override.aes = list(shape = 22, size = 6, alpha = 1, colour = "grey15")
  )
) +

############################################################
## Labels
############################################################

labs(
  x = NULL,
  y = expression(log("Coefficient of variation")),
  fill = "Range status"

) +

############################################################
## Theme
############################################################

theme_classic(base_size = 15) +

theme(
  legend.position = "top",
  legend.title = element_blank(),
  legend.text = element_text(size = 13),
  legend.key = element_rect(fill = "white", colour = NA),
  legend.key.size = unit(0.7, "cm"),

  axis.title = element_text(
    size = 20,
    face = "bold"
  ),

  axis.text = element_text(
    size = 15
  ),

  axis.ticks.x = element_blank(),
  panel.spacing = unit(1.15, "lines"),
  strip.background = element_blank(),
  strip.text = element_text(size = 16, face = "bold"),
  plot.margin = margin(8, 14, 12, 14)

)

############################################################
## Save
############################################################

ggsave(

  file.path(
    output_dir,
    "logCV_by_range_status.pdf"
  ),

  p1,

  width = 11,
  height = 9,
  dpi = 600

)












############################################################
## Local population variance analysis
## Section 1 - Read data
############################################################

rm(list = ls())
set.seed(12345)

suppressPackageStartupMessages({

  library(tidyverse)

  library(sf)

  library(ggplot2)

  library(rnaturalearth)
  library(rnaturalearthdata)

})

############################################################
## Paths
############################################################

data_path <- "../data/Bulbul_AVONET_data_integrated_trimmed.csv"

output_dir <- "../results/Sister_pair_analysis_8_50km"

if(!dir.exists(output_dir)){
  dir.create(output_dir,
             recursive = TRUE)
}

############################################################
## Read data
############################################################

bulbuls <- read.csv(
  data_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

############################################################
## Sister pairs
############################################################

sister_pairs <- tribble(

  ~pair_id, ~restricted_species, ~widespread_species,

  "pair_01", "Pycnonotus taivanus", "Pycnonotus sinensis",

  "pair_02", "Chlorocichla falkensteini", "Chlorocichla flaviventris",

  "pair_03", "Arizelocichla tephrolaema", "Arizelocichla masukuensis",

  "pair_04", "Criniger olivaceus", "Criniger calurus",

  "pair_05", "Phyllastrephus viridiceps", "Phyllastrephus albigularis",

  "pair_06", "Phyllastrephus fischeri", "Phyllastrephus cabanisi",

  "pair_07", "Alophoixus frater", "Alophoixus tephrogenys",

  "pair_08", "Iole palawanensis", "Iole viridescens",

  "pair_09", "Hypsipetes mindorensis", "Hypsipetes philippinus",

  "pair_10", "Hypsipetes thompsoni", "Hypsipetes leucocephalus",

  "pair_11", "Microtarsus priocephalus", "Microtarsus melanocephalos",

  "pair_12", "Rubigula gularis", "Rubigula flaviventris",

  "pair_13", "Rubigula montis", "Rubigula dispar",

  "pair_14", "Pycnonotus davisoni", "Pycnonotus finlaysoni",

  "pair_15", "Pycnonotus capensis", "Pycnonotus barbatus"

)

############################################################
## SECTION 2
## Explore spatial locality structure
############################################################

library(tidyverse)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)

############################################################
## Keep only focal sister species
############################################################

species_keep <- unique(c(
  sister_pairs$restricted_species,
  sister_pairs$widespread_species
))

############################################################
## Create lookup table
############################################################

pair_lookup <-

  sister_pairs %>%

  mutate(

    facet_label = paste0(
      restricted_species,
      "\nvs\n",
      widespread_species
    )

  ) %>%

  pivot_longer(

    cols = c(
      restricted_species,
      widespread_species
    ),

    names_to = "status",

    values_to = "AviList"

  ) %>%

  mutate(

    range_status = ifelse(
      status == "restricted_species",
      "Restricted",
      "Widespread"
    )

  ) %>%

  select(
    AviList,
    pair_id,
    range_status,
    facet_label
  )

############################################################
## Join pair information
############################################################

plot_dat <-

  bulbuls %>%

  filter(

    AviList %in% species_keep,

    !is.na(Lat),

    !is.na(Lon)

  ) %>%

  left_join(
    pair_lookup,
    by = "AviList"
  )

############################################################
## Count duplicate coordinates
############################################################

plot_points <-

  plot_dat %>%

  group_by(

    pair_id,

    facet_label,

    AviList,

    range_status,

    Lon,

    Lat

  ) %>%

  summarise(

    n_specimens = n(),

    .groups = "drop"

  )

############################################################
## Load world map
############################################################

world <- ne_countries(
  scale = "medium",
  returnclass = "sf"
)

############################################################
## Multi-page PDF
############################################################

pdf(
  file.path(
    output_dir,
    "Sister_pair_locality_maps.pdf"
  ),
  width = 8,
  height = 6
)

for(current_pair in unique(plot_points$pair_id)){

  dat <-

    plot_points %>%

    filter(pair_id == current_pair)

  xmin <- min(dat$Lon) - 3
  xmax <- max(dat$Lon) + 3

  ymin <- min(dat$Lat) - 3
  ymax <- max(dat$Lat) + 3

  p <-

    ggplot() +

    geom_sf(

      data = world,

      fill = "grey97",

      colour = "grey75",

      linewidth = 0.25

    ) +

    geom_point(

      data = dat,

      aes(

        Lon,

        Lat,

        colour = range_status,

        size = n_specimens,

        alpha = n_specimens

      )

    ) +

    scale_colour_manual(

      values = c(

        Restricted = "#2C7BB6",

        Widespread = "#D95F02"

      )

    ) +

    scale_size_continuous(
      range = c(2,8)
    ) +

    scale_alpha_continuous(
      range = c(0.3,1)
    ) +

    coord_sf(

      xlim = c(xmin,xmax),

      ylim = c(ymin,ymax),

      expand = FALSE

    ) +

    # scale bar to show geographic distance
    ggspatial::annotation_scale(
      location = "bl",
      width_hint = 0.25
    ) +

    labs(

      title = unique(dat$facet_label),

      x = "Longitude",

      y = "Latitude",

      colour = "",

      size = "Specimens\nper coordinate",

      alpha = "Specimens\nper coordinate"

    ) +

    theme_bw(base_size = 12) +

    theme(

      legend.position = "right",

      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      )

    )

  print(p)

}

dev.off()

############################################################
## Histogram of duplicate coordinates
############################################################

coord_summary <-

  plot_points %>%

  count(n_specimens)

p_hist <-

ggplot(

  coord_summary,

  aes(

    n_specimens,

    n

  )

) +

geom_col() +

theme_bw(base_size = 12) +

labs(

  x = "Specimens sharing identical coordinates",

  y = "Number of unique coordinates"

)

ggsave(

  filename = file.path(

    output_dir,

    "Specimens_per_coordinate_histogram.pdf"

  ),

  plot = p_hist,

  width = 7,

  height = 5

)

############################################################
## Save coordinate table
############################################################

write.csv(

  plot_points,

  file.path(

    output_dir,

    "Unique_coordinate_summary.csv"

  ),

  row.names = FALSE

)

############################################################
## SECTION 3
## Nearest-neighbour distances
############################################################

library(FNN)

############################################################
## Convert to sf object
############################################################

specimen_sf <-

  plot_dat %>%

  st_as_sf(

    coords = c("Lon","Lat"),

    crs = 4326

  ) %>%

  st_transform(6933)

############################################################
## Calculate nearest-neighbour distance
############################################################

nn_results <- list()

for(sp in unique(specimen_sf$AviList)){

  dat <-

    specimen_sf %>%

    filter(AviList == sp)

  ## Need at least two specimens

  if(nrow(dat) < 2)
    next

  coords <-

    st_coordinates(dat)

  ## Find the first two neighbours.
  ## The first neighbour is always the point itself.

  nn <- get.knn(
    coords,
    k = 2
  )

  dat$nearest_neighbour_km <-

    nn$nn.dist[,2] / 1000

  nn_results[[sp]] <- dat

}

nn_data <-

  bind_rows(nn_results) %>%

  st_drop_geometry()

write.csv(

  nn_data,

  file.path(

    output_dir,

    "Nearest_neighbour_distances.csv"

  ),

  row.names = FALSE

)

p1 <-

ggplot(

  nn_data,

  aes(nearest_neighbour_km)

)+

geom_histogram(

  bins = 40,

  colour = "black",

  fill = "grey70"

)+

theme_bw(base_size = 12)+

labs(

  x = "Nearest neighbour distance (km)",

  y = "Number of specimens"

)

ggsave(

file.path(

output_dir,

"Nearest_neighbour_histogram.pdf"

),

plot = p1,

width = 7,

height = 5

)

species_order <-

nn_data %>%

group_by(AviList) %>%

summarise(

median_nn = median(nearest_neighbour_km)

) %>%

arrange(median_nn)

p2 <-

nn_data %>%

left_join(
species_order,
by="AviList"
) %>%

mutate(

AviList = factor(
AviList,
levels = species_order$AviList
)

) %>%

ggplot(

aes(

AviList,

nearest_neighbour_km

)

)+

geom_boxplot(

outlier.size = 0.4

)+

coord_flip()+

theme_bw(base_size = 11)+

labs(

x = "",

y = "Nearest neighbour distance (km)"

)

ggsave(

file.path(

output_dir,

"Nearest_neighbour_boxplots.pdf"

),

plot = last_plot(),

width = 8,

height = 10

)

nn_summary <-

nn_data %>%

group_by(

pair_id,

AviList,

range_status

)%>%

summarise(

n_specimens = n(),

mean_nn = mean(nearest_neighbour_km),

median_nn = median(nearest_neighbour_km),

sd_nn = sd(nearest_neighbour_km),

p90_nn = quantile(
nearest_neighbour_km,
0.90
),

max_nn = max(nearest_neighbour_km),

.groups="drop"

)

write.csv(

nn_summary,

file.path(

output_dir,

"Nearest_neighbour_summary.csv"

),

row.names=FALSE

)
