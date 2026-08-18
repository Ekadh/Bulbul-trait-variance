library(ape)
library(ggtree)
library(dplyr)
library(ggplot2)
library(tidyr)

# This script can be run either with `code` or the project directory as the
# working directory.
project_root <- if (dir.exists("data")) "." else ".."
data_path <- file.path(project_root, "data")
results_path <- file.path(project_root, "results", "Phylogeny")
dir.create(results_path, recursive = TRUE, showWarnings = FALSE)

##############################################################
## Data and tree
##############################################################

bulbul <- read.csv(
  file.path(data_path, "Bulbul_AVONET_data_integrated_trimmed.csv"),
  stringsAsFactors = FALSE
)

tree <- read.tree(
  file.path(data_path, "Phylogeny", "mcc_dated_clements_McTavish_2025.nex")
)[[6]]
tree$tip.label <- gsub("_", " ", tree$tip.label)

# Retain all taxa represented by either the specimen data or the one species
# that is intentionally shown with zero specimens.
tree <- keep.tip(
  tree,
  intersect(tree$tip.label, c(unique(bulbul$AviList), "Hypsipetes moheliensis"))
)

traits <- c(
  "Beak.Length_Culmen", "Beak.Width", "Beak.Depth",
  "Tarsus.Length", "Wing.Length", "Tail.Length"
)

bulbul_counts <- bulbul %>%
  mutate(n_traits_measured = rowSums(!is.na(pick(all_of(traits))))) %>%
  filter(n_traits_measured >= 3) %>%
  count(AviList, name = "n")

bulbul_ESI <- read.csv(file.path(data_path, "BIRDBASE", "Bulbul_Birdbase_data.csv"))
bulbul_ESI$ESI <- 2 - bulbul_ESI$ESI
bulbul_range <- read.csv(file.path(data_path, "Bulbul_rangesize.csv"))

species_data <- tibble(AviList = tree$tip.label) %>%
  left_join(bulbul_counts, by = "AviList") %>%
  mutate(n = replace_na(n, 0L)) %>%
  left_join(
    bulbul_ESI %>% select(AviList.v1.2025, ESI),
    by = c("AviList" = "AviList.v1.2025")
  ) %>%
  left_join(bulbul_range %>% select(AviList, Range.Size), by = "AviList")

##############################################################
## Shared tree/ring geometry
##############################################################

# ggtree calculates the precise fan angle for every node.  Those angles—not a
# separately generated sequence—are used below for the trait rings.  This is
# what keeps each ring mark aligned with its corresponding tree tip.
fan_layout <- ggtree(tree, layout = "fan", open.angle = 180)$data
nodes <- fan_layout %>%
  transmute(
    node, parent, label, isTip, radius = x,
    layout_angle = angle * pi / 180
  )

# Put the terminal taxa on an evenly spaced semicircle.  Internal-node angles
# are then calculated from their descendants, rather than retained from the
# original layout; this keeps every branch connected smoothly to the new tip
# positions.
equal_tip_angles <- nodes %>%
  filter(isTip) %>%
  arrange(layout_angle) %>%
  transmute(node, tip_angle = seq(0, pi, length.out = n()))

nodes <- nodes %>%
  left_join(equal_tip_angles, by = "node")

children_by_parent <- nodes %>%
  filter(node != parent) %>%
  { split(.$node, .$parent) }
angle_by_node <- setNames(rep(NA_real_, nrow(nodes)), nodes$node)
tip_nodes <- nodes %>% filter(isTip)
angle_by_node[as.character(tip_nodes$node)] <- tip_nodes$tip_angle

node_angle <- function(node) {
  node_key <- as.character(node)
  if (!is.na(angle_by_node[node_key])) return(angle_by_node[node_key])

  child_nodes <- children_by_parent[[node_key]]
  angle_by_node[node_key] <<- mean(vapply(child_nodes, node_angle, numeric(1)))
  angle_by_node[node_key]
}

invisible(vapply(nodes$node, node_angle, numeric(1)))

nodes <- nodes %>%
  mutate(angle = unname(angle_by_node[as.character(node)])) %>%
  select(-layout_angle, -tip_angle)

tree_radius <- max(nodes$radius[nodes$isTip])

tip_data <- nodes %>%
  filter(isTip) %>%
  transmute(AviList = label, angle) %>%
  left_join(species_data, by = "AviList")

# Scale the three variables independently.  Missing trait/range values have a
# zero-length mark rather than shifting the alignment of any other species.
scale_to_radius <- function(x, upper) {
  out <- rep(0, length(x))
  observed <- !is.na(x) & is.finite(x)
  if (sum(observed) == 1L) out[observed] <- upper
  if (sum(observed) > 1L) out[observed] <- scales::rescale(x[observed], to = c(0, upper))
  out
}

tip_data <- tip_data %>%
  mutate(
    ESI_scaled = scale_to_radius(ESI, 0.9),
    Range_scaled = scale_to_radius(log10(if_else(Range.Size > 0, Range.Size, NA_real_)), 1.2),
    N_scaled = scale_to_radius(pmin(n, 100), 1.5)
  )

esi_inner_radius <- tree_radius + 0.40
range_inner_radius <- tree_radius + 1.50
n_inner_radius <- tree_radius + 3.25

make_ring_coordinates <- function(data, inner_radius, value_column, prefix) {
  value <- data[[value_column]]
  data %>%
    mutate(
      "{prefix}_x" := inner_radius * cos(angle),
      "{prefix}_y" := inner_radius * sin(angle),
      "{prefix}_xend" := (inner_radius + value) * cos(angle),
      "{prefix}_yend" := (inner_radius + value) * sin(angle)
    )
}

ring_data <- tip_data %>%
  make_ring_coordinates(esi_inner_radius, "ESI_scaled", "esi") %>%
  make_ring_coordinates(range_inner_radius, "Range_scaled", "range") %>%
  make_ring_coordinates(n_inner_radius, "N_scaled", "n")

# Pale background sectors for the focal genera.  A separate sector is drawn
# for each contiguous genus run, so a non-contiguous genus is not shaded across
# intervening taxa.
major_genera <- tibble(
  genus = c("Pycnonotus", "Hypsipetes", "Alophoixus", "Phyllastrephus", "Arizelocichla"),
  colour = c("#D98C8C", "#7FAEC8", "#8CBF9B", "#D0AD62", "#B48CC4")
)

tip_angle_half_step <- median(diff(sort(tip_data$angle))) / 2
genus_runs <- tip_data %>%
  transmute(genus = sub(" .*", "", AviList), angle) %>%
  arrange(angle) %>%
  mutate(run = cumsum(genus != lag(genus, default = first(genus)))) %>%
  filter(genus %in% major_genera$genus) %>%
  group_by(genus, run) %>%
  summarise(
    start_angle = max(0, min(angle) - tip_angle_half_step),
    end_angle = min(pi, max(angle) + tip_angle_half_step),
    .groups = "drop"
  ) %>%
  left_join(major_genera, by = "genus") %>%
  mutate(sector_id = row_number())

genus_sectors <- bind_rows(lapply(seq_len(nrow(genus_runs)), function(i) {
  sector <- genus_runs[i, ]
  sector_angle <- seq(sector$start_angle, sector$end_angle, length.out = 80)
  bind_rows(
    tibble(sector_id = sector$sector_id, genus = sector$genus, x = 0, y = 0),
    tibble(
      sector_id = sector$sector_id,
      genus = sector$genus,
      x = tree_radius * cos(sector_angle),
      y = tree_radius * sin(sector_angle)
    )
  )
}))

# Draw the fan tree in the same Cartesian coordinate system as the rings.
# A rectangular phylogram edge becomes a radial segment plus an arc; therefore
# the branch endpoints are exactly the coordinates used to place the rings.
edge_data <- nodes %>%
  filter(!is.na(parent), parent > 0, node != parent) %>%
  transmute(child = node, parent, child_radius = radius, child_angle = angle) %>%
  left_join(
    nodes %>% transmute(parent = node, parent_radius = radius, parent_angle = angle),
    by = "parent"
  ) %>%
  mutate(edge_id = row_number())

radial_edges <- edge_data %>%
  transmute(
    x = parent_radius * cos(child_angle),
    y = parent_radius * sin(child_angle),
    xend = child_radius * cos(child_angle),
    yend = child_radius * sin(child_angle)
  )

arc_edges <- bind_rows(lapply(seq_len(nrow(edge_data)), function(i) {
  edge <- edge_data[i, ]
  arc_angle <- seq(edge$parent_angle, edge$child_angle, length.out = 20)
  tibble(
    edge_id = edge$edge_id,
    x = edge$parent_radius * cos(arc_angle),
    y = edge$parent_radius * sin(arc_angle)
  )
}))

##############################################################
## Figure
##############################################################

col_ESI <- "#6A3D9A"
col_RANGE <- "#D4A017"
col_N <- "#2C7FB8"

max_radius <- max(
  esi_inner_radius + max(ring_data$ESI_scaled),
  range_inner_radius + max(ring_data$Range_scaled),
  n_inner_radius + max(ring_data$N_scaled)
)
padding <- 1.2

legend_df <- tibble(
  x1 = c(-max_radius * 0.80, -1.0, max_radius * 0.45),
  x2 = x1 + 1.25,
  y = -3.25,
  colour = c(col_ESI, col_RANGE, col_N),
  title = c("Ecological\nspecialisation", "Geographic\nrange size", "Specimen\ncount"),
  subtitle = c(
    paste0("Maximum = ", round(max(ring_data$ESI, na.rm = TRUE), 2)),
    paste0("Maximum = ", format(round(max(ring_data$Range.Size, na.rm = TRUE)), big.mark = ","), " km²"),
    "Values capped at 100"
  )
)

genus_legend <- major_genera %>%
  mutate(
    x1 = seq(-max_radius * 0.86, max_radius * 0.46, length.out = n()),
    x2 = x1 + 0.80,
    y = -5.30
  )

base_plot <- ggplot() +
  geom_polygon(
    data = genus_sectors,
    aes(x, y, group = sector_id, fill = genus),
    alpha = 0.22,
    colour = NA,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = setNames(major_genera$colour, major_genera$genus)) +
  geom_path(
    data = arc_edges,
    aes(x, y, group = edge_id),
    linewidth = 0.35,
    colour = "black"
  ) +
  geom_segment(
    data = radial_edges,
    aes(x, y, xend = xend, yend = yend),
    linewidth = 0.35,
    colour = "black"
  )

ring_layers <- list(
  geom_segment(
    data = ring_data,
    aes(x = esi_x, y = esi_y, xend = esi_xend, yend = esi_yend),
    colour = col_ESI, linewidth = 1.6, lineend = "round"
  ),
  geom_segment(
    data = ring_data,
    aes(x = range_x, y = range_y, xend = range_xend, yend = range_yend),
    colour = col_RANGE, linewidth = 1.6, lineend = "round"
  ),
  geom_segment(
    data = ring_data,
    aes(x = n_x, y = n_y, xend = n_xend, yend = n_yend),
    colour = col_N, linewidth = 1.6, lineend = "round"
  )
)

plot_limits <- coord_equal(
  xlim = c(-(max_radius + padding), max_radius + padding),
  ylim = c(-6.1, max_radius + padding),
  expand = FALSE
)

p_tree <- base_plot + plot_limits + theme_void()
p_rings <- ggplot() + ring_layers + plot_limits + theme_void()

p_combined <- base_plot +
  ring_layers +
  geom_segment(
    data = legend_df,
    aes(x = x1, y = y, xend = x2, yend = y, colour = colour),
    linewidth = 5, lineend = "round", show.legend = FALSE
  ) +
  scale_colour_identity() +
  geom_text(
    data = legend_df,
    aes(x = x2 + 0.45, y = y + 0.28, label = title),
    fontface = "bold", size = 6.2, hjust = 0, lineheight = 0.95
  ) +
  geom_text(
    data = legend_df,
    aes(x = x2 + 0.45, y = y - 0.95, label = subtitle),
    size = 5.4, hjust = 0
  ) +
  geom_text(
    aes(x = -max_radius * 0.86, y = -4.72, label = "Highlighted genera"),
    fontface = "bold", size = 5, hjust = 0
  ) +
  geom_segment(
    data = genus_legend,
    aes(x = x1, y = y, xend = x2, yend = y),
    linewidth = 5, lineend = "round", colour = genus_legend$colour
  ) +
  geom_text(
    data = genus_legend,
    aes(x = x2 + 0.30, y = y, label = genus),
    size = 4.4, hjust = 0
  ) +
  plot_limits +
  theme_void()

ggsave(
  file.path(results_path, "Bulbul_phylogeny_concentric_rings.png"),
  p_combined, width = 14, height = 14, dpi = 600, bg = "transparent"
)

# These two optional transparent exports share the exact same dimensions and
# coordinate limits as the combined figure, so they can also be overlaid.
ggsave(
  file.path(results_path, "Bulbul_phylogeny_fan.pdf"),
  p_tree, width = 14, height = 14, dpi = 600, bg = "transparent"
)
ggsave(
  file.path(results_path, "Bulbul_concentric_rings_publication.png"),
  p_rings, width = 14, height = 14, dpi = 600, bg = "transparent"
)
