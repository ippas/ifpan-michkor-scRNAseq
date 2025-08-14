# ##############################################################################
# ---- prepare data ----
# ##############################################################################
merged_pilotDexCortCtrl@meta.data %>% 
  group_by(bc1_well, treatment, seurat_clusters_res0.5) %>% 
  nest() %>%
  mutate(n_cells = map(data, ~nrow(.x))) %>% 
  unnest(n_cells) %>% 
  select(-data) %>% 
  as.data.frame() -> cell_counts_per_cluster

# ##############################################################################
# ---- prepare barplot, summary number of cells per cluster ----
# ##############################################################################

# Podsumowanie
counts_df <- cell_counts_per_cluster %>%
  group_by(seurat_clusters_res0.5) %>%
  summarise(n_cells = sum(n_cells), .groups = "drop")

total_cells <- sum(counts_df$n_cells)
n_clusters  <- nrow(counts_df)
max_cells   <- max(counts_df$n_cells)

# Rozmiary tekstów (pt) i przeliczenie do geoms
axis_text_pt   <- 12
axis_title_pt  <- 14
plot_title_pt  <- 16
legend_title_pt<- 14
legend_text_pt <- 12
bar_label_pt   <- 12
top_label_pt   <- 14
pt_to_geom     <- function(pt) pt / 2.845

library(ggplot2)

ggplot(
  counts_df,
  aes(x = seurat_clusters_res0.5, y = n_cells, fill = seurat_clusters_res0.5)
) +
  geom_col() +
  # liczby nad słupkami
  geom_text(
    aes(label = n_cells),
    angle = 45, hjust = -0.1, vjust = -0.5,
    size = pt_to_geom(bar_label_pt)
  ) +
  # znak sumy z indeksem dolnym i wynikiem
  annotate(
    "text",
    x = (n_clusters + 1) / 2,
    y = max_cells * 1.06,
    fontface = "bold",
    label = paste0("Sigma[cells] == ", total_cells),
    parse = TRUE,                        # włącza notację matematyczną
    size = pt_to_geom(18)
  ) +
  scale_fill_manual(values = cluster_colors) +
  coord_cartesian(ylim = c(0, max_cells * 1.12), clip = "off") +
  labs(
    x = "Cluster (res0.5)",
    y = "n cells",
    fill = "Cluster",
    title = "Cell counts per cluster"
  ) +
  theme_classic() +
  theme(
    axis.text.x  = element_text(size = axis_text_pt),
    axis.text.y  = element_text(size = axis_text_pt),
    axis.title.x = element_text(size = axis_title_pt),
    axis.title.y = element_text(size = axis_title_pt),
    plot.title   = element_text(size = plot_title_pt),
    legend.title = element_text(size = legend_title_pt),
    legend.text  = element_text(size = legend_text_pt)
  )


# ggplot(cell_counts_per_cluster %>%
#          group_by(seurat_clusters_res0.5, treatment) %>%
#          summarise(n_cells = sum(n_cells), .groups = "drop"),
#        aes(x = seurat_clusters_res0.5, y = n_cells, fill = seurat_clusters_res0.5)) +
#   geom_col() +
#   facet_grid(~ treatment) +
#   scale_fill_manual(values = setNames(cluster_colors, cluster_levels)) +
#   labs(x = "Cluster (res0.5)", y = "n cells",
#        fill = "Cluster",
#        title = "Cell counts per cluster, faceted by treatment") +
#   theme_classic()



# ##############################################################################
# ---- number per cells, facet by treatment ----
# ##############################################################################
merged_pilotDexCortCtrl@meta.data %>% 
  group_by(bc1_well, treatment, seurat_clusters_res0.5) %>% 
  nest() %>%
  mutate(n_cells = map(data, ~nrow(.x))) %>% 
  unnest(n_cells) %>% 
  select(-data) %>% 
  as.data.frame() -> cell_counts_per_cluster


# --- Kolejność facetów ---
desired_order <- c("CTRL", "DEX", "CORT")

# --- Dane pod faceting (liczebności) ---
counts_df_facet <- cell_counts_per_cluster %>%
  mutate(treatment = factor(as.character(treatment), levels = desired_order)) %>%
  group_by(seurat_clusters_res0.5, treatment) %>%
  summarise(n_cells = sum(n_cells), .groups = "drop")

# --- Sumy i maksima per treatment (do anotacji Σcells) ---
stats_df <- counts_df_facet %>%
  group_by(treatment) %>%
  summarise(
    total_cells = sum(n_cells),
    max_cells   = max(n_cells),
    n_clusters  = n_distinct(seurat_clusters_res0.5),
    .groups = "drop"
  ) %>%
  mutate(treatment = factor(as.character(treatment), levels = desired_order))

# --- Rozmiary tekstu ---
axis_text_pt    <- 12
axis_title_pt   <- 14
plot_title_pt   <- 16
legend_title_pt <- 14
legend_text_pt  <- 12
bar_label_pt    <- 12
pt_to_geom      <- function(pt) pt / 2.845

# --- Kolory klastrów (jeśli masz wektor cluster_colors) ---
cluster_levels <- sort(unique(counts_df_facet$seurat_clusters_res0.5))

ggplot(
  counts_df_facet,
  aes(x = seurat_clusters_res0.5, y = n_cells, fill = seurat_clusters_res0.5)
) +
  geom_col() +
  # liczby nad słupkami
  geom_text(
    aes(label = n_cells),
    angle = 45, hjust = -0.1, vjust = -0.5,
    size = pt_to_geom(bar_label_pt)
  ) +
  # Σcells w każdym facie
  geom_text(
    data = stats_df,
    aes(
      x = (n_clusters + 1) / 2,
      y = max_cells * 1.06,
      label = paste0("Sigma[cells] == ", total_cells)
    ),
    inherit.aes = FALSE,
    parse = TRUE,
    fontface = "bold",
    size = pt_to_geom(18)
  ) +
  facet_grid(treatment ~ ., scales = "free_y") +  # pionowo: CTRL, DEX, CORT
  coord_cartesian(clip = "off") +
  labs(
    x = "Cluster (res0.5)",
    y = "n cells",
    fill = "Cluster",
    title = "Cell counts per cluster, faceted by treatment"
  ) +
  theme_classic() +
  theme(
    axis.text.x  = element_text(size = axis_text_pt),
    axis.text.y  = element_text(size = axis_text_pt),
    axis.title.x = element_text(size = axis_title_pt),
    axis.title.y = element_text(size = axis_title_pt),
    plot.title   = element_text(size = plot_title_pt),
    legend.title = element_text(size = legend_title_pt),
    legend.text  = element_text(size = legend_text_pt)
  )



# ##############################################################################
# ---- normalize number cells  ----
# ##############################################################################

# --- 1) Zliczenia per klaster i treatment ---
counts_df_facet <- cell_counts_per_cluster %>%
  group_by(seurat_clusters_res0.5, treatment) %>%
  summarise(n_cells = sum(n_cells), .groups = "drop")

# --- 2) Wymuś kolejność facetów: CTRL, DEX, CORT ---
desired_order <- c("CTRL", "DEX", "CORT")
counts_df_facet <- counts_df_facet %>%
  mutate(treatment = factor(as.character(treatment), levels = desired_order))

# --- 3) Suma komórek per treatment + frakcje ---
totals_df <- counts_df_facet %>%
  group_by(treatment) %>%
  summarise(total_cells = sum(n_cells), .groups = "drop")

props_df <- counts_df_facet %>%
  left_join(totals_df, by = "treatment") %>%
  mutate(prop = n_cells / total_cells) %>%
  mutate(treatment = factor(as.character(treatment), levels = desired_order))  # ponownie dla pewności

# --- 4) Statystyki do anotacji (pozycja Σ i środek osi X) ---
stats_df <- props_df %>%
  group_by(treatment) %>%
  summarise(
    max_prop    = max(prop),
    n_clusters  = n_distinct(seurat_clusters_res0.5),
    total_cells = unique(total_cells),
    .groups = "drop"
  ) %>%
  mutate(
    sum_label = paste0("Sigma[cells] == ", total_cells),
    treatment = factor(as.character(treatment), levels = desired_order)
  )

# --- 5) Rozmiary tekstu ---
axis_text_pt    <- 12
axis_title_pt   <- 14
plot_title_pt   <- 16
legend_title_pt <- 14
legend_text_pt  <- 12
bar_label_pt    <- 12
pt_to_geom      <- function(pt) pt / 2.845

# --- 6) Kolory klastrów ---
cluster_levels <- sort(unique(props_df$seurat_clusters_res0.5))

# --- 7) Wykres ---
ggplot(
  props_df,
  aes(x = seurat_clusters_res0.5, y = prop, fill = seurat_clusters_res0.5)
) +
  geom_col() +
  # etykiety % nad słupkami
  geom_text(
    aes(label = percent(prop, accuracy = 0.1)),
    angle = 45, hjust = -0.1, vjust = -0.5,
    size = pt_to_geom(bar_label_pt)
  ) +
  # Σcells w każdym panelu
  geom_text(
    data = stats_df,
    aes(
      x = (n_clusters + 1) / 2,
      y = max_prop * 1.06,
      label = sum_label
    ),
    inherit.aes = FALSE,
    parse = TRUE,
    fontface = "bold",
    size = pt_to_geom(18)
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = setNames(cluster_colors, cluster_levels)) +
  facet_grid(rows = vars(treatment)) +   # respektuje levels() z desired_order
  coord_cartesian(ylim = c(0, max(props_df$prop) * 1.12), clip = "off") +
  labs(
    x = "Cluster (res0.5)",
    y = "fraction of cells",
    fill = "Cluster",
    title = "Fraction of cells per cluster, faceted by treatment"
  ) +
  theme_classic() +
  theme(
    axis.text.x  = element_text(size = axis_text_pt),
    axis.text.y  = element_text(size = axis_text_pt),
    axis.title.x = element_text(size = axis_title_pt),
    axis.title.y = element_text(size = axis_title_pt),
    plot.title   = element_text(size = plot_title_pt),
    legend.title = element_text(size = legend_title_pt),
    legend.text  = element_text(size = legend_text_pt)
  )

