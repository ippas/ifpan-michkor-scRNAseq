###############################################################
### 1) FILTRUJEMY TYLKO ASTROCYTY
###############################################################

df_long_ast <- df_long %>% 
  filter(celltype == "Astrocytes")

df_sum_ast <- df_summary %>%
  ungroup() %>%
  filter(celltype == "Astrocytes")

###############################################################
### 2) USTAWIAMY POZYCJE SŁUPKÓW
###############################################################

treat_levels <- c("CTRL","CORT","DEX")

offsets <- c(
  "CTRL" = -0.25,
  "CORT" =  0.00,
  "DEX"  =  0.25
)

global_max <- max(df_sum_ast$mean_CPM + df_sum_ast$se_CPM, na.rm = TRUE)


###############################################################
### 3) TWORZYMY BRACKET CTRL–CORT
###############################################################

local_df_cort <- df_sum_ast %>%
  filter(treatment %in% c("CTRL","CORT")) %>%
  arrange(match(treatment, c("CTRL","CORT")))

local_max_cort <- max(local_df_cort$mean_CPM + local_df_cort$se_CPM, na.rm = TRUE)

bracket_cort <- tibble(
  x_ctrl  = 1 + offsets["CTRL"],
  x_treat = 1 + offsets["CORT"],
  y       = local_max_cort + 0.10 * global_max,
  label   = df_sum_ast$p_label[df_sum_ast$treatment == "CORT"]
)


###############################################################
### 4) TWORZYMY BRACKET CTRL–DEX
###############################################################

local_df_dex <- df_sum_ast %>%
  filter(treatment %in% c("CTRL","DEX")) %>%
  arrange(match(treatment, c("CTRL","DEX")))

local_max_dex <- max(local_df_dex$mean_CPM + local_df_dex$se_CPM, na.rm = TRUE)

bracket_dex <- tibble(
  x_ctrl  = 1 + offsets["CTRL"],
  x_treat = 1 + offsets["DEX"],
  y       = local_max_dex + 0.20 * global_max,
  label   = df_sum_ast$p_label[df_sum_ast$treatment == "DEX"]
)


###############################################################
### 5) WSZYSTKIE BRACKETY W JEDNYM TIBBLE
###############################################################

df_brackets <- bind_rows(bracket_cort, bracket_dex)


###############################################################
### 6) OSTATECZNY WYKRES
###############################################################

ggplot() +
  
  # surowe punkty
  geom_jitter(
    data = df_long_ast,
    aes(x = 1, y = CPM, color = treatment),
    position = position_jitterdodge(jitter.width = 0.10, dodge.width = 0.75),
    size = 2.2,
    alpha = 0.9
  ) +
  
  # słupki (CTRL, CORT, DEX)
  geom_col(
    data = df_sum_ast,
    aes(x = 1 + offsets[treatment], y = mean_CPM, fill = treatment),
    width = 0.22,
    alpha = 0.75
  ) +
  
  # errorbary
  geom_errorbar(
    data = df_sum_ast,
    aes(
      x = 1 + offsets[treatment],
      ymin = mean_CPM - se_CPM,
      ymax = mean_CPM + se_CPM
    ),
    width = 0.05
  ) +
  
  # poziome kreski bracketów
  geom_segment(
    data = df_brackets,
    aes(x = x_ctrl, xend = x_treat, y = y, yend = y),
    linewidth = 0.8
  ) +
  
  # pionowy koniec bracketu przy CTRL
  geom_segment(
    data = df_brackets,
    aes(x = x_ctrl, xend = x_ctrl,
        y = y, yend = y - 0.02*global_max),
    linewidth = 0.8
  ) +
  
  # pionowy koniec bracketu przy CORT/DEX
  geom_segment(
    data = df_brackets,
    aes(x = x_treat, xend = x_treat,
        y = y, yend = y - 0.02*global_max),
    linewidth = 0.8
  ) +
  
  # tekst p-value nad bracketami
  geom_text(
    data = df_brackets,
    aes(
      x = (x_ctrl + x_treat) / 2,
      y = y + 0.03 * global_max,
      label = label
    ),
    size = 5,
    fontface = "bold"
  ) +
  
  # style
  scale_color_manual(values=c("CTRL"="grey40","CORT"="firebrick2","DEX"="steelblue3")) +
  scale_fill_manual(values=c("CTRL"="grey40","CORT"="firebrick2","DEX"="steelblue3")) +
  
  theme_classic(base_size = 16) +
  labs(
    title = paste0("Astrocytes — ", unique(df_long_ast$gene_symbol)),
    x = "",
    y = "CPM"
  ) +
  
  scale_x_continuous(breaks = 1, labels = "Astrocytes") +
  theme(legend.position = "bottom")
