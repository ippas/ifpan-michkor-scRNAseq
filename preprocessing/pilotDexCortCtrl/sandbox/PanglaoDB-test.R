

panglaoDB_all_markers <- read.table(
  "data/celltype_references/PanglaoDB/PanglaoDB_markers_27_Mar_2020.tsv",
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  # check.names = FALSE,
  row.names = NULL
)


panglaoDB_all_markers %>% filter(cell.type == "Pericytes")

panglaoDB_all_markers %>% 
  filter(organ == "Brain") %>% 
  filter(cell.type %in% c("Immature neurons", "Astrocytes", "Bergmann glia", "Choroid plexus cells", "Meningeal cells", "Microglia", "Neurons", "Oligodendrocyte progenitor cells",
                          "Pinealocytes", "Radial glia cells", "Satellite glial cells", "Tanycytes")) %>% 
  select(c(official.gene.symbol, cell.type)) %>% 
  mutate(cell.type = str_replace_all(cell.type, " ", "_")) %>% 
  group_by(cell.type) %>% 
  split(.$cell.type) %>% 
  lapply(function(x) unique(x$official.gene.symbol)) -> panglaoDB_list_selectedMarkers


panglaoDB_all_markers %>% 
  # filter(organ == "Brain") %>% 
  # filter(cell.type %in% c("Immature neurons", "Astrocytes", "Bergmann glia", "Choroid plexus cells", "Meningeal cells", "Microglia", "Neurons", "Oligodendrocyte progenitor cells",
  #                         "Pinealocytes", "Radial glia cells", "Satellite glial cells", "Tanycytes")) %>% 
  select(c(official.gene.symbol, cell.type)) %>% 
  mutate(cell.type = str_replace_all(cell.type, " ", "_")) %>% 
  group_by(cell.type) %>% 
  split(.$cell.type) %>% 
  lapply(function(x) unique(x$official.gene.symbol)) -> panglaoDB_list_allMarkers

panglaoDB_list_allMarkers$Pericytes


panglaoDB_all_markers %>% filter(cell.type == "Pericytes")
