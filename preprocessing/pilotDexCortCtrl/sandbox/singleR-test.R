BiocManager::install("ensembldb")



install

# Loading test data.
library(TENxPBMCData)
new.data <- TENxPBMCData("pbmc4k")

# Loading reference data with Ensembl annotations.
library(celldex)
ref.data <- HumanPrimaryCellAtlasData(ensembl=TRUE)

# Performing predictions.
library(SingleR)
predictions <- SingleR(test=new.data, assay.type.test=1, 
                       ref=ref.data, labels=ref.data$label.main)


table(predictions$labels) %>% sum


# Zakładamy, że masz obiekt seurat_obj
# Usuń wszystko po myślniku (czyli symbol)
gene_ids <- sub("-.*", "", rownames(seurat_obj))

# Sprawdź duplikaty
dup_ids <- gene_ids[duplicated(gene_ids)]

length(dup_ids)  # ile ich jest?
head(dup_ids)    # co to za geny?

# Usuń duplikaty – zachowaj tylko pierwszy występujący
unique_idx <- !duplicated(gene_ids)
seurat_obj <- seurat_obj[unique_idx, ]
rownames(seurat_obj) <- gene_ids[unique_idx]


sce_obj <- as.SingleCellExperiment(seurat_obj)

predictions <- SingleR(
  test = sce_obj,
  ref = ref.data,
  labels = ref.data$label.main
)


predictions@metadata


# test my data
