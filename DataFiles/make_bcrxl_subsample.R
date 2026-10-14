# Builds the mass cytometry (CyTOF) subsample used in Practicum 12.
# Run once by the instructor; students load the resulting CSV.
#
# Data: Bodenmiller et al. (2012) Nat Biotechnol 30:858. PBMCs from 8 healthy
# donors, each split into an unstimulated ("Reference") and a B-cell receptor
# crosslinked ("BCR-XL") sample. Population labels are the manual gates from
# Bruggner et al. (2014) PNAS, as distributed in the HDCytoData package
# (Weber & Soneson 2019).
#
# Requires Bioconductor: BiocManager::install("HDCytoData")

library(HDCytoData)
library(SummarizedExperiment)
library(tidyverse)

cytof_se <- Bodenmiller_BCR_XL_SE()

marker_info <- as.data.frame(colData(cytof_se))
keep_markers <- marker_info$marker_class %in% c("type", "state")

# Raw (untransformed) intensities; students apply the arcsinh transform
cytof <- as_tibble(assay(cytof_se, "exprs")[, keep_markers]) |>
  rename_with(\(x) str_replace_all(x, "-", "_")) |>   # HLA-DR -> HLA_DR
  mutate(across(everything(), \(x) round(x, 2))) |>
  bind_cols(as_tibble(rowData(cytof_se))) |>
  transmute(patient = patient_id,
            condition = if_else(group_id == "BCR-XL", "stimulated", "unstimulated"),
            population = population_id,
            pick(where(is.numeric)))

# Equal number of cells from each of the 16 samples
set.seed(5030)
cells_per_sample <- 625
cytof_sub <- cytof |>
  group_by(patient, condition) |>
  slice_sample(n = cells_per_sample) |>
  ungroup()

message(nrow(cytof_sub), " cells kept")
print(count(cytof_sub, population))

# Which markers are lineage ("type") vs signaling ("state") markers
marker_info |>
  filter(keep_markers) |>
  transmute(marker = str_replace_all(marker_name, "-", "_"), marker_class) |>
  print()

write_csv(cytof_sub, "BCR_XL_cytof_subsample.csv")
