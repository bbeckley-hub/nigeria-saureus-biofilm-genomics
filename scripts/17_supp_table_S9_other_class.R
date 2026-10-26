#!/usr/bin/env Rscript
# Supplementary Table S9: Determinants classified as "Other" after
# three-step AMR cleaning (non-classified resistance determinants)

library(dplyr)

amr <- read.csv('amr_clean_pass2_long.csv', stringsAsFactors = FALSE)

classify_amr <- function(gene) {
  g <- tolower(trimws(gene))
  if (grepl('^(mepa|mepb|lmrs|nora|norb|norc)$', g)) return('Intrinsic_efflux')
  if (grepl('^mec', g) || grepl('^bla', g) || grepl('^pc1', g)) return('Beta_lactam')
  if (grepl('^erm', g) || grepl('^msr', g) || grepl('^mph', g)) return('Macrolide')
  if (grepl('^aac', g) || grepl('^aph', g) || grepl('^ant', g) ||
      grepl('^aad', g) || grepl('^str$', g) || grepl('^sat', g)) return('Aminoglycoside')
  if (grepl('^tet', g)) return('Tetracycline')
  if (grepl('^dfr', g)) return('Trimethoprim')
  if (grepl('^cat', g)) return('Chloramphenicol')
  if (grepl('^lnu', g) || grepl('^lin', g)) return('Lincosamide')
  if (grepl('^qac', g)) return('Biocide')
  if (grepl('^fos', g)) return('Fosfomycin')
  if (grepl('^mup', g) || grepl('^iles', g)) return('Mupirocin')
  return('Other')
}

# Filter to "Other" class and summarise
other <- amr %>%
  mutate(Class = sapply(Gene, classify_amr)) %>%
  filter(Class == "Other") %>%
  count(Gene, name = "n_isolates") %>%
  arrange(desc(n_isolates))

# Add mechanism annotations where known
mechanism_lookup <- c(
  arsr = "Arsenic resistance (regulator)",
  arsb = "Arsenic resistance (efflux pump)",
  arsc = "Arsenic resistance (reductase)",
  cadc = "Cadmium resistance (regulator)",
  cadd = "Cadmium resistance (efflux)",
  mera = "Mercury resistance (mercuric reductase)",
  merb = "Mercury resistance (organomercurial lyase)",
  mert = "Mercury resistance (transport)",
  ble  = "Bleomycin resistance"
)

other$Mechanism <- ifelse(
  other$Gene %in% names(mechanism_lookup),
  mechanism_lookup[other$Gene],
  "Other non-classified determinant"
)

other$Target <- ifelse(
  grepl("^ars", other$Gene), "Arsenic",
  ifelse(grepl("^cad", other$Gene), "Cadmium",
  ifelse(grepl("^mer", other$Gene), "Mercury",
  ifelse(grepl("^ble", other$Gene), "Bleomycin",
         "Unknown / non-classified")))
)

write.csv(other, "Supplementary_Table_S9_AMR_other_class_definition.csv",
          row.names = FALSE)

cat("Saved Supplementary_Table_S9_AMR_other_class_definition.csv\n")
cat("Determinants in Other class:", nrow(other), "\n")
print(other)
