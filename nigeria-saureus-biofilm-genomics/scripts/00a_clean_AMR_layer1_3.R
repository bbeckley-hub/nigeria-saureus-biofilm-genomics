#!/usr/bin/env Rscript
# Three-step cleaning of amr_genes.csv: Layer 1 (category filter) + Layer 3 (canonicalisation)
library(dplyr)
library(tidyr)

amr <- read.csv('amr_genes.csv', stringsAsFactors = FALSE)

cat("BEFORE: ", nrow(amr), "rows | ", n_distinct(amr$Gene), "unique gene names\n\n")

# ============================================================
# STEP 1: Category filter (remove non-AMR entries)
# ============================================================
biofilm_genes <- c("icaC", "cna", "hld", "aur", "scn")

virulence_genes <- c(
  "hlgA","hlgB","hlgC","lukD","lukE","lukF-PV","lukS-PV",
  "sea","seb","sec2","sed","seh","sei","sej","sek","sel","sel26","sel27",
  "sel28","sel30","selX","selZ","sem","sen","seo","sep","seq","ser","ses",
  "set","seu","sey",
  "splA","splB","splE","tst","eta","etb","sak","ednB",
  "sasG","sdrM"
)

mutation_genes <- c(
  "gyrA_S84L","parC_S80F","parC_S80Y","parE_R444S",
  "folP_E208K","folP_F17L",
  "murA_D278E","murA_E291D","murA_G257D","murA_T396N",
  "glpT_A100V","glpT_L27F","glpT_V213I"
)

regulators <- c("arlR","arlS","mgrA","rsmA","mepR","qacR","blaR1","blaI","blaI_of_Z")

# Housekeeping genes captured by AMRFinderPlus --plus mode (not AMR determinants)
housekeeping_genes <- c("dha1", "dhap", "kdpd", "rlmh", "sepa")

exclude_regex <- paste0("^(", paste(c(virulence_genes,
                                     mutation_genes,
                                     regulators,
                                     housekeeping_genes), collapse = "|"), ")$")

amr_layer1 <- amr %>%
  filter(!Gene %in% biofilm_genes) %>%
  filter(!grepl(exclude_regex, Gene, ignore.case = TRUE))

cat("AFTER STEP 1 (remove biofilm/virulence/mutations/regulators/housekeeping): ",
    nrow(amr_layer1), "rows | ", n_distinct(amr_layer1$Gene), "unique names\n")
cat("  Removed:", nrow(amr) - nrow(amr_layer1), "rows\n\n")

# ============================================================
# STEP 3: Canonicalise names
# ============================================================
normalize_gene <- function(g) {
  g <- gsub("^\\([A-Za-z]+\\)", "", g)
  g <- gsub("_[0-9]+$", "", g)
  g <- gsub("_(bifunctional_protein|gen)$", "", g)
  g <- gsub("_beta-lactamase_\\(blaZ\\)", "", g)
  g <- gsub("_of_Z$", "", g)
  g <- gsub("Staphylococcus_aureus_", "", g)
  g <- gsub("Staphylococcus_intermedius_", "", g)
  g <- gsub("Enterococcus_faecalis_", "", g)
  g <- gsub("_conferring_resistance_to_mupirocin$", "", g)
  g <- gsub("_chloramphenicol_acetyltransferase$", "", g)
  g <- gsub("'", "", g)
  g <- gsub("-", "", g)
  g <- gsub("_", "", g)
  g <- gsub("\\.", "", g)
  g <- tolower(g)
  g <- gsub("^aac.*aph.*2.*ia.*$", "aac6aph2", g)
  g <- gsub("^aac6.*$", "aac6", g)
  g <- gsub("^aac3$", "aac3", g)
  g <- gsub("^aac6prime$", "aac6", g)
  g <- gsub("^aad6$", "aad6", g)
  g <- gsub("^aadd1$", "aadd", g)
  g <- gsub("^aadd$", "aadd", g)
  g <- gsub("^ant4.*$", "ant4", g)
  g <- gsub("^ant6.*$", "ant6", g)
  g <- gsub("^aph2.*$", "aph2", g)
  g <- gsub("^aph3.*$", "aph3", g)
  g <- gsub("^aph3prime$", "aph3", g)
  g <- gsub("^ant4prime$", "ant4", g)
  g <- gsub("^ant6prime$", "ant6", g)
  g <- gsub("^aph2dprime$", "aph2", g)
  g <- gsub("^bla.*$", "blaZ", g)
  g <- gsub("^pc1.*$", "blaZ", g)
  g <- gsub("^ble.*$", "ble", g)
  g <- gsub("^ermb.*$", "ermB", g)
  g <- gsub("^ermc.*$", "ermC", g)
  g <- gsub("^msra.*$", "msrA", g)
  g <- gsub("^mphc.*$", "mphC", g)
  g <- gsub("^lnua.*$", "lnuA", g)
  g <- gsub("^lina.*$", "linA", g)
  g <- gsub("^tetk.*$", "tetK", g)
  g <- gsub("^tetl.*$", "tetL", g)
  g <- gsub("^tetm.*$", "tetM", g)
  g <- gsub("^tet38$", "tet38", g)
  g <- gsub("^meca.*$", "mecA", g)
  g <- gsub("^dfrg.*$", "dfrG", g)
  g <- gsub("^dfrc.*$", "dfrC", g)
  g <- gsub("^dfrs1$", "dfrS1", g)
  g <- gsub("^cata.*$", "catA", g)
  g <- gsub("^cat$", "catA", g)
  g <- gsub("^catpc221.*$", "catA", g)
  g <- gsub("^fosb.*$", "fosB", g)
  g <- gsub("^mupa.*$", "mupA", g)
  g <- gsub("^qacab$", "qacAB", g)
  g <- gsub("^qaca.*$", "qacA", g)
  g <- gsub("^qacc.*$", "qacC", g)
  g <- gsub("^qacj.*$", "qacJ", g)
  g <- gsub("^nora.*$", "norA", g)
  g <- gsub("^norb.*$", "norB", g)
  g <- gsub("^norc.*$", "norC", g)
  g <- gsub("^mepa.*$", "mepA", g)
  g <- gsub("^mepb.*$", "mepB", g)
  g <- gsub("^lmrs.*$", "lmrS", g)
  g <- gsub("^sat.*$", "sat4", g)
  g <- gsub("^str.*$", "str", g)
  g <- gsub("^arsb.*$", "arsB", g)
  g <- gsub("^arsc.*$", "arsC", g)
  g <- gsub("^arsr.*$", "arsR", g)
  g <- gsub("^cadc.*$", "cadC", g)
  g <- gsub("^cadd.*$", "cadD", g)
  g <- gsub("^mera.*$", "merA", g)
  g <- gsub("^merb.*$", "merB", g)
  g <- gsub("^mert.*$", "merT", g)
  g
}

amr_layer3 <- amr_layer1 %>%
  mutate(Canonical = sapply(Gene, normalize_gene))

cat("AFTER STEP 3 (canonicalisation): ",
    n_distinct(amr_layer3$Canonical), "unique canonical genes\n")
cat("  Top 20 canonical genes:\n")
print(amr_layer3 %>% count(Canonical, sort = TRUE) %>% head(20))

# ============================================================
# STEP 2: Cross-database deduplication
# ============================================================
cat("\n========== DUPLICATE CHECK AFTER STEPS 1 + 3 ==========\n")
dup_pairs <- amr_layer3 %>%
  separate_rows(Genomes, sep = ";") %>%
  count(Genomes, Canonical) %>%
  filter(n > 1)
cat("Duplicate gene-isolate pairs remaining:", nrow(dup_pairs), "\n")
if (nrow(dup_pairs) > 0) {
  cat("Max duplicates for a single gene-isolate:", max(dup_pairs$n), "\n\n")
}

amr_clean <- amr_layer3 %>%
  separate_rows(Genomes, sep = ";") %>%
  distinct(Genomes, Canonical) %>%
  rename(Isolate = Genomes, Gene = Canonical)

cat("FINAL unique gene-isolate pairs:", nrow(amr_clean), "\n")
cat("Unique canonical genes:", n_distinct(amr_clean$Gene), "\n\n")

per_isolate <- amr_clean %>%
  count(Isolate, name = "AMR_gene_count") %>%
  arrange(desc(AMR_gene_count))

cat("========== AMR DETERMINANT COUNT PER ISOLATE ==========\n")
cat("Range:", min(per_isolate$AMR_gene_count), "-", max(per_isolate$AMR_gene_count), "\n")
cat("Median:", median(per_isolate$AMR_gene_count), "\n")
cat("Mean:", round(mean(per_isolate$AMR_gene_count), 1), "\n\n")
print(head(per_isolate, 10))
print(tail(per_isolate, 10))

write.csv(amr_clean,    "amr_clean_long.csv",        row.names = FALSE)
write.csv(per_isolate,  "amr_clean_per_isolate.csv", row.names = FALSE)
cat("\nSaved: amr_clean_long.csv, amr_clean_per_isolate.csv\n")
