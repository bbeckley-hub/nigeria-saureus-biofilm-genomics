#!/usr/bin/env Rscript
# Second-pass: consolidate remaining name variants in amr_clean_long.csv
library(dplyr)
library(tidyr)

amr <- read.csv('amr_clean_long.csv', stringsAsFactors = FALSE)

cat("BEFORE: ", nrow(amr), "pairs | ", n_distinct(amr$Gene), "unique genes\n\n")

final_normalize <- function(g) {
  g <- gsub("\\(", "", g); g <- gsub("\\)", "", g)
  g <- gsub("'", "", g); g <- gsub("-", "", g)
  g <- gsub("_", "", g); g <- gsub("\\.", "", g)
  g <- tolower(g)
  g <- gsub("^ermb$",  "ermB", g)
  g <- gsub("^ermc$",  "ermC", g)
  g <- gsub("^erma$",  "ermA", g)
  g <- gsub("^lnua$",  "lnuA", g)
  g <- gsub("^lina$",  "linA", g)
  g <- gsub("^tetk$",  "tetK", g)
  g <- gsub("^tetl$",  "tetL", g)
  g <- gsub("^tetm$",  "tetM", g)
  g <- gsub("^tet38$", "tet38", g)
  g <- gsub("^ant6ia$",   "ant6", g)
  g <- gsub("^ant4ia$",   "ant4", g)
  g <- gsub("^ant6$",     "ant6", g)
  g <- gsub("^ant4$",     "ant4", g)
  g <- gsub("^aph3iiia$", "aph3", g)
  g <- gsub("^aph3iii$",  "aph3", g)
  g <- gsub("^aph3$",     "aph3", g)
  g <- gsub("^aph2ia$",   "aph2", g)
  g <- gsub("^aph2$",     "aph2", g)
  g <- gsub("^aac6aph2$", "aac6aph2", g)
  g <- gsub("^aadd$",     "aadD", g)
  g <- gsub("^aad6$",     "aad6", g)
  g <- gsub("^chloramphenicolacetyltransferase$", "catA", g)
  g <- gsub("^cata$",     "catA", g)
  g <- gsub("^catb$",     "catB", g)
  g <- gsub("^sat4$",     "sat4", g)
  g <- gsub("^sat$",      "sat4", g)
  g
}

amr2 <- amr %>%
  mutate(Gene = sapply(Gene, final_normalize)) %>%
  distinct(Isolate, Gene)

cat("AFTER: ", nrow(amr2), "pairs | ", n_distinct(amr2$Gene), "unique genes\n\n")

per_isolate <- amr2 %>%
  count(Isolate, name = "AMR_gene_count")

cat("Range:", min(per_isolate$AMR_gene_count),
    "-", max(per_isolate$AMR_gene_count), "\n")
cat("Median:", median(per_isolate$AMR_gene_count), "\n")
cat("Mean:", round(mean(per_isolate$AMR_gene_count), 1), "\n\n")

before <- amr %>% count(Isolate, name = "count_before")
after  <- per_isolate %>% rename(count_after = AMR_gene_count)
comp   <- inner_join(before, after, by = "Isolate") %>%
  mutate(delta = count_before - count_after)
cat("Isolates affected by residual variant merging:",
    sum(comp$delta != 0), "of", nrow(comp), "\n")
cat("Max change:", max(comp$delta), "\n\n")

write.csv(amr2,        "amr_clean_pass2_long.csv",        row.names = FALSE)
write.csv(per_isolate, "amr_clean_pass2_per_isolate.csv", row.names = FALSE)
cat("Saved: amr_clean_pass2_long.csv, amr_clean_pass2_per_isolate.csv\n")
