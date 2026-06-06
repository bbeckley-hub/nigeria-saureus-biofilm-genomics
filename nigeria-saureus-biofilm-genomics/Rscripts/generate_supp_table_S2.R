# generate_supp_table_S2.R
library(dplyr)
library(rstatix)

df <- read.csv('biofilm_analysis_dataset.csv', stringsAsFactors = FALSE)
df_st <- df %>% filter(MLST != 'Not typed', MLST %in% c('ST772','ST152','ST30','ST508','ST8','ST121','ST789','ST1','ST15'))
df_st$MLST <- factor(df_st$MLST)
dunn <- df_st %>% dunn_test(biofilm_score ~ MLST, p.adjust.method = "bonferroni")
write.csv(dunn, "Supplementary_Table_S2_Dunn_test.csv", row.names = FALSE)
print("Supplementary_Table_S2_Dunn_test.csv created")
