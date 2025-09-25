##==Read in pheno file for adipose cohort=======================================
##==============================================================================
pheno_data_adipose <- read.delim("M:\\VKraus_Lab\\Melissa.Orenduff.orend001\\Belsky_Adipose\\New_checked_pheno\\Adipose_miRNA_pheno_complete_checked_ordered_125.txt", sep="\t", header=T)
dim(pheno_data_adipose)  #[1] 125  45
head(rownames(pheno_data_adipose))
rownames(pheno_data_adipose) <- pheno_data_adipose$NAME
colnames(pheno_data_adipose)
View(pheno_data_adipose)

##======read in count data frame for adipose cohort - saved as csv==============
##==============================================================================
count_adipose <- read_csv("M:\\Melissa.Orenduff.orend001\\Belsky_Adipose\\New_checked_pheno\\Adipose_combined_transposed_LM_version.csv")
dim(count_adipose)  #[1] 125 993
View(count_adipose)
head(count_adipose)
colnames(count_adipose)

##==Setting up 12 month data for linear model===================================
##==============================================================================
count_adipose[1:5, 990:993,]
colnames(count_adipose)
count_12m_bl_adipose <- dplyr::filter(count_adipose, Metabolon.Timepoint %in% c("Month 12", "Baseline")) %>%
  dplyr::group_by(DEIDNUM) %>%
  dplyr::mutate(n_tpts = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n_tpts > 1)
View(count_12m_bl_adipose)
unique_deidnums <- unique(count_12m_bl_adipose$DEIDNUM)
mirna_cols <- 2:991
log_12m_array <- matrix(NA, nrow = length(unique_deidnums), ncol = length(mirna_cols)) %>%
  magrittr::set_rownames(unique_deidnums) %>%
  magrittr::set_colnames(colnames(count_12m_bl_adipose)[mirna_cols])
View(log_12m_array)
#***
for(deid in unique_deidnums){
  print(deid)
  timepoint_col <- count_12m_bl_adipose$Metabolon.Timepoint[count_12m_bl_adipose$DEIDNUM == deid]
  patient_data <- count_12m_bl_adipose[count_12m_bl_adipose$DEIDNUM == deid, mirna_cols] %>% as.matrix() %>% apply(2, as.numeric)
  log_12m_array[paste0("", deid), ] <- log2(patient_data[timepoint_col == "Month 12", ]+1)


}

log_bl_array <- matrix(NA, nrow = length(unique_deidnums), ncol = length(mirna_cols)) %>%
  magrittr::set_rownames(unique_deidnums) %>%
  magrittr::set_colnames(colnames(count_12m_bl_adipose)[mirna_cols])
for(deid in unique_deidnums){
  print(deid)
  timepoint_col <- count_12m_bl_adipose$Metabolon.Timepoint[count_12m_bl_adipose$DEIDNUM == deid]
  patient_data <- count_12m_bl_adipose[count_12m_bl_adipose$DEIDNUM == deid, mirna_cols] %>% as.matrix() %>% apply(2, as.numeric)
  log_bl_array[paste0("", deid), ] <- log2(patient_data[timepoint_col == 'Baseline', ]+1)


}


##==========================================================================================
## linear regression model so that control variables can be added:
predictors <- pheno_data_adipose %>%
  dplyr::filter(DEIDNUM %in% row.names(log_12m_array), Metabolon.Timepoint == "Month 12") %>%
  dplyr::select(DEIDNUM, GENDER, AGEBL, DEIDSITE, RACE, BMI_strata_BL, TREATMENT, pctcr, pdeltawt) %>%
  dplyr::distinct() %>%
  dplyr::mutate(DEIDSITE = as.factor(DEIDSITE),
                TREATMENT = factor(TREATMENT, levels=c("AL", "CR")))
## here, predictors$DEIDNUM == row.names(log_12m_array) should be all TRUE
lm_adipose_results_12_tx_4 <- lapply(1:ncol(log_12m_array), function(j){
  x <- log_12m_array[, j]
  lm_df <- cbind(log_ratio=x, predictors, BASELINE=log_bl_array[, j])
  lm_res <- lm(log_ratio~GENDER+AGEBL+DEIDSITE+RACE+BMI_strata_BL+TREATMENT+pdeltawt+BASELINE, data=lm_df)
  lm_summary <- summary(lm_res)$coefficients
  lm_summary_df <- cbind(miRNA = colnames(log_12m_array)[j],
                         coefficient = row.names(lm_summary),
                         data.frame(lm_summary))
}) %>%
  data.table::rbindlist() %>%
  magrittr::set_colnames(c("miRNA", "coefficient", "Estimate", "Std.Error", "t.value", "p.value")) %>%
  dplyr::mutate(ci95_lwr = Estimate - 1.96*`Std.Error`, ci95_upr = Estimate + 1.96*`Std.Error`) %>%
  dplyr::mutate(adjusted.p.value = p.adjust(p.value, method = "BH"))


View(lm_results_12_pctcr)

write.table(lm_adipose_results_12_tx_4, file = "M:\\VKraus_Lab\\Melissa.Orenduff.orend001\\Belsky_Adipose\\Results_linear_model\\Adipose_12_months\\Adipose_treatment_12_model_4.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)

# to add BMI need to create new column with % change (12-BL)/BL



##==Setting up 24 month data for linear model===================================
##==============================================================================
count_adipose[1:5, 990:993,]
colnames(count_adipose)
count_24m_bl_adipose <- dplyr::filter(count_adipose, Metabolon.Timepoint %in% c("Month 24", "Baseline")) %>%
  dplyr::group_by(DEIDNUM) %>%
  dplyr::mutate(n_tpts = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n_tpts > 1)
View(count_24m_bl_adipose)
unique_deidnums <- unique(count_24m_bl_adipose$DEIDNUM)
mirna_cols <- 2:991
log_24m_array <- matrix(NA, nrow = length(unique_deidnums), ncol = length(mirna_cols)) %>%
  magrittr::set_rownames(unique_deidnums) %>%
  magrittr::set_colnames(colnames(count_24m_bl_adipose)[mirna_cols])
View(log_24m_array)
#***
for(deid in unique_deidnums){
  print(deid)
  timepoint_col <- count_24m_bl_adipose$Metabolon.Timepoint[count_24m_bl_adipose$DEIDNUM == deid]
  patient_data <- count_24m_bl_adipose[count_24m_bl_adipose$DEIDNUM == deid, mirna_cols] %>% as.matrix() %>% apply(2, as.numeric)
  log_24m_array[paste0("", deid), ] <- log2(patient_data[timepoint_col == "Month 24", ]+1)


}

log_bl_array <- matrix(NA, nrow = length(unique_deidnums), ncol = length(mirna_cols)) %>%
  magrittr::set_rownames(unique_deidnums) %>%
  magrittr::set_colnames(colnames(count_24m_bl_adipose)[mirna_cols])
for(deid in unique_deidnums){
  print(deid)
  timepoint_col <- count_24m_bl_adipose$Metabolon.Timepoint[count_24m_bl_adipose$DEIDNUM == deid]
  patient_data <- count_24m_bl_adipose[count_24m_bl_adipose$DEIDNUM == deid, mirna_cols] %>% as.matrix() %>% apply(2, as.numeric)
  log_bl_array[paste0("", deid), ] <- log2(patient_data[timepoint_col == 'Baseline', ]+1)


}


##==========================================================================================
## linear regression model so that control variables can be added - removed race b/c only white are included in this cohort and model would not work
predictors <- pheno_data_adipose %>%
  dplyr::filter(DEIDNUM %in% row.names(log_24m_array), Metabolon.Timepoint == "Month 24") %>%
  dplyr::select(DEIDNUM, GENDER, AGEBL, DEIDSITE, BMI_strata_BL, TREATMENT, pctcr, pdeltawt) %>%
  dplyr::distinct() %>%
  dplyr::mutate(DEIDSITE = as.factor(DEIDSITE),
                TREATMENT = factor(TREATMENT, levels=c("AL", "CR")))
## here, predictors$DEIDNUM == row.names(log_12m_array) should be all TRUE
lm_adipose_results_24_pctcr_model_4 <- lapply(1:ncol(log_24m_array), function(j){
  x <- log_24m_array[, j]
  lm_df <- cbind(log_ratio=x, predictors, BASELINE=log_bl_array[, j])
  lm_res <- lm(log_ratio~GENDER+AGEBL+DEIDSITE+BMI_strata_BL+pctcr+pdeltawt+BASELINE, data=lm_df)
  lm_summary <- summary(lm_res)$coefficients
  lm_summary_df <- cbind(miRNA = colnames(log_24m_array)[j],
                         coefficient = row.names(lm_summary),
                         data.frame(lm_summary))
}) %>%
  data.table::rbindlist() %>%
  magrittr::set_colnames(c("miRNA", "coefficient", "Estimate", "Std.Error", "t.value", "p.value")) %>%
  dplyr::mutate(ci95_lwr = Estimate - 1.96*`Std.Error`, ci95_upr = Estimate + 1.96*`Std.Error`) %>%
  dplyr::mutate(adjusted.p.value = p.adjust(p.value, method = "BH"))


View(lm_adipose_results_24_pctcr)

write.table(lm_adipose_results_24_pctcr_model_4, file = "M:\\VKraus_Lab\\Melissa.Orenduff.orend001\\Belsky_Adipose\\Results_linear_model\\Adipose_24_months\\Adipose_pctcr_24_model_4.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)
