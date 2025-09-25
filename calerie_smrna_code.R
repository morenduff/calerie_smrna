##==Read in pheno
pheno_data <- read.delim(file)

##==read in count data 
count <- read_csv(file)

##==Setting up data for linear model===================================
##==============================================================================
count_data <- dplyr::filter(count,variable %in% c("variable name", "variable name")) %>%
  dplyr::group_by(variable name) %>%
  dplyr::mutate(n_tpts = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n_tpts > 1)

unique_deidnums <- unique(file name$variable name)

log_array <- matrix(NA, nrow = length(unique_deidnums), ncol = length(mirna_cols)) %>%
  magrittr::set_rownames(unique_deidnums) %>%
  magrittr::set_colnames(colnames(count)[mirna_cols])

for(deid in unique_deidnums){
  print(deid)
  timepoint_col <- count$X[count$X == Y]
  patient_data <- count[name$variable name == deid, mirna_cols] %>% as.matrix() %>% apply(2, as.numeric)
  log_12m_array[paste0("", deid), ] <- log2(patient_data[timepoint_col == "variable name", ]+1)
}

log_array <- matrix(NA, nrow = length(unique_deidnums), ncol = length(mirna_cols)) %>%
  magrittr::set_rownames(unique_deidnums) %>%
  magrittr::set_colnames(colnames(count_12m_bl_adipose)[mirna_cols])
for(deid in unique_deidnums){
  print(deid)
  timepoint_col <- count[name$variable name == deid]
  patient_data <- count [name$variable name == deid, mirna_cols] %>% as.matrix() %>% apply(2, as.numeric)
  log_bl_array[paste0("", deid), ] <- log2(patient_data[timepoint_col == 'variable name', ]+1)
}


##==========================================================================================
## linear regression model so that control variables can be added:
predictors <- pheno %>%
  dplyr::filter(varible name %in% row.names(log_12m_array), varible == "variable name") %>%
  dplyr::select(variables) %>%
  dplyr::distinct() %>%
  dplyr::mutate(variable name = as.factor(variable name),
                TREATMENT = factor(variable name, levels=c("X", "X")))
## here, predictors$DEIDNUM == row.names(log_12m_array) should be all TRUE
lm_results <- lapply(1:ncol(log_12m_array), function(j){
  x <- log_12m_array[, j]
  lm_df <- cbind(log_ratio=x, predictors, BASELINE=log_bl_array[, j])
  lm_res <- lm(predictors;log_ratio~GENDER+AGEBL+DEIDSITE+RACE+BMI_strata_BL+TREATMENT+pdeltawt+BASELINE, data=lm_df)
  lm_summary <- summary(lm_res)$coefficients
  lm_summary_df <- cbind(miRNA = colnames(log_12m_array)[j],
                         coefficient = row.names(lm_summary),
                         data.frame(lm_summary))
}) %>%
  data.table::rbindlist() %>%
  magrittr::set_colnames(c("miRNA", "coefficient", "Estimate", "Std.Error", "t.value", "p.value")) %>%
  dplyr::mutate(ci95_lwr = Estimate - 1.96*`Std.Error`, ci95_upr = Estimate + 1.96*`Std.Error`) %>%
  dplyr::mutate(adjusted.p.value = p.adjust(p.value, method = "BH"))


