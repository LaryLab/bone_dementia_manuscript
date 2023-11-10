invisible(sapply(list.files(
  pattern = '(directories|modules|utils)',
  recursive = TRUE,
  full.names = TRUE),
  function(f)
    sys.source(f, env = attach(NULL))))

# DATA --------------------------------------------------------------------


param_grid <- param_grid %>%
  # TABLE 2, MODEL 1
  add_row(
    Table = '2',
    Model = '1.1',
    Outcome = 'DEMSTATUS10',
    BaselineBMD = 'BMDBASE',
    Age = 'AGECONTS',
    APOE = NA,
    BMI = 'BMIBASE',
    Estrogen = NA,
    Mobility = 'HELPTRANSFERRING',
    OtherBMD = NA,
    Sex = 'SEX',
    Smoke = NA
  ) %>%
  slice(rep(1:n(), each = 8)) %>%
  mutate(Row = row_number()) %>%
  mutate(Model = case_when(Row == 2 ~ '1.2',
                           Row == 3 ~ '1.3',
                           Row == 4 ~ '2.1',
                           Row == 5 ~ '2.1*',
                           Row == 6 ~ '2.2',
                           Row == 7 ~ '2.3',
                           Row == 8 ~ '3',
                           TRUE ~ Model),
         BaselineBMD = case_when(Row == 2 ~ 'BMDBASEQUARTILE',
                                 Row == 3 ~ 'BMDBASELOWESTQUARTILE',
                                 Row == 4 ~ 'BMDPERCENTPERYEAR',
                                 Row == 5 ~ 'BMDPERCENTPERYEAR',
                                 Row == 6 ~ 'BMDPCTQUARTILE',
                                 Row == 7 ~ 'BMDPCTQUARTILE4',
                                 Row == 8 ~ 'BMDBASE',
                                 TRUE ~ BaselineBMD),
         Age = case_when(Row == 5 ~ 'AGEGROUPS',
                         TRUE ~ Age),
         OtherBMD = case_when(Row == 8 ~ 'BMDPERCENTPERYEAR',
                              TRUE ~ OtherBMD)) %>%
  slice(rep(1:n(), times = 6)) %>%
  mutate(Row = row_number(),
         Table = case_when(Row %in% 9:24 ~ 'S2',
                           Row %in% 25:32 ~ 'S3',
                           Row %in% 33:40 ~ 'S4',
                           Row %in% 41:48 ~ 'S5',
                           TRUE ~ Table),
         Outcome = case_when(Table == 'S3' ~ 'ADSTATUS10',
                             TRUE ~ Outcome),
         APOE = case_when(Table != 'S4' ~ 'APOE4STAT',
                          TRUE ~ APOE),
         Estrogen = case_when(Table == 'S5' ~ 'ESTROGEN',
                              TRUE ~ Estrogen),
         Sex = case_when(Row %in% 9:16 ~ 'Male',
                         Row %in% 17:24 ~ 'Female',
                         TRUE ~ Sex),
         Smoke = case_when(Table == 'S5' ~ 'CURRENTSMOKER',
                           TRUE ~ Smoke))
# Add interaction rows
param_grid <- bind_rows(param_grid, slice(param_grid, rep(4, each = 14))) %>%
  mutate(Row = row_number()) %>%
  mutate(Table = case_when(Row %in% 49:58 ~ '3',
                           Row %in% 59:62 ~ 'S1',
                           TRUE ~ Table),
         BaselineBMD = case_when(Row == 49 ~ 'BMDPERCENTPERYEAR*APOE4STAT',
                                 Row == 52 ~ 'BMDPERCENTPERYEAR*AGECONTS',
                                 Row == 56 ~ 'BMDPERCENTPERYEAR*BMDBASELOWESTQUARTILE',
                                 Row == 59 ~ 'BMDDELTA',
                                 Row == 60 ~ 'BMDDELTA',
                                 Row == 61 ~ 'BMDDELTANorm',
                                 Row == 62 ~ 'BMDDELTANorm',
                                 TRUE ~ BaselineBMD),
         APOE = case_when(
           Row == 50 ~ '1',
           Row == 51 ~ '0',
           TRUE ~ APOE
         ),
         Age = case_when(
           Row == 53 ~ '<= 74',
           Row == 54 ~ '74 < age <= 79',
           Row == 55 ~ '> 79',
           TRUE ~ Age
         ),
         OtherBMD = case_when(Row == 57 ~ 'BMDBASELOWESTQUARTILENo',
                              Row == 58 ~ 'BMDBASELOWESTQUARTILEYes',
                              Row %in% c(60, 62) ~ 'BMDINIT',
                              TRUE ~ OtherBMD)) %>%
  mutate(Model = case_when(
    Table == 'S1' & is.na(OtherBMD) & BaselineBMD == 'BMDDELTA' ~ '1.1 version 2',
    Table == 'S1' & !is.na(OtherBMD) & BaselineBMD == 'BMDDELTA' ~ '1.1 version 3',
    Table == 'S1' & is.na(OtherBMD) & BaselineBMD == 'BMDDELTANorm' ~ '1.1 version 2*',
    Table == 'S1' & !is.na(OtherBMD) & BaselineBMD == 'BMDDELTANorm' ~ '1.1 version 3*',
    Sex == 'Male' ~ paste(Model, 'Male', sep = ':'),
    Sex == 'Female' ~ paste(Model, 'Female', sep = ':'),
    Row == 49 ~ '4',
    Table == '3' & APOE == 0 ~ '4:noE4',
    Table == '3' & APOE == 1 ~ '4:E4',
    Row == 52 ~ '5',
    Row == 53 ~ '5:<= 74',
    Row == 54 ~ '5:74 < age <= 79',
    Row == 55 ~ '5:> 79',
    Row == 56 ~ '6',
    Row == 57 ~ '6:q1-3',
    Row == 58 ~ '6:q4',
    TRUE ~ Model,
  )) %>%
  relocate(Row, Table, Model)
# View(param_grid)
# saveRDS(param_grid, paste0(data_subdir, 'BMDParamGrid.rds'))

# HELPER FUNCTIONS --------------------------------------------------------

summary_output <- function(model_object){

  summary_object <- finalfit::fit2df(
    model_object,
    digits = c(2, 2, 3),
    explanatory_name = 'Variable',
    estimate_name = 'Hazard Ratio',
    p_name = 'p'
  )
  summary_object$Variable <- as.character(summary_object$Variable)
  summary_object$`Hazard Ratio` <- as.character(summary_object$`Hazard Ratio`)

  return(summary_object)
}

# Fit model with specified formula
fit_model <- function(row_info, kinship_matrix, df, stratified = FALSE) {
  #' @param row_info Data frame row with information on model to be fitted.
  #' @param kinship_matrix Matrix with kinship coding.
  #' @param df Data frame with Framingham data.
  #' @param stratified Default FALSE. If TRUE, then search for stratifying variable.
  #' @return Summary of results from coxme model.
  
  outcome <- row_info %>% dplyr::select(Outcome) 
  
  predictors <- row_info %>%
    dplyr::select(-Table,-Model,-Outcome,-Row)

  if (stratified == FALSE) {
    
    # If not stratifying, then just concatenate the covariates into a model formula
     predictors <- predictors %>%
      dplyr::select(where( ~ !any(is.na(.)))) %>%
      paste(collapse = ' + ')
    predictors <- paste(predictors, '(1|RANID)', sep = ' + ')
    
  }
  
  if (stratified == TRUE) {
    
    # Effect modification variables are denoted [variable]:[group]
    strat_var <- names(predictors)[grep(':', predictors)]
    print(strat_var)
    tmp <- strsplit(predictors[, strat_var], split = ':')
    var_name <- tmp[[1]][1]
    var_group <- tmp[[1]][2]
    
    # Filter data by stratum of interest
    df <- df %>%
      filter(!!as.name(var_name) == !!var_group)
    
    predictors <- predictors %>%
      dplyr::select(-!!as.name(strat_var)) %>%
      dplyr::select(where( ~ !any(is.na(.)))) %>%
      paste(collapse = ' + ')
    print(predictors)
    predictors <- paste(predictors, '(1|RANID)', sep = ' + ')
  }

  # if(row_info$Sex == 'SEX') {
  # 
  #   # Case: No APOE
  #   if (!is.na(row_info$APOE) & row_info$APOE != 'APOE4STAT') {
  # 
  #     df <- df %>% filter(APOE4STAT == as.numeric(row_info$APOE))
  #     row_info$APOE <- NA
  #   }
  #   # Case: Age Groups
  #   # if (row_info$Age %notin% c('AGEGROUPS', 'AGECONTS')) {
  #   # 
  #   #   df <- df %>% filter(AGEGROUPS == row_info$Age)
  #   #   row_info$Age <- NA
  #   # }
  # 
  #   # Case: Age Tertiles
  #   # if (row_info$Age %in% c('(60,74]', '(74,79]', '(79,90]')) {
  #   #   
  #   #   df <- df %>% filter(AGETERTILES == row_info$Age)
  #   #   row_info$Age <- NA
  #   #   print(dim(df))
  #   # }
  #   
  #   # Case: Age Median
  #   if (row_info$Age %in% c('AGEMEDIAN >= Median')) {
  #     df <- df %>% filter(AGEMEDIAN == '>= Median')
  #     row_info$Age <- NA
  #   }
  #   
  #   if (row_info$Age %in% c('AGEMEDIAN < Median')) {
  #     df <- df %>% filter(AGEMEDIAN == '< Median')
  #     row_info$Age <- NA
  #   }
  #   
  #   # Case: Age Q1-3 vs. Q4
  #   if (row_info$Age %in% c('AGEQuartile: Q4')) {
  #     df <- df %>% filter(AGEQ4 == 'Q4')
  #     row_info$Age <- NA
  #   }
  #   
  #   if (row_info$Age %in% c('AGEQuartile: Q1-3')) {
  #     df <- df %>% filter(AGEQ4 == 'Q1-3')
  #     row_info$Age <- NA
  #   }
  #   
  #   # Case: Baseline BMD quartiles 1-3
  #   if (!is.na(row_info$OtherBMD) & row_info$OtherBMD == 'BMDBASELOWESTQUARTILENo') {
  # 
  #     df <- df %>% filter(BMDBASELOWESTQUARTILE == 'No')
  #     row_info$OtherBMD <- NA_real_
  #   }
  #   # Case: Baseline BMD quartile 4
  #   if (!is.na(row_info$OtherBMD) & row_info$OtherBMD == 'BMDBASELOWESTQUARTILEYes') {
  #     df <- df %>% filter(BMDBASELOWESTQUARTILE == 'Yes')
  #     row_info$OtherBMD <- NA_real_
  #   }
  # 
  #   predictors <- row_info %>%
  #     dplyr::select(-Table,-Model,-Outcome,-Row) %>%
  #     dplyr::select(where( ~ !any(is.na(.)))) %>%
  #     paste(collapse = ' + ')
  # }
  # 
  # if(row_info$Sex == 'Male') {
  #   df <- df %>% filter(SEX == 'Male') %>% droplevels()
  #   # print(c(levels(data$SEX), dim(data)))
  #   predictors <- row_info %>%
  #     dplyr::select(-Table,-Model,-Outcome,-Sex,-Row) %>%
  #     dplyr::select(where( ~ !any(is.na(.)))) %>%
  #     paste(collapse = ' + ')
  # }
  # 
  # if (row_info$Sex == 'Female') {
  #   df <- df %>% filter(SEX == 'Female') %>% droplevels()
  #   # print(c(levels(data$SEX), dim(data)))
  #   predictors <- row_info %>%
  #     dplyr::select(-Table,-Model,-Outcome,-Sex,-Row) %>%
  #     dplyr::select(where( ~ !any(is.na(.)))) %>%
  #     paste(collapse = ' + ')
  # }
  # 
  # predictors <- paste(predictors, '(1|RANID)', sep = ' + ')
  
  formula_to_fit <- paste0(sprintf('Surv(SURVTIME10, %s) ~ ', outcome), predictors)
  print(formula_to_fit)
  formula_to_fit <- as.formula(formula_to_fit)
  # if(!is.na(age_var) && !is.na(bmd_init_var) && !is.na(bmd_final_var)){
  #   formula_to_fit <- as.formula(sprintf('Surv(surv.time10, DEM_STATUS10) ~ APOE4stat + BMIFinal + Sex + help_transferring + %s + %s + %s + (1|ranid)', 
  #                                        age_var, bmd_init_var, bmd_final_var))
  # }
  # 
  # if(!is.na(age_var) && is.na(bmd_init_var) && !is.na(bmd_final_var)){
  #   formula_to_fit <- as.formula(sprintf('Surv(surv.time10, DEM_STATUS10) ~ APOE4stat + BMIFinal + Sex + help_transferring + %s + %s + (1|ranid)', 
  #                                        age_var, bmd_final_var))
  # }
  # 
  # if(is.na(age_var) && !is.na(bmd_init_var) && !is.na(bmd_final_var)){
  #   formula_to_fit <- as.formula(sprintf('Surv(surv.time10, DEM_STATUS10) ~ APOE4stat + BMIFinal + Sex + help_transferring + %s + %s + (1|ranid)', 
  #                                        bmd_init_var, bmd_final_var))
  # }
  # 
  # if(is.na(age_var) && is.na(bmd_init_var) && !is.na(bmd_final_var)){
  #   formula_to_fit <- as.formula(sprintf('Surv(surv.time10, DEM_STATUS10) ~ APOE4stat + BMIFinal + Sex + help_transferring + %s + (1|ranid)', 
  #                                        bmd_final_var))
  # }
  
  ids_to_select <- as.character(df$RANID)
  rows_to_select <- which(rownames(kinship_matrix) %in% ids_to_select)
  cols_to_select <- which(colnames(kinship_matrix) %in% ids_to_select)
  kmat_subset <- kinship_matrix[rows_to_select, cols_to_select]

  mdl <- coxme(formula_to_fit, data = df, varlist = list(kmat_subset))
  print(mdl)
  output <- summary_output(mdl)
  print(output)
  output$Table <- rep(row_info$Table, nrow(output))
  output$Model <- rep(row_info$Model, nrow(output))
  
  output <- output %>%
    relocate(Table, Model, Variable, `Hazard Ratio`)

  return(output)
}

# MARKDOWN FORMATS --------------------------------------------------------

relabelOutput <- function(df) {
  
  # WARNING: Major hardcoding ahead!
  df$Variable <- gsub('AGECONTS', 'Age (Years)', df$Variable)
  
  df$Variable <- gsub('AGEGROUPS', 'Age Groups', df$Variable)
  
  df$Variable <- gsub('APOE4STAT', 'ApoE4 Status', df$Variable)
  
  df$Variable <- gsub('BMDINIT', 'Initial BMD (g/cm\u00B2)', df$Variable)
  
  # df$Variable <- gsub(bmd_initmedian_var, 'Baseline BMD Below Median',
  #                     df$Variable)
  
  df$Variable <- gsub('BMIBASE', 'BMI', df$Variable)
  
  df$Variable <- gsub('HELPTRANSFERRINGYes', 'Need Help Transferring (Mobility)', df$Variable)
  
  df$Variable <- gsub('BMDBASELOWESTQUARTILE',
                      'Baseline BMD quartile 4 vs. quartiles 1-3', df$Variable)
  
  df$Variable <- gsub('BMDBASEQUARTILE', 'Baseline BMD Quartile', df$Variable)
  
  df$Variable <- gsub('BMDPCTQUARTILE4Yes', 'Bone Loss quartile 4 vs. quartiles 1-3', df$Variable)
  
  df$Variable <- gsub('BMDPCTQUARTILE', 'Bone Loss Quartile', df$Variable)
  
  df$Variable <- gsub('BMDPERCENTPERYEAR', 'Bone Loss', df$Variable)
  
  df$Variable <- gsub('BMDBASE', 'Baseline BMD (g/cm\u00B2)', df$Variable)
  
  df$Variable <- gsub('CURRENTSMOKER', 'Current Smoker', df$Variable)
  
  df$Variable <- gsub('SEXMale', 'Male', df$Variable)
  
  df$Variable <- gsub('ESTROGEN', 'Female Estrogen Use', df$Variable)
  
  # Interaction variables
  df$Variable <- gsub(':', ' x ', df$Variable)
  
  # word_path <- '/home/ckhoo@mmcf.mehealth.org/BMD/WordDocs/'
  # df %>%
  #   flextable::as_grouped_data(groups = NULL) %>%
  #   flextable::as_flextable() %>%
  #   flextable::padding() %>%
  #   flextable::autofit() %>%
  #   flextable::font(fontname = 'Arial', part = 'all') %>%
  #   save_as_docx(path = paste0(word_path, 'BMDTables-', Sys.Date(), '.docx'))
  
  return(df)
}

# datatable(output_df)

word_path <- '/home/ckhoo@mmcf.mehealth.org/BMD/WordDocs/'

initial_table <- function(df) {
  df %>%
    flextable::as_grouped_data(groups = NULL) %>%
    flextable::as_flextable() %>%
    flextable::padding() %>%
    flextable::autofit() %>%
    flextable::align_text_col(align = 'center') %>%
    flextable::font(fontname = 'Arial', part = 'all')
}

border_style <- fp_border(style = "solid", width=2)

# FIT MODELS --------------------------------------------------------------

param_grid <- readRDS(paste0(data_subdir, 'BMDParamGrid.rds'))
# # Kinship matrix
kinship_df <- readRDS(paste0(data_subdir, 'KinshipMatrix.rds'))
data <- readRDS(paste0(data_subdir, 'BMDKinship.rds'))

# NOTE: 03/18/22 - Test sex interaction effects for models in Table 2.
#   Once models are finalized, we would want to modify param_grid to reflect
#   actual models that went into the manuscript for reproducibility purposes.
# subset_grid <- param_grid %>%
#   filter(Table == 2 & !Model == '2.1*') %>%
#   mutate(BaselineBMD = paste(BaselineBMD, 'SEX', sep = '*'),
#          OtherBMD = case_when(
#            Model == 3 ~ 'BMDPERCENTPERYEAR*SEX',
#            TRUE ~ OtherBMD
#          ))
# results <- lapply(1:nrow(subset_grid), function(row) fit_model(subset_grid[row, ], kinship_df, data))
# results_df <- do.call(rbind, results)
# View(results_df)
# saveRDS(results_df, paste0(res_subdir, 'BMDKinshipResults-Tbl2SexInteraction.rds'))

# Save as Word doc
# outputDF <- relabelOutput(results_df %>% filter(grepl('BMD', Variable)))
# outputDF %>%
#   flextable::as_grouped_data(groups = NULL) %>%
#   flextable::as_flextable() %>%
#   flextable::padding() %>%
#   flextable::autofit() %>%
#   flextable::font(fontname = 'Arial', part = 'all') %>%
#   save_as_docx(path = paste0(word_path, 'BMDTables-SexInteraction.docx'))
# View(outputDF)


# MODEL FOR SURVIVAL CURVE ------------------------------------------------

# subset <- data %>%
#   dplyr::select(SURVTIME10,
#     DEMSTATUS10,
#     BMDPCTQUARTILE4,
#     AGECONTS,
#     APOE4STAT,
#     BMIBASE,
#     HELPTRANSFERRING,
#     SEX
#   ) %>%
#   dplyr::rename(Age = AGECONTS,
#                 ApoE4 = APOE4STAT,
#                 BMI = BMIBASE,
#                 Mobility = HELPTRANSFERRING,
#                 Sex = SEX) %>%
#   mutate(BMDPCTQUARTILE4 = recode_factor(
#     BMDPCTQUARTILE4,
#     No = 'Bone Loss Quartiles 1-3',
#     Yes = 'Bone Loss Quartile 4'
#   ))
# dd <- datadist(subset)
# options(datadist = 'dd')
# 
# get_surv_curv <- function(){
# 
#   row_info <- param_grid %>%
#     filter(Table == '2' & BaselineBMD == 'BMDPCTQUARTILE4')
# 
#   outcome <- row_info %>% dplyr::select(Outcome)
# 
#   predictors <- row_info %>%
#     dplyr::mutate(Age = case_when(Age == 'AGECONTS' ~ 'Age', TRUE ~ Age),
#                   APOE = case_when(APOE == 'APOE4STAT' ~ 'ApoE4', TRUE ~ APOE),
#                   BMI = case_when(BMI == 'BMIBASE' ~ 'BMI', TRUE ~ BMI),
#                   Mobility = case_when(Mobility == 'HELPTRANSFERRING' ~ 'Mobility',
#                                        TRUE ~ Mobility),
#                   Sex = case_when(Sex == 'SEX' ~ 'Sex', TRUE ~ Sex)) %>%
#     dplyr::select(-Table,-Model,-Outcome,-Row) %>%
#     dplyr::select(where( ~ !any(is.na(.)))) %>%
#     paste(collapse = ' + ')
# 
#   # data <- data %>%
#   #   mutate(BMDPCTQUARTILE4 = ifelse(BMDPCTQUARTILE4 == 'Yes', 'Bone Loss: Quartile 4',
#   #                                   'Bone Loss: Quartiles 1-3'))
# 
#   formula_to_fit <- paste0(sprintf('Surv(SURVTIME10, %s) ~ ', outcome), predictors)
#   # formula_to_fit <- gsub('BMDPCTQUARTILE4', 'strata(BMDPCTQUARTILE4)', formula_to_fit)
#   print(formula_to_fit)
#   formula_to_fit <- as.formula(formula_to_fit)
# 
#   mdl <- cph(formula_to_fit, data = subset,
#            x = TRUE,
#            y = TRUE,
#            surv = TRUE,
#            model = TRUE)
#   # mdl <- coxph(formula = formula_to_fit, data = data)
#   # mdl <- survfit(mdl, data = data)
#   return(mdl)
# }
# 
# mdl_object <- get_surv_curv()
# 
# # png('./Figures/SurvivalCurves.png', width = 720, height = 720, units = 'px')
# # survplot(
# #   mdl_object,
# #   BMDPCTQUARTILE4,
# #   label.curves = list(keys = 'lines',
# #                       keyloc = c(1500, .775),
# #                       cex = 1.15, transparent = TRUE),
# #   col = c('red', 'dodgerblue4'),
# #   lty = 1,
# #   lwd = 2,
# #   col.fill = c('lightpink', 'cadetblue1'),
# #   conf.int = .95,
# #   conf = 'bands',
# #   ylim = c(.75, 1)
# # )
# # dev.off()
# 

# TMP ---------------------------------------------------------------------

ageSummary <- summary(data$AGECONTS)
# data$AGETERTILES <- cut(data$AGECONTS, breaks = c(ageSummary[['Min.']],
#                                                   ageSummary[['Median']],
#                                                   ageSummary[['3rd Qu.']],
#                                                   ageSummary[['Max.']]))
# 
# with(data, table(DEMSTATUS10, AGETERTILES))
# with(data, prop.table(table(DEMSTATUS10, AGETERTILES)))
# data$AGEMEDIAN <- ifelse(data$AGECONTS >= ageSummary[['Median']],
#                          '>= Median', '< Median')
data$AGEQuartile <- ifelse(data$AGECONTS >= ageSummary[['3rd Qu.']],
                     'Q4', 'Q1-3')

# tbl <- param_grid %>%
#   filter(Table == 3 & grepl('5$', Model)) %>%
#   slice(rep(1:n(), times = 4)) %>%
#   mutate(Row = row_number(),
#          Model = case_when(Row == 1 ~ '5: >= Median',
#                            Row == 2 ~ '5: < Median',
#                            Row == 3 ~ '5: Age Q4',
#                            Row == 4 ~ '5: Age Q1-3'),
#          Age = case_when(Row == 1 ~ 'AGEMEDIAN >= Median',
#                          Row == 2 ~ 'AGEMEDIAN < Median',
#                          Row == 3 ~ 'AGEQuartile: Q4',
#                          Row == 4 ~ 'AGEQuartile: Q1-3'),
#          BaselineBMD = 'BMDPERCENTPERYEAR')
# #   mutate(
# #     BaselineBMD = case_when(
# #       BaselineBMD == 'BMDPERCENTPERYEAR*AGECONTS' ~ 'BMDPERCENTPERYEAR*AGETERTILES',
# #       TRUE ~ BaselineBMD
# #     ),
# #     Age = case_when(
# #       Age == 'AGECONTS' ~ 'AGETERTILES',
# #       Age == '<= 74' ~ '(60,74]',
# #       Age == '74 < age <= 79' ~ '(74,79]',
# #       Age == '> 79' ~ '(79,90]',
# #       TRUE ~ Outcome
# #     ),
# #     Model = case_when(
# #       Model == '5:<= 74' ~ '5: (60,74]',
# #       Model == '5:74 < age <= 79' ~ '5: (74,79]',
# #       Model == '5:> 79' ~ '5: (79,90]',
# #       TRUE ~ Model
# #     )
# #   )
# # View(tbl)
# 
outputs <- vector('list', length=2)
for(row in 3:4){

  outputs[[row]] <- fit_model(tbl[row, ], kinship_df, data)
}
# View(outputs)
#
outputDF <- do.call(rbind, outputs)
outputDF <- relabelOutput(outputDF)
# outputDF %>%
#   flextable::as_grouped_data(groups = NULL) %>%
#   flextable::as_flextable() %>%
#   flextable::padding() %>%
#   flextable::autofit() %>%
#   flextable::font(fontname = 'Arial', part = 'all') %>%
#   save_as_docx(path = paste0(word_path, 'BMDTables-AgeMedian_and_Q4.docx'))
# View(outputDF)

# saveRDS(outputDF, 'BMDKinshipResults-AgeQ4.rds')

tbl2_wout_mobility <- tbl2 %>%
  mutate(Mobility = NA)

tbl3_wout_mobility <- tbl3 %>%
  mutate(Mobility = NA)

tbl4_wout_mobility_nostrat <- tbl4 %>%
  mutate(Mobility = NA) %>%
  filter(Row %in% c(1, 4, 7))

tbl4_wout_mobility_strat <- tbl4 %>%
  mutate(Mobility = NA) %>%
  filter(Row %in% c(2:3, 5:6))

tbl2res <- lapply(4:4, function(row) fit_model(tbl2_wout_mobility[row, ], kinship_df, data))

results <- lapply(1:nrow(tbl2), function(row) fit_model(tbl2_wout_mobility[row, ], kinship_df, data))
results_df <- do.call(rbind, results)
# View(results_df)

res3 <- lapply(1:nrow(tbl3), function(row) fit_model(tbl3_wout_mobility[row, ], 
                                                     kinship_df, data,
                                                     stratified = TRUE))
resdf3 <- do.call(rbind, res3)
# View(resdf3)

res4_p1 <-
  lapply(1:nrow(tbl4_wout_mobility_strat), function(row)
    fit_model(tbl4_wout_mobility_strat[row,],
              kinship_df, data,
              stratified = T))
res4_p2 <-
  lapply(1:nrow(tbl4_wout_mobility_nostrat), function(row)
    fit_model(tbl4_wout_mobility_nostrat[row,],
              kinship_df, data,
              stratified = F))
res4 <- rbind(do.call(rbind, res4_p1), do.call(rbind, res4_p2))

nomobility <- rbind(results_df, resdf3, res4) %>% arrange(Table, Model)

# nomobility %>%
#   flextable() %>%
#   save_as_docx(path = './WordDocs/Tables2-4_WithoutMobility.docx')
