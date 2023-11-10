# Run analyses using only the first BMD measurement, i.e.
#   Exam 20 for Original cohort; Exam 6 for Offspring.
#
# DIRECTORIES -------------------------------------------------------------

sys.source(paste(script_dir, 'utils.R', sep = '/'), envir = attach(NULL))
set.seed(42)

# LIBRARIES ---------------------------------------------------------------

librarian::shelf(coxme,
                 dplyr,
                 finalfit,
                 haven,
                 plyr,
                 survival,
                 survminer,
                 tidyverse)
sys.source(paste(home_dir, 'utils.R', sep = '/'), envir = attach(NULL))

# DATA & PARAMETERS -------------------------------------------------------

data <- readRDS(paste(data_dir, 'SingleBMD.rds', sep = '/'))
# NOTE: 
kinshipDf <- readRDS(paste(data_dir, 'Kinship_for_SingleBMD.rds', sep = '/'))

tbl2 <- tbl3 <- tbl4 <- tblS1 <- tblS2 <- tblS3 <- data.frame(
  Table = character(),
  Model = character(),
  Survival = character(),
  Outcome = character(),
  PrimaryBMD = character(),
  Age = character(),
  APOE = character(),
  BMI = character(),
  Estrogen = character(),
  Mobility = character(),
  SecondaryBMD = character(),
  Sex = character(),
  Smoke = character(),
  stringsAsFactors = FALSE
)

tbl2 <- tbl2 %>%
  add_row(
    Table = '2',
    Model = NA,
    Survival = 'SurvTime10',
    Outcome = 'DemStatus10',
    PrimaryBMD = 'BMDInit',
    Age = 'AgeInit',
    APOE = 'APOEStatus',
    BMI = 'BMIInit',
    Estrogen = NA,
    Mobility = NA,
    SecondaryBMD = NA,
    Sex = 'Sex',
    Smoke = NA
  ) %>%
  dplyr::slice(rep(1:n(), each = 3)) %>%
  dplyr::mutate(Row = row_number(),
         Model = case_when(
           Row == 1 ~ '1.1',
           Row == 2 ~ '1.2',
           Row == 3 ~ '1.3'
         ),
         PrimaryBMD = case_when(
           Row == 2 ~ 'BMDInitLowestQuartileLabels',
           Row == 3 ~ 'BMDInitLowestQuartile',
           TRUE ~ PrimaryBMD
         )) #%>%
  # # Make a second set to model with age quartiles.
  # dplyr::slice(rep(1:n(), times = 2)) %>%
  # # Relabel row numbers
  # dplyr::mutate(Row = row_number(),
  #               Age = case_when(
  #                 Row %in% 4:6 ~ 'AgeInitQuartiles',
  #                 TRUE ~ Age
  #               ))

# Table 2 for single BMD measurements.
tbl2_results <- lapply(1:nrow(tbl2), function(row) {
  row_info <- tbl2[row, ]
  mdl <- fit_model(row_info, kinship_matrix = kinshipDf, df = data)
  tidy_output(mdl, row_info)
}) 
# saveRDS(do.call(rbind, tbl2_results), paste0(res_dir, '/results/SingleBMD_Table2.rds'))

# Table 3: Sex-stratified models.
tbl3 <- tbl2 %>%
  slice(rep(1:nrow(tbl2), each = 2)) %>%
  dplyr::mutate(
    Row = row_number(),
    Table = '3',
    Sex = case_when(Row %% 2 == 1 ~ 'Sex:Female',
                    TRUE ~ 'Sex:Male'),
    Model = case_when(Row %% 2 == 1 ~ paste(Model, 'Female'),
                      TRUE ~ paste(Model, 'Male'))
  )
tbl3_results <- lapply(1:nrow(tbl3), function(row) {
  row_info <- tbl3[row, ]
  mdl <- fit_model(row_info, kinship_matrix = kinshipDf, df = data, stratified = TRUE)
  tidy_output(mdl, row_info)
})
saveRDS(do.call(rbind, tbl3_results), paste0(res_dir, '/SingleBMD_Table3.rds'))

# Table 4: Effect modification models. 
# NOTE 06/28/22: Removed Mobility.
tbl4 <- tbl4 %>%
  add_row(
    Table = '4',
    Model = NA,
    Survival = 'SurvTime10',
    Outcome = 'DemStatus10',
    PrimaryBMD = 'BMDInit',
    Age = 'AgeInit',
    APOE = 'APOEStatus',
    BMI = 'BMIInit',
    Estrogen = NA,
    Mobility = NA,
    SecondaryBMD = NA,
    Sex = 'Sex',
    Smoke = NA
  ) %>%
  slice(rep(1:n(), times = 7)) %>%
  dplyr::mutate(Row = row_number(),
                PrimaryBMD = case_when(
                  Row == 1 ~ 'BMDInit:APOEStatus',
                  Row == 4 ~ 'BMDInit:AgeInit',
                  Row == nrow(.) ~ 'BMDInit:BMDInitLowestQuartile',
                  TRUE ~ PrimaryBMD
                ),
                APOE = case_when(
                  Row == 2 ~ 'APOEStatus:0',
                  Row == 3 ~ 'APOEStatus:1',
                  TRUE ~ APOE
                ),
                Age = case_when(
                  Row == 5 ~ 'AgeCat:Under71',
                  Row == 6 ~ 'AgeCat:AtLeast71',
                  TRUE ~ Age
                ),
                SecondaryBMD = case_when(
                  Row == 7 ~ 'BMDInitLowestQuartile',
                  TRUE ~ SecondaryBMD
                ),
                Model = c('4', '4: no E4', '4: E4', '5', '5: age < 71', '5: age >= 71', '6'))
# Bone Loss interaction models.
tmp1 <- tbl4 %>%
  dplyr::filter(grepl(':', PrimaryBMD))
# APOE- or age-stratified models.
tmp2 <- tbl4 %>%
  dplyr::filter(!grepl(':', PrimaryBMD))
tmp1_results <- lapply(1:nrow(tmp1), function(row) {
  row_info <- tmp1[row, ]
  mdl <- fit_model(row_info, kinship_matrix = kinshipDf, df = data)
  tidy_output(mdl, row_info)
})
tmp2_results <- lapply(1:nrow(tmp2), function(row) {
  row_info <- tmp2[row, ]
  mdl <- fit_model(row_info, kinship_matrix = kinshipDf, df = data, stratified = TRUE)
  tidy_output(mdl, row_info)
})
tbl4_results <- bind_rows(do.call(rbind, tmp1_results),
                          do.call(rbind, tmp2_results))
tbl4_results <- tbl4_results[order(tbl4_results$Model), ]
saveRDS(tbl4_results, paste0(res_dir, '/SingleBMD_Table4.rds'))

# PROPORTIONAL HAZARDS ----------------------------------------------------

# Check for Model 1.1
ph_mdl1p1 <- coxph(Surv(SurvTime10, DemStatus10) ~ BMDInit + AgeInit + APOEStatus + 
                     BMIInit + Sex, data = data)
cox.zph(ph_mdl1p1)
