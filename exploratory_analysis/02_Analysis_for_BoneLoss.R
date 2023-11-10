# Run analyses using only the first BMD measurement, i.e.
#   Exam 20 for Original cohort; Exam 6 for Offspring.
# SETUP -------------------------------------------------------------------

sys.source(paste(script_dir, 'utils.R', sep = '/'), envir = attach(NULL))
set.seed(42)

# DATA & PARAMETERS -------------------------------------------------------

data <- readRDS(paste(data_dir, 'BoneLoss.rds', sep = '/'))
kinshipDf <- readRDS(paste(data_dir, 'Kinship_for_BoneLoss.rds', sep = '/'))

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

# NOTE 06/28/22: Mobility is removed from Table 2.
tbl2 <- tbl2 %>%
  add_row(
    Table = '2',
    Model = NA,
    Survival = 'SurvTime10',
    Outcome = 'DemStatus10',
    PrimaryBMD = 'BMDBase',
    Age = 'AgeBase',
    APOE = 'APOEStatus',
    BMI = 'BMIBase',
    Estrogen = NA,
    Mobility = NA,
    SecondaryBMD = NA,
    Sex = 'Sex',
    Smoke = NA
  ) %>%
  slice(rep(1:nrow(.), each = 7)) %>%
  dplyr::mutate(
    Row = row_number(),
    Model = case_when(
      Row == 1 ~ '1.1',
      Row == 2 ~ '1.2',
      Row == 3 ~ '1.3',
      Row == 4 ~ '2.1',
      Row == 5 ~ '2.2',
      Row == 6 ~ '2.3',
      Row == 7 ~ '3'
    ),
    PrimaryBMD = case_when(
      Row == 2 ~ 'BMDBaseHighestQuartileLabels',
      Row == 3 ~ 'BMDBaseQ4',
      Row == 4 ~ 'BoneLoss',
      Row == 5 ~ 'BoneLossLowestQuartileLabels',
      Row == 6 ~ 'BoneLossQ4',
      TRUE ~ PrimaryBMD
    ),
    SecondaryBMD = case_when(Row == 7 ~ 'BoneLoss',
                             TRUE ~ SecondaryBMD)
  ) %>%
  dplyr::slice(rep(1:nrow(.), times = 3)) %>%
  dplyr::mutate(
    Row = row_number(),
    Model = case_when(
      Row %in% 1:7 ~ paste(Model, '- 10-year FU'),
      Row %in% 8:14 ~ paste(Model, '- 5-year FU'),
      Row %in% 15:nrow(.) ~ paste(Model, '- 3-year FU')
    ),
    Survival = case_when(Row %in% 8:14 ~ 'SurvTime5',
                         Row %in% 15:nrow(.) ~ 'SurvTime5',
                         TRUE ~ Survival),
    Outcome = case_when(Row %in% 8:14 ~ 'DemStatus3',
                        Row %in% 15:nrow(.) ~ 'DemStatus3',
                        TRUE ~ Outcome)
  )

tbl2_results <- lapply(1:nrow(tbl2), function(row) {
  row_info <- tbl2[row, ]
  mdl <- fit_model(row_info, kinship_matrix = kinshipDf, df = data)
  tidy_output(mdl, row_info)
})
tbl2_results <- do.call(rbind, tbl2_results)
saveRDS(tbl2_results, paste0(res_dir, '/BoneLoss_Table2.rds'))

# Exactly the same as above, but saving models without tidying to compute SEs for meta-analysis.
# tbl2_results <- lapply(1:nrow(tbl2), function(row) {
#   row_info <- tbl2[row, ]
#   mdl <- fit_model(row_info, kinship_matrix = kinshipDf, df = data)
# })
# saveRDS(tbl2_results, paste0(res_dir, '/BoneLoss_Table2_Models.rds'))

# Table 2 with age quartiles
# tbl2_age_quartiles <- tbl2 %>%
#   filter(Row %in% 1:7) %>%
#   mutate(Age = case_when(
#     Age == 'AgeBase' ~ 'AgeBaseQuartiles',
#     TRUE ~ Age
#   ))
# tbl2_age_tertiles <- tbl2_age_quartiles %>%
#   mutate(Age = 'AgeBaseTertiles')
# tbl2_age_tertiles_results <-
#   lapply(1:nrow(tbl2_age_tertiles), function(row)
#     fit_model(
#       tbl2_age_tertiles[row,],
#       kinship_matrix = kinshipDf,
#       df = data
#     ))
# saveRDS(do.call(rbind, tbl2_age_tertiles_results), paste0(home_dir, '/results/BoneLoss_Table2_AgeTertiles.rds'))

# Table 3: Sex-stratified models.
tbl3 <- tbl2 %>%
  slice(rep(1:7, each = 2)) %>%
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
saveRDS(do.call(rbind, tbl3_results), paste0(res_dir, '/BoneLoss_Table3.rds'))

# Table 4: Effect modification models. 
# NOTE 06/28/22: Removed Mobility.
tbl4 <- tbl4 %>%
  add_row(
    Table = '4',
    Model = NA,
    Survival = 'SurvTime10',
    Outcome = 'DemStatus10',
    PrimaryBMD = 'BoneLoss',
    Age = 'AgeBase',
    APOE = 'APOEStatus',
    BMI = 'BMIBase',
    Estrogen = NA,
    Mobility = NA,
    SecondaryBMD = NA,
    Sex = 'Sex',
    Smoke = NA
  ) %>%
  slice(rep(1:n(), times = 7)) %>%
  dplyr::mutate(Row = row_number(),
         PrimaryBMD = case_when(
           Row == 1 ~ 'BoneLoss*APOEStatus',
           Row == 4 ~ 'BoneLoss*AgeBase',
           Row == nrow(.) ~ 'BoneLoss*BMDBaseQ4',
           TRUE ~ PrimaryBMD
         ),
         APOE = case_when(
           Row == 2 ~ 'APOEStatus:0',
           Row == 3 ~ 'APOEStatus:1',
           TRUE ~ APOE
         ),
         Age = case_when(
           Row == 5 ~ 'AgeCat:Under79',
           Row == 6 ~ 'AgeCat:AtLeast79',
           TRUE ~ Age
         ),
         SecondaryBMD = case_when(
           Row == 7 ~ 'BMDBaseQ4',
           TRUE ~ SecondaryBMD
         ),
         Model = c('4', '4: no E4', '4: E4', '5', '5: age < 79', '5: age >= 79', '6'))
# Bone Loss interaction models.
tmp1 <- tbl4 %>%
  dplyr::filter(grepl('\\*', PrimaryBMD))
# APOE- or age-stratified models.
tmp2 <- tbl4 %>%
  dplyr::filter(!grepl('\\*', PrimaryBMD))
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
saveRDS(tbl4_results, paste0(home_dir, '/results/BoneLoss_Table4.rds'))

coxme(Surv(SurvTime10, DemStatus10) ~ BoneLoss:APOEStatus + AgeBase + 
        APOEStatus + BMIBase + Sex + Mobility + (1|RANID), data = data, varlist = list(kinshipDf))

coxme(Surv(SurvTime10, DemStatus10) ~ BoneLoss*APOEStatus + AgeBase + 
        APOEStatus + BMIBase + Sex + Mobility + (1|RANID), data = data, varlist = list(kinshipDf))

# Table S3: Table 2 with AD as outcome.
tblS3 <- tbl2 %>%
  mutate(Outcome = 'ADStatus10')
tblS3_results <- lapply(1:nrow(tblS3), function(row) {
  row_info <- tblS3[row, ]
  mdl <- fit_model(row_info, kinship_matrix = kinshipDf, df = data)
  tidy_output(mdl, row_info)
})
saveRDS(do.call(rbind, tblS3_results), paste0(res_dir, '/BoneLoss_TableS3.rds'))

# Table S4: Alternate parametrization of BMD.
tblS4 <- tbl2[1:2,] %>%
  mutate(
    PrimaryBMD = case_when(Row == 1 ~ 'DeltaBMDStd',
                           Row == 2 ~ 'DeltaBMDStd'),
    SecondaryBMD = case_when(Row == 2 ~ 'BMDInit',
                             TRUE ~ SecondaryBMD)
  )
tblS4_results <- lapply(1:nrow(tblS4), function(row) {
  row_info <- tblS4[row, ]
  mdl <- fit_model(row_info, kinship_matrix = kinshipDf, df = data)
  tidy_output(mdl, row_info)
})
saveRDS(do.call(rbind, tblS4_results), paste0(res_dir, '/BoneLoss_TableS4.rds'))

# Table S5: Table 2 with mobility and estrogen.
tblS5 <- tbl2 %>%
  mutate(Sex = NA,
         Estrogen = 'EstrogenUse',
         Mobility = 'Mobility',
         Smoke = 'CurrentSmoker') %>%
  dplyr::filter(grepl('(1.1|2.1).*10-year', Model))
tblS5_results <- lapply(1:nrow(tblS5), function(row) {
  row_info <- tblS5[row, ]
  mdl <- fit_model(row_info, kinship_matrix = kinshipDf, df = data)
  tidy_output(mdl, row_info)
})
saveRDS(do.call(rbind, tblS5_results), paste0(res_dir, '/BoneLoss_TableS5.rds'))

# MISCELLANEOUS MODELS ----------------------------------------------------

tmpfit <- coxme(Surv(SurvTime10, DemStatus10) ~ BoneLoss:APOEStatus + AgeBase + 
                  APOEStatus + BMIBase + Sex + (1 | RANID), 
                data = data, varlist = list(kinshipDf))

misc <- tbl2[1, ] %>%
  mutate(PrimaryBMD = 'BMDBase*APOEStatus')
misc_results <- lapply(1:nrow(misc), function(row) {
  row_info <- misc[row, ]
  mdl <- fit_model(row_info, kinship_matrix = kinshipDf, df = data)
  tidy_output(mdl, row_info)
})
misc_results

# Table 4 with baseline BMD instead of bone loss.
misc <- tbl4 %>% 
  dplyr::filter(grepl('4', Model)) %>%
  dplyr::mutate(PrimaryBMD = gsub('BoneLoss', 'BMDBase', PrimaryBMD),
                PrimaryBMD = gsub(':', '*', PrimaryBMD))
misc_results <- lapply(2:3, function(row) {
  row_info <- misc[row, ]
  mdl <- fit_model(row_info, kinship_matrix = kinshipDf, df = data, stratified = TRUE)
  tidy_output(mdl, row_info)
})
View(do.call(rbind, misc_results))

p1 <- ggplot(data[!is.na(data$APOEStatus),], aes(AgeBase, BoneLoss, color = factor(DemStatus10))) + 
  geom_point() + 
  theme(legend.position = "top") + 
  labs(color = "Dementia Status (10-year)") + 
  facet_grid(APOEStatus ~ Sex, labeller = label_both)
             
p2 <- ggplot(data[!is.na(data$APOEStatus),], aes(AgeBase, BMDBase, color = factor(DemStatus10))) + 
  geom_point() + 
  theme(legend.position = "top") +
  labs(color = "Dementia Status (10-year)") +
  facet_grid(APOEStatus ~ Sex, labeller = label_both)

p3 <- ggplot(data[!is.na(data$APOEStatus),], aes(SurvTime10, BoneLoss, color = factor(DemStatus10))) + 
  geom_point() + 
  theme(legend.position = "top") + 
  labs(color = "Dementia Status (10-year)") + 
  facet_grid(APOEStatus ~ Sex, labeller = label_both)

p4 <- ggplot(data[!is.na(data$APOEStatus),], aes(SurvTime10, BMDBase, color = factor(DemStatus10))) + 
  geom_point() + 
  theme(legend.position = "top") + 
  labs(color = "Dementia Status (10-year)") + 
  facet_grid(APOEStatus ~ Sex, labeller = label_both)

png("plots.png", res = 300, height = 15, width = 20, units = 'in')
ggarrange(plotlist = list(p1, p2))
dev.off()
