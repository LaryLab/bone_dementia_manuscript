# Run analyses using only the first BMD measurement, i.e.
#   Exam 20 for Original cohort; Exam 6 for Offspring.
# SETUP -------------------------------------------------------------------

#script_dir <- "/cluster/tufts/larylabdbgap/ckhoo01/cobre-dementia/scripts"
#data_dir <- "/cluster/tufts/larylabdbgap/ckhoo01/cobre-dementia/data"
sys.source(paste(script_dir, 'utils.R', sep = '/'), envir = attach(NULL))
#library(coxme)

library(meta)
# 4/8/2023 on Discovery now
# installing meta package

# checking all numbers 4/1/2023
# rotterdam bone loss (education adjusted))
beta_rot = -0.008	# matches "Rotterdam education adjusted.xlsx" last modifed 10/21/2022)
se_rot = 0.039 # matches above p=0.837
exp(beta_rot)
# fhs bone loss (education adjusted)
beta_fhs = 0.26127424 # beta = 0.26127424 SE = 0.07678422 
se_fhs = 0.07678422
exp(beta_fhs)
# meta-analysis
TE = c(beta_fhs,beta_rot)
seTE = c(se_fhs,se_rot)
metabl = metagen(TE,seTE)
summary(metabl)
# FHS results
# 0.2613 [ 0.1108; 0.4118] 
exp(0.2613)
exp(0.1108)
exp(0.4118)
# [1] 1.298617
# [1] 1.117171
# [1] 1.509532
# RS results
# -0.0080 [-0.0844; 0.0684] 
exp(-0.0080)
exp(-0.0844)
exp(.0684)
# [1] 0.9920319
# [1] 0.9190636
# [1] 1.070794
# common effect model results p=0.17 (p=0.1745 to more precision)
result = c(0.0472, -0.0209, 0.1154)
exp(result)
# 1.0483317 0.9793169 1.1223223
# random effect model results p=0.38 (p=0.3779 to more precision)
result = c( 0.1185,-0.1449,0.3819)
exp(result)
# 1.1258069 0.8651088 1.4650656
#png(paste(out_dir, paste0('boneloss_educ_corrected', Sys.Date(), '.png'), sep = '/'), width = 10, height = 5,  units = 'in', res = 300)
forest(metabl, digits = 3, colgap.forest = unit(0, 'mm'), fontfamily = "sans", colgap.studlab = unit(-5, "mm"))
#dev.off()

# baseline BMD dementia meta analysis
rott_beta <- -0.446
rott_se <- 0.573 # p=0.437
fram_beta <- -0.75420793 # beta = -0.75420793 SE = 0.65218041
fram_se <- 0.65218041
# 8/4/2023 update MAP with <60 exclusion
# 0.3259 ( 0.0916, 1.1590) -1.1211(0.6473, 0.0833)
# map_beta = -1.1211
# map_se = 0.6473
#map_beta <- -1.1183 # est and SE for "bmd044" in "R44.models.rtf" in first table # p=0.0841
#map_se <- 0.6473
map_beta <- -1.1211
map_se <- 0.6473
out <-metagen(TE=c(fram_beta,rott_beta,map_beta),seTE=c(fram_se,rott_se,map_se))
summary(out)
forest(out)
# common = c(-0.7452,-1.4478, -0.0427) # p=0.0376
# updated 8/4/2023 common effect model  -0.7461 [-1.4486; -0.0436] p = 0.0374
common = c(-0.7461,-1.4486,-0.0436)
exp(common)
# [1] 0.4746394 0.2350869 0.9581988 # this is previous
# [1] 0.4742124 0.2348989 0.9573368 this is updated
# random = c(-0.7452,-1.4478, -0.0427) # p=0.0376
random = c(-0.7461,-1.4486,-0.0436) # updated
exp(random)
# [1] 0.4746394 0.2350869 0.9581988
fhs = c(-0.7542, -2.0325, 0.5240)
exp(fhs)
# [1] 0.4703868 0.1310076 1.6887692
rott = c(-0.4460,-1.5691, 0.6771)
exp(rott)
# [1] 0.6401838 0.2082325 1.9681618
map = c(-1.1183,-2.3870,0.1504)
exp(map)
# 0.32683494 0.09190499 1.16229907

# baseline BMD AD meta analysis
fram_beta <- -0.52440861 # beta = -0.52440861 SE = 0.71074219
fram_se <- 0.71074219
# from Sam: beta: -0.37943, SE: 0.67465, P-value: 0.57284
# from Sam: HR: 0.68 (0.18-2.6, p=0.57)
# checking RS numbers 0.68 (0.18-2.6, p=0.57) # added by Sam to Table 4 in yellow highlight in "Bone Loss Dementia Manuscript 11_3_2022_SG.docx sent via e-mail on 11/9/2022 
# Beta: -0.37943, SE: 0.67465, P-value: 0.57284 from e-mail from Sam dated 11/9/2022
rott_beta <- -0.37943
rott_se <- 0.67465
#map_beta <- -1.0893 # # est and SE for "bmd044" in "R44.models.rtf" in second table # p=0.1008
# update 8/4/2023 with exlcusion of those <60
# 0.3355 ( 0.0913, 1.2323) -1.0921(0.6638, 0.0999)
map_beta <- -1.0921
#map_se <- 0.6638
map_se <- 0.6638
out <-metagen(TE=c(fram_beta,rott_beta,map_beta),seTE=c(fram_se,rott_se,map_se))
summary(out)
# -0.6739 [-1.4458; 0.0981] p=0.0871 # fixed
# update 8/4 with exclusion of those <60
# -0.7318 [-1.4397; -0.0239] -2.03  0.0427
exp(-0.7318)
exp(-1.4397)
exp(-0.0239)
# exp(-0.6739)
# exp(-1.4458)
# exp(0.0981)
# [1] 0.5097168
# [1] 0.2355576
# [1] 1.103073
# new below
# > exp(-0.7318)
# [1] 0.4810423
# > exp(-1.4397)
# [1] 0.2369988
# > exp(-0.0239)
# [1] 0.9763833
# -0.6739 [-1.4458; 0.0981] -1.71  0.0871 # random
# update random 8/4/2023 -0.7318 [-1.4397; -0.0239] -2.03  0.0427
# FHS
fhs = c(-0.5244,-1.9174,0.8686)
rs = c(-0.3794,-1.7017,0.9429)
#map = c(-1.0893,-2.3903,0.2117)
exp(fhs)
exp(rs)
#exp(map)
# [1] 0.5919104 0.1469886 2.3835715
# [1] 0.6842718 0.1823732 2.5674161
# [1] 0.3364519 0.0916022 1.2357771
# new map 0.3355 ( 0.0913, 1.2323 )

forest(out, digits = 4, colgap.forest = unit(0, 'mm'), fontfamily = "sans", colgap.studlab = unit(-5, "mm"))

set.seed(42)
# DATA & PARAMETERS -------------------------------------------------------

data <- readRDS(paste(data_dir, 'BoneLoss.rds', sep = '/'))
kinshipDf <- readRDS(paste(data_dir, 'Kinship_for_BoneLoss.rds', sep = '/'))

data2 <- readRDS(paste(data_dir, 'BoneLoss2.rds', sep = '/'))
kinshipDf2 <- readRDS(paste(data_dir, 'Kinship_for_BoneLoss2.rds', sep = '/'))
# these are primary models for FHS
library(coxme)
mod <- coxme(Surv(SurvTime10, DemStatus10) ~ BoneLoss + AgeBase + EducationGroup +
        APOEStatus + BMIBase + Sex + (1|RANID), data = data, varlist = list(kinshipDf))
summary(mod)
exp(confint(mod))
# reran 4/1/2023 
# coef exp(coef)   se(coef)     z       p
# BoneLoss               0.26127424 1.2985837 0.07678422  3.40 0.00067
# 2.5 %   97.5 %
#   BoneLoss              1.1171490 1.509485

# 7/27/2023 also correct for baseline BMD
mod <- coxme(Surv(SurvTime10, DemStatus10) ~ BMDBase + BoneLoss + AgeBase + EducationGroup +
               APOEStatus + BMIBase + Sex + (1|RANID), data = data, varlist = list(kinshipDf))
summary(mod)
exp(confint(mod))
# 1.3074785 1.1095617 1.540698, p=0.0014

mod <- coxme(Surv(SurvTime10, DemStatus10) ~ BMDBase + AgeBase + EducationGroup +
               APOEStatus + BMIBase + Sex + (1|RANID), data = data, varlist = list(kinshipDf))
summary(mod)
exp(confint(mod))
# coef exp(coef)   se(coef)     z      p
# BMDBase               -0.75420793 0.4703830 0.65218041 -1.16 0.2500
# 2.5 %   97.5 %
#   BMDBase               0.1310131 1.688840

mod <- coxme(Surv(SurvTime10, ADStatus10) ~ BMDBase + AgeBase + EducationGroup +
               APOEStatus + BMIBase + Sex + (1|RANID), data = data, varlist = list(kinshipDf))
summary(mod)
exp(confint(mod))
# coef exp(coef)   se(coef)     z      p
# BMDBase               -0.52440861 0.5919053 0.71074219 -0.74 0.4600
# 2.5 %    97.5 %
#   BMDBase               0.1469831 2.3836204

# 7/27/2023
# power calculation
# BMD
HR = 0.4703830
sd = sd(data$BMDBase) # 0.1550656
# bone loss
HR = 1.30
sd = sd(data$BoneLoss) # 0.9399553
psi = 207/1643 # just FHS
psi = (207+210+189)/(1643+2138+655)
psi2 = (207+210)/(1643+2138)
library(powerSurvEpi)
# 4436
powerEpiCont.default(n=4437,theta = 0.47,sigma2 = 0.155^2, psi=psi,rho2=0) # 82% power
powerEpiCont.default(n=3781,theta = 1.30,sigma2 = 0.940^2, psi=psi2,rho2=0) # 99.99% power

# Table 1 saved in Table1_2022_11_21.docx
# code below not working do to package issue
# manual check

table(data$Sex,exclude=F)
100*table(data$Sex,exclude=F)/nrow(data)
median(data$AgeBase)
quantile(data$AgeBase,p=c(0.25,0.75))
table(data$EstrogenUse,exclude=F)
100*table(data$EstrogenUse,exclude=F)/nrow(data)
median(data$BoneLoss)
quantile(data$BoneLoss,p=c(0.25,0.75))
median(data$BMDBase)
quantile(data$BMDBase,p=c(0.25,0.75))
median(data$BMIBase,na.rm = T)
quantile(data$BMIBase,p=c(0.25,0.75),na.rm=T)
table(data$APOEStatus,exclude=F)
100*table(data$APOEStatus,exclude=F)/(nrow(data)-19)
table(data$EducationGroup,exclude=F)
100*table(data$EducationGroup,exclude=F)/(nrow(data))
table(data$DemStatus10,exclude=F)
100*table(data$DemStatus10,exclude=F)/(nrow(data))
table(data$ADStatus10,exclude=F)
100*table(data$ADStatus10,exclude=F)/(nrow(data))

# Table 1 Rotterdam: see what Sam filled in for BMD Dementia Manuscript 11_3_2022_SG.docx in Table S1
# note bone loss median updated in e-mail from Sam dated 11/9/2022
# Table 1 MAP checked with "R44.descriptives.rtf

# Table 1
tmp <- data %>%
  dplyr::select(
    Sex,
    AgeBase,
    EstrogenUse,
    Cohort,
    BoneLoss,
    BoneLossLowestQuartileLabels,
    BMDBase,
    BMD_tscore,
    CurrentSmoker,
    BMIBase,
    APOEStatus,
    #Mobility,
    EducationGroup,
    matches('(Dem|AD)Status10')
  ) %>%
  gtsummary::tbl_summary(
    by = 'Sex',
    missing_text = 'Missing',
    label = list(
      AgeBase ~ 'Age at Baseline Visit (years)',
      EstrogenUse ~ 'Estrogen',
      BoneLoss ~ 'Bone Loss',
      BoneLossLowestQuartileLabels ~ 'Bone Loss Quartiles',
      BMDBase ~ 'Baseline BMD (g/cm\u00B2)',
      BMD_tscore ~ 'Baseline BMD T-score',
      CurrentSmoker ~ 'Current Smoker',
      BMIBase ~ 'BMI at Baseline Visit',
      APOEStatus ~ '\u22651 APOE \u03B54 allele',
      EducationGroup ~ 'Education',
      # DemStatus3 ~ 'Dementia Status at Year 3',
      # DemStatus5 ~ 'Dementia Status at Year 5',
      DemStatus10 ~ 'Dementia Status at Year 10',
      # ADStatus3 ~ 'AD Status at Year 3',
      # ADStatus5 ~ 'AD Status at Year 5',
      ADStatus10 ~ 'AD Status at Year 10'
    )
  )
tmp

 # tmp %>%
#   flextable::as_flextable() %>%
#   flextable::font(fontname = 'Arial', part = 'all') %>%
#   flextable::save_as_docx(path = paste0(home_dir, '/worddocs/Table1_', gsub('-', '_', Sys.Date()), '.docx'))



tmp %>%
  # rrtable::df2flextable2(
  #   .,
  #   vanilla = T,
  #   font = 'Arial',
  #   fontsize = 11,
  #   colorheader = F,
  #   odd_header = 'transparent',
  #   even_body = 'transparent',
  #   vlines = F
  # ) %>%
  # flextable::font(fontname = 'Arial', part = 'all') %>% 
  flextable::save_as_docx(path = paste0(home_dir, '/worddocs/Table1_test_', gsub('-', '_', Sys.Date()), '.docx'))



# end Table 1

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
  Education = character(),
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
    Education = 'EducationGroup',
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
      Row == 2 ~ 'BMD_tscore',
      Row == 3 ~ 'BoneLoss',
      Row == 4 ~ 'BoneLossLowestQuartileLabels',
      Row == 5 ~ 'BoneLossQ4',
      Row == 6 ~ 'as.numeric(BoneLossLowestQuartileLabels)',
      TRUE ~ PrimaryBMD
    ),
    SecondaryBMD = case_when(Row == 7 ~ 'BoneLoss',
                             TRUE ~ SecondaryBMD)
  ) #%>%
  # dplyr::slice(rep(1:nrow(.), times = 3)) %>%
  # dplyr::mutate(
  #   Row = row_number(),
  #   Model = case_when(
  #     Row %in% 1:7 ~ paste(Model, '- 10-year FU'),
  #     Row %in% 8:14 ~ paste(Model, '- 5-year FU'),
  #     Row %in% 15:nrow(.) ~ paste(Model, '- 3-year FU')
  #   ),
  #   Survival = case_when(Row %in% 8:14 ~ 'SurvTime5',
  #                        Row %in% 15:nrow(.) ~ 'SurvTime5',
  #                        TRUE ~ Survival),
  #   Outcome = case_when(Row %in% 8:14 ~ 'DemStatus3',
  #                       Row %in% 15:nrow(.) ~ 'DemStatus3',
  #                       TRUE ~ Outcome)
  # )

tbl2_results <- lapply(1:nrow(tbl2), function(row) {
  row_info <- tbl2[row, ]
  mdl <- fit_model(row_info, kinship_matrix = kinshipDf, df = data)
  tidy_output(mdl, row_info)
})
tbl2_results <- do.call(rbind, tbl2_results)
saveRDS(tbl2_results, paste0(res_dir, '/BoneLoss_Table2.rds'))
write.csv(tbl2_results, paste0(res_dir, '/BoneLoss_Table2.csv'))

tbl2_results %>%
  rrtable::df2flextable2(
    .,
    vanilla = T,
    font = 'Arial',
    fontsize = 11,
    colorheader = F,
    odd_header = 'transparent',
    even_body = 'transparent',
    vlines = F
  ) %>%
  flextable::font(fontname = 'Arial', part = 'all') %>% 
  flextable::save_as_docx(path = paste0(home_dir, '/worddocs/Table2_test_', gsub('-', '_', Sys.Date()), '.docx'))




#Exactly the same as above, but saving models without tidying to compute SEs for meta-analysis.
tbl2_results <- lapply(1:nrow(tbl2), function(row) {
  row_info <- tbl2[row, ]
  mdl <- fit_model(row_info, kinship_matrix = kinshipDf, df = data)
})
saveRDS(tbl2_results, paste0(res_dir, '/BoneLoss_Table2_Models.rds'))

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
tbl3_results <- do.call(rbind, tbl3_results)
saveRDS(tbl3_results, paste0(res_dir, '/BoneLoss_Table3.rds'))
write.csv(tbl3_results, paste0(res_dir, '/BoneLoss_Table3.csv'))

tbl3_results %>%
  rrtable::df2flextable2(
    .,
    vanilla = T,
    font = 'Arial',
    fontsize = 11,
    colorheader = F,
    odd_header = 'transparent',
    even_body = 'transparent',
    vlines = F
  ) %>%
  flextable::font(fontname = 'Arial', part = 'all') %>% 
  flextable::save_as_docx(path = paste0(home_dir, '/worddocs/Table3_test_', gsub('-', '_', Sys.Date()), '.docx'))


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
    Education = 'EducationGroup',
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
write.csv(tbl4_results, paste0(res_dir, '/BoneLoss_Table4.csv'))

coxme(Surv(SurvTime10, DemStatus10) ~ BoneLoss + AgeBase + 
        APOEStatus + BMIBase + Sex + (1|RANID), data = data, varlist = list(kinshipDf))

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
write.csv(do.call(rbind,tblS3_results), paste0(res_dir, '/BoneLoss_TableS3.csv'))

coxme(Surv(SurvTime10, ADStatus10) ~ BMDBase + AgeBase + 
        APOEStatus + BMIBase + Sex + (1|RANID), data = data, varlist = list(kinshipDf))

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
write.csv(do.call(rbind,tblS4_results), paste0(res_dir, '/BoneLoss_TableS4.csv'))

# Table S5: Table 2 with smoking and estrogen.
tblS5 <- tbl2 %>%
  mutate(Sex = NA,
         Estrogen = 'EstrogenUse',
         Smoke = 'CurrentSmoker') %>%
  dplyr::filter(grepl('(1.1|2.1|2.2|2.3|3).*10-year', Model))
tblS5_results <- lapply(1:nrow(tblS5), function(row) {
  row_info <- tblS5[row, ]
  mdl <- fit_model(row_info, kinship_matrix = kinshipDf, df = data)
  tidy_output(mdl, row_info)
})
saveRDS(do.call(rbind, tblS5_results), paste0(res_dir, '/BoneLoss_TableS5.rds'))
write.csv(do.call(rbind,tblS5_results), paste0(res_dir, '/BoneLoss_TableS5.csv'))

# Table S6: Table 2 with smaller bone loss cohort
tblS6 <- tbl2[1,] %>%
  slice(rep(1:nrow(.), each = 2)) %>%
  dplyr::mutate(Row=row_number(),PrimaryBMD = case_when(Row == 1 ~ 'BoneLoss',Row == 2 ~ 'BoneLoss2', TRUE ~ PrimaryBMD), Table = 'S6')

tblS6_results <- lapply(1:nrow(tblS6), function(row) {
  row_info <- tblS6[row, ]
  mdl <- fit_model(row_info, kinship_matrix = kinshipDf2, df = data2)
  tidy_output(mdl, row_info)
})
tblS6_results <- do.call(rbind, tblS6_results)
saveRDS(tblS6_results, paste0(res_dir, '/BoneLoss_TableS6.rds'))
write.csv(tblS6_results, paste0(res_dir, '/BoneLoss_TableS6.csv'))

# Exactly the same as above, but saving models without tidying to compute SEs for meta-analysis.
tblS6_results2 <- lapply(1:nrow(tblS6), function(row) {
  row_info <- tblS6[row, ]
  mdl <- fit_model(row_info, kinship_matrix = kinshipDf2, df = data2)
})
saveRDS(tblS6_results2 , paste0(res_dir, '/BoneLoss_TableS6_Models.rds'))

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

