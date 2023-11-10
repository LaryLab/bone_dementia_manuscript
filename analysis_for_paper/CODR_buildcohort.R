# LIBRARIES ---------------------------------------------------------------
#rm(list=ls())
# CWL checked 4_1_2023 and regenerated files BoneLoss.rds and Kinship_for_BoneLoss.rds in data directory

# 4/8/2023 now on discovery
# including Carmen's .rprofile items
##### from .rprofile
home_dir <- getwd()
data_dir <- paste0(home_dir, '/data')
fig_dir <- paste0(home_dir, '/figures')
out_dir <- paste0(home_dir, '/output')
script_dir <- paste0(home_dir, '/scripts')
res_dir <- paste0(home_dir, '/results')

librarian::shelf(coxme,
                 dplyr,
                 finalfit,
                 flextable,
                 ggplot2,
                 haven,
                 survival,
                 survminer,
                 tidyr,
                 tidyverse,
                 quiet = T,
                 cran_repo = 'https://cran.r-project.org')
######################

if(!require('librarian')) install.packages('librarian')
library(librarian)
librarian::shelf(dplyr,
                 haven,
                 kinship2,
                 plyr,
                 tidyverse)
sys.source(paste(script_dir, 'utils.R', sep = '/'), envir = attach(NULL))

# FILES -------------------------------------------------------------------
#
# Recurring modifications to data files:
# - Convert to data frame class only, instead of leaving it as tibble.
#     This prevents the two-level column names that are common to SAS files;
#     in SAS, you can specify a variable name and a label for that column.
# - Change all column names to uppercase.
#
# *************************************************************************

# Columns that are used to combine datasets.
common_cols <- c('RANID', 'IDTYPE')

# Files that contain information for both cohorts.
common_files <- list(
  'fhs_master_apoe_19.sas7bdat',
  'vr_dates_2019_a_1175s_19.sas7bdat',
  'vr_demsurv_2018_a_1281s_19.sas7bdat'
)
common_data <- lapply(setNames(common_files, common_files), function(f) {
  output <- read_file_by_ext(f = f, path = data_dir) %>%
    as.matrix() %>%
    as.data.frame(stringsAsFactors = F) %>%
    rename_with(., toupper) %>%
    mutate(RANID = as.numeric(RANID),
           IDTYPE = as.numeric(IDTYPE))
  output
})

# Filenames for original cohort.
orig_files <- list('bmd0_2009s_19.sas7bdat',
                   'ex0_24s_v2_19.sas7bdat',
                   'vr_wkthru_ex32_0_0997s_19.sas7bdat')
orig_data <- lapply(setNames(orig_files, orig_files), read_file_by_ext, path = data_dir)

off_files <- list(
  't_bmd_ex07_1_0104s_v1_19.sas7bdat',
  't_bmdhs_2008_1_0748s_19.sas7bdat',
  'ex1_8s_19.sas7bdat',
  'meno1_8s_19.sas7bdat',
  'vr_wkthru_ex09_1_1001s_19.sas7bdat'
)
off_data <- lapply(setNames(off_files, off_files), read_file_by_ext, path = data_dir)

# APOE --------------------------------------------------------------------
#
# We are primarily interested in whether a subject has at least one copy of
#   the e4 allele (risk factor for dementia, especially AD).
#
# *************************************************************************

common_data[['APOE']] <-
  common_data[['fhs_master_apoe_19.sas7bdat']] %>%
  # Subset only Original (IDTYPE = 0) and Offspring (IDTYPE = 1)
  dplyr::filter(IDTYPE %in% 0:1) %>%
  # Status = 1 if at least one e4 allele
  mutate(APOEStatus = case_when(grepl('4', APOE) ~ 1,
                                TRUE ~ 0))

# DATES -------------------------------------------------------------------

common_data[['Dates']] <-
  common_data[['vr_dates_2019_a_1175s_19.sas7bdat']] %>%
  dplyr::filter(IDTYPE %in% 0:1) %>%
  # Create age columns to denote age at initial visit for both cohorts.
  mutate(
    AgeInit = case_when(
      IDTYPE == 0 ~ AGE20,
      IDTYPE == 1 ~ AGE6,
      TRUE ~ NA_real_
    ),
    AgeBase = case_when(
      IDTYPE == 0 ~ AGE24,
      IDTYPE == 1 ~ AGE8,
      TRUE ~ NA_real_
    ),
    AgeMid = case_when(
      IDTYPE == 0 ~ AGE22,
      IDTYPE == 1 ~ AGE7,
      TRUE ~ NA_real_
    ),
    InitExamDate = case_when(
      IDTYPE == 0 ~ DATE20,
      # IDTYPE == 1 ~ DATE6,
      TRUE ~ NA_real_
    ),
    BaseExamDate = case_when(
      IDTYPE == 0 ~ DATE24,
      # IDTYPE == 1 ~ DATE8,
      TRUE ~ NA_real_
    ),
    MidExamDate = case_when(
      IDTYPE == 0 ~ DATE22,
      # IDTYPE == 1 ~ DATE8,
      TRUE ~ NA_real_
    )
  ) %>%
  select(RANID, IDTYPE, SEX, matches('(Init|Base|Mid)'), matches('[A-Z](6|7|8|20|22|24)$'))

# DEMENTIA ----------------------------------------------------------------

common_data[['Dementia']] <-
  common_data[['vr_demsurv_2018_a_1281s_19.sas7bdat']] %>%
  dplyr::filter(IDTYPE %in% 0:1)

# BMD files
orig_data[['BMD']] <- orig_data[['bmd0_2009s_19.sas7bdat']] %>%
  as.matrix() %>%
  as.data.frame(stringsAsFactors = F) %>%
  rename_with(., toupper) %>%
  dplyr::filter(IDTYPE == 0) %>%
  dplyr::rename(BMDInit = F20NBMD,
                BMDBase = F24NBMD,
                BMDMid = F22NBMD) %>%
  dplyr::select(RANID, IDTYPE, BMDInit, BMDBase, BMDMid) %>%
  mutate(RANID = as.numeric(RANID),
         IDTYPE = as.numeric(IDTYPE)) %>%
  # Remove is missing BMD measurements at one or both exams.
  dplyr::filter(!is.na(BMDInit) & !is.na(BMDBase))
orig_data[['BMD']] <- orig_data[['BMD']] %>%
  left_join(., common_data[["Dates"]] %>% 
              dplyr::filter(IDTYPE == 0) %>%
              dplyr::select(RANID, IDTYPE, InitExamDate, MidExamDate, BaseExamDate), 
            by = common_cols)

off_data[['BMD6']] <- off_data[['t_bmd_ex07_1_0104s_v1_19.sas7bdat']] %>%
  as.matrix() %>%
  as.data.frame(stringsAsFactors = F) %>%
  rename_with(., toupper) %>%
  dplyr::filter(IDTYPE == 1) %>%
  dplyr::rename(BMDInit = F6_7NBMD,
                InitExamDate = F6_7SCDT) %>%
  mutate(RANID = as.numeric(RANID),
         IDTYPE = as.numeric(IDTYPE)) %>%
  dplyr::select(RANID, IDTYPE, InitExamDate, BMDInit)

off_data[['BMD8']] <- off_data[['t_bmdhs_2008_1_0748s_19.sas7bdat']] %>%
  as.matrix() %>%
  as.data.frame(stringsAsFactors = F) %>%
  rename_with(., toupper) %>%
  dplyr::filter(IDTYPE == 1) %>%
  dplyr::rename(BMDBase = F8CBNBMD,
                BaseExamDate = F8CBSCDT) %>%
  dplyr::rename(BMDMid = F7_5NBMD,
                MidExamDate = F7_5SCDT) %>%
  mutate(RANID = as.numeric(RANID),
         IDTYPE = as.numeric(IDTYPE)) %>%
  dplyr::select(RANID, IDTYPE, BaseExamDate, BMDBase, MidExamDate, BMDMid)

off_data[['BMD']] <- off_data[['BMD6']] %>%
  left_join(off_data[['BMD8']], by = common_cols) %>%
  # dplyr::filter(!is.na(BMDInit) & !is.na(BMDBase)) %>%
  dplyr::filter(!is.na(BMDBase)) %>%
  mutate(InitExamDate = as.numeric(InitExamDate),
         BaseExamDate = as.numeric(BaseExamDate),
         MidExamDate = as.numeric(MidExamDate))

common_data[['BMD']] <- bind_rows(orig_data[['BMD']], off_data[['BMD']]) %>%
  mutate(RANID = as.numeric(RANID),
         IDTYPE = as.numeric(IDTYPE),
         BMDInit = as.numeric(BMDInit),
         BMDMid = as.numeric(BMDMid),
         BMDBase = as.numeric(BMDBase),
         InitExamDate = as.numeric(InitExamDate),
         MidExamDate = as.numeric(MidExamDate),
         BaseExamDate = as.numeric(BaseExamDate))

# Clinical
orig_data[['Clinical']] <- orig_data[['ex0_24s_v2_19.sas7bdat']] %>%
  as.matrix() %>%
  as.data.frame(stringsAsFactors = F) %>%
  rename_with(., toupper) %>%
  dplyr::mutate(
    Estrogen = FQ186,
    HelpTransferring = FQ057) %>%
    # FQ057 = as.numeric(trimws(FQ057)),
    # HelpTransferring = case_when(
    #   FQ057 %in% 1:4 ~ 'Yes Help',
    #   FQ057 %in% 0 ~ 'No Help',
    #   TRUE ~ NA_character_)) %>%
  dplyr::select(RANID, IDTYPE, Estrogen, HelpTransferring) %>%
  dplyr::mutate(RANID = as.numeric(RANID),
                IDTYPE = as.numeric(IDTYPE))

off_data[['Clinical1']] <- off_data[['ex1_8s_19.sas7bdat']] %>%
  as.matrix() %>%
  as.data.frame(stringsAsFactors = F) %>%
  rename_with(., toupper) %>%
  mutate(HelpTransferring = H477) %>%
  # mutate(HelpTransferring = case_when(
  #   H477 %in% 1:4 ~ 'Yes Help',
  #   H477 %in% 0 ~ 'No Help',
  #   TRUE ~ NA_character_)) %>%
  dplyr::select(RANID, IDTYPE, HelpTransferring) %>%
  dplyr::mutate(RANID = as.numeric(RANID),
                IDTYPE = as.numeric(IDTYPE))

off_data[['Clinical2']] <- off_data[['meno1_8s_19.sas7bdat']] %>%
  as.matrix() %>%
  as.data.frame(stringsAsFactors = F) %>%
  rename_with(., toupper) %>%
  mutate(Estrogen = EST8) %>%
  dplyr::select(RANID, IDTYPE, Estrogen) %>%
  dplyr::mutate(RANID = as.numeric(RANID),
                IDTYPE = as.numeric(IDTYPE),
                Estrogen = as.character(Estrogen))

off_data[['Clinical']] <- off_data[['Clinical1']] %>%
  full_join(off_data[['Clinical2']], by = common_cols)

common_data[['Clinical']] <-
  bind_rows(orig_data[['Clinical']], off_data[['Clinical']]) %>%
  mutate(RANID = as.numeric(RANID),
         IDTYPE = as.numeric(IDTYPE),
         HelpTransferring = as.numeric(trimws(HelpTransferring)),
         Mobility = case_when(
           HelpTransferring == 0 ~ 'No Help',
           HelpTransferring %in% 1:4 ~ 'Yes Help',
           TRUE ~ NA_character_
         ))

# Workthru
orig_data[['Workthru']] <- orig_data[['vr_wkthru_ex32_0_0997s_19.sas7bdat']] %>%
  # as.matrix() %>%
  # as.data.frame(stringsAsFactors = F) %>%
  rename_with(., toupper) %>%dplyr::rename(
    BMIInit = BMI20,
    BMIBase = BMI24,
    CurrentSmoker = CURRSMK24) %>%
  # mutate(
  #   CurrentSmoker = case_when(CURRSMK24 == 0 ~ 'No',
  #                             CURRSMK24 == 1 ~ 'Yes',
  #                             TRUE ~ NA_character_)
  # ) %>%
  dplyr::select(RANID, IDTYPE, matches('(Init|Base)'), CurrentSmoker)

off_data[['Workthru']] <- off_data[['vr_wkthru_ex09_1_1001s_19.sas7bdat']] %>%
  # as.matrix() %>%
  # as.data.frame(stringsAsFactors = F) %>%
  rename_with(., toupper) %>%
  dplyr::rename(
    # AgeInit = AGE6,
    # AgeBase = AGE8,
    BMIInit = BMI6,
    BMIBase = BMI8,
    CurrentSmoker = CURRSMK8
  ) %>%
  # mutate(CurrentSmoker = case_when(CURRSMK8 == 0 ~ 'No',
  #                                  CURRSMK8 == 1 ~ 'Yes',
  #                                  TRUE ~ NA_character_)) %>%
  dplyr::select(RANID, IDTYPE, matches('(Init|Base)'), CurrentSmoker)

common_data[['Workthru']] <- bind_rows(orig_data[['Workthru']], off_data[['Workthru']]) %>%
  mutate(RANID = as.numeric(RANID),
         IDTYPE = as.numeric(IDTYPE))

# MERGE -------------------------------------------------------------------

# t score formula
# pheno$FN_Tscore <- ((-0.023 + 0.939 * pheno$f8cbnbmd  - 0.019)/1.087 - 0.858) / 0.120
# reference: McManus DD, Rong J, Huan T, Lacey S, Tanriverdi K, Munson PJ, Larson MG, Joehanes R, Murthy V, Shah R, Freedman JE, Levy D. Messenger RNA and MicroRNA transcriptomic signatures of cardiometabolic risk factors. BMC Genomics. 2017 Dec;18(1):139
one_year_cutoff <- 365.25
three_year_cutoff <- one_year_cutoff*3
five_year_cutoff <- one_year_cutoff*5
ten_year_cutoff <- one_year_cutoff*10
bmd_cutoffs <- c(0,.25, .5, .75, 1)

df <- common_data[['BMD']] %>%
  inner_join(common_data[['Dates']] %>% 
               dplyr::select(-InitExamDate, -BaseExamDate, -MidExamDate), 
             by = common_cols) %>%
  # mutate(InitExamDate = gsub(' NA', '', gsub('NA ', '', paste(InitExamDate.x, InitExamDate.y))),
  #        BaseExamDate = gsub(' NA', '', gsub('NA ', '', paste(BaseExamDate.x, BaseExamDate.y)))) %>%
  # mutate(InitExamDate = gsub('NA', NA_character_, InitExamDate),
  #        BaseExamDate = gsub('NA', NA_character_, BaseExamDate)) %>%
  # mutate(InitExamDate = as.numeric(InitExamDate),
  #        BaseExamDate = as.numeric(BaseExamDate)) %>%
  # dplyr::select(-InitExamDate.x, -InitExamDate.y, -BaseExamDate.x, -BaseExamDate.y) %>%
  left_join(common_data[['Dementia']], by = common_cols) %>%
  left_join(common_data[['APOE']], by = common_cols) %>%
  left_join(common_data[['Workthru']], by = common_cols) %>%
  left_join(common_data[['Clinical']], by = common_cols) %>%
  # filter((IDTYPE == 0 &
  #           ATT24 == 1 &
  #           BaseExamDate < DEM_SURVDATE) |
  #          (IDTYPE == 1 & ATT8 == 1 & BaseExamDate < DEM_SURVDATE)
  # ) %>%
  dplyr::filter((IDTYPE == 0 &
            ATT24 == 1 & 
            BaseExamDate < DEM_SURVDATE) |
           (IDTYPE == 1 & ATT8 == 1 & BaseExamDate < DEM_SURVDATE)
  ) %>%
  dplyr::filter(!is.na(BMDInit) & !is.na(BMDBase)) %>%
  dplyr::filter(AgeBase >= 60) %>% # Cutoff based on Alexa's recommendation.
  dplyr::filter(!is.na(InitExamDate)) %>%
  mutate(
    Cohort = case_when(IDTYPE == 0 ~ 'Original',
                       IDTYPE == 1 ~ 'Offspring',
                       
                       TRUE ~ NA_character_),
    # AgeInitQuartiles = cut(AgeInit, include.lowest = T,
    #                        breaks = quantile(AgeInit, probs = bmd_cutoffs)),
    # AgeBaseQuartiles = cut(AgeBase, include.lowest = T,
    #                        breaks = quantile(AgeBase, probs = bmd_cutoffs)),
    AgeBaseCats = case_when(
      AgeBase < 65 ~ '[60,65)',
      AgeBase < 70 ~ '[65,70)',
      AgeBase < 75 ~ '[70,75)',
      AgeBase < 80 ~ '[75,80)',
      AgeBase < 85 ~ '[80,85)',
      AgeBase <= 95 ~ '[85,95)',
    ),
    # Set youngest as reference category.
    # AgeBaseTertiles = relevel(factor(AgeBaseTertiles), ref = '[60,74]'),
    # 79 y/o is the third quartile.
    AgeCat = ifelse(AgeBase < 79, 'Under79', 'AtLeast79'),
    Sex = case_when(SEX == 1 ~ 'Male',
                    SEX == 2 ~ 'Female',
                    TRUE ~ NA_character_),
    AgeSexCats = paste(Sex,AgeBaseCats,sep="_"),
    # Estrogen = as.factor(Estrogen),
    EstrogenUse = ifelse(Sex == 'Male', 'Male',
                         ifelse(Sex == 'Female' & Estrogen == '1','Female Estrogen Use', 
                                ifelse(Sex == 'Female' & Estrogen %in% c('0','2'), 'Female No Estrogen Use', NA))),
    # Find survival time.
    SurvTime = DEM_SURVDATE - BaseExamDate,
    SurvTime3 = ifelse(SurvTime > three_year_cutoff, three_year_cutoff, SurvTime),
    SurvTime5 = ifelse(SurvTime > five_year_cutoff, five_year_cutoff, SurvTime),
    SurvTime10 = ifelse(SurvTime > ten_year_cutoff, ten_year_cutoff, SurvTime),
    # Find dementia/AD status at years 3, 5, and 10.
    DemStatus3 = ifelse(SurvTime > three_year_cutoff, 0, DEM_STATUS),
    DemStatus5 = ifelse(SurvTime > five_year_cutoff, 0, DEM_STATUS),
    DemStatus10 = ifelse(SurvTime > ten_year_cutoff, 0, DEM_STATUS),
    ADStatus3 = ifelse(SurvTime > three_year_cutoff, 0, AD_STATUS),
    ADStatus5 = ifelse(SurvTime > five_year_cutoff, 0, AD_STATUS),
    ADStatus10 = ifelse(SurvTime > ten_year_cutoff, 0, AD_STATUS),
    # Compute bone loss (see Documentation.Rmd for equation).
    BoneLoss = (100 * ((BMDInit - BMDBase) / BMDInit)) / ((BaseExamDate - InitExamDate) /
                                                            one_year_cutoff),
    BoneLoss2 = (100 * ((BMDMid - BMDBase) / BMDMid)) / ((BaseExamDate - MidExamDate) /
                                                            one_year_cutoff),
    # Create quartiles: Q4 = Highest bone loss quartile OR Q4 = Lowest baseline BMD quartile.
    BMD_tscore = ((-0.023 + 0.939 * BMDBase - 0.019)/1.087 - 0.858) / 0.120,
    BoneLossLowestQuartileLabels = cut(
      BoneLoss,
      include.lowest = T,
      breaks = quantile(BoneLoss, probs = bmd_cutoffs)
    ),
    BMDBaseHighestQuartileLabels = cut(
      BMDBase,
      include.lowest = T,
      breaks = quantile(BMDBase, probs = bmd_cutoffs)
    ),
    # BMDBaseHighestQuartile = ifelse((
    #   5 - as.numeric(BMDBaseHighestQuartileLabels)
    # ) == 1, 1, 0),
    # Assign reference quartiles. NOTE: Does not seem to work. See below.
    # BMDBaseHighestQuartileLabels = relevel(BMDBaseHighestQuartileLabels, ref = '(0.959,1.6]'),
    # BoneLossLowestQuartileLabels = relevel(BoneLossLowestQuartileLabels, ref = '[-5.7,-0.259]'),
    DeltaBMD = BMDInit - BMDBase,
    DeltaBMDStd = (DeltaBMD - mean(DeltaBMD)) / sd(DeltaBMD)
  ) %>%
  dplyr::select(-IDTYPE, -SEX)

# Assign Q1 to highest BMD quartile.
levels(df$BMDBaseHighestQuartileLabels) <- paste(paste0('Q', 4:1), levels(df$BMDBaseHighestQuartileLabels))
# Assign Q1 to lowest bone loss quartile.
levels(df$BoneLossLowestQuartileLabels) <- paste(paste0('Q', 1:4), levels(df$BoneLossLowestQuartileLabels))
# Assign Q1 as reference quartile.
df$BMDBaseHighestQuartileLabels <- relevel(df$BMDBaseHighestQuartileLabels, ref = 'Q1 (0.956,1.6]')
df$BoneLossLowestQuartileLabels <- relevel(df$BoneLossLowestQuartileLabels, ref = 'Q1 [-5.7,-0.299]')
# Assign Q4 as reference quartile.
df$BMDBaseQ4 <- ifelse(df$BMDBaseHighestQuartileLabels == 'Q4 [0.345,0.753]', 1, 0)
df$BoneLossQ4 <- ifelse(df$BoneLossLowestQuartileLabels == 'Q4 (0.779,4.86]', 1, 0)
# saveRDS(df, paste(data_dir, 'BoneLoss.rds', sep = '/'))

# add education
educ <- read_sas(paste(data_dir, 'vr_np_2018_a_1185s_19.sas7bdat', sep = '/')) %>%
  dplyr::filter(IDTYPE %in% 0:1)
names(educ) <- toupper(names(educ))
educ <- educ %>%
  distinct(RANID, .keep_all = T) %>%
  dplyr::mutate(Education = case_when(
    EDUCG == 0 ~ '0 - Did not graduate HS',
    EDUCG == 1 ~ '1 - High school',
    EDUCG == 2 ~ '2 - Some college',
    EDUCG == 3 ~ '3 - College graduate',
    TRUE ~ NA_character_
  )) %>% dplyr::select(RANID, EDUCG)

df <- df %>%
  dplyr::left_join(., educ, by = 'RANID') %>%
  dplyr::mutate(
    EducationGroup = NA_character_,
    EducationGroup = case_when(EDUCG %in% 0:1 ~ 'HS',
                               EDUCG %in% 2:3 ~ 'College',
                               TRUE ~ 'Missing')
  )

saveRDS(df, paste(data_dir, 'BoneLoss.rds', sep = '/'))
# df2 <- df  %>%
#    filter(!is.na(BoneLoss2))
df2 <- df[!is.na(df$BoneLoss2),]
saveRDS(df2, paste(data_dir, 'BoneLoss2.rds', sep = '/'))
# KINSHIP CODING ----------------------------------------------------------

# Kinship coding
kinship_df <-
  read_file_by_ext(f = list.files(pattern = 'ped', path = data_dir),
                   path = data_dir) %>%
  dplyr::rename_with(toupper) %>%
  mutate(Sex = ifelse(SEX == 1, 'Male', 'Female')) %>%
  dplyr::select(-SEX)
maxped = max(kinship_df$PEDNO) # 1542
missing_id.loss <- as.numeric(setdiff(df$RANID, kinship_df$RANID))
nmiss = length(missing_id.loss)
missing_sex.loss <- df$Sex[df$RANID %in% missing_id.loss]
missing_ped.loss <- data.frame(PEDNO = (maxped+1):(maxped+nmiss),
                               RANID = missing_id.loss,
                               FATHER = NA,
                               MOTHER = NA,
                               Sex = missing_sex.loss,
                               ITWIN = NA)
ped_all.l <- rbind(kinship_df, missing_ped.loss)
ped.l <- with(ped_all.l, kinship2::pedigree(id = RANID, FATHER, MOTHER, sex = Sex, famid = PEDNO))
kmat.l <- kinship2::kinship(ped.l)
kmat1.l <- as.matrix(kmat.l)
ids <- colnames(kmat1.l) %in% df$RANID
kmat.test <- kmat.l[ids, ids]
data = merge(df, ped_all.l, by = 'RANID', suffixes = c('', '.ped')) %>%
  dplyr::select(-Sex.ped)
saveRDS(kmat.test, paste(data_dir, 'Kinship_for_BoneLoss.rds', sep = '/'))
ids2 <- colnames(kmat1.l) %in% df2$RANID
kmat.test2 <- kmat.l[ids2, ids2]
data2 = merge(df2, ped_all.l, by = 'RANID', suffixes = c('', '.ped')) %>%
  dplyr::select(-Sex.ped)
saveRDS(kmat.test2, paste(data_dir, 'Kinship_for_BoneLoss2.rds', sep = '/'))
