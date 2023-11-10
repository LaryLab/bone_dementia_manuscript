# Build the data frame necessary for single time-point BMD modeling.
#   In particular, we are looking at the first BMD measurement in Framingham,
#   which corresponds to Exam 20 for the Original cohort and Exam 6 for Offspring.
#
# DIRECTORIES -------------------------------------------------------------

home_dir <-
  dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
data_dir <- paste(home_dir, 'data', sep = '/')

# LIBRARIES ---------------------------------------------------------------

librarian::shelf(dplyr,
                 haven,
                 kinship2,
                 plyr,
                 tidyverse)
sys.source(paste(home_dir, 'utils.R', sep = '/'), envir = attach(NULL))

# Table to keep track of data frames.
docTbl <- data.frame(
  Table = character(),
  File = character(),
  Criterion = character(),
  InitialN = integer(),
  OutputN = integer(),
  Notes = character()
)

# FILES -------------------------------------------------------------------
#
# Recurring modifications to data files:
# - Convert to data frame class only, instead of leaving it as tibble.
#     This prevents the two-level column names that are common to SAS files;
#     in SAS, you can specify a variable name and a label for that column.
# - Change all column names to uppercase.
#
# *************************************************************************

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

# Files for Original cohort only.
orig_files <- list(
  'bmd0_2009s_19.sas7bdat',
  'ex0_20s_v1_19.sas7bdat',
  'vr_wkthru_ex32_0_0997s_19.sas7bdat'
)
orig_data <- lapply(setNames(orig_files, orig_files), function(f) {
  output <- read_file_by_ext(f = f, path = data_dir) %>%
    as.matrix() %>%
    as.data.frame(stringsAsFactors = F) %>%
    rename_with(., toupper) %>%
    mutate(RANID = as.numeric(RANID),
           IDTYPE = as.numeric(IDTYPE))
  output
})

# Files for Offspring cohort only.
off_files <- list(
  't_bmd_ex07_1_0104s_v1_19.sas7bdat',
  'ex1_6s_v1_19.sas7bdat',
  'meno1_8s_19.sas7bdat',
  'vr_wkthru_ex09_1_1001s_19.sas7bdat'
)
off_data <- lapply(setNames(off_files, off_files), function(f) {
  output <- read_file_by_ext(f = f, path = data_dir) %>%
    as.matrix() %>%
    as.data.frame(stringsAsFactors = F) %>%
    rename_with(., toupper) %>%
    mutate(RANID = as.numeric(RANID),
           IDTYPE = as.numeric(IDTYPE))
  output
})

# Columns that are used to combine datasets.
common_cols <- c('RANID', 'IDTYPE')

# APOE --------------------------------------------------------------------
#
# We are primarily interested in whether a subject has at least one copy of
#   the e4 allele (risk factor for dementia, especially AD).
#
# *************************************************************************

common_data[['APOE']] <-
  common_data[['fhs_master_apoe_19.sas7bdat']] %>%
  # Subset only Original (IDTYPE = 0) and Offspring (IDTYPE = 1)
  filter(IDTYPE %in% 0:1) %>%
  # Status = 1 if at least one e4 allele
  mutate(APOEStatus = case_when(grepl('4', APOE) ~ 1,
                                TRUE ~ 0))

# DATES -------------------------------------------------------------------

common_data[['Dates']] <-
  common_data[['vr_dates_2019_a_1175s_19.sas7bdat']] %>%
  filter(IDTYPE %in% 0:1) %>%
  # Create age columns to denote age at initial visit for both cohorts.
  mutate(
    AgeInit = case_when(IDTYPE == 0 ~ AGE20,
                        IDTYPE == 1 ~ AGE6,
                        TRUE ~ NA_real_),
    InitExamDate = case_when(IDTYPE == 0 ~ DATE20,
                             IDTYPE == 1 ~ DATE6,
                             TRUE ~ NA_real_)
  ) %>%
  select(RANID,
         IDTYPE,
         SEX,
         matches('(Init|Base)'),
         matches('[A-Z](6|20)$'))

# DEMENTIA ----------------------------------------------------------------

common_data[['Dementia']] <-
  common_data[['vr_demsurv_2018_a_1281s_19.sas7bdat']] %>%
  filter(IDTYPE %in% 0:1)

# BMD ---------------------------------------------------------------------
#
# For each cohort, we want to keep only subjects with complete records
#   for BMD measurements.
# Note: A new copy of each data file is created instead of writing over
#   the raw file solely for convenience when mistakes arise and we 
#   need to start over from scratch.
#
# *************************************************************************

# orig_data[['BMD']] is going to be a modified copy of 
#   orig_data[['bmd0_2009s_19.sas7bdat']].
orig_data[['BMD']] <- 
  orig_data[['bmd0_2009s_19.sas7bdat']] %>%
  dplyr::rename(BMDInit = F20NBMD) %>%
  dplyr::select(RANID, IDTYPE, BMDInit)
 
off_data[['BMD']] <- 
  off_data[['t_bmd_ex07_1_0104s_v1_19.sas7bdat']] %>%
  dplyr::rename(BMDInit = F6_7NBMD) %>%
  dplyr::select(RANID, IDTYPE, BMDInit)

# Merge BMD data for both cohorts.
common_data[['BMD']] <- 
  bind_rows(orig_data[['BMD']], off_data[['BMD']])  %>%
  mutate(across(where(is.character), as.numeric)) %>%
  filter(!is.na(BMDInit))

docTbl <- docTbl %>%
  add_row(
    Table = 'Combined BMD',
    File = 'N/A',
    Criterion = 'Exclude if BMD is missing.',
    InitialN = nrow(bind_rows(orig_data[['BMD']], off_data[['BMD']])),
    OutputN = nrow(common_data[['BMD']]),
    Notes = 'Converted RANID, IDTYPE, and BMD columns from character to numeric.'
  )

# CLINICAL ----------------------------------------------------------------

orig_data[['Clinical']] <- orig_data[['ex0_20s_v1_19.sas7bdat']] %>%
  dplyr::mutate(Estrogen = FM202,
                Mobility = FM53) %>%
  dplyr::select(RANID, IDTYPE, Estrogen, Mobility) %>%
  mutate(Estrogen = as.numeric(Estrogen),
         Mobility = as.numeric(Mobility))

off_data[['Clinical1']] <- off_data[['ex1_6s_v1_19.sas7bdat']] %>%
  mutate(Mobility = F094) %>%
  dplyr::select(RANID, IDTYPE, Mobility) %>%
  mutate(Mobility = as.numeric(Mobility))

off_data[['Clinical2']] <- off_data[['meno1_8s_19.sas7bdat']] %>%
  mutate(Estrogen = EST8) %>%
  dplyr::select(RANID, IDTYPE, Estrogen) %>%
  mutate(Estrogen = as.numeric(Estrogen))

off_data[['Clinical']] <- off_data[['Clinical1']] %>%
  inner_join(off_data[['Clinical2']], by = common_cols)

common_data[['Clinical']] <-
  bind_rows(orig_data[['Clinical']], off_data[['Clinical']])

docTbl <- docTbl %>%
  add_row(
    Table = 'Clinical',
    File = 'TBD.',
    Criterion = 'N/A',
    InitialN = NA_real_,
    OutputN = nrow(common_data[['Clinical']]),
    Notes = 'Notes on Estrogen and Mobility TBD.'
  )

# WORKTHRU ----------------------------------------------------------------

orig_data[['Workthru']] <-
  orig_data[['vr_wkthru_ex32_0_0997s_19.sas7bdat']] %>%
  dplyr::rename(BMIInit = BMI20) %>%
  mutate(CurrentSmoker = case_when(CURRSMK20 == 0 ~ 'No',
                                   CURRSMK20 == 1 ~ 'Yes',
                                   TRUE ~ 'Missing')) %>%
  dplyr::select(RANID, IDTYPE, matches('(Init|Base)'), CurrentSmoker)

off_data[['Workthru']] <-
  off_data[['vr_wkthru_ex09_1_1001s_19.sas7bdat']] %>%
  dplyr::rename(BMIInit = BMI6) %>%
  mutate(CurrentSmoker = case_when(CURRSMK6 == 0 ~ 'No',
                                   CURRSMK6 == 1 ~ 'Yes',
                                   TRUE ~ 'Missing')) %>%
  dplyr::select(RANID, IDTYPE, matches('(Init|Base)'), CurrentSmoker)

common_data[['Workthru']] <- bind_rows(orig_data[['Workthru']], off_data[['Workthru']])

docTbl <- docTbl %>%
  add_row(
    Table = 'Workthru',
    File = 'TBD.',
    Criterion = 'N/A',
    InitialN = NA_real_,
    OutputN = nrow(common_data[['Workthru']]),
    Notes = 'TBD.'
  )

# MERGE -------------------------------------------------------------------

one_year_cutoff <- 365.25
three_year_cutoff <- one_year_cutoff*3
five_year_cutoff <- one_year_cutoff*5
ten_year_cutoff <- one_year_cutoff*10
bmd_cutoffs <- c(0,.25, .5, .75, 1)

df <- common_data[['BMD']] %>%
  inner_join(common_data[['Dates']], by = common_cols) %>%
  inner_join(common_data[['Dementia']], by = common_cols) %>%
  filter((IDTYPE == 0 &
            ATT20 == 1 &
            DATE20 < DEM_SURVDATE) |
           (IDTYPE == 1 & ATT6 == 1 & DATE6 < DEM_SURVDATE)
  ) %>%
  left_join(common_data[['APOE']], by = common_cols) %>%
  left_join(common_data[['Clinical']], by = common_cols) %>%
  left_join(common_data[['Workthru']], by = common_cols) %>%
  filter(AgeInit >= 60) %>%
  mutate(
    Cohort = case_when(IDTYPE == 0 ~ 'Original',
                       IDTYPE == 1 ~ 'Offspring',
                       TRUE ~ 'Missing'),
    Sex = case_when(SEX == 1 ~ 'Male',
                    SEX == 2 ~ 'Female',
                    TRUE ~ 'Missing'),
    AgeInitQuartiles = cut(AgeInit, include.lowest = T,
                           breaks = quantile(AgeInit, probs = bmd_cutoffs)),
    # 71 y/o is the median.
    AgeCat = ifelse(AgeBase < 71, 'Under71', 'AtLeast71'),
    EstrogenUse = ifelse(
      Sex == 'Male',
      'Male',
      ifelse(
        Sex == 'Female' & Estrogen == 1,
        'Female Estrogen Use',
        ifelse(Sex == 'Female' &
                 Estrogen %in% c(0,2), 'Female No Estrogen Use', NA)
      )
    ),
    Mobility = case_when(
      Mobility == 0 ~ 'No Help',
      Mobility %in% 1:2 ~ 'Yes Help',
      TRUE ~ NA_character_
    ),
    # Find survival time.
    SurvTime = DEM_SURVDATE - InitExamDate,
    # SurvTime3 = ifelse(SurvTime > three_year_cutoff, three_year_cutoff, SurvTime),
    # SurvTime5 = ifelse(SurvTime > five_year_cutoff, five_year_cutoff, SurvTime),
    SurvTime10 = ifelse(SurvTime > ten_year_cutoff, ten_year_cutoff, SurvTime),
    # Find dementia/AD status at years 3, 5, and 10.
    DemStatus3 = ifelse(SurvTime > three_year_cutoff, 0, DEM_STATUS),
    DemStatus5 = ifelse(SurvTime > five_year_cutoff, 0, DEM_STATUS),
    DemStatus10 = ifelse(SurvTime > ten_year_cutoff, 0, DEM_STATUS),
    ADStatus3 = ifelse(SurvTime > three_year_cutoff, 0, AD_STATUS),
    ADStatus5 = ifelse(SurvTime > five_year_cutoff, 0, AD_STATUS),
    ADStatus10 = ifelse(SurvTime > ten_year_cutoff, 0, AD_STATUS),
    BMDInit = as.numeric(BMDInit),
    BMDInitLowestQuartileLabels = cut(
      BMDInit,
      include.lowest = T,
      breaks = quantile(BMDInit, probs = bmd_cutoffs)
    ),
    BMDInitLowestQuartileLabels = relevel(BMDInitLowestQuartileLabels, ref = 'Q1 (0.925,1.56]'),
    BMDInitLowestQuartile = ifelse((
      5 - as.numeric(BMDInitLowestQuartileLabels)
    ) == 4, 'Yes', 'No')
  ) %>%
  dplyr::select(-IDTYPE,-SEX)

# Relabel such that Q1 is the reference quartile (highest BMD).
levels(df$BMDInitLowestQuartileLabels) <- paste(paste0('Q', 4:1), levels(df$BMDInitLowestQuartileLabels))
df$BMDInitLowestQuartileLabels <- relevel(df$BMDInitLowestQuartileLabels, ref = 'Q1 (0.925,1.56]')

# df %>%
#   dplyr::select(-RANID, -matches('(ATT|AGE|\\_STATUS|^APOE$|Date|DATE|Surv)')) %>%
#   tbl_summary(by = 'Sex', missing_text = 'Missing')
saveRDS(df, paste(data_dir, 'SingleBMD.rds', sep = '/'))
# NOTE 06/03/22: Need to double check patient counts, etc.

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
ped.l <- with(ped_all.l, pedigree(id = RANID, FATHER, MOTHER, sex = Sex, famid = PEDNO))
kmat.l <- kinship(ped.l)
kmat1.l <- as.matrix(kmat.l)
ids <- colnames(kmat1.l) %in% df$RANID
kmat.test <- kmat.l[ids, ids]
data = merge(df, ped_all.l, by = 'RANID', suffixes = c('', '.ped')) %>%
  dplyr::select(-Sex.ped)
saveRDS(kmat.test, paste(data_dir, 'Kinship_for_SingleBMD.rds', sep = '/'))
# rm(list = ls())
