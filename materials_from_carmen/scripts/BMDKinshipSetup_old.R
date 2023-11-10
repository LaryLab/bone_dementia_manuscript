# Slight modification (e.g. variable names) using Meghan's data; see `surv.w.kinship cwl manuscript.R`
invisible(sapply(list.files(
  pattern = '(directories|modules|utils)',
  recursive = TRUE,
  full.names = TRUE),
  function(f)
    sys.source(f, env = attach(NULL))))

# LOAD DATA ---------------------------------------------------------------

# Source: /home/nlittlefie@mmcf.mehealth.org/Meghan Nick/Meghan/Dem & BMD Analysis/Usable Data
original <- readRDS(paste0(data_subdir, 'Original Survival Analysis DF.rds'))
offspring <- readRDS(paste0(data_subdir, 'Offspring Survival Analysis DF.rds'))
# Identify and rename exams to be used for each cohort
original <- original %>%
  # Drop duplicated columns
  dplyr::select(-idtype, -age20, -age24) %>%
  dplyr::rename_with(stringr::str_replace,
                     pattern = '\\.',
                     replacement = '') %>%
  dplyr::rename_with(stringr::str_replace,
                     pattern = '_',
                     replacement = '') %>%
  dplyr::rename(
                # AGECONTS = AGE24,
                BMDINIT = BMD20,
                BMDBASE = BMD24,
                # BMIBASE = BMI24
                ) %>%
  dplyr::rename_with(toupper)
  
offspring <- offspring %>%
# Drop duplicated columns
  dplyr::select(-age6, -age8) %>%
  dplyr::rename_with(stringr::str_replace,
                     pattern = '\\.',
                     replacement = '') %>%
  dplyr::rename_with(stringr::str_replace,
                     pattern = '_',
                     replacement = '') %>%
  dplyr::rename(BMDINIT = BMDbase,
                BMDBASE = BMDafter,
                # BMIBASE = BMI8
                ) %>%
  # dplyr::rename(AGECONTS = AGE8) %>%
  dplyr::rename_with(toupper)
  
# Manipulate additional variables
data <- bind_rows(original, offspring)
data$AGECONTS <- rowSums(data[, c('AGE8', 'AGE24')], na.rm = TRUE)
data$BMIBASE <- rowSums(data[, c('BMI8', 'BMI24')], na.rm = TRUE)

data <- data %>%
  mutate(HELPTRANSFERRING = ifelse(TRANSFERRING %in% c("1", "2"), "Yes",
                                   ifelse(TRANSFERRING == "0", "No", NA))) %>%
  mutate(DEMSTATUS10 = ifelse(SURVTIME > 3652.5, 0, DEMSTATUS)) %>%
  mutate(ADSTATUS10 = ifelse(SURVTIME > 3652.5, 0, ADSTATUS)) %>%
  mutate(SURVTIME10 = ifelse(SURVTIME > 3652.5, 3652.5, SURVTIME)) %>%
  mutate(DEMSTATUS5 = ifelse(SURVTIME > 1826.25, 0, DEMSTATUS)) %>%
  mutate(ADSTATUS5 = ifelse(SURVTIME > 1826.25, 0, ADSTATUS)) %>%
  mutate(SURVTIME5 = ifelse(SURVTIME > 1826.25, 1826.25, SURVTIME)) %>%
  mutate(DEMSTATUS2 = ifelse(SURVTIME > 730.5, 0, DEMSTATUS)) %>%
  mutate(SURVTIME2 = ifelse(SURVTIME > 730.5, 730.5, SURVTIME)) %>%
  mutate(DEMSTATUS1 = ifelse(SURVTIME > 365.25, 0, DEMSTATUS)) %>%
  mutate(SURVTIME1 = ifelse(SURVTIME > 365.25, 365.25, SURVTIME)) %>%
  mutate(SEX = ifelse(is.na(SEX), SEX_EXAMDATES, SEX)) %>%
  mutate(Cohort = ifelse(IDTYPE == 0, "Original", "Offspring")) %>%
  # Combines estrogen usage in women with sex
  mutate(SEXEstrogen = ifelse(SEX == 1, "Male", 
                              ifelse(SEX == 2 & ESTROGEN == 1, "Female Estrogen Use", 
                                     ifelse(SEX == 2 & ESTROGEN == 0, "Female No Estrogen Use", NA)))) %>%
  mutate(AGESQUARED = AGECONTS * AGECONTS) %>%
  mutate(AGEGTE70 = ifelse(AGECONTS >= 70, 1, 0)) %>%
  # Bone loss quartile
  mutate(BMDPCTQUARTILE = as.factor(cut(BMDPERCENTPERYEAR,
                                        include.lowest = T, 
                                        breaks = quantile(data$BMDPERCENTPERYEAR, 
                                                          probs = c(0, .25, .5, .75, 1))))) %>%
  mutate(BMDPCTQUARTILE4 = ifelse(BMDPCTQUARTILE == levels(BMDPCTQUARTILE)[4], "Yes", "No")) %>%
  # Baseline BMD quartile (final BMD before dementia follow-up)
  mutate(BMDBASEQUARTILE = as.factor(cut(BMDBASE, 
                                         include.lowest = T,
                                         breaks = quantile(BMDBASE, probs = c(0, .25, .5, .75, 1))))) %>%
  # Baseline BMD quartile reordered as highest BMD to lowest BMD
  # mutate(quartilefinalltmp = as.factor(cut(
  #   BMDBASE,
  #   include.lowest = T,
  #   breaks =
  #     quantile(BMDBASE, probs = c(0, .25, .5, .75, 1)),
  #   labels = F
  # ))) %>%
  mutate(BMDBASEQUARTILECAT = as.factor(5 - as.numeric(BMDBASEQUARTILE))) %>%
  # mutate(BMDBaseQuartileCattmp = as.factor(5 - as.numeric(quartilefinalltmp))) %>%
  # baseline BMD quartile with "yes" as lowest quartile and "no" as remaining 3 quartiles
  mutate(BMDBASELOWESTQUARTILE = ifelse(BMDBASEQUARTILE == levels(BMDBASEQUARTILE)[1], "Yes", "No")) %>%
  # mutate(BMDBASEQUARTILE4tmp = ifelse(quartilefinalltmp == 4, "Yes", "No")) %>%
  mutate(BONELOSSGROUPS = as.factor(
    case_when(
      BMDPERCENTPERYEAR < 0.5 ~ 0,
      BMDPERCENTPERYEAR >= 1.5 ~ 3,
      BMDPERCENTPERYEAR >= 1 ~ 2,
      BMDPERCENTPERYEAR >= 0.5 ~ 1
    )
  )) %>%
  mutate(ANNUALBONELOSS = as.factor(ifelse(BMDPERCENTPERYEAR >= 1, "\u22651%", "<1%"))) %>%
  mutate(BONELOSSGTE1p5 = as.factor(ifelse(BMDPERCENTPERYEAR >= 1.5, "Yes", "No"))) %>%
  mutate(APOE = factor(APOE, levels = c("24", "22", "23", "33", "34", "44"))) %>%
  mutate(
    AGEGROUPS = case_when(
      AGECONTS <= 74 ~ "<= 74",
      (AGECONTS > 74) & (AGECONTS <= 79) ~ "74 < age <= 79",
      AGECONTS > 79 ~ "> 79"
    )
  ) %>%
  mutate(BMDDELTANorm = scale(BMDDELTA)) %>%
  filter(AGECONTS >= 60)

# Kinship coding
kinship_df <- read_file_by_ext(f = list.files(pattern = 'ped', path = data_subdir),
                               path = data_subdir) %>%
  dplyr::rename_with(toupper)
maxped = max(kinship_df$PEDNO) # 1542
missing_id.loss <- as.numeric(setdiff(data$RANID, kinship_df$RANID))
nmiss = length(missing_id.loss)
missing_sex.loss <- data$SEX[data$RANID %in% missing_id.loss]
missing_ped.loss <- data.frame(PEDNO = (maxped+1):(maxped+nmiss), 
                               RANID = missing_id.loss, 
                               FATHER = NA, 
                               MOTHER = NA, 
                               SEX = missing_sex.loss, 
                               ITWIN = NA)
ped_all.l <- rbind(kinship_df, missing_ped.loss)
ped.l <- with(ped_all.l, pedigree(id = RANID, FATHER, MOTHER, sex = SEX, famid = PEDNO))
kmat.l <- kinship(ped.l)
kmat1.l <-as.matrix(kmat.l)
ids <- colnames(kmat1.l) %in% data$RANID
kmat.test <- kmat.l[ids, ids]
data = merge(data, ped_all.l, by = 'RANID', suffixes = c('', 'ped'))

# Refactor selected variables
data <- data %>%
  dplyr::mutate(BMDBASEQUARTILE = relevel(BMDBASEQUARTILE, ref = '(0.974,1.6]'),
                SEX = as.factor(ifelse(SEX == 1, 'Male', 'Female')),
                APOE4STAT = ifelse(APOE4STAT == 0, 'No E4', 'E4')) 

# saveRDS(data, paste0(data_subdir, 'BMDKinship.rds'))

# PARAMETERS --------------------------------------------------------------
# 
# Given how specific each model is, most of these rows will be manually coded.
#   It is easiest (though not necessarily most expedient) to create separate
#   data frames for each table and then merge them into one data frame.
#   The output is saved to ./Data/ParamGrid.rds
#
# .........................................................................

tbl2 <- tbl3 <- tbl4 <- tblS1 <- tblS2 <- tblS3 <- data.frame(
    Table = character(),
    Model = character(),
    Outcome = character(),
    BaselineBMD = character(),
    Age = character(),
    APOE = character(),
    BMI = character(),
    Estrogen = character(),
    Mobility = character(),
    OtherBMD = character(),
    Sex = character(),
    Smoke = character(),
    stringsAsFactors = FALSE
  )

tbl2 <- tbl2 %>%
  add_row(
    Table = '2',
    Model = NA,
    Outcome = 'DEMSTATUS10',
    BaselineBMD = 'BMDBASE',
    Age = 'AGECONTS',
    APOE = 'APOE4STAT',
    BMI = 'BMIBASE',
    Estrogen = NA,
    Mobility = 'HELPTRANSFERRING',
    OtherBMD = NA,
    Sex = 'SEX',
    Smoke = NA
  ) %>%
  slice(rep(1:n(), each = 7)) %>%
  mutate(Row = row_number(),
         Model = case_when(
           Row == 1 ~ '1.1',
           Row == 2 ~ '1.2',
           Row == 3 ~ '1.3',
           Row == 4 ~ '2.1',
           Row == 5 ~ '2.2',
           Row == 6 ~ '2.3',
           Row == 7 ~ '3'
         ),
         BaselineBMD = case_when(
           Row == 2 ~ 'BMDBASEQUARTILE',
           Row == 3 ~ 'BMDBASELOWESTQUARTILE',
           Row == 4 ~ 'BMDPERCENTPERYEAR',
           Row == 5 ~ 'BMDPCTQUARTILE',
           Row == 6 ~ 'BMDPCTQUARTILE4',
           TRUE ~ BaselineBMD
         ),
         OtherBMD = case_when(
           Row == 7 ~ 'BMDPERCENTPERYEAR',
           TRUE ~ OtherBMD
         ))

tbl3 <- tbl2 %>%
  slice(rep(1:n(), each = 2)) %>%
  mutate(Row = row_number(),
         Table = 3,
         Model = case_when(
           Row %% 2 == 1 ~ paste(Model, 'Female', sep=':'),
           Row %% 2 == 0 ~ paste(Model, 'Male', sep=':')
         ),
         Sex = case_when(
           Row %% 2 == 1 ~ 'SEX:Female',
           Row %% 2 == 0 ~ 'SEX:Male'
         ))

tbl4 <- tbl4 %>%
  add_row(
    Table = '4',
    Model = NA,
    Outcome = 'DEMSTATUS10',
    BaselineBMD = 'BMDPERCENTPERYEAR',
    Age = NA,
    APOE = NA,
    BMI = 'BMIBASE',
    Estrogen = NA,
    Mobility = 'HELPTRANSFERRING',
    OtherBMD = NA,
    Sex = 'SEX',
    Smoke = NA
  ) %>%
  slice(rep(1:n(), each = 7)) %>%
  mutate(Row = row_number(),
         Model = c('4', '4:E4', '4:No E4', '5', '5:age < 79', '5:age >= 79', '6'),
         BaselineBMD = case_when(
           Row == 1 ~ 'BMDPERCENTPERYEAR*APOE4STAT',
           Row == 4 ~ 'BMDPERCENTPERYEAR*AGECONTS',
           Row == 7 ~ 'BMDPERCENTPERYEAR*BMDBASELOWESTQUARTILE',
           TRUE ~ BaselineBMD
         ),
         Age = case_when(
           Row %in% c(1:4, 7) ~ 'AGECONTS',
           Row == 5 ~ 'AGEQuartile:Q1-3',
           Row == 6 ~ 'AGEQuartile:Q4'
         ),
         APOE = case_when(
           Row %in% c(1, 4:7) ~ 'APOE4STAT',
           Row == 2 ~ 'APOE4STAT:E4',
           Row == 3 ~ 'APOE4STAT:No E4'
         ))

tblS1 <- tblS1 %>%
  add_row(
    Table = 'S1',
    Model = NA,
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
  slice(rep(1:n(), each = 7))

tblS2 <- tblS2 %>%
  add_row(
    Table = 'S2',
    Model = NA,
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
  slice(rep(1:n(), each = 2))

tblS3 <- tblS3 %>%
  add_row(
    Table = 'S3',
    Model = NA,
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
  slice(rep(1:n(), each = 2))

# HELPER FUNCTIONS --------------------------------------------------------

summary_output <- function(model_object){
  #' @param model_object Output from coxme models.
  #' @return A data frame with the hazard ratio, 95% CI, and p-value.
  
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
