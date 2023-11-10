sys.source(paste(script_dir, 'utils.R', sep = '/'), envir = attach(NULL))

data <- readRDS(paste0(data_dir, 'BoneLoss.rds'))
kinshipDf <- readRDS(paste0(data_dir, 'Kinship_for_BoneLoss.rds'))

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
  )) 
  
data <- data %>%
  dplyr::left_join(., educ, by = 'RANID') %>%
  dplyr::mutate(
    EducationGroup = NA_character_,
    EducationGroup = case_when(EDUCG %in% 0:1 ~ 'HS',
                               EDUCG %in% 2:3 ~ 'College',
                               TRUE ~ 'Missing')
  )

table(data$EducationGroup, useNA = "ifany")

tbl2 <- data.frame(
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
  # Added column 
  Education = character(),
  stringsAsFactors = FALSE
)

# Table 2 with Education
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
    Smoke = NA,
    Education = 'EducationGroup'
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
  )

# tbl2_results <- lapply(1:nrow(tbl2), function(row) {
#   row_info <- tbl2[row, ]
#   mdl <- fit_model(row_info, kinship_matrix = kinshipDf, df = data)
#   tidy_output(mdl, row_info)
# })
# tbl2_results <- do.call(rbind, tbl2_results)
tbl2_results <- lapply(1:nrow(tbl2), function(row) {
  row_info <- tbl2[row, ]
  mdl <- fit_model(row_info, kinship_matrix = kinshipDf, df = data)
})
saveRDS(tbl2_results, paste0(res_dir, "BoneLoss_Table2_Education_Models.rds"))
View(tbl2_results)

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
  flextable::save_as_docx(path = paste0(home_dir, '/worddocs/Table2_w_Education_', gsub('-', '_', Sys.Date()), '.docx'))
