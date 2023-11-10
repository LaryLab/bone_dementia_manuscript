# Results from Samuel Ghatan (07/11/22);
#   see ~/results/Rotterdam Table 2 and Sex Stratified Table (corrected quantiles).docx

if(!require('textreadr')) install.packages('textreadr')
library(textreadr)

# read_docx converts the document into a character vector, where each element is typically a line from
#   the document or a cell from a table.
results <- read_docx(grep('Rotterdam.*Table 2.*docx', list.files(path = res_dir, full.names = T), value = T))

# This particular document contains two tables.
# The first table runs from 1:94, with the first element being the caption.
tbl2 <- list(
  results[9:15],
  results[16:22],
  results[23:29],
  results[30:36],
  results[37:43],
  results[45:51],
  results[52:58],
  results[59:65],
  results[66:72],
  results[73:79],
  results[81:87],
  results[88:94]
)
tbl2 <- lapply(1:length(tbl2), function(i) matrix(tbl2[[i]], nrow = 1))
tbl2 <- do.call(rbind, tbl2) %>%
  as.data.frame() %>%
  setNames(., c('Parametrization', results[2:7]))
# Format to match Framingham/manuscript.
tbl2 <- tbl2 %>%
  mutate(Exposure = c(rep('Baseline BMD', 5), rep('Bone Loss', 5), rep('Baseline BMD and Bone Loss', 2)),
         Model = c('1.1', rep('1.2', 3), '1.3', '2.1', rep('2.2', 3), '2.3', rep('3', 2)),
         Parametrization = case_when(
           Parametrization == 'Baseline BMD' ~ 'continuous (g/cm\u00B2)',
           Parametrization %in% c('BMD_quart2', 'BMD_change_yr_q2') ~ 'Q2 vs. Q1',
           Parametrization %in% c('BMD_quart3', 'BMD_change_yr_q3') ~ 'Q3 vs. Q1',
           Parametrization %in% c('BMD_quart4', 'BMD_change_yr_q4') ~ 'Q4 vs. Q1',
           Parametrization %in% c('BMD_quantile4th vs 1-3', 'BMD_change_yr quantile4th vs 1-3') ~ 'Q4 vs. Q1-3',
           Parametrization == 'BMD_change_yr' & Model == '2.1' ~ 'continuous (% decline/year)',
           Parametrization == 'BMD_change_yr' & Model == '3' ~ 'Bone Loss',
           Parametrization == 'sec_fnbmd' ~ 'Baseline BMD (g/cm\u00B2)',
           TRUE ~ Parametrization
         )) %>%
  relocate(Exposure, Model, Parametrization)
# Format results to match Framingham: OR (Lower CI, Upper CI, p=pvalue)
format_OR <- formatC(round(as.numeric(tbl2$OR), digits = 2), digits = 2, format = 'f')
format_pval <- formatC(round(as.numeric(tbl2$`p-value`), digits = 3), digits = 3, format = 'f')
format_lowerCI <- formatC(round(as.numeric(tbl2$`Lower 95% CI`), digits = 2), digits = 2, format = 'f')
format_upperCI <- formatC(round(as.numeric(tbl2$`Upper 95% CI`), digits = 2), digits = 2, format = 'f')
tbl2 <- tbl2 %>%
  mutate(`Rotterdam: Hazard Ratio (95% CI, p-value)` = 
           paste0(format_OR, ' (', format_lowerCI, '-', format_upperCI, ', p=', format_pval, ')'))
# View(tbl2)
saveRDS(tbl2, paste(res_dir, 'Rotterdam_BoneLoss_Table2.rds', sep = '/'))

tbl3 <- list(
  results[110:116],
  results[119:125],
  results[128:134],
  results[135:141],
  results[142:148],
  results[151:157],
  results[158:164],
  results[165:171],
  results[174:180],
  results[183:189],
  results[193:199],
  results[202:208],
  results[211:217],
  results[218:224],
  results[225:231],
  results[234:240],
  results[241:247],
  results[248:254],
  results[257:263],
  results[266:272],
  results[276:282],
  results[283:289],
  results[292:298],
  results[299:305]
)
tbl3 <- lapply(1:length(tbl3), function(i) matrix(tbl3[[i]], nrow = 1))
tbl3 <- do.call(rbind, tbl3) %>%
  as.data.frame() %>%
  setNames(., c('Parametrization', results[2:7])) 
tbl3 <- tbl3 %>%
  mutate(Exposure = c(rep('Baseline BMD', 10), rep('Bone Loss', 10), rep('Baseline BMD and Bone Loss', 4)), 
         Model = rep(c('1.1', rep('1.2', 3), '1.3', '2.1', rep('2.2', 3), '2.3', rep('3', 2)), each = 2),
         Sex = c(rep(c('Female', 'Male', rep('Female', 3), rep('Male', 3), 'Female', 'Male'), times = 2),
         rep('Female', 2), rep('Male', 2))) %>%
  relocate(Exposure, Sex, Model, Parametrization)
# Format results to match Framingham: OR (Lower CI, Upper CI, p=pvalue)
format_OR <- formatC(round(as.numeric(tbl3$OR), digits = 2), digits = 2, format = 'f')
format_pval <- formatC(round(as.numeric(tbl3$`p-value`), digits = 3), digits = 3, format = 'f')
format_lowerCI <- formatC(round(as.numeric(tbl3$`Lower 95% CI`), digits = 2), digits = 2, format = 'f')
format_upperCI <- formatC(round(as.numeric(tbl3$`Upper 95% CI`), digits = 2), digits = 2, format = 'f')
tbl3 <- tbl3 %>%
  mutate(`Rotterdam: Hazard Ratio (95% CI, p-value` = 
           paste0(format_OR, ' (', format_lowerCI, '-', format_upperCI, ', p=', format_pval, ')'))
saveRDS(tbl3, paste(res_dir, 'Rotterdam_BoneLoss_Table3.rds', sep = '/'))

# ADJUSTED FOR EDUCATION --------------------------------------------------

rotterdam_edu <- read_docx(paste0(res_dir, 'Rotterdam_EducationAdjusted.docx'))

# Only need values for Model 1.1 (Baseline BMD)
matrix(rotterdam_edu[11:12], nrow = 1) %>%
  as.data.frame() %>%
  dplyr::mutate(across(everything(), as.numeric)) %>%
  dplyr::rename(Beta = V1, 
                SE = V2) %>%
  saveRDS(., paste0(res_dir, 'Rotterdam_BaselineBMD_EducationAdjusted.rds'))
