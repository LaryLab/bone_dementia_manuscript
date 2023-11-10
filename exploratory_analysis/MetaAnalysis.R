setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
librarian::shelf(coxme,
                 dplyr,
                 grid,
                 meta)

# DATA --------------------------------------------------------------------
#
# 6/21/22: Latest Rotterdam results from Beta_SE.docx.
#
# *************************************************************************

# se_calc <- function(lower_ci, upper_ci) 
#   (upper_ci - lower_ci) / (2*1.96)

# Want only Models 1.1 and 2.1
# framingham <- readRDS(paste(res_dir, 'BoneLoss_Table2_Models.rds', sep = '/'))[c(1,4)]
# rotterdam <- readRDS(paste(res_dir, 'Rotterdam_BoneLoss_Table2.rds', sep = '/')) %>%
#   dplyr::filter(grepl('(1.1|2.1)', Model))

# Baseline BMD education-adjusted results
#framingham <- readRDS(paste0(res_dir, "BoneLoss_Table2_Education_Models.rds"))
#framingham <- readRDS(paste0(res_dir, "BoneLoss_Table2.rds"))
framingham2 <- readRDS(paste0(res_dir, "BoneLoss_TableS6_Models.rds"))
rotterdam <- readRDS(paste0(res_dir, 'Rotterdam_BaselineBMD_EducationAdjusted.rds'))
#map_results <- read_rtf(paste0(res_dir, 'MAP_EducationAdjusted.rtf')) # redone

# Currently only need the first model instead of all Table 2 models.
meta_analysis_df <- lapply(1:1, function(x) {
  # Index for BMD/Bone Loss variables.
  idx <- grep('(BMDBase|BoneLoss)$', names(framingham[[x]]$coef), value = T)
  # Beta
  fram_beta <- framingham[[x]]$coef[idx]
  rott_beta <- rotterdam[x, 'Beta']
  map_beta <- -0.8973
  # SE
  fram_se <- sqrt(diag(vcov(framingham[[x]])))[idx]
  rott_se <- rotterdam[x, 'SE']
  map_se <- 0.5947
  output <- data.frame('beta_F' = fram_beta,
                       'se_F' = fram_se,
                       'beta_R' = rott_beta,
                       'se_R' = rott_se,
                       'beta_M' = map_beta,
                       'se_M' = map_se)
  print(output)
})

meta_analysis_df <- lapply(1:length(meta_analysis_df), function(i) 
  meta_analysis_df[[i]] <- meta_analysis_df[[i]] %>%
    dplyr::mutate(across(everything(), as.numeric)))

# Model 1.1
mdl_1p1 <- metagen(
  TE = c(meta_analysis_df[[1]]$beta_F, meta_analysis_df[[1]]$beta_R, meta_analysis_df[[1]]$beta_M),
  seTE = c(meta_analysis_df[[1]]$se_F,  meta_analysis_df[[1]]$se_R, meta_analysis_df[[1]]$se_M)
)
summary(mdl_1p1)
png(paste(out_dir, paste0('Model1point1-', Sys.Date(), '.png'), sep = '/'), width = 10, height = 5,  units = 'in', res = 300)
# forest(mdl_1p1, digits = 4, colgap.forest = unit(0, 'mm'), fontfamily = "sans", colgap.studlab = unit(-5, "mm"))
forest(mdl_1p1, digits = 4)
dev.off()

# Model 2.1
# mdl_2p1 <- metagen(
#   TE = c(meta_analysis_df[[2]]$beta_F, meta_analysis_df[[2]]$beta_R, meta_analysis_df[[2]]$beta_M),
#   seTE = c(meta_analysis_df[[2]]$se_F,  meta_analysis_df[[2]]$se_R, meta_analysis_df[[2]]$se_M)
# )
# summary(mdl_2p1)
# png(paste(out_dir, paste0('Model2point1-', Sys.Date(), '.png'), sep = '/'), width = 10, height = 5,  units = 'in', res = 300)
# forest(mdl_2p1, digits = 4, colgap.forest = unit(0, 'mm'), fontfamily = "sans", colgap.studlab = unit(-5, "mm"))
# dev.off()


