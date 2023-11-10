library(rms)
library(survminer)

data <- readRDS(paste(data_dir, 'BoneLoss.rds', sep = '/'))
kinshipDf <- readRDS(paste(data_dir, 'Kinship_for_BoneLoss.rds', sep = '/'))

# Transform into longitudinal data for plotting.

long_data <- data %>%
  dplyr::select(RANID, Cohort, BMDBase, BMDInit, BoneLossLowestQuartileLabels) %>%
  pivot_longer(
    cols = c("BMDBase", "BMDInit"),
    names_to = "Timepoint",
    values_to = "BMD"
  ) %>%
  dplyr::rename(`Bone Loss Quartiles` = BoneLossLowestQuartileLabels) %>%
  dplyr::mutate(
    Timepoint = case_when(Timepoint == "BMDInit" ~ "t1",
                          Timepoint == "BMDBase" ~ "t2"),
    `Bone Loss Quartiles` = as.character(`Bone Loss Quartiles`),
    `Bone Loss Quartiles` = case_when(
      grepl('Q1', `Bone Loss Quartiles`) ~ 'Q1',
      grepl('Q2', `Bone Loss Quartiles`) ~ 'Q2',
      grepl('Q3', `Bone Loss Quartiles`) ~ 'Q3',
      grepl('Q4', `Bone Loss Quartiles`) ~ 'Q4'
    )
  )
View(long_data)

png(paste(out_dir, "BMDovertime.png", sep = "/"), height = 600, width = 600, units = "px")
ggplot(long_data, aes(x = Timepoint, y = BMD)) + 
  geom_point(aes(color = `Bone Loss Quartiles`), show.legend = F) +
  geom_line(aes(color = `Bone Loss Quartiles`, group = RANID), show.legend = F) +
  facet_wrap(~ `Bone Loss Quartiles`) + 
  labs(title = "BMD Over Time", y = expression("BMD (g/cm\u00B2)")) + 
  theme_bw() + 
  theme(text = element_text(size = 25),
        strip.text.x = element_text(size = 25),
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 25)) + 
  scale_color_manual(values = c("black", "gray47", "gray80", "brown1"))
dev.off()

# Survival plots
subset <- data %>%
  dplyr::select(SurvTime10,
                DemStatus10,
                BoneLossQ4,
                AgeBase,
                APOEStatus,
                BMIBase,
                Sex) %>%
  dplyr::mutate(BoneLossQ4 = ifelse(BoneLossQ4 == 1, "Bone Loss Q4", "Bone Loss Q1-3"),
                SurvTime10 = SurvTime10/365.25)
dd <- rms::datadist(subset)
options(datadist = 'dd')
mdl_formula <- "Surv(SurvTime10, DemStatus10) ~ BoneLossQ4 + AgeBase + APOEStatus + BMIBase + Sex"

get_surv_curv <- function(formula_to_fit, df = subset){

  formula_to_fit <- as.formula(formula_to_fit)

  mdl <- rms::cph(formula_to_fit, data = df,
           x = TRUE,
           y = TRUE,
           surv = TRUE,
           model = TRUE)
  # mdl <- coxph(formula = formula_to_fit, data = df)
  # mdl <- survminer::surv_fit(formula_to_fit, data = df)
  return(mdl)
}

mdl_object <- get_surv_curv(mdl_formula, df = subset)
 
png(paste(out_dir, 'SurvivalCurves.png', sep = "/"), width = 720, height = 720, units = 'px')
par(mar=c(5,6,1,1), cex.axis = 2)
rms::survplot(
  mdl_object,
  BoneLossQ4,
  label.curves = list(keys = 'lines',
                      keyloc = c(4, .775),
                      cex = 3.5, 
                      transparent = TRUE),
  col = c('gray47', 'brown1'),
  lty = 1,
  lwd = 4,
  col.fill = c('gray70', 'gray87'),
  conf.int = .95,
  conf = 'bands',
  ylim = c(.75, 1),
  adj.subtitle = FALSE,
  xlab = "Years",
  cex.xlab = 3,
  cex.ylab = 3
)
dev.off()
