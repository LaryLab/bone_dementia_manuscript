sys.source(paste(script_dir, 'utils.R', sep = '/'), envir = attach(NULL))

data <- readRDS(paste0(data_dir, '/BoneLoss.rds'))
kinshipDf <- readRDS(paste0(data_dir, '/Kinship_for_BoneLoss.rds'))
set.seed(42)
model <- coxme(Surv(SurvTime10,DemStatus10)~AgeBase+BoneLoss+APOEStatus+BMIBase+EducationGroup+Sex+(1|colnames(kinshipDf)),varlist = list(kinshipDf),data=data)
summary(model)
model <- coxme(Surv(SurvTime10,DemStatus10)~AgeBase+BMDBase+APOEStatus+BMIBase+EducationGroup+Sex+(1|colnames(kinshipDf)),varlist = list(kinshipDf),data=data)
# BMDBase               -0.72784681 0.4829477 0.62370322 -1.17 0.2400
model <- coxme(Surv(SurvTime10,DemStatus10)~AgeBase+BMD_tscore+APOEStatus+BMIBase+EducationGroup+Sex+(1|colnames(kinshipDf)),varlist = list(kinshipDf),data=data)
# BMD_tscore            -0.10110792 0.9038355 0.08664095 -1.17 0.2400
model <- coxme(Surv(SurvTime10,DemStatus10)~AgeBase+BMDBaseHighestQuartileLabels+APOEStatus+BMIBase+EducationGroup+Sex+(1|colnames(kinshipDf)),varlist = list(kinshipDf),data=data)
#data$BMDBaseQuartOrdered = factor(data$BMDBaseHighestQuartileLabels,ordered = T)
table(data$BMDBaseHighestQuartileLabel)
table(as.numeric(data$BMDBaseHighestQuartileLabel))
model <- coxme(Surv(SurvTime10,DemStatus10)~AgeBase+as.numeric(BMDBaseHighestQuartileLabels)+APOEStatus+BMIBase+EducationGroup+Sex+(1|colnames(kinshipDf)),varlist = list(kinshipDf),data=data)
summary(model)
# as.numeric(BMDBaseHighestQuartileLabels) -0.08650342 0.9171324 0.07847421 -1.10 0.2700
model <- coxme(Surv(SurvTime10,DemStatus10)~AgeBase+BoneLossLowestQuartileLabels+APOEStatus+BMIBase+EducationGroup+Sex+(1|colnames(kinshipDf)),varlist = list(kinshipDf),data=data)
summary(model)
#                     Chisq    df p    AIC    BIC
# Integrated loglik 274.23 10.00 0 254.23 220.95
# Penalized loglik 428.62 82.61 0 263.41 -11.51
model <- coxme(Surv(SurvTime10,DemStatus10)~AgeBase+as.numeric(BoneLossLowestQuartileLabels)+APOEStatus+BMIBase+EducationGroup+Sex+(1|colnames(kinshipDf)),varlist = list(kinshipDf),data=data)
summary(model)
# as.numeric(BoneLossLowestQuartileLabels)  0.26575050 1.3044096 0.06673066  3.98 6.8e-05
# Chisq    df p    AIC    BIC
# Integrated loglik 272.48  8.00 0 256.48 229.86
# Penalized loglik 427.02 80.89 0 265.24  -3.96