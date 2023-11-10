librarian::shelf(dplyr,
                 haven,
                 finalfit)

# BASIC EXPRESSIONS -------------------------------------------------------

`%notin%` <- Negate(`%in%`)

# BASIC FUNCTIONS ---------------------------------------------------------

read_file_by_ext <- function(f, path, skip_lines = 10) {
  #' @description Read file by extension type. Currently caters to .gz, .txt, and .xml.
  #' @param f str. Filename.
  #' @param path str. Path to directory.
  #' @param skip_lines int. Number of lines to skip, usually 10 for FHS files.
  
  filename <- paste(path, f, sep = '/')
  
  if (grepl('gz', f)) {
    data <- read.delim(gzfile(filename), skip = skip_lines)
  }
  
  if (grepl('txt', f) && !grepl('gz', f)) {
    data <- read.delim(filename, skip = skip_lines)
  }
  
  if (grepl('xml', f)) {
    # If using .xml files, may want to consider xmlToDataFrame(),
    #   though not all .xml files can be read that way. Using
    #   xmlParse() for simplicity
    data <- xmlParse(filename)
  }
  
  if (grepl('csv', f)) {
    data <- read.csv(filename, stringsAsFactors = F)
  }
  
  if (grepl('sas', f)) {
    data <- haven::read_sas(paste(path, f, sep = '/'))
  }
  
  return(data)
}

sas_from_dir <- function(f, path) {
  #' @description Use haven::read_sas to read sas7bdat file from specified directory.
  #' @param f str.
  #' @param path str.
  
  return(haven::read_sas(paste0(path, f)))
}

load_standard <- function(files, path) {
  #' @description Load and standardize FHS data files by changing variable names to uppercase for standardization purposes.
  #' @param files vector containing full path of files to read.
  
  # Read files into a list (or a data frame if just one)
  data_list <- lapply(files, read_file_by_ext, path = path)
  
  # Convert all column names to uppercase for consistency
  data_list <- lapply(data_list, function(d) {
    names(d) <- toupper(names(d))
    return(d)
  })
  
  # Will leave the list of data frames unnamed for now since
  #   this can be very irregular.
  
  return(data_list)
  
}

# BMD: MORE SPECIFIC FUNCTIONS --------------------------------------------
#
# Functions mainly used for BMD project.
#
# *************************************************************************

summary_output <- function(model_object, ...) {
  #' @description Format coxme models to return hazard ratio, 95% confidence
  #' interval, and p-value using \link[finalfit]{fit2df}. For purposes of
  #' formatting the output for the manuscript, the default number of
  #' significant digits for each of those estimates are c(2,2,3) respectively.
  #' If modifying, specify optional arguments: hazard_digits, ci_digits, pval_digits.
  #' @param model_object Fitted coxme models.
  #' @return A data frame with the hazard ratio, 95% CI, and p-value.
  #' @importFrom finalfit fit2df
  
  # Set default number of significant digits.
  hazard_digits <- 2
  ci_digits <- 2
  pval_digits <- 3
  
  args <- list(...)
  
  if (length(args) != 0) {
    if ('hazard_digits' %in% names(args)) {
      hazard_digits <- args$hazard_digits
    }
    if ('ci_digits' %in% names(args)) {
      ci_digits <- args$ci_digits
    }
    if ('pval_digits' %in% names(args)) {
      pval_digits <- args$pval_digits
    }
  }
  
  summary_object <- finalfit::fit2df(
    model_object,
    digits = c(hazard_digits, ci_digits, pval_digits),
    explanatory_name = 'Variable',
    estimate_name = 'Hazard Ratio',
    p_name = 'p'
  )
  summary_object$Variable <- as.character(summary_object$Variable)
  summary_object$`Hazard Ratio` <-
    as.character(summary_object$`Hazard Ratio`)
  
  return(summary_object)
}

create_param_grid <- function(input_grid,
                              table,
                              model,
                              survival_time,
                              outcome,
                              primary_bmd,
                              age,
                              apoe,
                              bmi,
                              estrogen,
                              mobility,
                              secondary_bmd,
                              sex,
                              smoke) {
  #' @description Add to an existing data frame of parameters that identifies the model,
  #' i.e. in accordance with manuscript, and fits the models according to specified covariates.
  #' @param input_grid A data frame that contains the columns Table, Model, Survival, Outcome,
  #' PrimaryBMD, Age, APOE, BMI, Estrogen, Mobility, SecondaryBMD, Sex, Smoke.
  #' @return A data frame containing all the columns in `input_grid`.
  #' @import dplyr
  
  output_grid <- input_grid %>%
    add_row(
      Table = table,
      Model = model,
      Survival = survival_time,
      Outcome = outcome,
      PrimaryBMD = primary_bmd,
      Age = age,
      APOE = apoe,
      BMI = bmi,
      Estrogen = estrogen,
      Mobility = mobility,
      SecondaryBMD = secondary_bmd,
      Sex = sex,
      Smoke = smoke,
      Row = nrow(input_grid + 1)
    )
  
  return(output_grid)
}

format_model_formula <- function(row_info,
                                 stratified = FALSE,
                                 id_col = 'RANID') {
  #' @description Using information provided in `row_info`,
  #' generate formula to be passed into coxme models.
  #' @param row_info Data frame (see `create_param_grid`).
  #' @param stratified Boolean to indicate whether fitting effect modification models.
  #' @return Object of class formula.
  #' @import dplyr
  
  survival <- row_info %>% dplyr::pull(Survival)
  outcome <- row_info %>% dplyr::select(Outcome)
  
  predictors <- row_info %>%
    dplyr::select(-Table,-Model,-Survival,-Outcome,-Row)
  
  if (!stratified) {
    # If not stratifying, then just concatenate the covariates into a model formula
    predictors <- predictors %>%
      dplyr::select(where( ~ !any(is.na(.)))) %>%
      paste(collapse = ' + ')
  }
  
  if (stratified) {
    # Effect modification variables are denoted [variable]:[group]
    strat_var <- names(predictors)[grep(':', predictors)]
    # print(strat_var)
    # tmp <- strsplit(predictors[, strat_var], split = ':')
    # var_name <- tmp[[1]][1]
    # var_group <- tmp[[1]][2]
    
    predictors <- predictors %>%
      dplyr::select(-!!as.name(strat_var)) %>%
      dplyr::select(where( ~ !any(is.na(.)))) %>%
      paste(collapse = ' + ')
  }
  
  predictors <-
    paste(predictors, paste0('(1|', id_col, ')'), sep = ' + ')
  formula_to_fit <-
    paste0(sprintf('Surv(%s, %s) ~ ', survival, outcome), predictors)
  formula_to_fit <- formula(formula_to_fit)
  
  return(formula_to_fit)
}

test_model <- function(model_object) {
  #' @description Test model for proportional hazards assumption.
  #' @param model_object Fitted model object from coxme
  
  # Test for PH assumption.
  ph_test <- cox.zph(model_object)
  # Get diagnostic plots.
  diag_plots <- ggcoxzph(ph_test)
  
  return(list('ph_test' = ph_test,
              'diag_plots' = diag_plots))
}

fit_model <-
  function(row_info,
           kinship_matrix,
           df,
           stratified = FALSE,
           id_col = 'RANID') {
    #' @description Fit coxme models.
    #' @param row_info Data frame (see `create_param_grid`).
    #' @param kinship_matrix Matrix with kinship coding.
    #' @param df Data frame with Framingham data.
    #' @param stratified Default FALSE. If TRUE, then search for stratifying variable.
    #' @return Summary of results from coxme model.
    
    
    if (stratified) {
      predictors <- row_info %>%
        dplyr::select(-Table,-Model,-Survival,-Outcome,-Row)
      # Effect modification variables are denoted [variable]:[group]
      strat_var <- names(predictors)[grep(':', predictors)]
      # print(strat_var)
      tmp <- strsplit(predictors[, strat_var], split = ':')
      var_name <- tmp[[1]][1]
      var_group <- tmp[[1]][2]
      if(grepl('^[[:digit:]]$', var_group)) {
        var_group <- as.numeric(var_group)
      }

      df <- df[which(df[, var_name] == var_group), ]
      print(var_group)
    }
    
    formula_to_fit <- format_model_formula(row_info, stratified = stratified)
    
    print(formula_to_fit)
    
    ids_to_select <- df %>% dplyr::pull(!!id_col) %>% as.character()
    kmat_subset <-
      kinship_matrix[which(rownames(kinship_matrix) %in% ids_to_select),
                     which(colnames(kinship_matrix) %in% ids_to_select)]
    
    mdl <-
      coxme(formula_to_fit,
            df,
            varlist = list(kmat_subset))
    # print(mdl)
    # print(class(mdl))
    
    # Model diagnostics.
    # mdl_diag <- cox.zph(mdl)
    # print(mdl_diag)
    # output <- summary_output(mdl)
    # print(output)
    # output$Table <- rep(row_info$Table, nrow(output))
    # output$Model <- rep(row_info$Model, nrow(output))

    # output <- output %>%
    #   relocate(Table, Model, Variable, `Hazard Ratio`)
    
    return(mdl)
    # return(output)
    # return(list('mdl' = mdl,
    #             'mdl_diag' = mdl_diag,
    #             'output' = output))
  }

tidy_output <- function(mdl, row_info) {
  #' @description Given coxme model, clean and format output.
  #' @param mdl Output from coxme model.
  #' @param row_info Data frame row with model parameters.

  output <- summary_output(mdl)
  output$Table <- rep(row_info$Table, nrow(output))
  output$Model <- rep(row_info$Model, nrow(output))
  output <- output %>%
    relocate(Table, Model, Variable, `Hazard Ratio`)
  
  col <- 'Variable'
  
  # Age variables.
  output[, col] <- gsub('AgeInit', 'Age at Initial Visit', output[, col])
  output[, col] <- gsub('AgeBase', 'Age at Baseline Visit', output[, col])

  # BMD/Bone Loss variables.
  output[, col] <- gsub('BoneLoss', 'Bone Loss', output[, col])
  output[, col] <- gsub('BMDInit', 'Initial BMD', output[, col])
  output[, col] <- gsub('BMDBase', 'Baseline BMD', output[, col])
  output[, col] <- gsub('DeltaBMDStd', '\u2206BMD (Standardized)', output[, col])
  output[, col] <- gsub('HighestQuartileLabels', ' ', output[, col])
  output[, col] <- gsub('LowestQuartileLabels', ' ', output[, col])
 
  # BMI variables.
  output[, col] <- gsub('BMIInit', 'BMI at Initial Visit', output[, col])
  output[, col] <- gsub('BMIBase', 'BMI at Baseline Visit', output[, col])
  
  # Other variables.
  output[, col] <- gsub('APOEStatus', '\u22651 APOE \u03B54 allele', output[, col])
  output[, col] <- gsub('MobilityYes Help', 'Need Help Transferring', output[, col])
  output[, col] <- gsub('SexMale', 'Male', output[, col])
  
  output[, col] <- gsub('EstrogenUseFemale No Estrogen Use', 'Female No Estrogen Use', output[, col])
  output[, col] <- gsub('EstrogenUseMale', 'Male', output[, col])
  
  # output[, col] <- gsub('', '', output[, col])
  # output[, col] <- gsub('', '', output[, col])
  # output[, col] <- gsub('', '', output[, col])
  # output[, col] <- gsub('', '', output[, col])
  # output[, col] <- gsub('', '', input_df[, col])
    
  return(output)
}