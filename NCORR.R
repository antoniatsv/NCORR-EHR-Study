
########################################################################################################################



library(here)
library(readxl)
library(dplyr)
library(ggplot2)
library(pROC)
library(furniture)
library(tidyverse)
library(mice)
library(haven)
library(mice)
library(purrr)
library(predtools)

set.seed(123)

# Read SPSS data
thoracic_raw <- read_sav("6600_anonymous.sav")

colnames(thoracic_raw)


# Set BMI values less than 14 and over 64 to missing
thoracic_raw$BMI[thoracic_raw$BMI < 14 | thoracic_raw$BMI > 64] <- NA

# Set CreatinineumolL values less than 11.30 and over 654 to missing
thoracic_raw$CreatinineumolL[thoracic_raw$CreatinineumolL < 11.3 | thoracic_raw$CreatinineumolL > 654] <- NA

# Set CreatinineumolL values less than 18 and over 187 to missing
thoracic_raw$DLCOPredicted[thoracic_raw$DLCOPredicted < 18 | thoracic_raw$DLCOPredicted > 187] <- NA

# Replace the blanks in ECOG with NA (NOTE: the blank is a character "") check levels again 
thoracic_raw$ECOG <- as.factor(ifelse(thoracic_raw$ECOG == "", NA, as.character(thoracic_raw$ECOG)))


# Replace the blanks in ASA with NA (NOTE: the blank is a character "") check levels again 
thoracic_raw$ASA <- as.factor(ifelse(thoracic_raw$ASA == "", NA, as.character(thoracic_raw$ASA)))


# Replace the blanks in Dyspnoea with NA (NOTE: the blank is a character "") check levels again 
thoracic_raw$Dyspnoea <- as.factor(ifelse(thoracic_raw$Dyspnoea == "", NA, as.character(thoracic_raw$Dyspnoea)))

# Make sure the variables are assigned properly - i.e. continuous or categorical - ?? 
thoracic_raw$Age <- as.numeric(thoracic_raw$Age)

thoracic_raw$MaleSex <- as.factor(thoracic_raw$MaleSex)

thoracic_raw$Anaemia <- as.factor(thoracic_raw$Anaemia)

thoracic_raw$Arrhythmia <- as.factor(thoracic_raw$Arrhythmia)

thoracic_raw$Right <- as.factor(thoracic_raw$Right)

thoracic_raw$ResectedSegments <- as.numeric(thoracic_raw$ResectedSegments)

thoracic_raw$Thoracotomy <- as.factor(thoracic_raw$Thoracotomy)

thoracic_raw$Malignant <- as.factor(thoracic_raw$Malignant)

thoracic_raw$Deadat90days <- as.factor(thoracic_raw$Dead31to90Days)

thoracic_raw$DyspnoeaGroups <- as.factor(thoracic_raw$DyspnoeaGroups)

# Recode the levels to 0 for "Elective" and 1 for "Urgent"
thoracic_raw$Urgency <- ifelse(thoracic_raw$Urgency == "Elective", 0, 1)

thoracic_raw$Urgency <- as.factor(thoracic_raw$Urgency)

thoracic_raw$Pneumonectomy <- as.factor(thoracic_raw$Pneumonectomy)

thoracic_raw$ComorbidityScoreGroups <- as.factor(thoracic_raw$ComorbidityScoreGroups)

thoracic_raw$ComorbidityScoreGroups <- as.factor(thoracic_raw$ComorbidityScore1and2)

thoracic_raw$ComorbidityScoreGroups <- as.factor(thoracic_raw$ComorbidityScore3andAbove)


thoracic_raw$DeadatDischarge <- as.factor(thoracic_raw$DeadatDischarge)

########################################################################################################################

# Select only relevant variables to be included in the dataset

resect_vars <- c("Age", # 
                 "MaleSex", # 
                 "ECOG", # 
                 "DLCOPredicted",
                 "FEV1Predicted",
                 "BMI",
                 "CreatinineumolL",
                 "Anaemia",
                 "Arrhythmia",
                 "Right",
                 "ResectedSegments",
                 "Thoracotomy",
                 "Malignant",
                 "Deadat90days") # 


thoracoscore_vars <- c("Age55to65",
                       "AgeOver65", 
                       "MaleSex", # 
                       "ASA",
                       "ECOG", # 
                       "Dyspnoea",
                       "Urgency",
                       "Pneumonectomy",
                       "Malignant", 
                       "ComorbidityScore1and2",
                       "ComorbidityScore3andAbove",
                       "DeadatDischarge") #is this the equivalent to "in hospital mortality"?



df <- select(thoracic_raw, resect_vars, thoracoscore_vars) %>% 
  mutate(ID = 1:nrow(.))

col_names <- colnames(df) #assign col_names 

id_col_index <- which(col_names == "ID");  #identify the indexes of the ID column 

df <- df[, c(id_col_index, setdiff(seq_along(col_names), id_col_index))] 

df$ID <- sample(df$ID)


df_val <- df %>%
  filter(ID %in% sample(ID, nrow(df) / 2))

df_val_resect <- df_val %>% 
  select(-ASA, -Dyspnoea, -Urgency, -Pneumonectomy, -ComorbidityScore1and2, -ComorbidityScore3andAbove, -DeadatDischarge)

df_val_thoracoscore <- df_val %>% 
  select(-DLCOPredicted, -BMI, -CreatinineumolL, -Anaemia, -Arrhythmia, -Right, -ResectedSegments, -Thoracotomy, -Deadat90days)

df_imp <- df %>%
  filter(!(ID %in% df_val$ID))

df_imp_resect <- df_imp %>% 
  select(-ASA, -Dyspnoea, -Urgency, -Pneumonectomy, -ComorbidityScore1and2, -ComorbidityScore3andAbove, -DeadatDischarge)

df_imp_thoracoscore <- df_imp %>% 
  select(-DLCOPredicted, -BMI, -CreatinineumolL, -Anaemia, -Arrhythmia, -Right, -ResectedSegments, -Thoracotomy, -Deadat90days)

#### CCA 
CCA_function <- function(df) {
  df[complete.cases(df), ]
}

### MEAN ZERO IMP
mean_zero_imputation <- function(df) {
  for (col in names(df)) {
    if (is.numeric(df[[col]])) {
      df[[col]][is.na(df[[col]])] <- mean(df[[col]], na.rm = TRUE)
    } else if (is.factor(df[[col]]) || is.character(df[[col]])) {
      levels <- levels(df[[col]])
      lowest_level <- min(levels, na.rm = TRUE)
      df[[col]][is.na(df[[col]])] <- lowest_level
    }
  }
  return(df)
}


### MICE 
mice_function <- function(df, m = m, outcome_var = NULL, include_outcome) {
  
  print(is.data.frame(df))
  dummy_run <- mice(df, m = 1, maxit = 0)
  predmat <- dummy_run$predictorMatrix
  
  if (include_outcome == FALSE) {
    predmat[outcome_var, ] <- 0
    predmat[, outcome_var] <- 0   
  } 
  
  predmat[,"ID"] <- 0  
  predmat["ID",] <- 0
  
  print(predmat)
  print(dummy_run$method)
  
  method <- mice(df, predmat = predmat, m = m, print = FALSE)
  
  return(method)
}


### 'MASTER' IMPUTATION FUNCTION
imputation_function <- function(df = df, m = m, outcome_var, include_outcome) {
  
  MI_noY_val_resect <- mice_function(df = df_val_resect, m = m, include_outcome = FALSE)
  MI_withY_val_resect <- mice_function(df = df_val_resect, m = m, include_outcome = TRUE, outcome_var = "Deadat90days")
  MI_noY_imp_resect <- mice_function(df = df_imp_resect, m = m, include_outcome = FALSE)
  
  MI_noY_val_thoracoscore <- mice_function(df = df_val_thoracoscore, m = m, include_outcome = FALSE)
  MI_withY_val_thoracoscore <- mice_function(df = df_val_thoracoscore, m = m, include_outcome = TRUE, outcome_var = "DeadatatDischarge")
  MI_noY_imp_thoracoscore <- mice_function(df = df_imp_thoracoscore, m = m, include_outcome = FALSE)
  
  CCA_val_resect <- CCA_function(df = df_val_resect)
  CCA_imp_resect <- CCA_function(df = df_imp_resect)
  
  CCA_val_thoracoscore <- CCA_function(df = df_val_thoracoscore)
  CCA_imp_thoracoscore <- CCA_function(df = df_imp_thoracoscore)
  
  mean_zero_val_resect <- mean_zero_imputation(df = df_val_resect)
  mean_zero_imp_resect <- mean_zero_imputation(df = df_imp_resect)
  
  mean_zero_val_thoracoscore <- mean_zero_imputation(df = df_val_thoracoscore)
  mean_zero_imp_thoracoscore <- mean_zero_imputation(df = df_imp_thoracoscore)
  
  
  return(list(
    "MI_noY_val_resect" = MI_noY_val_resect,
    "MI_withY_val_resect" = MI_withY_val_resect,
    "MI_noY_imp_resect" = MI_noY_imp_resect, 
    "MI_noY_val_thoracoscore" = MI_noY_val_thoracoscore, 
    "MI_withY_val_thoracoscore" = MI_withY_val_thoracoscore, 
    "MI_noY_imp_thoracoscore" = MI_noY_imp_thoracoscore,
    "CCA_val_resect" = CCA_val_resect,
    "CCA_imp_resect" = CCA_imp_resect, 
    "CCA_val_thoracoscore" = CCA_val_thoracoscore,
    "CCA_imp_thoracoscore" = CCA_imp_thoracoscore,
    "mean_zero_val_resect" = mean_zero_val_resect, 
    "mean_zero_imp_resect" = mean_zero_imp_resect,
    "mean_zero_val_thoracoscore" = mean_zero_val_thoracoscore,
    "mean_zero_imp_thoracoscore" = mean_zero_imp_thoracoscore
  ))
  
}

imputed_datasets <- imputation_function(df, m = 5, outcome_var, include_outcome)


imputed_datasets <- imputed_datasets %>% 
  map_if(grepl("MI", names(.)), mice::complete, action = "long")

########################################################################################################################
################################ External validation of Resect and Thoracoscore ########################################


################################ Calculate LP and Pi  ########################################

resect_datasets_LP <- imputed_datasets[grepl("resect", names(imputed_datasets))] 
thoracoscore_datasets_LP <- imputed_datasets[grepl("thoracoscore", names(imputed_datasets))] 


for (dataset_name in names(resect_datasets_LP)) {
  dataset <- resect_datasets_LP[[dataset_name]]
  
  # Calculate LP
  LP <- -6.036 +
    (as.numeric(dataset$Age) * 0.041) +
    (as.numeric(dataset$MaleSex) * 0.493) +
    (as.numeric(dataset$ECOG) * 0.183) +
    (as.numeric(dataset$DLCOPredicted) * 0.029) -
    (as.numeric(dataset$BMI) * 0.056) +
    (as.numeric(dataset$CreatinineumolL) * 0.005) +
    (as.numeric(dataset$Anaemia) * 0.242) +
    (as.numeric(dataset$Arrhythmia) * 0.608) +
    (as.numeric(dataset$Right) * 0.379) +
    (as.numeric(dataset$ResectedSegments) * 0.179) +
    (as.numeric(dataset$Thoracotomy) * 0.634) +
    (as.numeric(dataset$Malignant) * 0.769)
  
  
  # Calculate Pi
  Pi <- exp(LP) / (1 + exp(LP))
  
  # Add LP and Pi as new columns
  dataset$LP <- LP
  dataset$Pi <- Pi
  resect_datasets_LP[[dataset_name]] <- dataset
  
}

# Calculate LP and Pi for datasets with 'thoracoscore' in their name
for (dataset_name in names(thoracoscore_datasets_LP)) {
  dataset <- thoracoscore_datasets_LP[[dataset_name]]
  
  # Calculate LP
  LP <- -7.3737 +
    (as.numeric(dataset$Age55to65) * 0.7679) + 
    (as.numeric(dataset$AgeOver65) * 1.0073) + 
    (as.numeric(dataset$MaleSex) * 0.4505) +
    (as.numeric(dataset$ASA) * 0.6057) +
    (as.numeric(dataset$ECOG) * 0.6890) +
    (as.numeric(dataset$Dyspnoea) * 0.9075) +
    (as.numeric(dataset$Urgency) * 0.8443) +
    (as.numeric(dataset$Pneumonectomy) * 1.2176) +
    (as.numeric(dataset$Malignant) * 1.2423) + 
    (as.numeric(dataset$ComorbidityScore1and2) * 0.7447) + 
    (as.numeric(dataset$ComorbidityScore3andAbove) * 0.9065)
  
  # Calculate Pi
  Pi <- exp(LP) / (1 + exp(LP))
  
  # Add LP and Pi as new columns
  dataset$LP <- LP
  dataset$Pi <- Pi
  
  thoracoscore_datasets_LP[[dataset_name]] <- dataset
}



# Take average LP for each patient across the m imputations for the MI resect datasets 
#LP_mean <- P %>% 
#  group_by(ID) %>% 
#  summarise(LP_mean = mean(LP))



################################ Calculate target measures  ########################################
calc_target_measures_function <- function(datasets, outcome_var = NULL, Pi) {
  library(pROC)
  
  outcome_var = c("Deadat90days", "DeadatDischarge")
  
  ### Brier ### 
  ####----------------------------------------------------
  Brier_individuals <- (Pi - outcome_var)^2
  Brier <- mean(Brier_individuals)
  
  ## Calibration intercept (i.e. calibration-in-the-large)
  ####-----------------------------------------------------
  LP <- log(Pi/ (1 - Pi)) 
  Cal_Int <- glm(Y ~ offset(LP), family = binomial(link = "logit"))
  
  ## Calibration slope
  ####-------------------------------------------------------
  Cal_Slope <- glm(outcome_var ~ LP, family = binomial(link = "logit"))
  
  ## Calibration plot for each dataset
  ####-------------------------------------------------------
  calibration_plot(data = datasets, 
                   obs = "Pi", 
                   pred = "outcome_var", 
                   title = "Calibration plot")
  
  ## Discrimination 
  ####--------------------------------------------------------
  AUC <- roc(response = outcome_var, 
             predictor = as.vector(Pi),
             direction = "<",
             levels = c(0,1))$auc
  
  
  
  ## Store performance results in a data.frame and return
  ####------------------------------------------------------------------------
  target_measures <- data.frame("Cal_Int" = as.numeric(coef(Cal_Int)),
                                "Cal_Slope" = as.numeric(coef(Cal_Slope)[2]),
                                "AUC" = as.numeric(AUC),
                                "Brier" = as.numeric(Brier)
  )
  
  
  return(target_measures)
  
}






