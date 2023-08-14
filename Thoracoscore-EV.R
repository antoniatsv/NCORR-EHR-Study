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
library(tidyverse)
library(magrittr)
library(dplyr)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(ggh4x)

set.seed(123)

#iter <- 100

#all_results <- data.frame()

#for (i in 1:iter) {

# Read SPSS data
thoracic_raw <- read_sav("6600_anonymous.sav")

colnames(thoracic_raw)


# Set BMI values less than 14 and over 64 to missing
thoracic_raw$BMI[thoracic_raw$BMI < 14 | thoracic_raw$BMI > 64] <- NA

# Set CreatinineumolL values less than 11.30 and over 654 to missing
thoracic_raw$CreatinineumolL[thoracic_raw$CreatinineumolL < 11.3 | thoracic_raw$CreatinineumolL > 654] <- NA

# Set DLCOP values less than 18 and over 187 to missing
thoracic_raw$DLCOPredicted[thoracic_raw$DLCOPredicted < 18 | thoracic_raw$DLCOPredicted > 187] <- NA

# Replace the blanks in ECOG with NA (NOTE: the blank is a character "") check levels again 
thoracic_raw$ECOG <- as.factor(ifelse(thoracic_raw$ECOG == "", NA, as.character(thoracic_raw$ECOG)))


# Replace the blanks in ASA with NA (NOTE: the blank is a character "") check levels again 
thoracic_raw$ASA <- as.factor(ifelse(thoracic_raw$ASA == "", NA, as.character(thoracic_raw$ASA)))


# Replace the blanks in Dyspnoea with NA (NOTE: the blank is a character "") check levels again 
thoracic_raw$Dyspnoea <- as.factor(ifelse(thoracic_raw$Dyspnoea == "", NA, as.character(thoracic_raw$Dyspnoea)))

# Make sure the variables are assigned properly - i.e. continuous or categorical - ?? 

# Recode the levels to 0 for "Elective" and 1 for "Urgent"
thoracic_raw$Urgency <- ifelse(thoracic_raw$Urgency == "Elective", 0, 1)

########################################################################################################################

# Select only relevant variables to be included in the dataset

resect_vars <- c("Age", # 
                 "MaleSex", # 
                 "ECOG", # 
                 "DLCOPredicted",
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
                       "ECOG3orAbove", #performance score
                       "NYHA3or4", #NYHA score
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


#####
df_val_resect <- df_val %>% 
  select(-ASA, -NYHA3or4, -Pneumonectomy, -ComorbidityScore1and2, -ComorbidityScore3andAbove, -Urgency, -DeadatDischarge)

df_val_thoracoscore <- df_val %>% 
  select(-DLCOPredicted, -BMI, -CreatinineumolL, -Anaemia, -Arrhythmia, -Right, -ResectedSegments, -Thoracotomy, -Deadat90days)

df_imp <- df %>%
  filter(!(ID %in% df_val$ID))

df_imp_resect <- df_imp %>% 
  select(-ASA, -NYHA3or4, -Pneumonectomy, -ComorbidityScore1and2, -ComorbidityScore3andAbove, -Urgency, -DeadatDischarge)

df_imp_thoracoscore <- df_imp %>% 
  select(-DLCOPredicted, -BMI, -CreatinineumolL, -Anaemia, -Arrhythmia, -Right, -ResectedSegments, -Thoracotomy, -Deadat90days)


rm(df_imp)
rm(df_val)

##############################################################################################################################
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
mice_function <- function(df, m = m, outcome_var, include_outcome) {
  
  print(is.data.frame(df))
  dummy_run <- mice(df, m = m, maxit = 0)
  predmat <- dummy_run$predictorMatrix
  
  if (include_outcome == FALSE) {
    predmat[outcome_var, ] <- 0
    predmat[, outcome_var] <- 0   
  } 
  
  predmat[,"ID"] <- 0  
  predmat["ID",] <- 0
  
  print(predmat)
  print(dummy_run$method)
  
  method <- mice(df, method = dummy_run$method, predictorMatrix = predmat, m = m, print = FALSE)
  
  return(method)
}


### 'MASTER' IMPUTATION FUNCTION
imputation_function <- function(df = df, m = m) {
  
  MI_noY_val_thoracoscore <- mice_function(df = df_val_thoracoscore, m = m, include_outcome = FALSE, outcome_var = "DeadatDischarge")
  MI_withY_val_thoracoscore <- mice_function(df = df_val_thoracoscore, m = m, include_outcome = TRUE, outcome_var = "DeadatDischarge")
  MI_noY_imp_thoracoscore <- mice_function(df = df_imp_thoracoscore, m = m, include_outcome = FALSE, outcome_var = "DeadatDischarge")
  MI_withY_imp_thoracoscore <- mice_function(df = df_imp_thoracoscore, m = m, include_outcome = TRUE, outcome_var = "DeadatDischarge")
  
  CCA_val_resect <- CCA_function(df = df_val_resect)

  
  CCA_val_thoracoscore <- CCA_function(df = df_val_thoracoscore)

  
  mean_zero_val_thoracoscore <- mean_zero_imputation(df = df_val_thoracoscore)
  mean_zero_imp_thoracoscore <- mean_zero_imputation(df = df_imp_thoracoscore)
  
  
  return(list(

    "MI_noY_val_thoracoscore" = MI_noY_val_thoracoscore, 
    "MI_withY_val_thoracoscore" = MI_withY_val_thoracoscore, 
    "MI_noY_imp_thoracoscore" = MI_noY_imp_thoracoscore,
    "MI_withY_imp_thoracoscore" = MI_withY_imp_thoracoscore,
    "CCA_val_thoracoscore" = CCA_val_thoracoscore,
    "mean_zero_val_thoracoscore" = mean_zero_val_thoracoscore,
    "mean_zero_imp_thoracoscore" = mean_zero_imp_thoracoscore
  ))
  
}

imputed_datasets <- imputation_function(df, m = 5)


T_datasets <- imputed_datasets[grepl("thoracoscore", names(imputed_datasets))]


#################################################################################
#################################################################################
for (dataset_name in names(T_datasets)) {
  datasetT <- T_datasets[[dataset_name]]
  
  if (is.mids(datasetT)) {
    datasetT <- mice::complete(datasetT, action = "long")

    
    # Calculate LP
    LP <- -7.3737 +
      (as.numeric(datasetT$Age55to65) * 0.7679) + 
      (as.numeric(datasetT$AgeOver65) * 1.0073) + 
      (as.numeric(datasetT$MaleSex) * 0.4505) +
      (as.numeric(datasetT$ASA) * 0.6057) +
      (as.numeric(datasetT$ECOG3orAbove) * 0.689) +
      (as.numeric(datasetT$NYHA3or4) * 0.9075) +
      (as.numeric(datasetT$Urgency) * 0.8443) +
      (as.numeric(datasetT$Pneumonectomy) * 1.2176) +
      (as.numeric(datasetT$Malignant) * 1.2423) + 
      (as.numeric(datasetT$ComorbidityScore1and2) * 0.7447) + 
      (as.numeric(datasetT$ComorbidityScore3andAbove) * 0.9065)
    
    Pi <- exp(LP) / (1 + exp(LP))
    
    datasetT$LP <- LP
    datasetT$Pi <- Pi
    
  }  else {
    
    # Calculate LP
    LP <- -7.3737 +
      (as.numeric(datasetT$Age55to65) * 0.7679) + 
      (as.numeric(datasetT$AgeOver65) * 1.0073) + 
      (as.numeric(datasetT$MaleSex) * 0.4505) +
      (as.numeric(datasetT$ASA) * 0.6057) +
      (as.numeric(datasetT$ECOG3orAbove) * 0.689) +
      (as.numeric(datasetT$NYHA3or4) * 0.9075) +
      (as.numeric(datasetT$Urgency) * 0.8443) +
      (as.numeric(datasetT$Pneumonectomy) * 1.2176) +
      (as.numeric(datasetT$Malignant) * 1.2423) + 
      (as.numeric(datasetT$ComorbidityScore1and2) * 0.7447) + 
      (as.numeric(datasetT$ComorbidityScore3andAbove) * 0.9065)
    
    # Calculate Pi
    Pi <- exp(LP) / (1 + exp(LP))
    
  }
  
  # Add LP and Pi as new columns
  datasetT$LP <- LP
  datasetT$Pi <- Pi
  
  T_datasets[[dataset_name]] <- datasetT
  
}





################################
# Create an empty df to store the results
target_measuresT <- data.frame()

for (dataset_name in names(T_datasets)) {
  datasetT <- T_datasets[[dataset_name]]
  
  # Specify the outcome variable
  outcome_var <- as.numeric(datasetT$DeadatDischarge) #set it up as.numeric at the start, as I was getting errors but I think something's gone wrong now
  
  # Specify the predicted probabilities
  Pi <- datasetT$Pi
  
  # Calculate Brier Score
  Brier_individuals <- (Pi - outcome_var)^2
  Brier <- mean(Brier_individuals)
  Brier_var <- var(Brier_individuals)/length(Pi)
  
  # Calculate Calibration Intercept
  LP <- log(Pi/ (1 - Pi))
  Cal_Int <- glm(outcome_var ~ offset(LP), family = binomial(link = "logit"))
  Cal_Int_var <- vcov(Cal_Int)[1,1]
  
  # Calculate Calibration Slope
  Cal_Slope <- glm(outcome_var ~ LP, family = binomial(link = "logit"))
  Cal_Slope_var <- vcov(Cal_Slope)[2,2]
  
  # Calculate AUC
  AUC <- roc(response = outcome_var, 
             predictor = as.vector(Pi),
             direction = "<",
             levels = c(0,1))$auc
  AUC_var <- var(AUC, method = "delong") #approximation method used for AUC to calculate the variance
  
  # Create a data frame with target measures for the current dataset
  measures <- data.frame("Dataset" = dataset_name,
                         "Cal_Int" = as.numeric(coef(Cal_Int)),
                         "Cal_Int_var" = Cal_Int_var,
                         "Cal_Slope" = as.numeric(coef(Cal_Slope)[2]),
                         "Cal_Slope_var" = as.numeric(Cal_Slope_var),
                         "AUC" = as.numeric(AUC),
                         "AUC_var" = as.numeric(AUC_var),
                         "Brier" = as.numeric(Brier),
                         "Brier_var" = as.numeric(Brier_var)
  )
  
  # Append the measures to the overall results data frame
  target_measuresT <- bind_rows(target_measuresT, measures)
}









