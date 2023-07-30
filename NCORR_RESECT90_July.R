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

setwd("~/Desktop/Antonia/Data")

set.seed(123)

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
  MI_withY_imp_resect <- mice_function(df = df_imp_resect, m = m, include_outcome = TRUE, outcome_var = "Deadat90days")
  
  MI_noY_val_thoracoscore <- mice_function(df = df_val_thoracoscore, m = m, include_outcome = FALSE)
  MI_withY_val_thoracoscore <- mice_function(df = df_val_thoracoscore, m = m, include_outcome = TRUE, outcome_var = "DeadatatDischarge")
  MI_noY_imp_thoracoscore <- mice_function(df = df_imp_thoracoscore, m = m, include_outcome = FALSE)
  MI_withY_imp_thoracoscore <- mice_function(df = df_imp_thoracoscore, m = m, include_outcome = TRUE, outcome_var = "DeadatatDischarge")
    
  CCA_val_resect <- CCA_function(df = df_val_resect)
  #CCA_imp_resect <- CCA_function(df = df_imp_resect)
    
  CCA_val_thoracoscore <- CCA_function(df = df_val_thoracoscore)
  #CCA_imp_thoracoscore <- CCA_function(df = df_imp_thoracoscore)
  
  mean_zero_val_resect <- mean_zero_imputation(df = df_val_resect)
  mean_zero_imp_resect <- mean_zero_imputation(df = df_imp_resect)
    
  mean_zero_val_thoracoscore <- mean_zero_imputation(df = df_val_thoracoscore)
  mean_zero_imp_thoracoscore <- mean_zero_imputation(df = df_imp_thoracoscore)
  

  return(list(
    "MI_noY_val_resect" = MI_noY_val_resect,
    "MI_withY_val_resect" = MI_withY_val_resect,
    "MI_noY_imp_resect" = MI_noY_imp_resect, 
    "MI_withY_imp_resect" = MI_withY_imp_resect,
    "MI_noY_val_thoracoscore" = MI_noY_val_thoracoscore, 
    "MI_withY_val_thoracoscore" = MI_withY_val_thoracoscore, 
    "MI_noY_imp_thoracoscore" = MI_noY_imp_thoracoscore,
    "MI_withY_imp_thoracoscore" = MI_withY_imp_thoracoscore,
    "CCA_val_resect" = CCA_val_resect,
    #"CCA_imp_resect" = CCA_imp_resect, 
    "CCA_val_thoracoscore" = CCA_val_thoracoscore,
    #"CCA_imp_thoracoscore" = CCA_imp_thoracoscore,
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
  datasetR <- resect_datasets_LP[[dataset_name]]
  
  # Calculate LP
  LP <- -6.036 +
    (as.numeric(datasetR$Age) * 0.041) +
    (as.numeric(datasetR$MaleSex) * 0.493) +
    (as.numeric(datasetR$ECOG) * 0.183) -
    (as.numeric(datasetR$DLCOPredicted) * 0.029) -
    (as.numeric(datasetR$BMI) * 0.056) +
    (as.numeric(datasetR$CreatinineumolL) * 0.005) +
    (as.numeric(datasetR$Anaemia) * 0.242) +
    (as.numeric(datasetR$Arrhythmia) * 0.608) +
    (as.numeric(datasetR$Right) * 0.379) +
    (as.numeric(datasetR$ResectedSegments) * 0.179) +
    (as.numeric(datasetR$Thoracotomy) * 0.634) +
    (as.numeric(datasetR$Malignant) * 0.769)

  
  # Calculate Pi
  Pi <- exp(LP) / (1 + exp(LP))
  
  # Add LP and Pi as new columns
  datasetR$LP <- LP
  datasetR$Pi <- Pi
  resect_datasets_LP[[dataset_name]] <- datasetR
  
}



# Calculate LP and Pi for datasets with 'thoracoscore' in their name
for (dataset_name in names(thoracoscore_datasets_LP)) {
  datasetT <- thoracoscore_datasets_LP[[dataset_name]]
  
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
  
  # Add LP and Pi as new columns
  datasetT$LP <- LP
  datasetT$Pi <- Pi
  
  thoracoscore_datasets_LP[[dataset_name]] <- datasetT
}




# Take average LP for each patient across the m imputations for the MI resect datasets 
#LP_mean <- LP %>% 
#  group_by(ID) %>% 
#  summarise(LP_mean = mean(LP))



################################ Calculate target measures  ########################################

library(pROC)


########################################## RESECT-90 ########################################
# Create an empty df to store the results
target_measuresR <- data.frame()

for (dataset_name in names(resect_datasets_LP)) {
  datasetR <- resect_datasets_LP[[dataset_name]]
  
  # Specify the outcome variable
  outcome_var <- as.numeric(datasetR$Deadat90days) #set it up as.numeric at the start, as I was getting errors but I think something's gone wrong now
  
  # Specify the predicted probabilities
  Pi <- datasetR$Pi
  
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
  target_measuresR <- rbind(target_measuresR, measures)
}



########################################## Thoracoscore ########################################
# Create an empty df to store the results
target_measuresT <- data.frame()

for (dataset_name in names(thoracoscore_datasets_LP)) {
  datasetT <- thoracoscore_datasets_LP[[dataset_name]]
  
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
  target_measuresT <- rbind(target_measuresT, measures)
}



# Function to create and print calibration plots for each dataset
create_calibration_plots <- function(dataset_list, dataset_names, outcome_var) {
  for (i in seq_along(dataset_list)) {
    dataset_name <- names(dataset_list)[i]
    plot_title <- dataset_names[i]
    
    # Extract necessary data from the dataset
    outcome_var_values <- dataset_list[[dataset_name]][[outcome_var]]
    Pi <- dataset_list[[dataset_name]]$Pi
    
    # Fit the spline model
    spline_model <- stats::glm(outcome_var_values ~ splines::ns(LP, df = 3),
                               data = dataset_list[[dataset_name]],
                               family = stats::binomial(link = "logit"))
    
    # Predict with spline model
    spline_preds <- stats::predict(spline_model, type = "response", se = TRUE)
    
    # Create the data frame for the calibration plot
    plot_data <- data.frame("p" = Pi,
                            "o" = spline_preds$fit)
    
    # Create the calibration plot using ggplot2
    calibration_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = p, y = o)) +
      ggplot2::geom_line(ggplot2::aes(linetype = "Calibration Curve", colour = "Calibration Curve")) +
      ggplot2::geom_abline(ggplot2::aes(intercept = 0, slope = 1, linetype = "Reference", colour = "Reference"), show.legend = FALSE) +
      ggplot2::geom_point(alpha = 0) +
      ggplot2::coord_fixed() +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::labs(color = "Guide name", linetype = "Guide name") +
      ggplot2::scale_linetype_manual(values = c("dashed", "solid"), breaks = c("Reference", "Calibration Curve"), labels = c("Reference", "Calibration Curve")) +
      ggplot2::scale_colour_manual(values = c("black", "blue"), breaks = c("Reference", "Calibration Curve")) +
      ggplot2::theme(legend.title = ggplot2::element_blank()) +
      ggplot2::ggtitle(plot_title)
    
    # Print the calibration plot
    print(calibration_plot)
  }
}

# List of dataset names and their new names for Resect datasets
resect_datasets <- list(MI_noY_val_resect = resect_datasets_LP$MI_noY_val_resect, 
                        MI_withY_val_resect = resect_datasets_LP$MI_withY_val_resect,
                        MI_noY_imp_resect = resect_datasets_LP$MI_noY_imp_resect,
                        MI_withY_imp_resect = resect_datasets_LP$MI_withY_imp_resect,
                        CCA_val_resect = resect_datasets_LP$CCA_val_resect,
                        mean_zero_val_resect = resect_datasets_LP$mean_zero_val_resect,
                        mean_zero_imp_resect = resect_datasets_LP$mean_zero_imp_resect)

new_namesR <- c("Resect: MI no Y at validation", "Resect: MI with Y at validation",
                "Resect: MI no Y at implementation", "Resect: MI with Y at implementation", "Resect: CCA at validation",
                "Resect: Mean + Risk Factor Absent at validation",
                "Resect: Mean + Risk Factor Absent at implementation")

# List of dataset names and their new names for Thoracoscore datasets
thoracoscore_datasets <- list(MI_noY_val_thoracoscore = thoracoscore_datasets_LP$MI_noY_val_thoracoscore, 
                              MI_withY_val_thoracoscore = thoracoscore_datasets_LP$MI_withY_val_thoracoscore,
                              MI_noY_imp_thoracoscore = thoracoscore_datasets_LP$MI_noY_imp_thoracoscore,
                              MI_withY_imp_thoracoscore = thoracoscore_datasets_LP$MI_withY_imp_thoracoscore,
                              CCA_val_thoracoscore = thoracoscore_datasets_LP$CCA_val_thoracoscore, 
                              mean_zero_val_thoracoscore = thoracoscore_datasets_LP$mean_zero_val_thoracoscore,
                              mean_zero_imp_thoracoscore = thoracoscore_datasets_LP$mean_zero_imp_thoracoscore)

new_namesT <- c("Thoracoscore: MI no Y at validation", "Thoracoscore: MI with Y at validation",
                "Thoracoscore: MI no Y at implementation", "Thoracoscore: MI with Y at implementation", "Thoracoscore: CCA at validation", 
                "Thoracoscore: Mean + Risk Factor Absent at validation",
                "Thoracoscore: Mean + Risk Factor Absent at implementation")

# Example usage for Resect datasets:
create_calibration_plots(dataset_list = resect_datasets, 
                         dataset_names = new_namesR,
                         outcome_var = "Deadat90days")

# Example usage for Thoracoscore datasets:
create_calibration_plots(dataset_list = thoracoscore_datasets,
                         dataset_names = new_namesT,
                         outcome_var = "DeadatDischarge")


### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### 

### data manipulation before bias calc 

#turn data into long format and rename TM
target_measuresR_long <- target_measuresR %>%
  pivot_longer(cols = c("Cal_Int", "Cal_Slope", "AUC", "Brier"), names_to = "target_measures", values_to = "estimates") %>%
  select(Dataset, target_measures, estimates) %>%
  mutate(target_measures = recode(target_measures,
                                  Cal_Int = 'Calibration Intercept',
                                  Cal_Slope = 'Calibration Slope',
                                  AUC = 'AUC',
                                  Brier = 'Brier Score')) %>% 
  mutate(Dataset = case_when(
    Dataset == "MI_noY_val_resect" ~ "MI no Y at validation",
    Dataset == "MI_withY_val_resect" ~ "MI with Y at validation",
    Dataset == "MI_noY_imp_resect" ~ "MI no Y at implementation",
    Dataset == "MI_withY_imp_resect" ~ "MI with Y at implementation",
    Dataset == "CCA_val_resect" ~ "CCA at validation",
    Dataset == "mean_zero_val_resect" ~ "Mean + Risk Factor Absent at validation",
    Dataset == "mean_zero_imp_resect" ~ "Mean + Risk Factor Absent at implementation",
    TRUE ~ Dataset
  ))



#split the data into val and imp 
R_val <- target_measuresR_long[target_measuresR_long$Dataset %like% "validation", ]
R_imp <- target_measuresR_long[target_measuresR_long$Dataset %like% "implementation", ]



#############################################################################
#############################################################################

#turn data into long format and rename TM
target_measuresT_long <- target_measuresT %>%
  pivot_longer(cols = c("Cal_Int", "Cal_Slope", "AUC", "Brier"), names_to = "target_measures", values_to = "estimates") %>%
  select(Dataset, target_measures, estimates) %>%
  mutate(target_measures = recode(target_measures,
                                  Cal_Int = 'Calibration Intercept',
                                  Cal_Slope = 'Calibration Slope',
                                  AUC = 'AUC',
                                  Brier = 'Brier Score')) %>% 
  mutate(Dataset = case_when(
    Dataset == "MI_noY_val_thoracoscore" ~ "MI no Y at validation",
    Dataset == "MI_withY_val_thoracoscore" ~ "MI with Y at validation",
    Dataset == "MI_noY_imp_thoracoscore" ~ "MI no Y at implementation",
    Dataset == "MI_withY_imp_thoracoscore" ~ "MI with Y at implementation",
    Dataset == "CCA_val_thoracoscore" ~ "CCA at validation",
    Dataset == "mean_zero_val_thoracoscore" ~ "Mean + Risk Factor Absent at validation",
    Dataset == "mean_zero_imp_thoracoscore" ~ "Mean + Risk Factor Absent at implementation",
    TRUE ~ Dataset
  ))



#split the data into val and imp 
T_val <- target_measuresT_long[target_measuresT_long$Dataset %like% "validation", ]
T_imp <- target_measuresT_long[target_measuresT_long$Dataset %like% "implementation", ]



####### BIAS CALCULATIONS #########
### RESECT 

## MI no Y at implementation 
MInoY_bias_R <- subset(R_imp, Dataset == "MI no Y at implementation") %>%
  rename_at(vars("estimates"), function(x) paste0("true_", x)) %>%
  left_join(R_val, MInoY_bias_R,
            multiple = "all",
            by = "target_measures") %>%
            mutate(bias = estimates - true_estimates) 

# MI with Y at implementation 
MIwithY_bias_R <- subset(R_imp, Dataset == "MI with Y at implementation") %>%
  rename_at(vars("estimates"), function(x) paste0("true_", x)) %>%
  left_join(R_val, MIwithY_bias_R,
            multiple = "all",
            by = "target_measures") %>%
  mutate(bias = estimates - true_estimates) 

# Mean + Risk Factor Absent at Implementation
Mean_Zero_R <- subset(R_imp, Dataset == "Mean + Risk Factor Absent at implementation") %>%
  rename_at(vars("estimates"), function(x) paste0("true_", x)) %>%
  left_join(R_val, Mean_Zero_R,
            multiple = "all",
            by = "target_measures") %>%
  mutate(bias = estimates - true_estimates) 


all_bias_R <- MInoY_bias_R %>% 
  bind_rows(MIwithY_bias_R) %>%
  bind_rows(Mean_Zero_R)


### Thoracoscore

## MI no Y at implementation 
MInoY_bias_T <- subset(T_imp, Dataset == "MI no Y at implementation") %>%
  rename_at(vars("estimates"), function(x) paste0("true_", x)) %>%
  left_join(T_val, MInoY_bias_T,
            multiple = "all",
            by = "target_measures") %>%
  mutate(bias = estimates - true_estimates) 

# MI with Y at implementation 
MIwithY_bias_T <- subset(T_imp, Dataset == "MI with Y at implementation") %>%
  rename_at(vars("estimates"), function(x) paste0("true_", x)) %>%
  left_join(T_val, MIwithY_bias_T,
            multiple = "all",
            by = "target_measures") %>%
  mutate(bias = estimates - true_estimates) 

# Mean + Risk Factor Absent at Implementation
Mean_Zero_T <- subset(T_imp, Dataset == "Mean + Risk Factor Absent at implementation") %>%
  rename_at(vars("estimates"), function(x) paste0("true_", x)) %>%
  left_join(T_val, Mean_Zero_T,
            multiple = "all",
            by = "target_measures") %>%
  mutate(bias = estimates - true_estimates) 

all_bias_T <- MInoY_bias_T %>% 
  bind_rows(MIwithY_bias_T) %>%
  bind_rows(Mean_Zero_T)


#### plot bias 
all_bias_R$target_measures <- factor(all_bias_R$target_measures, levels = c("AUC", "Calibration Intercept", "Calibration Slope", "Brier Score"))
all_bias_T$target_measures <- factor(all_bias_T$target_measures, levels = c("AUC", "Calibration Intercept", "Calibration Slope", "Brier Score"))

plot_R <- ggplot(data = all_bias_R, aes(x = bias, y = Dataset.y, color = factor(target_measures),
                                      shape = factor(target_measures))) +
  geom_point(size = 3, stroke = 0.5) +
  guides(color = guide_legend(reverse = TRUE)) +
  scale_shape_manual(values = c(8, 17, 16, 15)) +
  scale_color_brewer(palette = "Set1") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  xlab("Bias") +
  ylab("Validation Data Imputation Methods") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    strip.text = element_text(size = 16),
    panel.background = element_rect(fill = "gray90"),
    panel.spacing.x = unit(0.5, "lines")
  ) +
  ggh4x::facet_grid2(target_measures ~ Dataset.x, scales = "free_x", independent = "x") +
  scale_x_continuous(limits = function(x) c(-max(abs(x)), max(abs(x)))) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1.5),
    strip.text = element_text(size = 14, hjust = 0.5),
    strip.placement = "outside"
  ) +
  ggtitle("Missingness mechanisms at model implementation") +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))

plot_R <- plot_R + theme(panel.grid.major = element_line(size = 1.5))

print(plot_R)


#####
plot_T <- ggplot(data = all_bias_T, aes(x = bias, y = Dataset.y, color = factor(target_measures),
                                        shape = factor(target_measures))) +
  geom_point(size = 3, stroke = 0.5) +
  guides(color = guide_legend(reverse = TRUE)) +
  scale_shape_manual(values = c(8, 17, 16, 15)) +
  scale_color_brewer(palette = "Set1") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  xlab("Bias") +
  ylab("Validation Data Imputation Methods") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    strip.text = element_text(size = 16),
    panel.background = element_rect(fill = "gray90"),
    panel.spacing.x = unit(0.5, "lines")
  ) +
  ggh4x::facet_grid2(target_measures ~ Dataset.x, scales = "free_x", independent = "x") +
  scale_x_continuous(limits = function(x) c(-max(abs(x)), max(abs(x)))) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1.5),
    strip.text = element_text(size = 14, hjust = 0.5),
    strip.placement = "outside"
  ) +
  ggtitle("Missingness mechanisms at model implementation") +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))

plot_T <- plot_T + theme(panel.grid.major = element_line(size = 1.5))

print(plot_T)

##

