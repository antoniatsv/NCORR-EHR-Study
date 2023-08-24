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

########################################################################################################################

# Read SPSS data
thoracic_raw <- read_sav("6600_anonymous.sav")

colnames(thoracic_raw)

########################################################################################################################

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

# If there's a # then variables are not used by either R or T models, but still included in Table 1 as patients' characteristics
CPMs_vars <- c("Age",
               "DLCOPredicted",
               "BMI",
               "CreatinineumolL",
               "ResectedSegments",
               "PPOFEV1", 
               "MaleSex",
               "ECOG",
               "Anaemia",
               "Arrhythmia",
               "Right",
               "Thoracotomy",
               "Malignant",
               "Age55to65",
               "AgeOver65",
               "ASA",
               "ECOG3orAbove",
               "NYHA3or4",
               "Urgency",
               "Pneumonectomy",
               "ComorbidityScore1and2",
               "ComorbidityScore3andAbove",
               "Deadat90days",
               "DeadatDischarge",
               "Diabetes", #
               "Hypertension", #
               "Smoking", #
               "IHD", #
               "COPD", #
               "CerebrovascularDisease", #
               "PVD" #
)


### Sample from the thoracic_raw data to create a clean df with only relevant variables that will be used in the analysis 
df <- select(thoracic_raw, CPMs_vars) %>% 
  mutate(ID = 1:nrow(.))

col_names <- colnames(df) #assign col_names 

id_col_index <- which(col_names == "ID");  #identify the indexes of the ID column 

df <- df[, c(id_col_index, setdiff(seq_along(col_names), id_col_index))] 

df$ID <- sample(df$ID)


########################################################################################################################
# Table 1 of the Results' section

# Calculate the % of missingness overall in the whole df
missingness_df <- mean(is.na(df)) * 100 

# Count missing values in each variable
missing_counts <- colSums(is.na(df))

# Calculate proportion of missing values in each variable
missing_prop <- missing_counts / nrow(df)

# Convert proportions to percentages
missing_percent <- missing_prop * 100

# Return results as a data frame
result_missingness <- data.frame(missing_percent = missing_percent)


md.pattern(df, rotate.names = TRUE)

########################################################################################################################
#### summary continuous variables ######

summary_stats <- function(df) {
  # Define continuous variables
  cont_vars <- c("Age",
                 "DLCOPredicted",
                 "BMI",
                 "CreatinineumolL",
                 "ResectedSegments",
                 "PPOFEV1" #
  )
  
  # Define categorical variables
  cat_vars <- c("MaleSex",
                "ECOG",
                "Anaemia",
                "Arrhythmia",
                "Right",
                "Thoracotomy",
                "Malignant",
                "Age55to65",
                "AgeOver65",
                "ASA",
                "ECOG3orAbove",
                "NYHA3or4",
                "Urgency",
                "Pneumonectomy",
                "ComorbidityScore1and2",
                "ComorbidityScore3andAbove",
                "Deadat90days",
                "DeadatDischarge",
                "Diabetes", #
                "Hypertension", #
                "Smoking", #
                "IHD", #
                "COPD", #
                "CerebrovascularDisease", #
                "PVD" #
  )
  
  
  
  ########################################################################################################################
  
  # Calculate summary statistics for continuous variables
  cont_summary <- sapply(df[cont_vars], function(x) c(mean(x, na.rm = TRUE),
                                                      sd(x, na.rm = TRUE),
                                                      median(x, na.rm = TRUE),
                                                      min(x, na.rm = TRUE),
                                                      max(x, na.rm = TRUE),
                                                      length(x)))
  
  # Calculate percentage of missing values for continuous variables
  cont_missing_percent <- sapply(df[cont_vars], function(x) round(sum(is.na(x)) / length(x) * 100, 2))
  
  # Combine summary statistics and missingness information for continuous variables
  cont_summary_df <- data.frame(t(cont_summary), missing_percent = cont_missing_percent)
  colnames(cont_summary_df) <- c("Mean", "SD", "Median", "Min", "Max", "N", "Missing_Percent")
  rownames(cont_summary_df) <- cont_vars
  
  # Convert categorical variables to factors
  df[cat_vars] <- lapply(df[cat_vars], factor)
  
  # Calculate summary statistics for categorical variables
  cat_summary <- sapply(df[cat_vars], function(x) c(levels(x)[which.max(table(x))],
                                                    paste(levels(x), collapse = ", "),
                                                    length(x)))
  
  # Calculate percentage of missing values for categorical variables
  cat_missing_percent <- sapply(df[cat_vars], function(x) round(sum(is.na(x)) / length(x) * 100, 2))
  
  # Combine summary statistics and missingness information for categorical variables
  cat_summary_df <- data.frame(t(cat_summary), missing_percent = cat_missing_percent)
  colnames(cat_summary_df) <- c("Mode", "Levels", "N", "Missing_Percent")
  rownames(cat_summary_df) <- cat_vars
  
  # Plot bar charts for each categorical variable
  for (var in cat_vars) {
    p <- ggplot(data = df, aes_string(x = var)) +
      geom_bar(stat = "count", fill = "steelblue") +
      labs(title = var, x = var, y = "Frequency") + 
      geom_text(stat = "count", aes(label = scales::percent(..count../sum(..count..))),
                vjust = -0.5, size = 3)
    print(p)
  }
  
  return(list(cont_summary_df = cont_summary_df, cat_summary_df = cat_summary_df))
}

# Call the summary_stats function with the resect_df dataframe
summary_tables <- summary_stats(df)

# Access the summary tables for continuous variables and categorical variables
summary_table_continuous <- summary_tables$cont_summary_df
summary_table_categorical <- summary_tables$cat_summary_df


##############################################################################################################################
# DATA IMPUTATION 

############################################
### CCA 
CCA_function <- function(df) {
  df[complete.cases(df), ]
}

############################################
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

############################################
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

############################################
### 'MASTER' IMPUTATION FUNCTION
imputation_function <- function(df, m, ID_val) {
  
  # split data into df_val and df_imp by fitering IDs
  df_val <- df %>%
    filter(ID %in% ID_val)
  
  df_imp <- df %>%
    filter(!(ID %in% df_val$ID))
  
  # split df_val data into resect and thoracoscore df_vals
  df_val_resect <- df_val %>% 
    select(-ASA, -NYHA3or4, -Pneumonectomy, -ComorbidityScore1and2, -ComorbidityScore3andAbove, -Urgency, -DeadatDischarge)
  
  df_val_thoracoscore <- df_val %>% 
    select(-DLCOPredicted, -BMI, -CreatinineumolL, -Anaemia, -Arrhythmia, -Right, -ResectedSegments, -Thoracotomy, -Deadat90days)
  
  
  # split df_imp data into resect and thoracoscore df_imps
  df_imp_resect <- df_imp %>% 
    select(-ASA, -NYHA3or4, -Pneumonectomy, -ComorbidityScore1and2, -ComorbidityScore3andAbove, -Urgency, -DeadatDischarge)
  
  df_imp_thoracoscore <- df_imp %>% 
    select(-DLCOPredicted, -BMI, -CreatinineumolL, -Anaemia, -Arrhythmia, -Right, -ResectedSegments, -Thoracotomy, -Deadat90days)
  
  
  MI_noY_val_resect <- mice_function(df = df_val_resect, m = m, include_outcome = FALSE, outcome_var = "Deadat90days")
  MI_withY_val_resect <- mice_function(df = df_val_resect, m = m, include_outcome = TRUE, outcome_var = "Deadat90days")
  MI_noY_imp_resect <- mice_function(df = df_imp_resect, m = m, include_outcome = FALSE, outcome_var = "Deadat90days")
  MI_withY_imp_resect <- mice_function(df = df_imp_resect, m = m, include_outcome = TRUE, outcome_var = "Deadat90days")
  
  MI_noY_val_thoracoscore <- mice_function(df = df_val_thoracoscore, m = m, include_outcome = FALSE, outcome_var = "DeadatDischarge")
  MI_withY_val_thoracoscore <- mice_function(df = df_val_thoracoscore, m = m, include_outcome = TRUE, outcome_var = "DeadatDischarge")
  MI_noY_imp_thoracoscore <- mice_function(df = df_imp_thoracoscore, m = m, include_outcome = FALSE, outcome_var = "DeadatDischarge")
  MI_withY_imp_thoracoscore <- mice_function(df = df_imp_thoracoscore, m = m, include_outcome = TRUE, outcome_var = "DeadatDischarge")
  
  CCA_val_resect <- CCA_function(df = df_val_resect)
  
  
  CCA_val_thoracoscore <- CCA_function(df = df_val_thoracoscore)
  
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
    "CCA_val_thoracoscore" = CCA_val_thoracoscore,
    "mean_zero_val_resect" = mean_zero_val_resect, 
    "mean_zero_imp_resect" = mean_zero_imp_resect,
    "mean_zero_val_thoracoscore" = mean_zero_val_thoracoscore,
    "mean_zero_imp_thoracoscore" = mean_zero_imp_thoracoscore
  ))
  
}

IDs_for_val <- sample(df$ID, nrow(df)/2)

imputed_datasets <- imputation_function(df = df, m = 5, ID_val = IDs_for_val)


########################################################################################################################
### External Validation of Resect-90 (lines 354 to 632)

R_datasets <- imputed_datasets[grepl("resect", names(imputed_datasets))] 

#################################################################################
for (dataset_name in names(R_datasets)) {
  datasetR <- R_datasets[[dataset_name]]
  
  if (is.mids(datasetR)) {
    datasetR <- mice::complete(datasetR, action = "long")
    
    # Calculate LP for each observation
    datasetR$LP <- -6.036 +
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
    
    LP_mean <- datasetR %>%
      dplyr::group_by(.id) %>%
      dplyr::summarise(LP_mean = mean(LP)) %>%
      pull(LP_mean)
    
    Pi <- exp(LP_mean) / (1 + exp(LP_mean))
    
    
    datasetR <- data.frame("LP" = LP_mean,
                           "Pi" = Pi,
                           "Deadat90days" = datasetR$Deadat90days[1:3300])
    
    
  }  else {
    
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
    
    
    Pi <- exp(LP) / (1 + exp(LP))
    
    # Add LP and Pi as new columns
    datasetR$LP <- LP
    datasetR$Pi <- Pi
    
  }
  
  R_datasets[[dataset_name]] <- datasetR
  
  
}


################################
# Create an empty df to store the results
target_measuresR <- data.frame()



for (dataset_name in names(R_datasets)) {
  datasetR <- R_datasets[[dataset_name]]
  
  #}
  # Specify the outcome variable
  outcome_var <- as.numeric(datasetR$Deadat90days) #set it up as.numeric at the start, as I was getting errors but I think something's gone wrong now
  
  # Specify the predicted probabilities
  Pi <- datasetR$Pi
  
  # Calculate Brier Score
  Brier_individuals <- (Pi - outcome_var)^2
  Brier <- mean(Brier_individuals)
  Brier_var <- var(Brier_individuals)/length(Pi)
  Brier_SE <- sqrt(Brier_var) 
  Brier_LCI <- Brier - (1.96*Brier_SE)
  Brier_UCI <- Brier + (1.96*Brier_SE)
  
  # Calculate Calibration Intercept
  LP <- log(Pi/ (1 - Pi))
  Cal_int_mod <- glm(outcome_var ~ offset(LP), family = binomial(link = "logit"))
  Cal_Int <- as.numeric(coef(Cal_int_mod)[1])
  Cal_Int_var <- vcov(Cal_int_mod)[1,1]
  Cal_Int_SE <- sqrt(Cal_Int_var) 
  Cal_Int_LCI <- Cal_Int - (1.96*Cal_Int_SE)
  Cal_Int_UCI <- Cal_Int + (1.96*Cal_Int_SE)
  
  
  # Calculate Calibration Slope
  Cal_Slope_mod <- glm(outcome_var ~ LP, family = binomial(link = "logit"))
  Cal_Slope <- as.numeric(coef(Cal_Slope_mod)[2])
  Cal_Slope_var <- vcov(Cal_Slope_mod)[2,2]
  Cal_Slope_SE<- sqrt(Cal_Slope_var) 
  Cal_Slope_LCI <- Cal_Slope - (1.96*Cal_Slope_SE)
  Cal_Slope_UCI <- Cal_Slope + (1.96*Cal_Slope_SE)
  
  # Calculate AUC
  
  AUC <- roc(response = outcome_var, 
             predictor = as.vector(Pi),
             direction = "<",
             levels = c(0,1))$auc
  AUC_var <- var(AUC, method = "delong") #approximation method used for AUC to calculate the variance
  AUC_SE <- sqrt(AUC_var) 
  AUC_LCI <- AUC - (1.96*AUC_SE)
  AUC_UCI <- AUC + (1.96*AUC_SE)
  
  
  # Create a data frame with target measures for the current dataset
  measures <- data.frame("Dataset" = dataset_name,
                         "Cal_Int" = Cal_Int,
                         "Cal_Int_var" = Cal_Int_var,
                         "Cal_Int_LCI" = as.numeric(Cal_Int_LCI),
                         "Cal_Int_UCI" = as.numeric(Cal_Int_UCI),
                         
                         "Cal_Slope" = Cal_Slope,
                         "Cal_Slope_var" = as.numeric(Cal_Slope_var),
                         "Cal_Slope_LCI" = as.numeric(Cal_Slope_LCI),
                         "Cal_Slope_UCI" = as.numeric(Cal_Slope_UCI),
                         
                         "AUC" = as.numeric(AUC),
                         "AUC_var" = as.numeric(AUC_var),
                         "AUC_LCI" = as.numeric(AUC_LCI),
                         "AUC_UCI" = as.numeric(AUC_UCI),
                         
                         
                         "Brier" = as.numeric(Brier),
                         "Brier_var" = as.numeric(Brier_var),
                         "Brier_LCI" = as.numeric(Brier_LCI),
                         "Brier_UCI" = as.numeric(Brier_UCI)
  )
  
  # Append the measures to the overall results data frame
  target_measuresR <- bind_rows(target_measuresR, measures)
}



#################

################################################################################

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
      ggplot2::xlim(c(0,1)) + ylim(c(0,1)) + 
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
resect_datasets <- list(MI_noY_val_resect = R_datasets$MI_noY_val_resect, 
                        MI_withY_val_resect = R_datasets$MI_withY_val_resect,
                        MI_noY_imp_resect = R_datasets$MI_noY_imp_resect,
                        MI_withY_imp_resect = R_datasets$MI_withY_imp_resect,
                        CCA_val_resect = R_datasets$CCA_val_resect,
                        mean_zero_val_resect = R_datasets$mean_zero_val_resect,
                        mean_zero_imp_resect = R_datasets$mean_zero_imp_resect)

new_namesR <- c("Resect: MI no Y at validation", 
                "Resect: MI with Y at validation",
                "Resect: MI no Y at implementation", 
                "Resect: MI with Y at implementation", 
                "Resect: CCA at validation",
                "Resect: Mean + Risk Factor Absent at validation",
                "Resect: Mean + Risk Factor Absent at implementation")


# Example usage for Resect datasets:
create_calibration_plots(dataset_list = resect_datasets, 
                         dataset_names = new_namesR,
                         outcome_var = "Deadat90days")



################################################################################
#turn data into long format and rename TM

# Pivot the first set of columns
target_measuresR_long <- target_measuresR %>%
  pivot_longer(cols = c("Cal_Int", "Cal_Slope", "AUC", "Brier"), 
               names_to = "target_measures", 
               values_to = "estimates") %>%
  select(Dataset, target_measures, estimates) %>%
  mutate(target_measures = recode(target_measures,
                                  Cal_Int = 'Calibration Intercept',
                                  Cal_Slope = 'Calibration Slope',
                                  AUC = 'AUC',
                                  Brier = 'Brier Score')) %>%
  left_join(target_measuresR %>%
              pivot_longer(cols = c("Cal_Int_LCI", "Cal_Slope_LCI", "AUC_LCI", "Brier_LCI"), 
                           names_to = "target_measures", 
                           values_to = "LCI_values") %>%
              select(Dataset, target_measures, LCI_values) %>%
              mutate(target_measures = recode(target_measures,
                                              Cal_Int_LCI = 'Calibration Intercept',
                                              Cal_Slope_LCI = 'Calibration Slope',
                                              AUC_LCI = 'AUC',
                                              Brier_LCI = 'Brier Score')),
            by = c("Dataset", "target_measures")) %>%
  left_join(target_measuresR %>%
              pivot_longer(cols = c("Cal_Int_UCI", "Cal_Slope_UCI", "AUC_UCI", "Brier_UCI"), 
                           names_to = "target_measures", 
                           values_to = "UCI_values") %>%
              select(Dataset, target_measures, UCI_values) %>%
              mutate(target_measures = recode(target_measures,
                                              Cal_Int_UCI = 'Calibration Intercept',
                                              Cal_Slope_UCI = 'Calibration Slope',
                                              AUC_UCI = 'AUC',
                                              Brier_UCI = 'Brier Score')),
            by = c("Dataset", "target_measures"))



#split the data into val and imp 
R_val <- target_measuresR_long[target_measuresR_long$Dataset %like% "val", ]
R_imp <- target_measuresR_long[target_measuresR_long$Dataset %like% "imp", ]

R_val <- R_val %>%
  mutate(Dataset = case_when(
    Dataset == "MI_noY_val_resect" ~ "MI no Y",
    Dataset == "MI_withY_val_resect" ~ "MI with Y",
    Dataset == "CCA_val_resect" ~ "CCA",
    Dataset == "mean_zero_val_resect" ~ "Mean + RFA",
    TRUE ~ Dataset
  ))


R_imp <- R_imp %>%
  mutate(Dataset = case_when(
    Dataset == "MI_noY_imp_resect" ~ "MI no Y",
    Dataset == "MI_withY_imp_resect" ~ "MI with Y",
    Dataset == "CCA_imp_resect" ~ "CCA",
    Dataset == "mean_zero_imp_resect" ~ "Mean + RFA",
    TRUE ~ Dataset
  ))


##############
## MI no Y at implementation 
MInoY_bias_R <- subset(R_imp, Dataset == "MI no Y") %>%
  rename_at(vars("estimates"), function(x) paste0("true_", x)) %>%
  rename_at(vars("LCI_values"), function(x) paste0("imp_", x)) %>%
  rename_at(vars("UCI_values"), function(x) paste0("imp_",x)) %>%
  left_join(R_val, MInoY_bias_R,
            multiple = "all",
            by = "target_measures") %>%
  mutate(bias = estimates - true_estimates,
         LCI_diff = LCI_values - imp_LCI_values,
         UCI_diff = UCI_values - imp_UCI_values) %>%
  rename("Dataset_imp" = "Dataset.x") %>%
  rename("Dataset_val" = "Dataset.y")

#####


## MI no Y at implementation 
MIwithY_bias_R <- subset(R_imp, Dataset == "MI with Y") %>%
  rename_at(vars("estimates"), function(x) paste0("true_", x)) %>%
  rename_at(vars("LCI_values"), function(x) paste0("imp_", x)) %>%
  rename_at(vars("UCI_values"), function(x) paste0("imp_",x)) %>%
  left_join(R_val, MIwithY_bias_R,
            multiple = "all",
            by = "target_measures") %>%
  mutate(bias = estimates - true_estimates,
         LCI_diff = LCI_values - imp_LCI_values,
         UCI_diff = UCI_values - imp_UCI_values) %>%
  rename("Dataset_imp" = "Dataset.x") %>%
  rename("Dataset_val" = "Dataset.y")


#####
# Mean + RFA at implementation
## MI no Y at implementation 
Mean_RFA_bias_R <- subset(R_imp, Dataset == "Mean + RFA") %>%
  rename_at(vars("estimates"), function(x) paste0("true_", x)) %>%
  rename_at(vars("LCI_values"), function(x) paste0("imp_", x)) %>%
  rename_at(vars("UCI_values"), function(x) paste0("imp_",x)) %>%
  left_join(R_val, Mean_RFA_bias_R,
            multiple = "all",
            by = "target_measures") %>%
  mutate(bias = estimates - true_estimates,
         LCI_diff = LCI_values - imp_LCI_values,
         UCI_diff = UCI_values - imp_UCI_values) %>%
  rename("Dataset_imp" = "Dataset.x") %>%
  rename("Dataset_val" = "Dataset.y")



all_bias_R <- MInoY_bias_R %>% 
  bind_rows(MIwithY_bias_R) %>%
  bind_rows(Mean_RFA_bias_R)



all_bias_R$target_measures <- factor(all_bias_R$target_measures, levels = c("AUC", "Calibration Intercept", "Calibration Slope", "Brier Score"))

#####
##### plot "bias" 
plot_R <- ggplot(data = all_bias_R, aes(x = bias, y = Dataset_val, color = factor(target_measures),
                                        shape = factor(target_measures))) +
  geom_errorbar(aes(xmin = LCI_diff, xmax = UCI_diff), width = .1) + ### hashtag needs removing once the UCI and LCI are calculated
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
  ggh4x::facet_grid2(target_measures ~ Dataset_imp, scales = "free_x", independent = "x") +
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


##############################################################################################################################

T_datasets <- imputed_datasets[grepl("thoracoscore", names(imputed_datasets))]


for (dataset_name in names(T_datasets)) {
  datasetT <- T_datasets[[dataset_name]]
  
  if (is.mids(datasetT)) {
    datasetT <- mice::complete(datasetT, action = "long")
    
    
    datasetT$LP <- -7.3737 +
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
    
    LP_mean <- datasetT %>%
      dplyr::group_by(.id) %>%
      dplyr::summarise(LP_mean = mean(LP)) %>%
      pull(LP_mean)
    
    Pi <- exp(LP_mean) / (1 + exp(LP_mean))
    
    
    datasetT <- data.frame("LP" = LP_mean,
                           "Pi" = Pi,
                           "DeadatDischarge" = datasetT$DeadatDischarge[1:3300])
    
  } else {
    
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
    
    datasetT$LP <- LP
    datasetT$Pi <- Pi
    
  }
  
  T_datasets[[dataset_name]] <- datasetT
  
}


################################
# Create an empty df to store the results
target_measuresT <- data.frame()


for (dataset_name in names(T_datasets)) {
datasetT <- T_datasets[[dataset_name]]

# Specify the outcome variable
outcome_var <- as.numeric(datasetT$DeadatDischarge)

# Specify the predicted probabilities
Pi <- datasetT$Pi

# Calculate Brier Score
Brier_individuals <- (Pi - outcome_var)^2
Brier <- mean(Brier_individuals)
Brier_var <- var(Brier_individuals)/length(Pi)
Brier_SE <- sqrt(Brier_var) 
Brier_LCI <- Brier - (1.96*Brier_SE)
Brier_UCI <- Brier + (1.96*Brier_SE)

# Calculate Calibration Intercept
LP <- log(Pi/ (1 - Pi))
Cal_int_mod <- glm(outcome_var ~ offset(LP), family = binomial(link = "logit"))
Cal_Int <- as.numeric(coef(Cal_int_mod)[1])
Cal_Int_var <- vcov(Cal_int_mod)[1,1]
Cal_Int_SE <- sqrt(Cal_Int_var) 
Cal_Int_LCI <- Cal_Int - (1.96*Cal_Int_SE)
Cal_Int_UCI <- Cal_Int + (1.96*Cal_Int_SE)


# Calculate Calibration Slope
Cal_Slope_mod <- glm(outcome_var ~ LP, family = binomial(link = "logit"))
Cal_Slope <- as.numeric(coef(Cal_Slope_mod)[2])
Cal_Slope_var <- vcov(Cal_Slope_mod)[2,2]
Cal_Slope_SE<- sqrt(Cal_Slope_var) 
Cal_Slope_LCI <- Cal_Slope - (1.96*Cal_Slope_SE)
Cal_Slope_UCI <- Cal_Slope + (1.96*Cal_Slope_SE)

# Calculate AUC

AUC <- roc(response = outcome_var, 
           predictor = as.vector(Pi),
           direction = "<",
           levels = c(0,1))$auc
AUC_var <- var(AUC, method = "delong") #approximation method used for AUC to calculate the variance
AUC_SE <- sqrt(AUC_var) 
AUC_LCI <- AUC - (1.96*AUC_SE)
AUC_UCI <- AUC + (1.96*AUC_SE)


# Create a data frame with target measures for the current dataset
measures <- data.frame("Dataset" = dataset_name,
                       "Cal_Int" = Cal_Int,
                       "Cal_Int_var" = Cal_Int_var,
                       "Cal_Int_LCI" = as.numeric(Cal_Int_LCI),
                       "Cal_Int_UCI" = as.numeric(Cal_Int_UCI),
                       
                       "Cal_Slope" = Cal_Slope,
                       "Cal_Slope_var" = as.numeric(Cal_Slope_var),
                       "Cal_Slope_LCI" = as.numeric(Cal_Slope_LCI),
                       "Cal_Slope_UCI" = as.numeric(Cal_Slope_UCI),
                       
                       "AUC" = as.numeric(AUC),
                       "AUC_var" = as.numeric(AUC_var),
                       "AUC_LCI" = as.numeric(AUC_LCI),
                       "AUC_UCI" = as.numeric(AUC_UCI),
                       
                       
                       "Brier" = as.numeric(Brier),
                       "Brier_var" = as.numeric(Brier_var),
                       "Brier_LCI" = as.numeric(Brier_LCI),
                       "Brier_UCI" = as.numeric(Brier_UCI)
)

# Append the measures to the overall results data frame
target_measuresT <- bind_rows(target_measuresT, measures)

}


################################################################################

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
      ggplot2::xlim(c(0,1)) + ylim(c(0,1)) + 
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

# List of dataset names and their new names for Thoracoscore datasets
thoracoscore_datasets <- list(MI_noY_val_thoracoscore = T_datasets$MI_noY_val_thoracoscore, 
                              MI_withY_val_thoracoscore = T_datasets$MI_withY_val_thoracoscore,
                              MI_noY_imp_thoracoscore = T_datasets$MI_noY_imp_thoracoscore,
                              MI_withY_imp_thoracoscore = T_datasets$MI_withY_imp_thoracoscore,
                              CCA_val_thoracoscore = T_datasets$CCA_val_thoracoscore,
                              mean_zero_val_thoracoscore = T_datasets$mean_zero_val_thoracoscore,
                              mean_zero_imp_thoracoscore = T_datasets$mean_zero_imp_thoracoscore)

new_namesT <- c("Thoracoscore: MI no Y at validation", 
                "Thoracoscore: MI with Y at validation",
                "Thoracoscore: MI no Y at implementation", 
                "Thoracoscore: MI with Y at implementation", 
                "Thoracoscore: CCA at validation",
                "Thoracoscore: Mean + Risk Factor Absent at validation",
                "Thoracoscore: Mean + Risk Factor Absent at implementation")

# Example usage for Resect datasets:
create_calibration_plots(dataset_list = thoracoscore_datasets, 
                         dataset_names = new_namesT,
                         outcome_var = "DeadatDischarge")



################################################################################
#turn data into long format and rename TM


# Pivot the first set of columns
target_measuresT_long <- target_measuresT %>%
  pivot_longer(cols = c("Cal_Int", "Cal_Slope", "AUC", "Brier"), 
               names_to = "target_measures", 
               values_to = "estimates") %>%
  select(Dataset, target_measures, estimates) %>%
  mutate(target_measures = recode(target_measures,
                                  Cal_Int = 'Calibration Intercept',
                                  Cal_Slope = 'Calibration Slope',
                                  AUC = 'AUC',
                                  Brier = 'Brier Score')) %>%
  left_join(target_measuresT %>%
              pivot_longer(cols = c("Cal_Int_LCI", "Cal_Slope_LCI", "AUC_LCI", "Brier_LCI"), 
                           names_to = "target_measures", 
                           values_to = "LCI_values") %>%
              select(Dataset, target_measures, LCI_values) %>%
              mutate(target_measures = recode(target_measures,
                                              Cal_Int_LCI = 'Calibration Intercept',
                                              Cal_Slope_LCI = 'Calibration Slope',
                                              AUC_LCI = 'AUC',
                                              Brier_LCI = 'Brier Score')),
            by = c("Dataset", "target_measures")) %>%
  left_join(target_measuresT %>%
              pivot_longer(cols = c("Cal_Int_UCI", "Cal_Slope_UCI", "AUC_UCI", "Brier_UCI"), 
                           names_to = "target_measures", 
                           values_to = "UCI_values") %>%
              select(Dataset, target_measures, UCI_values) %>%
              mutate(target_measures = recode(target_measures,
                                              Cal_Int_UCI = 'Calibration Intercept',
                                              Cal_Slope_UCI = 'Calibration Slope',
                                              AUC_UCI = 'AUC',
                                              Brier_UCI = 'Brier Score')),
            by = c("Dataset", "target_measures"))



#split the data into val and imp 
T_val <- target_measuresT_long[target_measuresT_long$Dataset %like% "val", ]
T_imp <- target_measuresT_long[target_measuresT_long$Dataset %like% "imp", ]

T_val <- T_val %>%
  mutate(Dataset = case_when(
    Dataset == "MI_noY_val_thoracoscore" ~ "MI no Y",
    Dataset == "MI_withY_val_thoracoscore" ~ "MI with Y",
    Dataset == "CCA_val_thoracoscore" ~ "CCA",
    Dataset == "mean_zero_val_thoracoscore" ~ "Mean + RFA",
    TRUE ~ Dataset
  ))


T_imp <- T_imp %>%
  mutate(Dataset = case_when(
    Dataset == "MI_noY_imp_thoracoscore" ~ "MI no Y",
    Dataset == "MI_withY_imp_thoracoscore" ~ "MI with Y",
    Dataset == "CCA_imp_thoracoscore" ~ "CCA",
    Dataset == "mean_zero_imp_thoracoscore" ~ "Mean + RFA",
    TRUE ~ Dataset
  ))


##############
## MI no Y at implementation 
MInoY_bias_T <- subset(T_imp, Dataset == "MI no Y") %>%
  rename_at(vars("estimates"), function(x) paste0("true_", x)) %>%
  rename_at(vars("LCI_values"), function(x) paste0("imp_", x)) %>%
  rename_at(vars("UCI_values"), function(x) paste0("imp_",x)) %>%
  left_join(T_val, MInoY_bias_T,
            multiple = "all",
            by = "target_measures") %>%
  mutate(bias = estimates - true_estimates,
         LCI_diff = LCI_values - imp_LCI_values,
         UCI_diff = UCI_values - imp_UCI_values) %>%
  rename("Dataset_imp" = "Dataset.x") %>%
  rename("Dataset_val" = "Dataset.y")

#####


## MI no Y at implementation 
MIwithY_bias_T <- subset(T_imp, Dataset == "MI with Y") %>%
  rename_at(vars("estimates"), function(x) paste0("true_", x)) %>%
  rename_at(vars("LCI_values"), function(x) paste0("imp_", x)) %>%
  rename_at(vars("UCI_values"), function(x) paste0("imp_",x)) %>%
  left_join(T_val, MIwithY_bias_T,
            multiple = "all",
            by = "target_measures") %>%
  mutate(bias = estimates - true_estimates,
         LCI_diff = LCI_values - imp_LCI_values,
         UCI_diff = UCI_values - imp_UCI_values) %>%
  rename("Dataset_imp" = "Dataset.x") %>%
  rename("Dataset_val" = "Dataset.y")


#####
# Mean + RFA at implementation
## MI no Y at implementation 
Mean_RFA_bias_T <- subset(T_imp, Dataset == "Mean + RFA") %>%
  rename_at(vars("estimates"), function(x) paste0("true_", x)) %>%
  rename_at(vars("LCI_values"), function(x) paste0("imp_", x)) %>%
  rename_at(vars("UCI_values"), function(x) paste0("imp_",x)) %>%
  left_join(T_val, Mean_RFA_bias_T,
            multiple = "all",
            by = "target_measures") %>%
  mutate(bias = estimates - true_estimates,
         LCI_diff = LCI_values - imp_LCI_values,
         UCI_diff = UCI_values - imp_UCI_values) %>%
  rename("Dataset_imp" = "Dataset.x") %>%
  rename("Dataset_val" = "Dataset.y")



all_bias_T <- MInoY_bias_T %>% 
  bind_rows(MIwithY_bias_T) %>%
  bind_rows(Mean_RFA_bias_T)



all_bias_T$target_measures <- factor(all_bias_T$target_measures, levels = c("AUC", "Calibration Intercept", "Calibration Slope", "Brier Score"))

####
##### plot "bias" 
plot_T <- ggplot(data = all_bias_T, aes(x = bias, y = Dataset_val, color = factor(target_measures),
                                        shape = factor(target_measures))) +
  geom_errorbar(aes(xmin = LCI_diff, xmax = UCI_diff), width = .1) + ### hashtag needs removing once the UCI and LCI are calculated
  geom_point(size = 3, stroke = 0.5) +
  guides(color = guide_legend(reverse = TRUE)) +
  scale_shape_manual(values = c(8, 17, 16, 15)) +
  scale_color_brewer(palette = "Set1") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  xlab("Difference in the performance metrics under different combinations of imputation methods") +
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
  ggh4x::facet_grid2(target_measures ~ Dataset_imp, scales = "free_x", independent = "x") +
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












































