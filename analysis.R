## ----setup, include=FALSE--------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, dpi=600)


## ----libraries, message=FALSE, warning=FALSE-------------------------------------------------------------------------------------------
# load libraries
library(readr)
library(ggplot2)
library(MVN)
library(GGally)
library(car)
library(biotools)
library(performance)
library(DescTools)
library(knitr)
library(kableExtra)
library(effectsize)
library(emmeans)
library(psych)
library(tidyr)
library(VIM)
library(mice)
library(naniar)
library(purrr)
library(Rmisc)
library(dplyr)
library(mitools)
library(heplots)


## ----assign-classrooms-----------------------------------------------------------------------------------------------------------------
# Load the classroom data into a data frame
classroom_data <- read.csv("classroom_data.csv", stringsAsFactors = FALSE)

# Split the data frame by grade level
grade6 <- classroom_data[classroom_data$grade == 6, ]
grade8 <- classroom_data[classroom_data$grade == 8, ]

# For grade 6
n_grade6 <- nrow(grade6)
# Randomly decide if group A or B gets the extra classroom when odd
if (n_grade6 %% 2 == 1) {
  extra_to_a <- sample(c(TRUE, FALSE), 1)
  n_group_a_grade6 <- floor(n_grade6 / 2) + ifelse(extra_to_a, 1, 0)
  n_group_b_grade6 <- floor(n_grade6 / 2) + ifelse(extra_to_a, 0, 1)
} else {
  n_group_a_grade6 <- n_grade6 / 2
  n_group_b_grade6 <- n_grade6 / 2
}
grade6$group <- c(rep("A", n_group_a_grade6), rep("B", n_group_b_grade6))
grade6$group <- sample(grade6$group)  # Randomize the assignments

# For grade 8
n_grade8 <- nrow(grade8)
# Randomly decide if group A or B gets the extra classroom when odd
if (n_grade8 %% 2 == 1) {
  extra_to_a <- sample(c(TRUE, FALSE), 1)
  n_group_a_grade8 <- floor(n_grade8 / 2) + ifelse(extra_to_a, 1, 0)
  n_group_b_grade8 <- floor(n_grade8 / 2) + ifelse(extra_to_a, 0, 1)
} else {
  n_group_a_grade8 <- n_grade8 / 2
  n_group_b_grade8 <- n_grade8 / 2
}
grade8$group <- c(rep("A", n_group_a_grade8), rep("B", n_group_b_grade8))
grade8$group <- sample(grade8$group)  # Randomize the assignments

# Recombine the data
classroom_data <- rbind(grade6, grade8)

# Export the final data frame to a new CSV file
write.csv(classroom_data, "classroom_data_with_groups.csv", row.names = FALSE)


## ----load-data-------------------------------------------------------------------------------------------------------------------------
# load data frame
raw_pbl_study_data <- read_csv('scd-s_responses.csv',
                           na = c("NA", " ", ""),
                           col_types = cols(
                               group = col_character(),
                               grade_reported = col_factor(),
                               gender = col_factor(),
                               race_reported = col_factor(),
                               race = col_factor(),
                               .default = col_double()
                           ))

# Convert group to a factor
raw_pbl_study_data$group <- as.factor(raw_pbl_study_data$group)

# Convert participant_id to text instead of number
raw_pbl_study_data$participant_id <- as.character(raw_pbl_study_data$participant_id)


## ----group-codes-----------------------------------------------------------------------------------------------------------------------
raw_pbl_study_data$group <- factor(raw_pbl_study_data$group,
                            levels = c("A", "B", "C"),
                            labels = c("CR-PBL", "S-PBL", "Control"))


## ----gender-terms----------------------------------------------------------------------------------------------------------------------
raw_pbl_study_data$gender <- factor(raw_pbl_study_data$gender,
                            levels = c("Boy", "Girl", "Other"),
                            labels = c("Male", "Female", "Other"))


## ----grade-levels----------------------------------------------------------------------------------------------------------------------
raw_pbl_study_data$grade <- with(raw_pbl_study_data, case_when(
        grade_reported == "Grade 6" ~ "6",
        grade_reported == "Grade 7" ~ "6",
        grade_reported == "Grade 8" ~ "8"
    )) |> as.factor()


## ----set-factor-orders-----------------------------------------------------------------------------------------------------------------
# Modify factor levels to establish reference level and order
raw_pbl_study_data <- raw_pbl_study_data |> mutate(
        group = factor(group, levels = c("Control", "S-PBL", "CR-PBL")),
        grade = factor(grade, levels = c("6", "8")),
        gender = factor(gender, levels = sort(as.character(unique(gender)))),
        race = factor(race, levels = sort(as.character(unique(race))))
    )


## ----compute-subscales-----------------------------------------------------------------------------------------------------------------
# Compute SCD-S Subscales from individual items. Since each survey item label begins with the subscale code, R can search for and add all items that start with the same code.

raw_pbl_study_data <- raw_pbl_study_data |>
  mutate(
    QI = rowSums(across(starts_with("QI_"))),
    PCC = rowSums(across(starts_with("PCC_"))),
    CuSo = rowSums(across(starts_with("CuSo_"))),
    CrCs = rowSums(across(starts_with("CrCs_")))
  )


## ----rename-variables------------------------------------------------------------------------------------------------------------------
# Rename the variables for convenience of reference
school_climate_data <- raw_pbl_study_data |> rename(
    covariate_ela = PSSA,
    dv1_qual_interaction = QI,
    dv2_prom_cult_competence = PCC,
    dv3_cult_soc = CuSo,
    dv4_crit_con_soc = CrCs
)


## ----remove-unused-columns-------------------------------------------------------------------------------------------------------------
analysis_vars <- c("participant_id", "group", "grade", "gender", "race", 
                   "covariate_ela", "dv1_qual_interaction", "dv2_prom_cult_competence", 
                   "dv3_cult_soc", "dv4_crit_con_soc")

school_climate_data <- school_climate_data[analysis_vars]


## ----missing-patterns------------------------------------------------------------------------------------------------------------------
missing_pattern <- md.pattern(school_climate_data[, c("covariate_ela", 
                                                      "dv1_qual_interaction",
                                                      "dv2_prom_cult_competence", 
                                                      "dv3_cult_soc", 
                                                      "dv4_crit_con_soc")])


## ----missingness-by-group--------------------------------------------------------------------------------------------------------------
missing_by_group <- school_climate_data |>
  group_by(group) |>
  summarise(
    n = n(),
    missing_ela = sum(is.na(covariate_ela)),
    pct_missing_ela = (missing_ela/n)*100,
    missing_dv1 = sum(is.na(dv1_qual_interaction)),
    pct_missing_dv1 = (missing_dv1/n)*100,
    .groups = 'drop'
  )

print("Missing values by treatment group:")
print(missing_by_group)


## ----mcar-test-------------------------------------------------------------------------------------------------------------------------
mcar_test_result <- mcar_test(school_climate_data[, c("group", "covariate_ela", 
                                                      "dv1_qual_interaction",
                                                      "dv2_prom_cult_competence", 
                                                      "dv3_cult_soc", 
                                                      "dv4_crit_con_soc")])
print(mcar_test_result)


## ----summary-raw-----------------------------------------------------------------------------------------------------------------------
print("Data summary before imputation:")
print(summary(school_climate_data[, c("covariate_ela", 
                                      "dv1_qual_interaction",
                                      "dv2_prom_cult_competence", 
                                      "dv3_cult_soc", 
                                      "dv4_crit_con_soc")]))


## ----imputation------------------------------------------------------------------------------------------------------------------------
set.seed(123)  # For reproducibility
imputed_data <- mice(school_climate_data, 
                     m = 5,  # Number of imputations; 5 is commonly recommended when a small percentage of data is missing
                     maxit = 50,  # Number of iterations
                     printFlag = FALSE)

# display any logged events in case there are issues
imputed_data$loggedEvents


## ----check-convergence-----------------------------------------------------------------------------------------------------------------
plot(imputed_data, c("covariate_ela"))


## ----density-plot----------------------------------------------------------------------------------------------------------------------
# Density plot of the entire dataset
densityplot(imputed_data, ~ covariate_ela)

# Density plot of just the Control group, since that is where all missing data occurred
densityplot(imputed_data, ~ covariate_ela, subset = group == "Control")



## ----complete-imputations--------------------------------------------------------------------------------------------------------------
imputed_climate_data <- complete(imputed_data, action = "all")


## ----linearity-crpbl, fig.width=10, fig.height=8---------------------------------------------------------------------------------------
# Use the first imputed dataset for linearity checks
first_imputed <- imputed_climate_data[[1]]

# Data with grade, gender, and race removed; includes group, covariate and all dependent variables
numeric_data <- first_imputed[, c(
    "group",
    "covariate_ela",
    "dv1_qual_interaction",
    "dv2_prom_cult_competence",
    "dv3_cult_soc",
    "dv4_crit_con_soc"
)]

# Create subgroups with only the data for one group
crpbl_data <- numeric_data |> filter(group == "CR-PBL") |> select(-group)  # removes group column
spbl_data <- numeric_data |> filter(group == "S-PBL") |> select(-group)
control_data <- numeric_data |> filter(group == "Control") |> select(-group)

# Create pairs plot with smoothed lines in lower triangle
pairs(crpbl_data,
     upper.panel = NULL,  # Blank upper panel
     diag.panel = NULL,   # Blank diagonal
     lower.panel = function(x, y, ...) {
       points(x, y, col = rgb(0,0,0,0.5), pch = 16, cex = 0.6)
       lines(lowess(x, y, f = 0.8), col = "black", lwd = 2)
     })


## ----linearity-spbl, fig.width=10, fig.height=8----------------------------------------------------------------------------------------
pairs(spbl_data,
     upper.panel = NULL,  # Blank upper panel
     diag.panel = NULL,   # Blank diagonal
     lower.panel = function(x, y, ...) {
       points(x, y, col = rgb(0,0,0,0.5), pch = 16, cex = 0.6)
       lines(lowess(x, y, f = 0.8), col = "black", lwd = 2)
     })


## ----linearity-control, fig.width=10, fig.height=8-------------------------------------------------------------------------------------
pairs(control_data,
     upper.panel = NULL,  # Blank upper panel
     diag.panel = NULL,   # Blank diagonal
     lower.panel = function(x, y, ...) {
       points(x, y, col = rgb(0,0,0,0.5), pch = 16, cex = 0.6)
       lines(lowess(x, y, f = 0.8), col = "black", lwd = 2)
     })


## ----covariate-transformation----------------------------------------------------------------------------------------------------------
# Add square of covariate to all 5 imputed datasets only
# covariate is centered to avoid multicollinearity issues in later assumption tests

# Calculate a single centering constant using pooled mean across all imputations
all_values <- unlist(lapply(imputed_climate_data, function(df) df$covariate_ela))
global_centering_constant <- mean(all_values, na.rm = TRUE)

cat("Using global centering constant (pooled mean):", global_centering_constant, "\n")

# Apply the SAME centering constant to all imputations
imputed_climate_data <- map(imputed_climate_data, function(df) {
  df |> mutate(
    covariate_ela_centered = covariate_ela - global_centering_constant,
    covariate_ela_squared = covariate_ela_centered^2
  )
})

# Verify that centering worked consistently
cat("Verification - Centered covariate means across imputations:\n")
centered_means <- sapply(imputed_climate_data, function(df) mean(df$covariate_ela_centered, na.rm = TRUE))
print(round(centered_means, 6))  # Should all be close to zero
print(round(mean(centered_means), 6))  # Should be zero


## ----linearity-transformed-------------------------------------------------------------------------------------------------------------
# Check polynomial covariate adequacy - model fit results only
check_polynomial_fit <- function(data) {
  
  dvs <- c("dv1_qual_interaction", "dv2_prom_cult_competence", 
           "dv3_cult_soc", "dv4_crit_con_soc")
  
  results <- data.frame(
    DV = character(),
    Linear_R2 = numeric(),
    Quadratic_R2 = numeric(),
    Improvement = numeric(),
    stringsAsFactors = FALSE
  )
  
for(dv in dvs) {
    # Create clean data for modeling inline
    complete_cases <- complete.cases(data[[dv]], data$covariate_ela_centered)
    clean_y <- data[[dv]][complete_cases]
    clean_x <- data$covariate_ela_centered[complete_cases]
    
    # Fit models directly
    linear_fit <- lm(clean_y ~ clean_x)
    quad_fit <- lm(clean_y ~ clean_x + I(clean_x^2))
  
    # Extract R-squared values
    linear_r2 <- summary(linear_fit)$r.squared
    quad_r2 <- summary(quad_fit)$r.squared
    improvement <- quad_r2 - linear_r2

    # Store results
    results <- rbind(results, data.frame(
      DV = dv,
      Linear_R2 = round(linear_r2, 3),
      Quadratic_R2 = round(quad_r2, 3),
      Improvement = round(improvement, 3)
    ))
  }
  
  return(results)
}

# Analyze first imputed dataset
first_imputed <- imputed_climate_data[[1]]
fit_results <- check_polynomial_fit(first_imputed)

# Display results
print(fit_results)


## ----hors-results----------------------------------------------------------------------------------------------------------------------
hors_results <- map(imputed_climate_data, ~{
    hors_model <- lm(cbind(dv1_qual_interaction, dv2_prom_cult_competence, dv3_cult_soc, dv4_crit_con_soc) ~
                           group * (covariate_ela_centered + covariate_ela_squared),
                           data = .x)
    hors_result <- Manova(hors_model, type = "III")
    summary(hors_result, multivariate = TRUE)
})

# Extract both interaction sections from the full summary
iwalk(hors_results, ~{
  cat("Imputation", .y, ":\n")
  
  # Capture full output and extract just the interaction section
  full_output <- capture.output(print(.x))
  
  # Find both interaction section
  linear_interaction <- grep("Term: group:covariate_ela_centered", full_output)
  squared_interaction <- grep("Term: group:covariate_ela_squared", full_output)
  
  # Function to extract a section
  extract_section <- function(start_line) {
    if(length(start_line) > 0) {
      # Find the next "Term:" or end of output
      next_term <- grep("^Term:", full_output)
      next_term <- next_term[next_term > start_line][1]
      
      if(is.na(next_term)) {
        end_line <- length(full_output)
      } else {
        end_line <- next_term - 1
      }
      
      return(full_output[start_line:end_line])
    }
    return(NULL)
  }
  
  # Print linear covariate interaction
  if(length(linear_interaction) > 0) {
    cat("LINEAR COVARIATE INTERACTION:\n")
    linear_section <- extract_section(linear_interaction)
    cat(linear_section, sep = "\n")
    cat("\n")
  }
  
  # Print squared covariate interaction  
  if(length(squared_interaction) > 0) {
    cat("QUADRATIC COVARIATE INTERACTION:\n")
    squared_section <- extract_section(squared_interaction)
    cat(squared_section, sep = "\n")
    cat("\n")
  }
  
  # If no interactions found, show what terms are available
  if(length(linear_interaction) == 0 && length(squared_interaction) == 0) {
    term_lines <- grep("^Term:", full_output, value = TRUE)
    cat("Available terms:\n")
    cat(term_lines, sep = "\n")
    cat("\n")
  }
  
  cat(strrep("=", 60), "\n\n")
})


## ----boxs-m----------------------------------------------------------------------------------------------------------------------------
box_test_result <- biotools::boxM(data = imputed_climate_data[[1]][, c("dv1_qual_interaction",
                                    "dv2_prom_cult_competence",
                                    "dv3_cult_soc",
                                    "dv4_crit_con_soc")],
                     group = imputed_climate_data[[1]]$group)

print(box_test_result)


## ----log-determinants------------------------------------------------------------------------------------------------------------------
box_test_result$logDet


## ----correlation-matrix-2--------------------------------------------------------------------------------------------------------------
print(round(cov2cor(box_test_result$pooled), 3))


## ----winsorize-function----------------------------------------------------------------------------------------------------------------
# Function to Winsorize covariate in one dataset
winsorize_covariate <- function(data) {
    # First, rename the existing covariate column to indicate it is the original raw data
    data <- data |> rename(covariate_ela_raw = covariate_ela)
    
    # Compute probabilities for ±3 SD
    mean_ela <- mean(data$covariate_ela_raw)
    sd_ela <- sd(data$covariate_ela_raw)
    lower_bound <- mean_ela - 3*sd_ela
    upper_bound <- mean_ela + 3*sd_ela
    
    # Winsorize the ELA scores
    data$covariate_ela <- Winsorize(data$covariate_ela_raw, val = c(lower_bound, upper_bound))
    
    # Create summary of winsorization
    changed_ela_scores <- which(data$covariate_ela != data$covariate_ela_raw)
    
    winsor_summary <- list(
      changes = changed_ela_scores,
      original_values = data$covariate_ela_raw[changed_ela_scores],
      new_values = data$covariate_ela[changed_ela_scores],
      bounds = c(lower_bound, upper_bound)
    )
    
    # Return both the modified data and the summary
    return(list(data = data, summary = winsor_summary))
}


## ----winsorize-covariate---------------------------------------------------------------------------------------------------------------
# Winsorize each dataset
winsorized_results <- map(imputed_climate_data, winsorize_covariate)

# Extract the datasets and summaries
imputed_climate_data <- map(winsorized_results, ~.x$data)
winsorization_summaries <- map(winsorized_results, ~.x$summary)

# View summaries to see what was changed in each dataset
winsorization_summaries


## ----display-winsorization-------------------------------------------------------------------------------------------------------------
# Display winsorization tables for each dataset
iwalk(winsorization_summaries, ~{
    cat("Dataset", .y, ": ")
    
    if(length(.x$changes) > 0) {
        cat("")
        print(
            data.frame(
            Original = .x$original_values,
            Winsorized = .x$new_values
          ) |> 
            knitr::kable(
              caption = sprintf("Dataset $d: Values Winsorized (n = %d)", 
                           .y, length(.x$changes)),
              digits = 2
            )
        )
        cat("\n")
    } else {
      cat("No values required winsorization\n")
    }
})


## ----standardized-residuals------------------------------------------------------------------------------------------------------------
# Compute standardized residuals
first_dataset_with_residuals <- within(imputed_climate_data[[1]], {
    dv1_resid <- rstandard(lm(dv1_qual_interaction ~ group))
    dv2_resid <- rstandard(lm(dv2_prom_cult_competence ~ group))
    dv3_resid <- rstandard(lm(dv3_cult_soc ~ group))
    dv4_resid <- rstandard(lm(dv4_crit_con_soc ~ group))
})


## ----collect-outliers------------------------------------------------------------------------------------------------------------------
# Gather all outliers in a separate list for later analysis if necessary
school_climate_outliers <- subset(first_dataset_with_residuals,
                                  abs(dv1_resid) > 3 |
                                  abs(dv2_resid) > 3 |
                                  abs(dv3_resid) > 3 |
                                  abs(dv4_resid) > 3)


## ----remove-outliers-------------------------------------------------------------------------------------------------------------------
# save the original imputed datasets with outliers
imputed_climate_with_outliers <- imputed_climate_data

# Get row IDs of outlier cases
outlier_row_ids <- as.numeric(rownames(school_climate_outliers))

# Now remove outlier cases from working data frame in all imputed datasets
imputed_climate_data <- map(imputed_climate_data, ~{
    .x[-outlier_row_ids, ]  # removes outlier rows
})


## ----qq-before-transform, fig.width=10, fig.height=8-----------------------------------------------------------------------------------
# Extract DV columns
dv_cols <- c("dv1_qual_interaction", "dv2_prom_cult_competence", "dv3_cult_soc", "dv4_crit_con_soc")

# Check normality across all datasets
normality_results_list <- map(imputed_climate_data,
                              ~mvn(as.matrix(.x[dv_cols]),
                                   mvn_test = "hz",
                                   univariate_test = "SW")
                              )

normality_results_with_outliers <- map(imputed_climate_with_outliers,
                              ~mvn(as.matrix(.x[dv_cols]),
                                   mvn_test = "hz",
                                   univariate_test = "SW")
                              )
# Get the DV matrix and create the QQ plot
dv_matrix_first <- imputed_climate_data[[1]][dv_cols]
multivariate_diagnostic_plot(dv_matrix_first, type = "qq")


## ----normality-before-transform--------------------------------------------------------------------------------------------------------
normality_result <- mvn(as.matrix(imputed_climate_data[[1]][dv_cols]),
                       mvn_test = "hz",
                       univariate_test = "SW")

# Table format
kable(normality_result$multivariate_normality, 
      caption = "Multivariate Normality Test (After DV Transformations)",
      digits = 4)


## ----skewness-check--------------------------------------------------------------------------------------------------------------------
# Check skewness for each DV
dv_cols <- c("dv1_qual_interaction", "dv2_prom_cult_competence", 
             "dv3_cult_soc", "dv4_crit_con_soc")

# Use just the first dataset
data <- imputed_climate_data[[1]]

skewness_results <- data.frame(
  DV = dv_cols,
  Skewness = map_dbl(dv_cols, ~skew(data[[.x]], na.rm = TRUE)),
  Kurtosis = map_dbl(dv_cols, ~kurtosi(data[[.x]], na.rm = TRUE))
)

print(skewness_results)

# Add interpretation
cat("\nInterpretation:\n")
skewness_results |>
  mutate(
    Skew_Severity = case_when(
      abs(Skewness) < 0.5 ~ "Mild",
      abs(Skewness) < 1.0 ~ "Moderate", 
      TRUE ~ "Severe"
    ),
    Skew_Direction = ifelse(Skewness > 0, "Right-skewed", "Left-skewed")
  ) |>
  select(DV, Skewness, Skew_Direction, Skew_Severity, Kurtosis) |>
  print()


## ----transform-dvs---------------------------------------------------------------------------------------------------------------------
# Calculate the reflect-sqrt transformation ranges for consistency
dv1_temp_ref = sqrt(16 - imputed_climate_data[[1]]$dv1_qual_interaction)
dv2_temp_ref = sqrt(31 - imputed_climate_data[[1]]$dv2_prom_cult_competence)
dv3_temp_ref = sqrt(16 - imputed_climate_data[[1]]$dv3_cult_soc)
dv4_temp_ref = sqrt(21 - imputed_climate_data[[1]]$dv4_crit_con_soc)

# Calculate re-reflection constants (same for all datasets)
dv1_reflect_constant <- max(dv1_temp_ref, na.rm = TRUE) + min(dv1_temp_ref, na.rm = TRUE)
dv2_reflect_constant <- max(dv2_temp_ref, na.rm = TRUE) + min(dv2_temp_ref, na.rm = TRUE)
dv3_reflect_constant <- max(dv3_temp_ref, na.rm = TRUE) + min(dv3_temp_ref, na.rm = TRUE)
dv4_reflect_constant <- max(dv4_temp_ref, na.rm = TRUE) + min(dv4_temp_ref, na.rm = TRUE)

# Apply transformation to ALL datasets using the same constants
imputed_climate_data <- map(imputed_climate_data, function(df) {
  df |> mutate(
    # Save original values for all DVs
    dv1_qual_interaction_raw = dv1_qual_interaction,
    dv2_prom_cult_competence_raw = dv2_prom_cult_competence, 
    dv3_cult_soc_raw = dv3_cult_soc,
    dv4_crit_con_soc_raw = dv4_crit_con_soc,
    
    # Step 1: Apply reflect and square root for each DV
    dv1_temp = sqrt(16 - dv1_qual_interaction_raw),
    dv2_temp = sqrt(31 - dv2_prom_cult_competence_raw),
    dv3_temp = sqrt(16 - dv3_cult_soc_raw),
    dv4_temp = sqrt(21 - dv4_crit_con_soc_raw),
    
    # Step 2: Re-reflect using consistent constants
    dv1_qual_interaction = dv1_reflect_constant - dv1_temp,
    dv2_prom_cult_competence = dv2_reflect_constant - dv2_temp,
    dv3_cult_soc = dv3_reflect_constant - dv3_temp,
    dv4_crit_con_soc = dv4_reflect_constant - dv4_temp
  ) |>
  select(-dv1_temp, -dv2_temp, -dv3_temp, -dv4_temp)  # Remove temporary variables
})

# Show the constants for verification
cat("Re-reflection constants (consistent across all datasets):\n")
cat("DV1:", round(dv1_reflect_constant, 3), "\n")
cat("DV2:", round(dv2_reflect_constant, 3), "\n") 
cat("DV3:", round(dv3_reflect_constant, 3), "\n")
cat("DV4:", round(dv4_reflect_constant, 3), "\n\n")

# Check transformations for all DVs
cat("=== TRANSFORMATION SUMMARY FOR ALL DVs ===\n\n")
data_check <- imputed_climate_data[[1]]

# Define DV information
dv_info <- data.frame(
  dv_name = c("dv1_qual_interaction", "dv2_prom_cult_competence", "dv3_cult_soc", "dv4_crit_con_soc"),
  dv_label = c("Quality Interaction", "Cultural Competence", "Cultural Socialization", "Critical Consciousness"),
  stringsAsFactors = FALSE
)

# Loop through each DV and report results
for(i in 1:nrow(dv_info)) {
  dv_name <- dv_info$dv_name[i]
  dv_label <- dv_info$dv_label[i]
  raw_name <- paste0(dv_name, "_raw")
  
  cat(dv_label, ":\n")
  
  # Original range
  orig_min <- min(data_check[[raw_name]], na.rm = TRUE)
  orig_max <- max(data_check[[raw_name]], na.rm = TRUE)
  cat("  Original range:", orig_min, "to", orig_max, "\n")
  
  # Transformed range  
  trans_min <- round(min(data_check[[dv_name]], na.rm = TRUE), 3)
  trans_max <- round(max(data_check[[dv_name]], na.rm = TRUE), 3)
  cat("  Transformed range:", trans_min, "to", trans_max, "\n")
  
  # Skewness comparison
  original_skew <- skew(data_check[[raw_name]], na.rm = TRUE)
  final_skew <- skew(data_check[[dv_name]], na.rm = TRUE)
  improvement <- abs(original_skew) - abs(final_skew)
  
  cat("  Original skewness:", round(original_skew, 3), "\n")
  cat("  Final skewness:", round(final_skew, 3), "\n")
  cat("  Improvement:", round(improvement, 3), "\n")
  
  # Directionality check
  direction_cor <- cor(data_check[[dv_name]], data_check[[raw_name]], use = "complete.obs")
  cat("  Correlation with original:", round(direction_cor, 3), "\n")
  
  cat("\n")
}

# Summary table of skewness improvements
cat("=== SKEWNESS SUMMARY TABLE ===\n")
skewness_summary <- data.frame(
  DV = dv_info$dv_label,
  Original_Skewness = numeric(4),
  Final_Skewness = numeric(4), 
  Improvement = numeric(4),
  Assessment = character(4),
  stringsAsFactors = FALSE
)

for(i in 1:nrow(dv_info)) {
  dv_name <- dv_info$dv_name[i]
  raw_name <- paste0(dv_name, "_raw")
  
  orig_skew <- skew(data_check[[raw_name]], na.rm = TRUE)
  final_skew <- skew(data_check[[dv_name]], na.rm = TRUE)
  improvement <- abs(orig_skew) - abs(final_skew)
  
  skewness_summary$Original_Skewness[i] <- round(orig_skew, 3)
  skewness_summary$Final_Skewness[i] <- round(final_skew, 3)
  skewness_summary$Improvement[i] <- round(improvement, 3)
  
  # Assessment
  if(abs(final_skew) < 0.5) {
    skewness_summary$Assessment[i] <- "Good"
  } else if(abs(final_skew) < 1.0) {
    skewness_summary$Assessment[i] <- "Acceptable" 
  } else {
    skewness_summary$Assessment[i] <- "Still concerning"
  }
}

print(skewness_summary)


## ----qq-after-transform, fig.width=10, fig.height=8------------------------------------------------------------------------------------
# Extract DV columns
dv_cols <- c("dv1_qual_interaction", "dv2_prom_cult_competence", "dv3_cult_soc", "dv4_crit_con_soc")

# Check normality across all datasets
normality_results_list <- map(imputed_climate_data,
                              ~mvn(as.matrix(.x[dv_cols]),
                                   mvn_test = "hz",
                                   univariate_test = "SW")
                              )

# Get the DV matrix and create the QQ plot
multivariate_diagnostic_plot(imputed_climate_data[[1]][dv_cols], type = "qq")


## ----normality-after-transform---------------------------------------------------------------------------------------------------------
normality_result <- mvn(as.matrix(imputed_climate_data[[1]][dv_cols]),
                       mvn_test = "hz",
                       univariate_test = "SW")

# Table format
kable(normality_result$multivariate_normality, 
      caption = "Multivariate Normality Test (After DV Transformations)",
      digits = 4)


## ----shapiro-wilk----------------------------------------------------------------------------------------------------------------------
# Display univariate normality results - single dataset only
mvn_residuals <- normality_result$univariate_normality

kable(mvn_residuals,
      caption = "Univariate Normality Tests (After DV Transformations)",
      digits = 3,
      format = "pipe",
      align = 'l')


## ----contrast-coding-------------------------------------------------------------------------------------------------------------------
for (i in 1:5) {
    contrasts(imputed_climate_data[[i]]$group) <- cbind(
        "PBLvsControl" = c(-2, 1, 1),
        "CRvsSPBL" = c(0, -1, 1)
    )
}


## ----descriptives-raw------------------------------------------------------------------------------------------------------------------
# Define the dependent variables and covariate
dvs <- c("dv1_qual_interaction", "dv2_prom_cult_competence", 
         "dv3_cult_soc", "dv4_crit_con_soc")
covariate <- "covariate_ela"
vars_to_summarize <- c(covariate, dvs)

# Compute descriptive statistics for all dependent variables and covariate
overall_desc <- describe(school_climate_data[, vars_to_summarize])

# Add confidence intervals to overall descriptives
overall_desc$df <- overall_desc$n - 1
overall_desc$t_critical <- qt(0.975, overall_desc$df)  # 95% CI, two-tailed
overall_desc$ci_lower <- overall_desc$mean - (overall_desc$se * overall_desc$t_critical)
overall_desc$ci_upper <- overall_desc$mean + (overall_desc$se * overall_desc$t_critical)
overall_desc$ci_half_width <- overall_desc$se * overall_desc$t_critical

print("Overall Descriptive Statistics with 95% Confidence Intervals:")
print(overall_desc[, c("n", "mean", "sd", "min", "max", "ci_lower", "ci_upper")])


## ----descriptives-raw-group------------------------------------------------------------------------------------------------------------
# Use summarySE for each variable
print("Descriptive Statistics by Treatment Group with 95% Confidence Intervals:")

group_desc_with_ci <- lapply(vars_to_summarize, function(var) {
  result <- summarySE(data = school_climate_data, 
                     measurevar = var,
                     groupvars = "group", 
                     conf.interval = 0.95,
                     na.rm = TRUE)
  
  # Add explicit CI bounds for clarity
  result$ci_lower <- result[[var]] - result$ci
  result$ci_upper <- result[[var]] + result$ci
  
  # Add min/max manually
  min_vals <- aggregate(school_climate_data[[var]], 
                       by = list(school_climate_data$group), 
                       FUN = min, na.rm = TRUE)
  max_vals <- aggregate(school_climate_data[[var]], 
                       by = list(school_climate_data$group), 
                       FUN = max, na.rm = TRUE)
  
  result$min <- min_vals$x[match(result$group, min_vals$Group.1)]
  result$max <- max_vals$x[match(result$group, max_vals$Group.1)]
  
  cat("\n", var, ":\n")
  print(result[, c("group", "N", var, "sd", "min", "max", "ci_lower", "ci_upper")])
  
  return(result)
})

# Name the list elements for easy access
names(group_desc_with_ci) <- vars_to_summarize

# Create a combined long-format data frame for further analysis
group_desc_combined <- do.call(rbind, lapply(names(group_desc_with_ci), function(var) {
  result <- group_desc_with_ci[[var]]
  result$variable <- var
  result$value <- result[[var]]  # Standardize the variable name
  return(result[, c("variable", "group", "N", "value", "sd", "min", "max", "ci_lower", "ci_upper")])
}))

print("\nCombined Group Statistics (Organized by Group):")

# Get unique groups and print each group's data
unique_groups <- unique(group_desc_combined$group)

for(group_name in unique_groups) {
  cat("\n", rep("=", 50), "\n")
  cat("GROUP:", group_name, "\n")
  cat(rep("=", 50), "\n")
  
  # Filter data for this group
  group_data <- group_desc_combined[group_desc_combined$group == group_name, ]
  
  # Print without the group column since it's redundant
  print(group_data[, c("variable", "N", "value", "sd", "min", "max", "ci_lower", "ci_upper")])
}



## ----descriptives-dvs------------------------------------------------------------------------------------------------------------------
# DVs are identical across imputations, so use first imputation
dv_desc <- describe(imputed_climate_data[[1]][, dvs])

# Add confidence intervals to overall descriptives
dv_desc$df <- dv_desc$n - 1
dv_desc$t_critical <- qt(0.975, dv_desc$df)  # 95% CI, two-tailed
dv_desc$ci_lower <- dv_desc$mean - (dv_desc$se * dv_desc$t_critical)
dv_desc$ci_upper <- dv_desc$mean + (dv_desc$se * dv_desc$t_critical)
dv_desc$ci_half_width <- dv_desc$se * dv_desc$t_critical

print("Descriptive Statistics - DVs (Imputed Data) with 95% Confidence Intervals:")
print(dv_desc[, c("n", "mean", "sd", "se", "min", "max", "ci_lower", "ci_upper")])

# Descriptive Statistics by Group with CIs
# Use summarySE for each DV (industry standard for grouped CIs)
print("Descriptive Statistics - DVs (By Group) with 95% Confidence Intervals:")

dv_by_group_ci <- lapply(dvs, function(var) {
  result <- summarySE(data = imputed_climate_data[[1]], 
                     measurevar = var,
                     groupvars = "group", 
                     conf.interval = 0.95)
  
  # Add explicit CI bounds for clarity
  result$ci_lower <- result[[var]] - result$ci
  result$ci_upper <- result[[var]] + result$ci
  
  # Add min/max manually
  min_vals <- aggregate(imputed_climate_data[[1]][[var]], 
                       by = list(imputed_climate_data[[1]]$group), 
                       FUN = min, na.rm = TRUE)
  max_vals <- aggregate(imputed_climate_data[[1]][[var]], 
                       by = list(imputed_climate_data[[1]]$group), 
                       FUN = max, na.rm = TRUE)
  
  result$min <- min_vals$x[match(result$group, min_vals$Group.1)]
  result$max <- max_vals$x[match(result$group, max_vals$Group.1)]
  
  cat("\n", var, ":\n")
  print(result[, c("group", "N", var, "sd", "se", "min", "max", "ci_lower", "ci_upper")])
  
  return(result)
})

# Name the list elements for easy access
names(dv_by_group_ci) <- dvs

# Optional: Create a combined long-format data frame for further analysis
dv_by_group_combined <- do.call(rbind, lapply(names(dv_by_group_ci), function(var) {
  result <- dv_by_group_ci[[var]]
  result$variable <- var
  result$value <- result[[var]]  # Standardize the variable name
  return(result[, c("variable", "group", "N", "value", "sd", "se", "min", "max", "ci_lower", "ci_upper")])
}))

print("\nCombined Group Statistics - DVs (Organized by Group):")

# Get unique groups and print each group's data
unique_groups <- unique(dv_by_group_combined$group)

for(group_name in unique_groups) {
  cat("\n", rep("=", 50), "\n")
  cat("GROUP:", group_name, "\n")
  cat(rep("=", 50), "\n")
  
  # Filter data for this group
  group_data <- dv_by_group_combined[dv_by_group_combined$group == group_name, ]
  
  # Print without the group column since it's redundant
  print(group_data[, c("variable", "N", "value", "sd", "se", "min", "max", "ci_lower", "ci_upper")])
}



## ----descriptives-covariate------------------------------------------------------------------------------------------------------------
# Define the transformed covariate variables
transformed_covariates <- c("covariate_ela_centered", "covariate_ela_squared")

# Function to calculate descriptives across all imputations
calculate_descriptives_across_imputations <- function(var_name) {
  
  # Calculate descriptives for each imputation
  imputation_results <- data.frame(
    dataset = paste("Imputation", 1:5),
    n = sapply(1:5, function(i) length(imputed_climate_data[[i]][[var_name]])),
    mean = sapply(1:5, function(i) mean(imputed_climate_data[[i]][[var_name]])),
    sd = sapply(1:5, function(i) sd(imputed_climate_data[[i]][[var_name]])),
    se = sapply(1:5, function(i) {
      x <- imputed_climate_data[[i]][[var_name]]
      sd(x)/sqrt(length(x))
    }),
    min = sapply(1:5, function(i) min(imputed_climate_data[[i]][[var_name]])),
    max = sapply(1:5, function(i) max(imputed_climate_data[[i]][[var_name]]))
  )
  
  # Add confidence intervals for individual imputations
  imputation_results$df <- imputation_results$n - 1
  imputation_results$t_critical <- qt(0.975, imputation_results$df)
  imputation_results$ci_lower <- imputation_results$mean - imputation_results$t_critical * imputation_results$se
  imputation_results$ci_upper <- imputation_results$mean + imputation_results$t_critical * imputation_results$se
  
  # Print individual results
  cat("\n", rep("=", 60), "\n")
  cat("VARIABLE:", toupper(var_name), "\n")
  cat(rep("=", 60), "\n")
  
  print("Individual Imputation Results:")
  print(imputation_results[, c("dataset", "n", "mean", "sd", "min", "max", "ci_lower", "ci_upper")])
  
  # Pooled descriptives with imputation uncertainty adjustment
  all_values <- unlist(lapply(1:5, function(i) imputed_climate_data[[i]][[var_name]]))
  
  pooled_n <- length(all_values)
  pooled_mean <- mean(all_values)
  pooled_sd <- sd(all_values)
  pooled_min <- min(all_values)
  pooled_max <- max(all_values)
  
  # Adjust for imputation uncertainty
  between_imputation_var <- var(imputation_results$mean)
  within_imputation_var <- mean(imputation_results$se^2)
  m <- 5
  adjusted_se <- sqrt(within_imputation_var + between_imputation_var + (between_imputation_var / m))
  
  # Conservative CI
  conservative_df <- min(imputation_results$df)
  t_critical <- qt(0.975, conservative_df)
  ci_lower <- pooled_mean - t_critical * adjusted_se
  ci_upper <- pooled_mean + t_critical * adjusted_se
  
  cat("\nPooled Statistics (All 5 Imputations):\n")
  cat("Total N:", pooled_n, "\n")
  cat("Mean:", sprintf("%.3f", pooled_mean), "\n")
  cat("SD:", sprintf("%.3f", pooled_sd), "\n")
  cat("Range:", sprintf("%.3f to %.3f", pooled_min, pooled_max), "\n")
  cat("95% CI:", sprintf("%.3f to %.3f", ci_lower, ci_upper), "\n")
  
  return(list(
    individual = imputation_results,
    pooled = list(n = pooled_n, mean = pooled_mean, sd = pooled_sd, 
                  min = pooled_min, max = pooled_max, 
                  ci_lower = ci_lower, ci_upper = ci_upper)
  ))
}

# Apply to both transformed covariates
print("COVARIATE DESCRIPTIVE STATISTICS ACROSS ALL IMPUTATIONS")

centered_results <- calculate_descriptives_across_imputations("covariate_ela_centered")
squared_results <- calculate_descriptives_across_imputations("covariate_ela_squared")

# Summary table
cat("\n", rep("=", 60), "\n")
cat("SUMMARY TABLE\n")
cat(rep("=", 60), "\n")

summary_table <- data.frame(
  Variable = c("Centered Covariate", "Squared Covariate"),
  Total_N = c(centered_results$pooled$n, squared_results$pooled$n),
  Mean = c(sprintf("%.3f", centered_results$pooled$mean), 
           sprintf("%.3f", squared_results$pooled$mean)),
  SD = c(sprintf("%.3f", centered_results$pooled$sd), 
         sprintf("%.3f", squared_results$pooled$sd)),
  Range = c(
    paste0(sprintf("%.3f", centered_results$pooled$min), " to ", 
           sprintf("%.3f", centered_results$pooled$max)),
    paste0(sprintf("%.3f", squared_results$pooled$min), " to ", 
           sprintf("%.3f", squared_results$pooled$max))
  ),
  CI_95 = c(
    paste0(sprintf("%.3f", centered_results$pooled$ci_lower), " to ", 
           sprintf("%.3f", centered_results$pooled$ci_upper)),
    paste0(sprintf("%.3f", squared_results$pooled$ci_lower), " to ", 
           sprintf("%.3f", squared_results$pooled$ci_upper))
  )
)

print(summary_table)


## ----dvs-by-group----------------------------------------------------------------------------------------------------------------------
# Unadjusted Means, Standard Deviations, and Confidence Intervals
group_desc <- imputed_climate_data[[1]] |>
  group_by(group) |>
  summarise(across(all_of(dvs), 
                   list(mean = ~mean(., na.rm = TRUE), 
                        sd = ~sd(., na.rm = TRUE),
                        n = ~sum(!is.na(.)),
                        se = ~sd(., na.rm = TRUE)/sqrt(sum(!is.na(.))))),
            .groups = "drop")

# Calculate confidence intervals using t-distribution
# Create helper columns for CI calculations
ci_data <- group_desc |>
  # Extract n values for each DV (they should all be the same since DVs have no missing data)
  mutate(
    # Get n from first DV (all DVs should have same n)
    n = get(paste0(dvs[1], "_n"))
  ) |>
  # Calculate df and t-critical
  mutate(
    df = n - 1,
    t_critical = qt(0.975, df)
  )

# Add CI calculations for each DV
for (dv in dvs) {
  mean_col <- paste0(dv, "_mean")
  se_col <- paste0(dv, "_se") 
  ci_lower_col <- paste0(dv, "_ci_lower")
  ci_upper_col <- paste0(dv, "_ci_upper")
  
  ci_data[[ci_lower_col]] <- ci_data[[mean_col]] - (ci_data[[se_col]] * ci_data$t_critical)
  ci_data[[ci_upper_col]] <- ci_data[[mean_col]] + (ci_data[[se_col]] * ci_data$t_critical)
}

print("Descriptive Statistics by Group (Unadjusted) with 95% Confidence Intervals:")

# Print results in a clean format for each DV
for (dv in dvs) {
  cat("\n", dv, ":\n")
  
  # Select relevant columns for this DV
  dv_results <- ci_data |>
    select(group, 
           n = all_of(paste0(dv, "_n")),
           mean = all_of(paste0(dv, "_mean")),
           sd = all_of(paste0(dv, "_sd")),
           se = all_of(paste0(dv, "_se")),
           ci_lower = all_of(paste0(dv, "_ci_lower")),
           ci_upper = all_of(paste0(dv, "_ci_upper")))
  
  print(dv_results, digits = 3)
}


## ----covariate-by-group----------------------------------------------------------------------------------------------------------------
# Covariate descriptives by treatment group across imputations
variables_to_analyze <- c("covariate_ela_centered", "covariate_ela_squared")

# Function to calculate group descriptives across imputations
calculate_group_descriptives <- function(var_name) {
  
  # Get unique groups (exclude overall)
  groups <- unique(as.character(imputed_climate_data[[1]]$group))
  
  # Store pooled results for summary table
  pooled_summary <- data.frame()
  
  cat("\n", rep("=", 60), "\n")
  cat("VARIABLE:", toupper(var_name), "\n")
  cat(rep("=", 60), "\n")
  
  # Process each group (treatment groups only)
  for(group_name in groups) {
    
    # Extract values for this group across all imputations
    group_values <- lapply(1:5, function(i) {
      subset_data <- imputed_climate_data[[i]][as.character(imputed_climate_data[[i]]$group) == group_name, ]
      subset_data[[var_name]]
    })
    
    # Individual imputation results for this group
    individual <- data.frame(
      dataset = paste("Imputation", 1:5),
      n = sapply(group_values, length),
      mean = sapply(group_values, mean),
      sd = sapply(group_values, sd),
      se = sapply(group_values, function(x) sd(x)/sqrt(length(x))),
      min = sapply(group_values, min),
      max = sapply(group_values, max)
    )
    
    # Add CIs for individual imputations
    individual$df <- individual$n - 1
    individual$t_critical <- qt(0.975, individual$df)
    individual$ci_lower <- individual$mean - individual$t_critical * individual$se
    individual$ci_upper <- individual$mean + individual$t_critical * individual$se
    
    # Print individual results
    cat("\nGroup:", group_name, "\n")
    print(individual[, c("dataset", "n", "mean", "sd", "min", "max", "ci_lower", "ci_upper")])
    
    # Pooled descriptives with imputation uncertainty adjustment
    all_values <- unlist(group_values)
    
    pooled_n <- length(all_values)
    pooled_mean <- mean(all_values)
    pooled_sd <- sd(all_values)
    pooled_min <- min(all_values)
    pooled_max <- max(all_values)
    
    # Adjust for imputation uncertainty
    between_imputation_var <- var(individual$mean)
    within_imputation_var <- mean(individual$se^2)
    m <- 5
    adjusted_se <- sqrt(within_imputation_var + between_imputation_var + (between_imputation_var / m))
    
    # Conservative CI
    conservative_df <- min(individual$df)
    t_critical <- qt(0.975, conservative_df)
    ci_lower <- pooled_mean - t_critical * adjusted_se
    ci_upper <- pooled_mean + t_critical * adjusted_se
    
    cat("\nPooled Statistics for", group_name, ":\n")
    cat("Total N:", pooled_n, "\n")
    cat("Mean:", sprintf("%.3f", pooled_mean), "\n")
    cat("SD:", sprintf("%.3f", pooled_sd), "\n")
    cat("Range:", sprintf("%.3f to %.3f", pooled_min, pooled_max), "\n")
    cat("95% CI:", sprintf("%.3f to %.3f", ci_lower, ci_upper), "\n")
    
    # Store for summary table
    pooled_summary <- rbind(pooled_summary, data.frame(
      Variable = var_name,
      Group = group_name,
      Total_N = pooled_n,
      Mean = sprintf("%.3f", pooled_mean),
      SD = sprintf("%.3f", pooled_sd),
      Range = paste0(sprintf("%.3f", pooled_min), " to ", sprintf("%.3f", pooled_max)),
      CI_95 = paste0(sprintf("%.3f", ci_lower), " to ", sprintf("%.3f", ci_upper))
    ))
  }
  
  return(pooled_summary)
}

# Process both variables
print("COVARIATE DESCRIPTIVE STATISTICS BY GROUP ACROSS ALL IMPUTATIONS")

all_pooled_results <- data.frame()
for(var_name in variables_to_analyze) {
  var_results <- calculate_group_descriptives(var_name)
  all_pooled_results <- rbind(all_pooled_results, var_results)
}

# Final summary table
cat("\n", rep("=", 80), "\n")
cat("COMPLETE SUMMARY TABLE - ALL VARIABLES AND GROUPS\n")
cat(rep("=", 80), "\n")

print(all_pooled_results)


## ----descriptives-by-gender------------------------------------------------------------------------------------------------------------
# 3. Descriptive Statistics by Gender
gender_desc <- imputed_climate_data[[1]] |>
  group_by(gender) |>
  summarise(across(c(covariate, all_of(dvs)), 
                   list(mean = ~mean(., na.rm = TRUE), 
                        sd = ~sd(., na.rm = TRUE))))
print("Descriptive Statistics by Gender (Unadjusted):")
print(gender_desc)


## ----descriptives-by-race--------------------------------------------------------------------------------------------------------------
# 4. Descriptive Statistics by Race/Ethnicity
race_desc <- imputed_climate_data[[1]] |>
  group_by(race) |>
  summarise(across(c(covariate, all_of(dvs)), 
                   list(mean = ~mean(., na.rm = TRUE), 
                        sd = ~sd(., na.rm = TRUE))))
print("Descriptive Statistics by Race/Ethnicity (Unadjusted):")
print(race_desc)


## ----adjusted-means--------------------------------------------------------------------------------------------------------------------
# Function to compute adjusted means with proper MI uncertainty
compute_adjusted_means_mi <- function(dv) {
  cat("Computing adjusted means for:", dv, "\n")
  
  # Fit models to each imputation directly
  formula_str <- paste(dv, "~ group + covariate_ela_centered + covariate_ela_squared")
  
  # Fit models manually to each dataset
  models_list <- map(imputed_climate_data, ~lm(as.formula(formula_str), data = .x))
  
  # Compute emmeans for each model
  emmeans_list <- map(models_list, ~{
    emm <- emmeans(.x, specs = "group")
    as.data.frame(emm)
  })
  
  # Pool using Rubin's rules manually for each group
  groups <- emmeans_list[[1]]$group
  pooled_results <- map_dfr(seq_along(groups), ~{
    # Extract estimates and variances for this group across imputations
    estimates <- map_dbl(emmeans_list, function(df) df$emmean[.x])
    variances <- map_dbl(emmeans_list, function(df) df$SE[.x]^2)
    
    # Apply Rubin's rules
    pooled_estimate <- mean(estimates)
    within_var <- mean(variances)
    between_var <- var(estimates)
    total_var <- within_var + between_var + (between_var/length(estimates))
    pooled_se <- sqrt(total_var)
    
    # Calculate degrees of freedom (Barnard-Rubin adjustment)
    r <- (between_var + between_var/length(estimates)) / within_var
    df_old <- (length(estimates) - 1) / r^2
    df_residual_mean <- mean(map_dbl(models_list, df.residual))
    df_adj <- df_residual_mean * (df_residual_mean + 1) / (df_residual_mean + 3) * (1 - r)
    
    # Ensure df is reasonable
    if(is.na(df_adj) || df_adj <= 0) {
      df_adj <- df_residual_mean
    }
    
    data.frame(
      group = groups[.x],
      emmean = pooled_estimate,
      SE = pooled_se,
      df = df_adj,
      ci_lower = pooled_estimate - qt(0.975, df_adj) * pooled_se,
      ci_upper = pooled_estimate + qt(0.975, df_adj) * pooled_se,
      # Keep the original column names for compatibility
      lower.CL = pooled_estimate - qt(0.975, df_adj) * pooled_se,
      upper.CL = pooled_estimate + qt(0.975, df_adj) * pooled_se
    )
  })
  
  return(pooled_results)
}

# Compute adjusted means for each dependent variable
adjusted_means_group <- lapply(dvs, compute_adjusted_means_mi)
names(adjusted_means_group) <- dvs

print("Adjusted Means by Group (Controlling for ELA Covariate, Pooled using Rubin's Rules):")
print("Note: CIs are 95% confidence intervals accounting for multiple imputation uncertainty")

# Print results with clear CI labeling and range information
for(dv_name in names(adjusted_means_group)) {
  cat("\n", dv_name, ":\n")
  
  # Select key columns for clean display
  display_results <- adjusted_means_group[[dv_name]][, c("group", "emmean", "SE", "df", "ci_lower", "ci_upper")]
  
  # Round for better display
  display_results$emmean <- round(display_results$emmean, 3)
  display_results$SE <- round(display_results$SE, 3)
  display_results$df <- round(display_results$df, 1)  # Round df to 1 decimal place
  display_results$ci_lower <- round(display_results$ci_lower, 3)
  display_results$ci_upper <- round(display_results$ci_upper, 3)
  
  print(display_results)
  
  # Add range information for adjusted means
  emmean_min <- min(display_results$emmean)
  emmean_max <- max(display_results$emmean)
  emmean_range <- emmean_max - emmean_min
  
  cat("Range of adjusted means:", sprintf("%.3f to %.3f", emmean_min, emmean_max), 
      " (width:", sprintf("%.3f", emmean_range), ")\n")
  
  # Show CI widths for comparison
  ci_widths <- display_results$ci_upper - display_results$ci_lower
  cat("95% CI widths:", paste(round(ci_widths, 3), collapse = ", "), "\n")
}


## ----correlation-matrix----------------------------------------------------------------------------------------------------------------
# 6. Correlation Matrix of Dependent Variables
cor_matrix <- cor(imputed_climate_data[[1]][, dvs], use = "complete.obs")
print("Correlation Matrix of Dependent Variables:")
print(round(cor_matrix, 3))


## ----missing-data-summary--------------------------------------------------------------------------------------------------------------
# 7. Missing Data Summary
missing_summary <- school_climate_data |>
  summarise(across(c(covariate, all_of(dvs)), 
                   ~sum(is.na(.)) / n() * 100))
print("Percentage of Missing Data:")
print(missing_summary)


## ----mancova-test----------------------------------------------------------------------------------------------------------------------
# Function to run both test statistics
run_mancova_both_tests <- function(dv_matrix, data) {
    manova_fit <- manova(dv_matrix ~ group + covariate_ela_centered + covariate_ela_squared, data = data)
    
    wilks_result <- Manova(manova_fit, type = 3, test.statistic = "Wilks")
    pillai_result <- Manova(manova_fit, type = 3, test.statistic = "Pillai")
    
    return(list(wilks = wilks_result, pillai = pillai_result))
}

# Create placeholders for data to pool results across imputations
mancova_results <- list()
mancova_results_with_outliers <- list()

# Iterate across all imputations and store results in the placeholders
for (i in 1:5) {
    # Create DV matrix for this imputation
    dv_matrix <- as.matrix(imputed_climate_data[[i]][, dvs])
    dv_matrix_with_outliers <- as.matrix(imputed_climate_with_outliers[[i]][, dvs])

    # Run both test statistics
    mancova_results[[i]] <- run_mancova_both_tests(dv_matrix, imputed_climate_data[[i]])
    mancova_results_with_outliers[[i]] <- run_mancova_both_tests(dv_matrix_with_outliers, imputed_climate_with_outliers[[i]])
    
    # Print results for comparison
    cat("\nImputation ", i, "\n==============\n")
    cat("\nClean data (n = 83) - Wilks' Lambda:\n")
    print(mancova_results[[i]]$wilks)
    cat("\nClean data (n = 83) - Pillai's Trace:\n")
    print(mancova_results[[i]]$pillai)
    
    cat("\nData with Outliers (n = 86) - Wilks' Lambda:\n")
    print(mancova_results_with_outliers[[i]]$wilks)
    cat("\nData with Outliers (n = 86) - Pillai's Trace:\n")
    print(mancova_results_with_outliers[[i]]$pillai)
}


## ----effect-sizes----------------------------------------------------------------------------------------------------------------------
# Function to extract group effect size from summary
calculate_group_effect_size <- function(manova_result, test_statistic) {
    # Get the summary which contains the test statistics
    summary_result <- summary(manova_result)
    
    # Find the group term in the summary
    # The summary contains multiple sections, we need the group section
    summary_text <- capture.output(summary_result)
    
    # Find the line with the group multivariate tests
    group_section_start <- grep("Term: group", summary_text)
    multivariate_start <- grep("Multivariate Tests: group", summary_text)
    
    # Extract the line with the test statistic we want
    if (test_statistic == "Wilks") {
        wilks_line <- summary_text[grep("Wilks", summary_text[multivariate_start:(multivariate_start + 10)])[1] + multivariate_start - 1]
        # Extract the test stat value (3rd column)
        lambda <- as.numeric(strsplit(trimws(wilks_line), "\\s+")[[1]][3])
        # Calculate partial eta squared: 1 - lambda^(1/s), where s = min(hypothesis df, error df)
        s <- min(2, manova_result$error.df)  # group has 2 df, error df is stored in object
        partial_eta_sq <- 1 - lambda^(1/s)
        
    } else if (test_statistic == "Pillai") {
        pillai_line <- summary_text[grep("Pillai", summary_text[multivariate_start:(multivariate_start + 10)])[1] + multivariate_start - 1]
        # Extract the test stat value (3rd column)
        pillai <- as.numeric(strsplit(trimws(pillai_line), "\\s+")[[1]][3])
        # Calculate partial eta squared: V / (s + V), where V = Pillai's trace
        s <- min(2, manova_result$error.df)  # group has 2 df
        partial_eta_sq <- pillai / (s + pillai)
    }
    
    return(partial_eta_sq)
}

# Calculate group effect sizes for both test statistics
group_eta_wilks <- numeric(5)
group_eta_pillai <- numeric(5)

for (i in 1:5) {
    dv_matrix <- as.matrix(imputed_climate_data[[i]][, dvs])
    manova_object <- manova(dv_matrix ~ group + covariate_ela_centered + covariate_ela_squared,
                            data = imputed_climate_data[[i]])
    
    manova_wilks <- Manova(manova_object, type = 3, test.statistic = "Wilks")
    manova_pillai <- Manova(manova_object, type = 3, test.statistic = "Pillai")
    
    group_eta_wilks[i] <- calculate_group_effect_size(manova_wilks, "Wilks")
    group_eta_pillai[i] <- calculate_group_effect_size(manova_pillai, "Pillai")
}


## ----effect-ranges---------------------------------------------------------------------------------------------------------------------
# Report ranges
cat("Group Effect Sizes (Partial Eta-Squared) Across Imputations:\n")
cat("Wilks' lambda: ", round(min(group_eta_wilks), 3), " to ", round(max(group_eta_wilks), 3), "\n")
cat("Pillai's trace: ", round(min(group_eta_pillai), 3), " to ", round(max(group_eta_pillai), 3), "\n")

# Show individual values
effect_size_table <- data.frame(
    Imputation = 1:5,
    Wilks = round(group_eta_wilks, 3),
    Pillai = round(group_eta_pillai, 3)
)
print(effect_size_table)


## ----contrast-functions----------------------------------------------------------------------------------------------------------------
climate_imputation_list <- imputationList(imputed_climate_data)

# Functions for planned contrasts
run_multivariate_contrast <- function(contrast_name, results_label) {
    cat("=== MULTIVARIATE PLANNED CONTRAST:", results_label, "===\n")
    coefficient_name <- paste0("group", contrast_name, " = 0")
    
    iwalk(imputed_climate_data, ~{
        cat("\nImputation", .y, ":\n")
        
        # Recreate the manova model
        dv_matrix <- as.matrix(.x[, dvs])
        manova_model <- manova(dv_matrix ~ group + covariate_ela_centered + covariate_ela_squared, data = .x)
        print(linearHypothesis(manova_model, coefficient_name))
    })
}

run_univariate_contrasts <- function(contrast_name) {
    contrast_coef_name <- paste0("group", contrast_name)
    contrast_weights_ss <- if(contrast_name == "PBLvsControl") sum(c(-2, 1, 1)^2) else sum(c(0, -1, 1)^2)
    
    map_dfr(dvs, ~{
        formula_str <- paste(.x, "~ group + covariate_ela_centered + covariate_ela_squared")
        models <- with(climate_imputation_list, lm(as.formula(formula_str)))
        pooled_model <- MIcombine(models)
        
        coeffs <- coef(pooled_model)
        vcov_matrix <- vcov(pooled_model)
        
        # Extract key information
        contrast_estimate <- coeffs[contrast_coef_name]
        contrast_se <- sqrt(vcov_matrix[contrast_coef_name, contrast_coef_name])
        df <- pooled_model$df[1]
        
        # Calculate statistics
        t_stat <- contrast_estimate / contrast_se
        p_value <- 2 * (1 - pt(abs(t_stat), df))
        
        # Cohen's d and CI
        cohens_d <- t_stat * sqrt(contrast_weights_ss / nrow(imputed_climate_data[[1]]))
        t_critical <- qt(0.975, df)
        d_se <- sqrt((nrow(imputed_climate_data[[1]]) * contrast_weights_ss + cohens_d^2 * contrast_weights_ss) / 
                    (nrow(imputed_climate_data[[1]]) * contrast_weights_ss))
        d_ci_lower <- cohens_d - t_critical * d_se
        d_ci_upper <- cohens_d + t_critical * d_se
        
        data.frame(
            DV = .x,
            Estimate = round(as.numeric(contrast_estimate), 3),
            t = round(as.numeric(t_stat), 3),
            df = round(as.numeric(df), 1),
            p = round(as.numeric(p_value), 4),
            Cohens_d = round(as.numeric(cohens_d), 3),
            d_CI_lower = round(as.numeric(d_ci_lower), 3),
            d_CI_upper = round(as.numeric(d_ci_upper), 3),
            stringsAsFactors = FALSE
        )
    })
}


## ----contrast-1-multi------------------------------------------------------------------------------------------------------------------
run_multivariate_contrast("PBLvsControl", "PjBL vs Control")


## ----contrast-1-dvs--------------------------------------------------------------------------------------------------------------------
contrast_summary_1 <- suppressWarnings(run_univariate_contrasts("PBLvsControl"))
print(contrast_summary_1)


## ----contrast-2-multi------------------------------------------------------------------------------------------------------------------
run_multivariate_contrast("CRvsSPBL", "CR-PjBL vs S-PjBL")


## ----contrast-2-dvs--------------------------------------------------------------------------------------------------------------------
contrast_summary_2 <- run_univariate_contrasts("CRvsSPBL")
print(contrast_summary_2)

