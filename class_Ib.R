# =====================================================
# Module I: Getting Started with R - Class Ib Assignment
# =====================================================

# Load necessary libraries
library(readr)
library(dplyr)

# Load dataset from raw_data folder
patient_data <- read_csv("raw_data/patient_info.csv")

# View the dataset in spreadsheet form
View(patient_data)

# Check structure and summary
str(patient_data)
summary(patient_data)

# Convert 'gender' to factor
patient_data$gender <- as.factor(patient_data$gender)

# Convert 'age' to numeric (if it's not already)
patient_data$age <- as.numeric(patient_data$age)

# FIXED: Create a binary smoking column: 1 = Yes, 0 = No
patient_data$smoking_binary <- ifelse(tolower(patient_data$smoker) == "yes", 1, 0)
patient_data$smoking_binary <- as.factor(patient_data$smoking_binary)

# Check updated structure
str(patient_data)

# Save cleaned dataset to clean_data folder
write_csv(patient_data, "clean_data/patient_info_clean.csv")

# Save R workspace
save.image("Aslam_Class_Ib_Assignment.RData")


