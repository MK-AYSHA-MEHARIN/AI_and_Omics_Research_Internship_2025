# ==============================
# Class Ic Assignment
# Author: Aysha Mehrin
# ==============================

# ------------------------------
# Exercise 1 — Checking cholesterol levels with if...else
# ------------------------------
cholesterol <- 200

if (cholesterol > 240) {
  print("High Cholesterol")
} else {
  print("Cholesterol is normal")
}

# ------------------------------
# Exercise 2 — Blood pressure check with if...else
# ------------------------------
Systolic_bp <- 130

if (Systolic_bp < 120) {
  print("Blood Pressure is normal")
} else {
  print("Blood Pressure is high")
}

# ------------------------------
# Exercise 3 — Converting character columns to factors
# ------------------------------

# ✅ Absolute path to Metadata.csv
metadata_file <- "C:/Users/Aslam/OneDrive/Desktop/AI_Omics_Internship_2025/Module_I/raw_data/Metadata.csv"

# Load metadata if file exists
if (file.exists(metadata_file)) {
  metadata <- read.csv(metadata_file)
  print("✅ Metadata loaded successfully!")
  
  # Convert character columns to factors
  char_cols <- sapply(metadata, is.character)
  metadata[char_cols] <- lapply(metadata[char_cols], as.factor)
  
  print("✅ Converted character columns to factors")
  
} else {
  stop(paste("❌ File not found:", metadata_file))
}

# ------------------------------
# Exercise 4 — Looping through patient IDs
# ------------------------------
for (i in 1:5) {
  print(paste("Patient ID:", i))
}

# ------------------------------
# Exercise 5 — Function for BMI calculation
# ------------------------------
calculate_BMI <- function(weight, height) {
  bmi <- weight / (height * height)
  return(bmi)
}

# Example call
BMI_value <- calculate_BMI(70, 1.75)
print(paste("BMI is:", BMI_value))
