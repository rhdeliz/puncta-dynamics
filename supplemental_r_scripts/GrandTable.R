args = commandArgs(trailingOnly=TRUE)

# Load libraries
pacman::p_load(data.table, arrow, kit)

# Set directory
setwd(args[1])

# Read table and give status
custom_fread <- function(x){
  df <- fread(x)
  id <- which(x == paths)
  status <- paste("Reading table", id, "of", NROW(paths))
  print(status)
  return(df)
}

# Essentials ----
# Find tables and remove calibration, cell paths
paths <- try(system("find . -name 'Essential.csv.gz'", intern = TRUE))
paths <- paths[!grepl("Calibrations", paths)]
paths <- paths[!grepl("/Cell_", paths)]
paths <- paths[!grepl("pending_processing", paths)]

table <- lapply(paths, custom_fread)
table <- rbindlist(table, fill = TRUE)

# Keep only ceirtain cohorts
Cohorts <- table[!fduplicated(table$COHORT), ]
Cohorts <- Cohorts$COHORT
Cohorts <- Cohorts[!grepl("half|inhibit|DMSO", Cohorts, ignore.case = T)]
table <- table[table$COHORT %in% Cohorts]

print("Writing table")
if(file.exists("Essential.gz.parquet")){
  file.remove("Essential.gz.parquet")
}
write_parquet(table, "Essential.gz.parquet")

if(file.exists("Essential.csv.gz")){
  file.remove("Essential.csv.gz")
}
fwrite(table, "Essential.csv.gz")


# Calibrations ----
# Find tables and remove calibration, cell paths
paths <- try(system("find . -name 'Essential.csv.gz'", intern = TRUE))
paths <- paths[grepl("Calibrations", paths)]
paths <- paths[!grepl("/Cell_", paths)]
paths <- paths[!grepl("pending_processing", paths)]

table <- lapply(paths, custom_fread)
table <- rbindlist(table, fill = TRUE)

print("Writing calibrations table")
if(file.exists("Calibrations.gz.parquet")){
  file.remove("Calibrations.gz.parquet")
}
write_parquet(table, "Calibrations.gz.parquet")

if(file.exists("Calibrations.csv.gz")){
  file.remove("Calibrations.csv.gz")
}
fwrite(table, "Calibrations.csv.gz")



# Parameters ----
# Find tables and remove calibration, cell paths
paths <- try(system("find . -name 'Parameters.csv.gz'", intern = TRUE))
paths <- paths[!grepl("Calibrations", paths)]
paths <- paths[!grepl("/Cell_", paths)]
paths <- paths[!grepl("pending_processing", paths)]

table <- lapply(paths, custom_fread)
table <- rbindlist(table, fill = TRUE)

# Keep only ceirtain cohorts
Cohorts <- table[!fduplicated(table$COHORT), ]
Cohorts <- Cohorts$COHORT
Cohorts <- Cohorts[!grepl("half|inhibit|DMSO", Cohorts, ignore.case = T)]
table <- table[table$COHORT %in% Cohorts]

print("Writing table")
if(file.exists("Parameters.gz.parquet")){
  file.remove("Parameters.gz.parquet")
}
write_parquet(table, "Parameters.gz.parquet")

if(file.exists("Parameters.csv.gz")){
  file.remove("Parameters.csv.gz")
}
fwrite(table, "Parameters.csv.gz")
