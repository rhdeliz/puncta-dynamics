setwd(OUTPUT_DIRECTORY)

# Loop in parallel to speed it up
TableFx <- function(TableX){
  tryCatch({
    # Temporary table
    TempTable <-
      TableX %>% 
      arrange(
        UNIVERSAL_TRACK_ID,
        FRAMES_ADJUSTED
      ) %>% 
      as.data.table()
    
    # Sort protein names
    TempTable$REFERENCE_PROTEIN = USE_REFERENCE_PROTEIN
    TempTable$QUERY_PROTEIN = factor(TempTable$QUERY_PROTEIN, levels = PROTEIN_ORDER)
    
    # Apply moving average
    if(MOVING_AVERAGE_WINDOW > 0){
      # TempTable$REFERENCE_TOTAL_INTENSITY_BACKUP = TempTable$REFERENCE_TOTAL_INTENSITY
      # TempTable$QUERY_TOTAL_INTENSITY_BACKUP = TempTable$QUERY_TOTAL_INTENSITY

      TempTable <- TempTable[, REFERENCE_TOTAL_INTENSITY := roll_median(REFERENCE_TOTAL_INTENSITY, n = 3, fill = NA, align = "center"), by = UNIVERSAL_TRACK_ID]
      TempTable <- TempTable[, QUERY_TOTAL_INTENSITY := roll_median(QUERY_TOTAL_INTENSITY, n = 3, fill = NA, align = "center"), by = UNIVERSAL_TRACK_ID]
      
      TempTable <- TempTable[, REFERENCE_TOTAL_INTENSITY := roll_median(REFERENCE_TOTAL_INTENSITY, n = MOVING_AVERAGE_WINDOW, fill = NA, align = "center"), by = UNIVERSAL_TRACK_ID]
      TempTable <- TempTable[, QUERY_TOTAL_INTENSITY := roll_median(QUERY_TOTAL_INTENSITY, n = MOVING_AVERAGE_WINDOW, fill = NA, align = "center"), by = UNIVERSAL_TRACK_ID]
      
      TempTable <- TempTable[!is.na(TempTable$REFERENCE_TOTAL_INTENSITY)]
    }
    
    # Calculate intensities
    TempTable <- TempTable[, LAG_REFERENCE_TOTAL_INTENSITY := shift(.(REFERENCE_TOTAL_INTENSITY), LEAD_LAG, type = "lag"), by = UNIVERSAL_TRACK_ID]
    TempTable <- TempTable[, LEAD_REFERENCE_TOTAL_INTENSITY := shift(.(REFERENCE_TOTAL_INTENSITY) ,LEAD_LAG, type = "lead"), by = UNIVERSAL_TRACK_ID]
    TempTable <- TempTable[, LAG_QUERY_TOTAL_INTENSITY := shift(.(QUERY_TOTAL_INTENSITY),LEAD_LAG, type = "lag"), by = UNIVERSAL_TRACK_ID]
    TempTable <- TempTable[, LEAD_QUERY_TOTAL_INTENSITY := shift(.(QUERY_TOTAL_INTENSITY),LEAD_LAG, type = "lead"), by = UNIVERSAL_TRACK_ID]
    # Calculate change in intensity
    TempTable <- TempTable[, DELTA_REFERENCE_TOTAL_INTENSITY := LEAD_REFERENCE_TOTAL_INTENSITY-LAG_REFERENCE_TOTAL_INTENSITY, by = UNIVERSAL_TRACK_ID]
    TempTable <- TempTable[, DELTA_QUERY_TOTAL_INTENSITY := LEAD_QUERY_TOTAL_INTENSITY-LAG_QUERY_TOTAL_INTENSITY, by = UNIVERSAL_TRACK_ID]

    # Calculate change in time
    TempTable <- TempTable[, LAG_TIME := shift(.(TIME),LEAD_LAG, type = "lag"), by = UNIVERSAL_TRACK_ID]
    TempTable <- TempTable[, LEAD_TIME := shift(.(TIME),LEAD_LAG, type = "lead"), by = UNIVERSAL_TRACK_ID]
    TempTable <- TempTable[, DELTA_TIME := LEAD_TIME-LAG_TIME, by = UNIVERSAL_TRACK_ID]
    TempTable <- TempTable[, DELTA_TIME := LEAD_TIME-LAG_TIME, by = UNIVERSAL_TRACK_ID]
    
    # Get reference data
    TempTable <- TempTable[, DELTA_REFERENCE_TOTAL_INTENSITY := DELTA_REFERENCE_TOTAL_INTENSITY/DELTA_TIME, by = UNIVERSAL_TRACK_ID]
    TempTable <- TempTable[, ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY := DELTA_REFERENCE_TOTAL_INTENSITY/REFERENCE_TOTAL_INTENSITY, by = UNIVERSAL_TRACK_ID]
    TempTable <- TempTable[, ROUNDED_REFERENCE_TOTAL_INTENSITY := round(REFERENCE_TOTAL_INTENSITY/STEP_SIZE)*STEP_SIZE]
    # Get query data
    TempTable <- TempTable[, DELTA_QUERY_TOTAL_INTENSITY := DELTA_QUERY_TOTAL_INTENSITY/DELTA_TIME, by = UNIVERSAL_TRACK_ID]
    TempTable <- TempTable[, ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY := DELTA_QUERY_TOTAL_INTENSITY/QUERY_TOTAL_INTENSITY, by = UNIVERSAL_TRACK_ID]
    TempTable <- TempTable[, ROUNDED_QUERY_TOTAL_INTENSITY := round(QUERY_TOTAL_INTENSITY/STEP_SIZE)*STEP_SIZE]
    
    # Adjust change in intensty of query protein
    TempTable <- TempTable[, ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY := DELTA_QUERY_TOTAL_INTENSITY/QUERY_TOTAL_INTENSITY, by = UNIVERSAL_TRACK_ID]
    
    # How long since it landed
    TempTable <- TempTable[, FRAMES_SINCE_LANDING := FRAME - min(TempTable$FRAME) + 1]
    TempTable <- TempTable[, MIN_FRAMES_SINCE_LANDING := min(FRAMES_SINCE_LANDING), by = UNIVERSAL_TRACK_ID]
    
    if(USE_MIN_FRAMES_SINCE_LANDING == F){

      TempTable <-
        TempTable %>% 
        filter(
          FRAMES_SINCE_LANDING <= FRAMES_SINCE_LANDING_THRESHOLD
        ) %>% 
        as.data.table()
      
      TempTable <- TempTable[, FRAMES_SINCE_LANDING_CAT := ceiling(FRAMES_SINCE_LANDING/50)*50]
      
    } else{
      
      TempTable <-
        TempTable %>% 
        filter(
          MIN_FRAMES_SINCE_LANDING <= FRAMES_SINCE_LANDING_THRESHOLD
        ) %>% 
        as.data.table()
      
      TempTable <- TempTable[, FRAMES_SINCE_LANDING_CAT := ceiling(MIN_FRAMES_SINCE_LANDING/50)*50]
    }
    
    # Columns to keep
    TempTable <-
      TempTable %>% 
      # Delete filtering variables
      select(
        COHORT,
        LIGAND_DENSITY_CAT,
        IMAGE,
        CELL,
        
        # Time
        FRAME,
        FRAMES_ADJUSTED,
        TIME_ADJUSTED,
        FRAMES_SINCE_LANDING,
        FRAMES_SINCE_LANDING_CAT,
        MIN_FRAMES_SINCE_LANDING,
        FPS,
        
        UNIVERSAL_TRACK_ID,
        SOURCE,
        
        REFERENCE_PROTEIN,
        REFERENCE_TOTAL_INTENSITY,
        ROUNDED_REFERENCE_TOTAL_INTENSITY,
        DELTA_REFERENCE_TOTAL_INTENSITY,
        ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY,

        QUERY_PROTEIN,
        QUERY_TOTAL_INTENSITY,
        ROUNDED_QUERY_TOTAL_INTENSITY,
        DELTA_QUERY_TOTAL_INTENSITY,
        ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY
      ) %>%
      filter(
        !is.na(DELTA_REFERENCE_TOTAL_INTENSITY)
      ) %>% 
      as.data.table()
    
    return(TempTable)
    # Error message
  }, error = function(e) {print("Error with CellFx")})}
# Run in parallel
StatTable <- mclapply(SplitTable, TableFx, mc.cores = detectCores(logical = F))
# Combine results from all cells
StatTable <- data.table::rbindlist(StatTable)

print(paste("StatTable Images=", NROW(unique(StatTable$IMAGE))))

# Table Name
SaveName <-
  paste0(
    "StatTable - ",
    "LeadLag ", LEAD_LAG, " - ",
    "StepSize ", round(STEP_SIZE, 2),
    ".gz.parquet"
  )
SaveName <- gsub(" - \\.", "\\.", SaveName)

file.remove(file.path(OUTPUT_DIRECTORY, SaveName))
write_parquet(StatTable, file.path(OUTPUT_DIRECTORY, SaveName))

remove(SplitTable)
