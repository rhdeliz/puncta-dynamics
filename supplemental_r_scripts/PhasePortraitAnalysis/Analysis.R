setwd(OUTPUT_DIRECTORY)

FPSTable <- rbindlist(SplitTable)

FPSTable <-
  FPSTable %>% 
  select(
    IMAGE,
    FRAME,
    TIME
  ) %>% 
  distinct() %>% 
  mutate(
    TIME_DELTA = TIME - lag(TIME)
  ) %>% 
  group_by(
    IMAGE
  ) %>% 
  summarize(
    FPS = median(TIME_DELTA, na.rm = T)
  ) %>% 
  mutate(
    FPS = 1/round(FPS)
  ) %>% 
  as.data.table()

SplitTable <- rbindlist(SplitTable)
SplitTable <- merge(SplitTable, FPSTable, by = "IMAGE")
SplitTable <- SplitTable %>% as_tibble() %>% group_split(IMAGE)

# Loop in parallel to speed it up
TableFx <- function(TableX){
  tryCatch({
    
    TempTable <- as.data.table(TableX)
    # Temporary table
    TempTable <-
      TempTable %>% 
      mutate(
        FRAMES_ADJUSTED = as.numeric(FRAMES_ADJUSTED)
      ) %>% 
      arrange(
        UNIVERSAL_TRACK_ID,
        FRAMES_ADJUSTED
      ) %>% 
      as.data.table()
    
    print(TempTable$UNIVERSAL_CELL_ID[1])
    # Sort protein names
    TempTable$QUERY_PROTEIN = factor(TempTable$QUERY_PROTEIN, levels = PROTEIN_ORDER)
    
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
    TempTable <- TempTable[, TIME_ADJUSTED := TIME-min(TIME), by = UNIVERSAL_TRACK_ID]
    
    # Get reference data
    TempTable <- TempTable[, DELTA_REFERENCE_TOTAL_INTENSITY := DELTA_REFERENCE_TOTAL_INTENSITY/DELTA_TIME, by = UNIVERSAL_TRACK_ID]
    TempTable <- TempTable[, ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY := DELTA_REFERENCE_TOTAL_INTENSITY/REFERENCE_TOTAL_INTENSITY, by = UNIVERSAL_TRACK_ID]
    # Get query data
    TempTable <- TempTable[, DELTA_QUERY_TOTAL_INTENSITY := DELTA_QUERY_TOTAL_INTENSITY/DELTA_TIME, by = UNIVERSAL_TRACK_ID]
    TempTable <- TempTable[, ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY := DELTA_QUERY_TOTAL_INTENSITY/QUERY_TOTAL_INTENSITY, by = UNIVERSAL_TRACK_ID]

    # Adjust change in intensty of query protein
    TempTable <- TempTable[, ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY := DELTA_QUERY_TOTAL_INTENSITY/QUERY_TOTAL_INTENSITY, by = UNIVERSAL_TRACK_ID]
    
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
        FPS,
        
        UNIVERSAL_TRACK_ID,
        
        REFERENCE_PROTEIN,
        ORIGINAL_REFERENCE_TOTAL_INTENSITY,
        REFERENCE_TOTAL_INTENSITY,
        DELTA_REFERENCE_TOTAL_INTENSITY,
        ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY,

        QUERY_PROTEIN,
        ORIGINAL_QUERY_TOTAL_INTENSITY,
        QUERY_TOTAL_INTENSITY,
        DELTA_QUERY_TOTAL_INTENSITY,
        ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY,
        
        MIN_REFERENCE,
        MIN_QUERY,
        MAX_REFERENCE,
        MAX_QUERY
      ) %>%
      filter(
        !is.na(DELTA_REFERENCE_TOTAL_INTENSITY)
      ) %>% 
      as.data.table()
    
    return(TempTable)
    # Error message
  }, error = function(e) {print("Error with CellFx")})
}
# Run in parallel
StatTable <- mclapply(SplitTable, TableFx, mc.cores = detectCores(logical = F))
# Combine results from all cells
StatTable <- rbindlist(StatTable)

Images <- StatTable[!fduplicated(StatTable$IMAGE), ]
Images <- Images$IMAGE
print(paste("StatTable Images=", NROW(Images)))
remove(Images)

# Table Name
SaveName <-
  paste0(
    "Normalized StatTable - ",
    "LeadLag ", LEAD_LAG,
    ".gz.parquet"
  )
SaveName <- gsub(" - \\.", "\\.", SaveName)

file.remove(file.path(OUTPUT_DIRECTORY, SaveName))
write_parquet(StatTable, file.path(OUTPUT_DIRECTORY, SaveName))

remove(SplitTable)
