setwd(OUTPUT_DIRECTORY)

# Import compressed table
if(file.exists(file.path(TABLE_PATH, "Essential.gz.parquet"))){
  Images <- read_parquet(file.path(TABLE_PATH, "Essential.gz.parquet"), col_select = "IMAGE")
  Images <- Images[!fduplicated(Images$IMAGE), ]
  Images <- Images$IMAGE
  print(paste("Table Images=", NROW(Images)))
  remove(Images)
  
  RemoveColumns <-
    c(
      "RELATIVE_PATH",
      "UNIVERSAL_SPOT_ID",
      "STARTING_NORMALIZED_INTENSITY",
      "MAX_NORMALIZED_INTENSITY",
      "START_TO_MAX_INTENSITY",
      "CELL_AREA",
      "ABSOLUTE_POSITION_X",
      "ABSOLUTE_POSITION_Y",
      "NEAREST_SPOT",
      "SPOTS_WITHIN_RADIUS",
      "ANALYSIS_TIME_STAMP",
      "FRAMES_SINCE_LANDING"
    )
  
  ColumnNames <- open_dataset(file.path(TABLE_PATH, "Essential.gz.parquet"))$schema$names
  ColumnNames <- ColumnNames[!grepl(paste(RemoveColumns, collapse = "|"), ColumnNames)]
  remove(RemoveColumns)
  
  Table <- read_parquet(
    file.path(TABLE_PATH, "Essential.gz.parquet"),
    col_select = all_of(ColumnNames)
  )
  
} else{
  Table <- fread(file.path(TABLE_PATH, "Essential.csv.gz"))
}

# Split by cell
SplitTable <-
  Table %>% 
  filter(
    !grepl("half", COHORT)
  ) %>%
  mutate(
    UNIVERSAL_CELL_ID = paste(IMAGE, CELL, sep = "...")
  ) %>% 
  as_tibble() %>% 
  group_split(
    UNIVERSAL_CELL_ID
  )

SplitTable <- SplitTable[order(sapply(SplitTable,nrow))]
SplitTable <- SplitTable[c(NROW(SplitTable):1)]

SplitFx <- function(TableX){
  tryCatch({
    # # Only distinct
    TempTable <- TableX %>% distinct() %>% as.data.table()
    
    # Limit frame number
    FilteredTable <- TempTable[TempTable$FRAMES_ADJUSTED <= MAX_FRAMES]
    
    # How long since it landed
    FilteredTable <- FilteredTable[, FRAMES_SINCE_LANDING := FRAME - min(FilteredTable$FRAME)]
    FilteredTable <- FilteredTable %>% filter(FRAMES_SINCE_LANDING <= FRAMES_SINCE_LANDING_THRESHOLD) %>% as.data.table()
    FilteredTable <- FilteredTable[, FRAMES_SINCE_LANDING_CAT := ceiling(FRAMES_SINCE_LANDING/50)*50]
    
    # Calculate number of frames
    FilteredTable <- FilteredTable[, LIFETIME := max(FRAMES_ADJUSTED, na.rm = T) - min(FRAMES_ADJUSTED, na.rm = T), by = UNIVERSAL_TRACK_ID]
    # Minimum frame count
    MIN_LIFE = LEAD_LAG*2+1
    FilteredTable <- FilteredTable[FilteredTable$LIFETIME >= MIN_LIFE]
    
    return(FilteredTable)
  }, error = function(e) {print("Error SplitFx")})
  
}
FilteredTable <- mclapply(SplitTable, SplitFx, mc.cores = detectCores(logical = F))
FilteredTable <- rbindlist(FilteredTable)

Images <- FilteredTable[!fduplicated(FilteredTable$IMAGE), ]
Images <- Images$IMAGE
print(paste("FilteredTable Images=", NROW(Images)))
remove(Images)

# Select columns
ComplementaryProteins <- names(FilteredTable)
ComplementaryProteins <- ComplementaryProteins[grepl("COMPLEMENTARY_PROTEIN_", ComplementaryProteins)]
ComplementaryProteins <- gsub("COMPLEMENTARY_PROTEIN_", "", ComplementaryProteins)

CorrectQuery <- function(ColumnX){
  
  ProteinName <- paste0("COMPLEMENTARY_PROTEIN_", ColumnX)
  ProteinIntensity <-  paste0("COMPLEMENTARY_NORMALIZED_INTENSITY_", ColumnX)
  
  FilteredTable <-
    FilteredTable %>% 
    filter(
      PROTEIN == USE_REFERENCE_PROTEIN
    ) %>% 
    mutate(
      REFERENCE_PROTEIN = USE_REFERENCE_PROTEIN,
      REFERENCE_TOTAL_INTENSITY = NORMALIZED_INTENSITY,
      DATE = substr(IMAGE, 0, 8),
      QUERY_PROTEIN = (!!as.name(ProteinName)),
      QUERY_TOTAL_INTENSITY = (!!as.name(ProteinIntensity)),
    ) %>% 
    as.data.table()
  
  return(FilteredTable)
}
FilteredTable <- mclapply(ComplementaryProteins, CorrectQuery, mc.cores = detectCores(logical = F))
FilteredTable <- rbindlist(FilteredTable)

Proteins <- FilteredTable[!kit::fduplicated(FilteredTable$QUERY_PROTEIN), ]
Proteins <- unique(Proteins$QUERY_PROTEIN)
Proteins <- Proteins[!Proteins %in% PROTEIN_ORDER]
PROTEIN_ORDER <- c(PROTEIN_ORDER, Proteins)

SplitFilteredTable <-
  FilteredTable %>%
  as_tibble() %>%
  group_split(
    IMAGE
  )

ScaleFx <- function(TableX){
  
  TempTable <- as.data.table(TableX)
  
  TempTable <-
    TempTable %>% 
    drop_na(
      REFERENCE_TOTAL_INTENSITY,
      QUERY_TOTAL_INTENSITY
    ) %>% 
    mutate(
      MIN_REFERENCE = quantile(REFERENCE_TOTAL_INTENSITY, LOWER_BOUNDARY),
      MIN_QUERY = quantile(QUERY_TOTAL_INTENSITY, LOWER_BOUNDARY),
      
      MAX_REFERENCE = quantile(REFERENCE_TOTAL_INTENSITY, UPPER_BOUNDARY),
      MAX_QUERY = quantile(QUERY_TOTAL_INTENSITY, UPPER_BOUNDARY)
    ) %>%
    mutate(
      ORIGINAL_REFERENCE_TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY,
      REFERENCE_TOTAL_INTENSITY = (REFERENCE_TOTAL_INTENSITY - MIN_REFERENCE)/(MAX_REFERENCE - MIN_REFERENCE),
      ORIGINAL_QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY,
      QUERY_TOTAL_INTENSITY = (QUERY_TOTAL_INTENSITY - MIN_QUERY)/(MAX_QUERY - MIN_QUERY)
    ) %>% 
    as.data.table()
  
  return(TempTable)
}
FilteredTable <- mclapply(SplitFilteredTable, ScaleFx, mc.cores = detectCores(logical = F))
FilteredTable <- rbindlist(FilteredTable)

Images <- FilteredTable[!fduplicated(FilteredTable$COHORT), ]
Images <- Images$IMAGE
print(paste("FilteredTable Images=", NROW(Images)))
remove(Images)

# Split table to run faster
SplitTable <-
  FilteredTable %>% 
  mutate(
    UNIVERSAL_CELL_ID = paste(IMAGE, CELL, sep = "...")
  ) %>% 
  ungroup() %>% 
  as_tibble() %>% 
  group_split(
    UNIVERSAL_CELL_ID
  )

# Keep only tibbles
SplitTable <- unclass(SplitTable)
SplitTable <- SplitTable[c(1:NROW(SplitTable))]

# Large to small
SplitTable <- SplitTable[order(sapply(SplitTable,nrow))]
SplitTable <- SplitTable[c(NROW(SplitTable):1)]

remove(FilteredTable)
