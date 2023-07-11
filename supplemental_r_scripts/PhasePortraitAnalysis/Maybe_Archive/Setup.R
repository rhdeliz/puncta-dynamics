setwd(OUTPUT_DIRECTORY)

# Import compressed table
# Change to Essentials once fixed
if(grepl("macOS", osVersion)){
  # Table <- fread(file.path(TABLE_PATH, "Grid.csv.gz"), fill=TRUE)
  Table <- fread(file.path(TABLE_PATH, "Essential.csv.gz"), fill=TRUE)
} else{
  library(R.utils)
  if(!file.exists(file.path(TABLE_PATH, "Essential.csv"))){
    gunzip(file.path(TABLE_PATH, "Essential.csv.gz"), remove = FALSE)
  }
  Table <- fread(file.path(TABLE_PATH, "Essential.csv"), fill=TRUE)
}

print(paste("Table Images=", NROW(unique(Table$IMAGE))))

# Fix table data
SplitTable <-
  Table %>% 
  as_tibble() %>% 
  group_split(
    IMAGE
  )

SplitTable <- SplitTable[order(sapply(SplitTable,nrow))]
SplitTable <- SplitTable[c(1:NROW(SplitTable))]

SplitFx <- function(TableX){
  # Only distinct
  TempTable <- TableX %>% distinct() %>% as.data.table()
  
  # Estimate frame rate using track with longest lifetime
  SelecTracks <-
    TempTable %>% 
    group_by(
      COHORT,
      LIGAND_DENSITY_CAT,
      IMAGE
    ) %>% 
    arrange(
      -LIFETIME
    ) %>% 
    filter(
      row_number()==1 
    ) %>% 
    ungroup() %>% 
    select(UNIVERSAL_TRACK_ID) %>% 
    as.data.table()
  
  # Select only track ID columns
  SelecTracks <- SelecTracks$UNIVERSAL_TRACK_ID
  
  # Calculate FPS
  FrameRateTable <-
    TempTable %>%
    filter(
      UNIVERSAL_TRACK_ID %in% SelecTracks
    ) %>% 
    arrange(
      FRAMES_ADJUSTED
    ) %>% 
    select(
      IMAGE,
      UNIVERSAL_TRACK_ID,
      TIME
    ) %>% 
    group_by(
      IMAGE
    ) %>% 
    distinct() %>% 
    mutate(
      FPS = TIME-lag(TIME)
    ) %>% 
    summarize(
      FPS = mean(FPS, na.rm = T)
    ) %>% 
    mutate(
      FPS = 1/round(FPS)
    ) %>% 
    mutate(
      FPS = ifelse(FPS == 0, ceiling(1/FPS), FPS)
    ) %>% 
    as.data.table()
  
  FilteredTable <- merge(TempTable, FrameRateTable, by = "IMAGE")
  
  # Limit frame number
  FilteredTable <- FilteredTable[FilteredTable$FRAMES_ADJUSTED <= MAX_FRAMES]
  # Calculate number of frames
  FilteredTable <- FilteredTable[, LIFETIME := max(FRAMES_ADJUSTED, na.rm = T) - min(FRAMES_ADJUSTED, na.rm = T), by = UNIVERSAL_TRACK_ID]
  # Minimum frame count
  MIN_LIFE = LEAD_LAG*2+1
  FilteredTable <- FilteredTable[FilteredTable$LIFETIME >= MIN_LIFE]
  
  return(FilteredTable)
}
FilteredTable <- mclapply(SplitTable, SplitFx, mc.cores = detectCores(logical = F))
FilteredTable <- rbindlist(FilteredTable)

print(paste("FilteredTable Images=", NROW(unique(FilteredTable$IMAGE))))

# Select columns
ComplementaryProteins <- names(FilteredTable)
ComplementaryProteins <- ComplementaryProteins[grepl("COMPLEMENTARY_PROTEIN_", ComplementaryProteins)]
ComplementaryProteins <- gsub("COMPLEMENTARY_PROTEIN_", "", ComplementaryProteins)

CorrectQuery <- function(ColumnX){
  
  ProteinName <- paste0("COMPLEMENTARY_PROTEIN_", ColumnX)
  ProteinIntensity <-  paste0("COMPLEMENTARY_NORMALIZED_INTENSITY_", ColumnX)
  
  ReferenceTable <-
    FilteredTable %>% 
    filter(
      PROTEIN == USE_REFERENCE_PROTEIN
    ) %>% 
    mutate(
      SOURCE = PROTEIN,
      REFERENCE_PROTEIN = USE_REFERENCE_PROTEIN,
      REFERENCE_TOTAL_INTENSITY = NORMALIZED_INTENSITY,
      DATE = substr(IMAGE, 0, 8),
      QUERY_PROTEIN = (!!as.name(ProteinName)),
      QUERY_TOTAL_INTENSITY = (!!as.name(ProteinIntensity)),
    ) %>% 
    filter(
      QUERY_PROTEIN %in% PROTEIN_ORDER
    ) %>% 
    as.data.table()
  
  QueryTable <-
    FilteredTable %>% 
    filter(
      PROTEIN != USE_REFERENCE_PROTEIN
    ) %>% 
    mutate(
      SOURCE = PROTEIN,
      REFERENCE_PROTEIN = USE_REFERENCE_PROTEIN,
      REFERENCE_TOTAL_INTENSITY = (!!as.name(ProteinIntensity)),
      DATE = substr(IMAGE, 0, 8),
      QUERY_PROTEIN = PROTEIN,
      QUERY_TOTAL_INTENSITY = NORMALIZED_INTENSITY,
    ) %>% 
    filter(
      QUERY_PROTEIN %in% PROTEIN_ORDER
    ) %>% 
    as.data.table()
  
  TempFilteredTable <- rbindlist(list(ReferenceTable, QueryTable))
  
  return(TempFilteredTable)
}
FilteredTable <- mclapply(ComplementaryProteins, CorrectQuery, mc.cores = detectCores(logical = F))
FilteredTable <- rbindlist(FilteredTable)
# remove(Table)

print(paste("FilteredTable Images=", NROW(unique(FilteredTable$IMAGE))))

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

# Small to large
SplitTable <- SplitTable[order(sapply(SplitTable,nrow))]
SplitTable <- SplitTable[c(1:NROW(SplitTable))]

remove(FilteredTable)
