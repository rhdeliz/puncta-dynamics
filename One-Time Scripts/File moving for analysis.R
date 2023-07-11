
# Ending of files
parameters_path = "/raven/u/deliz/new_pipeline/pending_processing/20200502/Input/parameter_tables"
intensity_ending = "_intensity.csv.gz"
colocalization_intensity_file = "colocalization.csv.gz"

if("pacman" %in% rownames(installed.packages()) == FALSE)
{install.packages("pacman")}

pacman::p_load(dplyr, stringr, parallel, tidyr, data.table, ff, dtplyr, compiler)
setDTthreads(parallel::detectCores(logical = F))
enableJIT(3)

# Get directories
directories_list = file.path(parameters_path, "directories.csv")
directories_list = fread(directories_list)
input_path = directories_list$path[directories_list$contains == "input"]
processing_path = directories_list$path[directories_list$contains == "processing"]
colocalization_path = file.path(processing_path, "06_Colocalization")
output_path = directories_list$path[directories_list$contains == "output"]
summary_path = file.path(input_path, "summary.csv")
# Get image list
file_list = fread(summary_path)

# Convert exposure units
ConvertUnits <- function(x){
  # Create array
  y = rep(NA, NROW(x))
  # Loop over x
  for(i in 1:NROW(x)){
    if(x[i] == "ms"){
      y[i] = 10^-3
    } else if(x[i] == "s"){
      y[i] = 1
    } else{
      y[i] = 1
    }
  }
  return(y)
}
# Prep up
file_list <-
  file_list %>% 
  mutate(
    image = dirname(dirname(protein_relative_path)),
    exposure_units = word(exposure, -1),
    exposure_measurement = as.numeric(word(exposure, 1)),
    direction = direction*pi/180,
    angle = angle*pi/180
  ) %>% 
  mutate(
    exposure = exposure_measurement*ConvertUnits(exposure_units),
    image = basename(image),
    date = substr(image, 0, 8)
  ) %>% 
  mutate(
    date =  as.Date(date, format = '%Y%m%d'),
  ) %>% 
  select(-c(
    exposure_units, exposure_measurement
  )) %>% 
  as.data.table()

# Separate calibrations and the rest
calibration_list = file_list[file_list$cohort=="Calibrations",]
Path <- calibration_list$protein_relative_path
Path <- dirname(Path)
Path <- file.path(colocalization_path, Path)
Path <- file.exists(Path)
calibration_list <- calibration_list[Path]

# Analyze calibrations
QuantifyCalibration <- function(FileX){
  # Get calibration image metadata
  temp_calibration_list = calibration_list[FileX]
  # Pull table
  Table <- file.path(colocalization_path, paste0(temp_calibration_list$protein_relative_path, intensity_ending))
  Table <- fread(Table)
  
  Table <-
    Table %>% 
    summarize(
      SD = mad(TOTAL_INTENSITY),
      SEM = mad(TOTAL_INTENSITY)/sqrt(n()),
      TOTAL_INTENSITY = median(TOTAL_INTENSITY),
      N = n()
    ) %>% 
    as.data.table()
  
  temp_calibration_list <- bind_cols(temp_calibration_list, Table)
  
  save_path = paste0(temp_calibration_list$protein_name, "_calibration.csv.gz")
  save_path = file.path(colocalization_path, dirname(dirname(temp_calibration_list$protein_relative_path)), save_path)
  # Remove old if it exists
  suppressWarnings({
    file.remove(save_path, showWarnings = FALSE)
  })
  data.table::fwrite(temp_calibration_list, save_path, row.names = F, na = "")
  
  return(temp_calibration_list)
}
calibration_list <- mclapply(1:NROW(calibration_list), QuantifyCalibration, mc.cores = detectCores(logical = F))
calibration_list <- rbindlist(calibration_list)

image_list = file_list[file_list$cohort!="Calibrations",]
Path <- image_list$protein_relative_path
Path <- dirname(Path)
Path <- file.path(colocalization_path, Path)
Path <- file.exists(Path)
image_list <- image_list[Path]

# Pairing
GetCalibrationImages <- function(ImageX){
  tryCatch({
    
    # Get image data
    temp_image_list = image_list[ImageX]
    temp_image_list$cell_path = dirname(temp_image_list$protein_relative_path)
    
    temp_calibration_list <-
      calibration_list %>% 
      # Get the correct channel
      filter(
        channel == temp_image_list$channel
      ) %>% 
      as.data.table()
    
    if(NROW(temp_calibration_list)>0){
      
      temp_calibration_list <-
        temp_calibration_list %>% 
        # Get the nearest power
        mutate(
          power_test = abs(power - temp_image_list$power)
        ) %>% 
        filter(
          power_test == min(power_test)
        ) %>% 
        # Get the nearest exposure
        mutate(
          exposure_test = abs(exposure - temp_image_list$exposure)
        ) %>% 
        filter(
          exposure_test == min(exposure_test)
        ) %>% 
        # Get the nearest direction
        mutate(
          direction_test = abs(direction - temp_image_list$direction)
        ) %>% 
        filter(
          direction_test == min(direction_test)
        ) %>% 
        # Get the nearest angle
        mutate(
          angle_test = abs(angle - temp_image_list$angle)
        ) %>% 
        filter(
          direction_test == min(direction_test)
        ) %>% 
        # Get the nearest date
        mutate(
          date_test = abs(date - temp_image_list$date)
        ) %>% 
        filter(
          date_test == min(date_test)
        ) %>% 
        # Lowest SEM
        filter(
          SEM == min(SEM)
        ) %>% 
        as.data.table()
      
      CalibrationImage = temp_calibration_list$image[1]
      CalibrationProtein = temp_calibration_list$protein_name[1]
      CalibrationIntensity = temp_calibration_list$TOTAL_INTENSITY[1]
      CalibrationError = temp_calibration_list$SD[1]
    } else{
      CalibrationImage = NA
      CalibrationProtein = temp_image_list$channel
      CalibrationIntensity = 1
      CalibrationError = NA
    }
    # Create export table
    ExportTable <- temp_image_list
    names(ExportTable) <- toupper(names(ExportTable))
    ExportTable$CALIBRATION_IMAGE = CalibrationImage
    ExportTable$FLUOROPHORE = CalibrationProtein
    ExportTable$CALIBRATION_TOTAL_INTENSITY = CalibrationIntensity
    ExportTable$CALIBRATION_STANDARD_DEVIATION = CalibrationError
    
    save_path = paste0(CalibrationProtein, "_calibration.csv.gz")
    save_path = file.path(colocalization_path, dirname(dirname(temp_image_list$protein_relative_path)), save_path)
    # Remove old if it exists
    suppressWarnings({
      file.remove(save_path, showWarnings = FALSE)
    })
    data.table::fwrite(temp_calibration_list, save_path, row.names = F, na = "")
    
    return(ExportTable)
  }, error = function(e){print(paste("ERROR with MoveNoColocalizationNeededFx ImageX =", ImageX))})
}
PairedList <- mclapply(1:NROW(image_list), GetCalibrationImages, mc.cores = detectCores(logical = F))
PairedList <- rbindlist(PairedList, fill = TRUE, use.names = TRUE)
PairedList <- PairedList %>% distinct() %>% as.data.table()

ANALYSIS_TIME_STAMP = Sys.time()
ANALYSIS_TIME_STAMP = as.POSIXlt(ANALYSIS_TIME_STAMP, tz = "UTC")
ANALYSIS_TIME_STAMP = as.character(ANALYSIS_TIME_STAMP)

Cells <- unique(PairedList$CELL_PATH)


CombineCellTables <- function(CellX){
  tryCatch({
    
    print(paste("CombineCellTables - CellX =", CellX))
    
    # Get tables
    CellPairedCalibrations <-
      PairedList %>%
      filter(
        CELL_PATH == Cells[CellX]
      ) %>% 
      mutate(
        PROTEIN = PROTEIN_NAME 
      ) %>% 
      select(
        PROTEIN,
        CALIBRATION_IMAGE,
        CALIBRATION_TOTAL_INTENSITY,
        CALIBRATION_STANDARD_DEVIATION
      ) %>% 
      as.data.table()
    
    # Get cell path
    CellPath <- file.path(colocalization_path, Cells[CellX])
    
    return(CellPath)
  }, error = function(e){print(paste("ERROR with CombineCellTables CellX =", CellX))})
}
CellAnalysis <- mclapply(1:NROW(Cells), CombineCellTables, mc.cores = detectCores(logical = F))


CellAnalysis <- unlist(CellAnalysis)
CellAnalysis <- CellAnalysis[file.exists(CellAnalysis)]

# Combine all cell tables
Images <- dirname(CellAnalysis)
Images <- unique(Images)

# Tables to merge
Outputs <- c("Parameters.csv.gz", "Analysis.csv.gz", "Essential.csv.gz")

nImages <- NROW(Images)
CombineImageTables <- function(ImageX){
  tryCatch({
    print(paste("ImageX =", ImageX, "of", nImages))
    # Get cell table paths
    CellsList <- CellAnalysis[dirname(CellAnalysis) == Images[ImageX]]
    # Function to combine tables
    TableByOutput <- function(OutputX){
      print(paste("Combining",basename(Images[ImageX]), "-", OutputX))
      # Get path
      Path <- file.path(CellsList, OutputX)
      Path <- Path[file.exists(Path)]
      # Pull tables
      Table <- lapply(Path, fread)
      Table <- rbindlist(Table, fill = TRUE, use.names = TRUE)
      # Save
      DestinationPath <- file.path(Images[ImageX], OutputX)
      file.remove(DestinationPath, showWarnings = F)
      fwrite(Table, DestinationPath, row.names = F, na = "")
      if(OutputX =="Essential.csv.gz"){
        return(DestinationPath)
      }
    }
    CellsList <- lapply(Outputs, TableByOutput)
    # CellsList <- rbindlist(CellsList)
    # # Make image summary for future analysis
    # CellSummary <-
    #   CellsList %>%
    #   filter(
    #     FRAMES_ADJUSTED == 0
    #   ) %>%
    #   group_by(
    #     LIGAND_DENSITY_CAT,
    #     COHORT,
    #     IMAGE,
    #     CELL,
    #     PROTEIN,
    #     CELL_AREA
    #   ) %>%
    #   summarize(
    #     SPOTS = n(),
    #     LIFETIME = mean(LIFETIME, na.rm = T),
    #     STARTING_NORMALIZED_INTENSITY = mean(STARTING_NORMALIZED_INTENSITY, na.rm = T),
    #     MAX_NORMALIZED_INTENSITY = mean(MAX_NORMALIZED_INTENSITY, na.rm = T),
    #     START_TO_MAX_INTENSITY = mean(START_TO_MAX_INTENSITY, na.rm = T),
    #     DURATION = max(TIME) - min(TIME)
    #   ) %>% 
    #   mutate(
    #     SPOTS_PER_AREA_PER_TIME_E6 = SPOTS/CELL_AREA/DURATION*10^6
    #   ) %>% 
    #   as_tibble()
    # 
    DestinationPath <- file.path(Images[ImageX], "CellSummary.csv.gz")
    # file.remove(DestinationPath, showWarnings = F)
    # fwrite(CellSummary, DestinationPath, row.names = F, na = "")
    # 
    # ProteinSummary <-
    #   CellSummary %>%
    #   group_by(
    #     LIGAND_DENSITY_CAT,
    #     COHORT,
    #     IMAGE,
    #     PROTEIN
    #   ) %>%
    #   summarize(
    #     CELLS = n(),
    #     SPOTS = sum(SPOTS),
    #     LIFETIME = mean(LIFETIME, na.rm = T),
    #     STARTING_NORMALIZED_INTENSITY = mean(STARTING_NORMALIZED_INTENSITY, na.rm = T),
    #     MAX_NORMALIZED_INTENSITY = mean(MAX_NORMALIZED_INTENSITY, na.rm = T),
    #     START_TO_MAX_INTENSITY = mean(START_TO_MAX_INTENSITY, na.rm = T)
    #   ) %>% 
    #   as.data.table()
    # 
    # DestinationPath <- file.path(Images[ImageX], "ProteinSummary.csv.gz")
    # file.remove(DestinationPath, showWarnings = F)
    # fwrite(ProteinSummary, DestinationPath, row.names = F, na = "")
    # 
    Progress = NROW(Images)/10
    Progress = round(Progress)
    if(Progress==0){
      Progress = 1
    }
    
    if(ImageX %% Progress == 0){
      Progress = ImageX/NROW(Images)
      Progress = Progress*100
      Progress = round(Progress)
      Progress = paste0("     ", Progress, "% complete")
      print(Progress)
    }
    
    return(Images[ImageX])
  }, error = function(e){print(paste("ERROR with CombineImageTables ImageX =", ImageX))})
}
ProcessedImages <- lapply(1:nImages, CombineImageTables)
ProcessedImages <- unlist(ProcessedImages)


