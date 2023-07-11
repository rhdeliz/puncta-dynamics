GFP <- fread("/Users/u_deliz/Desktop/Grid Data/GFP_intensity.csv.gz")
RFP <- fread("/Users/u_deliz/Desktop/Grid Data/mScarlet_intensity.csv.gz")

Paths <- c(
  "/Users/u_deliz/Desktop/Grid Data/GFP_intensity.csv.gz",
  "/Users/u_deliz/Desktop/Grid Data/mScarlet_intensity.csv.gz"
)

FluorophoreFx <- function(PathX){
  TempTable <- fread(PathX)
  TempTable <-
    TempTable %>% 
    filter(
      SPOT_AREA == 25
    ) %>% 
    mutate(
      DATE = substr(IMAGE, 1, 8)
    ) %>% 
    group_by(
      PROTEIN,
      DATE
    ) %>% 
    summarize(
      TOTAL_INTENSITY = median(TOTAL_INTENSITY, na.rm = T)
    ) %>% 
    as.data.table()
  return(TempTable)
}
Value <- lapply(Paths, FluorophoreFx)
Value <- rbindlist(Value)

GridTable <- fread("/Users/u_deliz/Desktop/Grid Data/Essential.csv.gz")

GridTable1 <-
  GridTable %>% 
  filter(
    PROTEIN == "MyD88"
  ) %>% 
  mutate(
    NORMALIZED_INTENSITY = NORMALIZED_INTENSITY/182.9642,
    COMPLEMENTARY_NORMALIZED_INTENSITY_1 = COMPLEMENTARY_NORMALIZED_INTENSITY_1/162.4642
  ) %>% 
  as.data.table()

GridTable2 <-
  GridTable %>% 
  filter(
    PROTEIN != "MyD88"
  ) %>% 
  mutate(
    NORMALIZED_INTENSITY = NORMALIZED_INTENSITY/162.4642,
    COMPLEMENTARY_NORMALIZED_INTENSITY_1 = COMPLEMENTARY_NORMALIZED_INTENSITY_1/182.9642
  ) %>% 
  as.data.table()

GridTable <- rbindlist(list(GridTable1, GridTable2))

fwrite(GridTable, "Grid.csv.gz")
