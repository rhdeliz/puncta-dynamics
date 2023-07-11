pacman::p_load(lemon, ggquiver, ggplot2, ggdark, scales, arrow, R.utils, metR, parallel,
               ggforce, data.table, viridis, RcppRoll, tidyr, dplyr, compiler, dtplyr, kit)

if(grepl("macOS", osVersion)){
  pacman::p_load(metR)
}

filter <- dplyr::filter
select <- dplyr::select

setwd("/Users/u_deliz/Desktop")

# Import tables
OriginalTable <- read_parquet("Essential.gz.parquet",
                      col_select = c("PROTEIN",
                                     "COHORT",
                                     "LIGAND_DENSITY_CAT",
                                     "IMAGE",
                                     "CELL",
                                     "UNIVERSAL_TRACK_ID",
                                     "FRAME",
                                     "FRAMES_ADJUSTED",
                                     "TIME",
                                     "NORMALIZED_INTENSITY",
                                     "COMPLEMENTARY_PROTEIN_1",
                                     "COMPLEMENTARY_NORMALIZED_INTENSITY_1"
                                     ))

Table <-
  OriginalTable %>% 
  filter(
    PROTEIN == "MyD88",
    FRAMES_ADJUSTED <= 199
  ) %>% 
  as_tibble() %>% 
  group_split(
    IMAGE
  )

Table <- Table[order(sapply(Table,nrow))]
Table <- Table[c(NROW(Table):1)]

AnalyzeFx <- function(df){
  
  df <- as.data.table(df)
  
  # Filter out data
  df <-
    df %>% 
    group_by(
      CELL
    ) %>% 
    mutate(
      FRAMES_SINCE_LANDING = FRAME - min(FRAME)
    ) %>% 
    filter(
      FRAMES_SINCE_LANDING <= 399
    ) %>% 
    as.data.table()
  
  return(df)
}
Table <- mclapply(Table, AnalyzeFx, mc.cores = detectCores(logical = F))

FPSFx <- function(df){
  df <-
    df %>% 
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
  return(df)
}
FPSTable <- mclapply(Table, FPSFx, mc.cores = detectCores(logical = F))
FPSTable <- rbindlist(FPSTable)

Table <- rbindlist(Table)

SummaryTable <-
  Table %>% 
  group_by(
    COHORT,
    LIGAND_DENSITY_CAT,
    IMAGE,
    FRAMES_ADJUSTED
  ) %>% 
  summarize(
    NORMALIZED_INTENSITY = median(NORMALIZED_INTENSITY),
    COMPLEMENTARY_NORMALIZED_INTENSITY_1 = median(COMPLEMENTARY_NORMALIZED_INTENSITY_1),
    N = n()
  ) %>% 
  filter(
    N >= 5
  ) %>% 
  as.data.table()

SummaryTable <- merge(SummaryTable, FPSTable, by = "IMAGE")

OldLigands <-
  SummaryTable %>% 
  select(IMAGE, LIGAND_DENSITY_CAT, COHORT, FPS) %>%
  distinct() %>% 
  as.data.table()

NewLigands <- read.csv("new_ligand.csv")
NewLigands <- NewLigands %>% select(IMAGE, LIGAND_DENSITY_CAT) %>% as.data.table()

VerificationNeeded <-
  OldLigands %>% 
  filter(
    !IMAGE %in% NewLigands$IMAGE
  ) %>% 
  as.data.table()

write.csv(VerificationNeeded, "verify_ligands.csv")

VerificationNeeded <-
  VerificationNeeded %>% 
  select(IMAGE, LIGAND_DENSITY_CAT) %>%
  as.data.table()

TestLigands <- rbindlist(list(NewLigands, VerificationNeeded))

SummaryTable$LIGAND_DENSITY_CAT = NULL
NewSummaryTable <- merge(SummaryTable, TestLigands, by = "IMAGE")

NewSummaryTable <-
  NewSummaryTable %>% 
  mutate(
    # LIGAND_DENSITY_CAT = stringr::str_pad(as.numeric(LIGAND_DENSITY_CAT), 3, pad = "0")
  ) %>% 
  arrange(
    COHORT, LIGAND_DENSITY_CAT, FPS, IMAGE
  ) %>% 
  group_by(
    COHORT, LIGAND_DENSITY_CAT, FPS, IMAGE
  ) %>%
  mutate(
    IMAGENUMBER = cur_group_id()
  ) %>%
  group_by(
    COHORT, LIGAND_DENSITY_CAT, FPS
  ) %>%
  mutate(
    IMAGENUMBER = IMAGENUMBER - min(IMAGENUMBER) + 1,
    COHORT = gsub("MyD88 ", "", COHORT)
  ) %>%
  mutate(
    # FACET = paste(COHORT, FPS, LIGAND_DENSITY_CAT, sep = "\n")
    FACET = paste(FPS, COHORT, sep = "\n")
  ) %>% 
  group_by(
    COHORT, FPS
  ) %>% 
  mutate(
    # NORMALIZED_INTENSITY = NORMALIZED_INTENSITY/max(NORMALIZED_INTENSITY)*100
  ) %>% 
  as.data.table()

NewSummaryTable %>% 
  mutate(
    TEST = !grepl("grid", COHORT, ignore.case = T)
  ) %>% 
  filter(
    TEST == T
  ) %>% 
  mutate(
    TEST = !grepl("KO", COHORT, ignore.case = T)
  ) %>% 
  filter(
    TEST == T,
    FPS == 0.25
    # COHORT == "TRAF6",
    # LIGAND_DENSITY_CAT == 10
  ) %>% 
  group_by(
    LIGAND_DENSITY_CAT, FPS, COHORT,
    IMAGENUMBER
  ) %>% 
  mutate(
    # NORMALIZED_INTENSITY = (NORMALIZED_INTENSITY-min(NORMALIZED_INTENSITY))/(max(NORMALIZED_INTENSITY)-min(NORMALIZED_INTENSITY))
  ) %>%
  as.data.table() %>% 
  ggplot(
    aes(
      # FRAMES_ADJUSTED,
      NORMALIZED_INTENSITY,
      COMPLEMENTARY_NORMALIZED_INTENSITY_1,
      color = as.factor(LIGAND_DENSITY_CAT),
      # color = as.factor(IMAGENUMBER),
      group = paste(LIGAND_DENSITY_CAT, IMAGENUMBER, COHORT)
      # group = paste(LIGAND_DENSITY_CAT, IMAGENUMBER)
      # group = IMAGENUMBER
    )
  ) +
  geom_path() +
  # scale_y_log10()+
  scale_color_viridis(
    option = "turbo",
    discrete = T
  ) +
  facet_wrap(
    ~FACET,
    scales = "free"
  ) +
  theme_classic()


Table <- read_parquet("Essential.gz.parquet")
write_parquet(Table, "BackupEssential.gz.parquet")

NewLigands <- read.csv("new_ligand.csv")

Table$LIGAND_DENSITY_CAT = NULL
NewLigands <- NewLigands %>% select(IMAGE, LIGAND_DENSITY_CAT) %>% as.data.table()
Table <- data.table::merge.data.table(Table, NewLigands, by = "IMAGE")
Table <- data.table::merge.data.table(Table, FPSTable, by = "IMAGE")
write_parquet(Table, "Essential.gz.parquet")
