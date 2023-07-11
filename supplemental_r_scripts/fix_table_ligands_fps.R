pacman::p_load(data.table, tidyr, dplyr, compiler, dtplyr, arrow)

Table <- read_parquet("Backup_Essential.gz.parquet")

# Add new ligand category ----
Table$LIGAND_DENSITY_CAT = NULL
NewLigands <- read.csv("new_ligand.csv")
NewLigands <- NewLigands %>% select(IMAGE, LIGAND_DENSITY_CAT) %>% as.data.table()
Table <- merge(Table, NewLigands, by = "IMAGE")

# Calculate frames per second ----
FPS_Table <-
  Table %>% 
  select(
    IMAGE,
    FRAME,
    TIME
  ) %>% 
  distinct() %>% 
  group_by(
    IMAGE
  ) %>% 
  mutate(
    FPS = TIME - lag(TIME)
  ) %>% 
  drop_na() %>% 
  mutate(
    FPS = round(FPS)
  ) %>% 
  summarize(
    FPS = median(FPS)
  ) %>% 
  mutate(
    FPS = 1/round(FPS)
  ) %>% 
  as.data.table()

Table$FPS = NULL
Table <- merge(Table, FPS_Table, by = "IMAGE")

write_parquet(Table, "Essential.gz.parquet")
