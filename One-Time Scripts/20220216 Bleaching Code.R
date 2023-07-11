# Load libraries
pacman::p_load(dtplyr, data.table, dplyr, ggplot2, ggdark)

# Import tables
setwd("/Users/u_deliz/Desktop")
Table <- fread("export.csv")

# Correct number of rows
Table <- Table[4:NROW(Table)]

# Convert to number
Table$TOTAL_INTENSITY_CH1 <- as.numeric(Table$TOTAL_INTENSITY_CH1)


ggplot(
  Table,
  aes(
    TOTAL_INTENSITY_CH1
  )
) +
  geom_density() +
  scale_x_log10()+
  dark_theme_classic()


LifetimeTable <-
  Table %>%
  filter(
    TRACK_ID != ""
  ) %>% 
  mutate(
    FRAME = as.numeric(FRAME)
  ) %>% 
  group_by(
    TRACK_ID
  ) %>% 
  summarize(
    LIFETIME = max(FRAME) - min(FRAME)
  ) %>% 
  as.data.table()


ggplot(
  LifetimeTable,
  aes(
    LIFETIME
  )
) +
  geom_density() +
  scale_x_log10()+
  dark_theme_classic()




BleachingRate <-
  Table %>%
  filter(
    TRACK_ID != ""
  ) %>% 
  mutate(
    FRAME = as.numeric(FRAME)
  ) %>% 
  group_by(
    FRAME
  ) %>% 
  summarize(
    N = n()
  ) %>% 
  mutate(
    N = N/max(N)
  ) %>% 
  as.data.table()

ggplot(
  BleachingRate,
  aes(
    FRAME,
    N/max(N)
  )
) +
  geom_point() +
  scale_y_log10() +
  dark_theme_classic()


HeatmapTable <-
  Table %>%
  filter(
    TRACK_ID != ""
  ) %>% 
  group_by(
    TRACK_ID
  ) %>% 
  mutate(
    FRAME = as.numeric(FRAME),
    LIFETIME = max(FRAME) - min(FRAME)
  ) %>% 
  mutate(
    POSITION_X = as.numeric(POSITION_X),
    POSITION_Y = as.numeric(POSITION_Y),
    
    POSITION_X = round(POSITION_X/88*3)*88/3,
    POSITION_Y = round(POSITION_Y/88*3)*88/3
  ) %>% 
  group_by(
    POSITION_X,
    POSITION_Y
  ) %>% 
  summarize(
    MED_TOTAL_INTENSITY_CH1 = median(TOTAL_INTENSITY_CH1),
    LIFETIME = median(LIFETIME)
  ) %>% 
  mutate(
    MED_TOTAL_INTENSITY_CH1 = MED_TOTAL_INTENSITY_CH1/median(Table$TOTAL_INTENSITY_CH1)
  ) %>% 
  as.data.table()
  
  
ggplot(
  HeatmapTable,
  aes(
    POSITION_X,
    POSITION_Y,
    fill = LIFETIME
  )
) + 
  geom_tile() +
  dark_theme_classic() +
  coord_fixed()



