Table <- fread("/Users/u_deliz/Desktop/Essential.csv.gz")
BIN_SIZE = 8

ImgTable <-
  Table %>%
  group_by(
    CELL
  ) %>% 
  mutate(
    FRAMES_SINCE_LANDING = FRAME - min(FRAME)
  ) %>% 
  filter(
    FRAMES_ADJUSTED == 0,
    FRAMES_SINCE_LANDING  <= 200
  ) %>% 
  group_by(
    COHORT,
    IMAGE
  ) %>% 
  summarize(
    N = n()
  ) %>% 
  arrange(
    -N
  ) %>% 
  as.data.table()

TempTable <-
  Table %>% 
  filter(
    IMAGE == '20210518 0.5nM 120-1E1_HOIL1_MyD88 grid_2.5 IL1 002',
    # CELL == 34
  ) %>% 
  group_by(
    CELL
  ) %>% 
  mutate(
    FRAMES_SINCE_LANDING = FRAME - min(FRAME),
    LIFETIME_FRAMES = max(FRAMES_ADJUSTED)
  ) %>% 
  filter(
    FRAMES_SINCE_LANDING  <= 300,
    LIFETIME_FRAMES >= 11
  ) %>% 
  group_by(
    CELL,
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    ROUNDED_POSITION_X = round(ABSOLUTE_POSITION_X/BIN_SIZE)*BIN_SIZE,
    ROUNDED_POSITION_Y = round(ABSOLUTE_POSITION_Y/BIN_SIZE)*BIN_SIZE,
    
    DELTA_X = ABSOLUTE_POSITION_X - lag(ABSOLUTE_POSITION_X),
    DELTA_Y = ABSOLUTE_POSITION_Y - lag(ABSOLUTE_POSITION_Y),
    
    TIME = TIME - lag(TIME),
    DELTA_X = DELTA_X/TIME,
    DELTA_Y = DELTA_Y/TIME
  ) %>% 
  as_tibble() %>% 
  group_by(
    IMAGE,
    ROUNDED_POSITION_X,
    ROUNDED_POSITION_Y
  ) %>% 
  summarize(
    DELTA_X = median(DELTA_X, na.rm = T)*1000,
    DELTA_Y = median(DELTA_Y, na.rm = T)*1000,
    N = n()
    # CELL = median(CELL)
  ) %>% 
  filter(
    N >= 5
  ) %>% 
  mutate(
    SPEED = sqrt(DELTA_Y^2 + DELTA_X^2)
  ) %>% 
  group_by(
    IMAGE
  ) %>% 
  filter(
    SPEED >= quantile(SPEED, .1, na.rm = T) &
      SPEED <= quantile(SPEED, .9, na.rm = T)
  ) %>% 
  drop_na() %>%
  mutate(
    # Get angle to account for negative magnitudes
    COLOR = atan2(-DELTA_X, -DELTA_Y),
    # Make angle
    COLOR = COLOR*180/pi+180,
    COLOR = COLOR/360,
    COLOR = hsv(COLOR, 1, 1)
  ) %>% 
  as.data.table()


StreamTempTable <-
  TempTable %>%
  as_tibble() %>% 
  complete(
    ROUNDED_POSITION_X = full_seq(c(min(ROUNDED_POSITION_X), max(ROUNDED_POSITION_X)), period = BIN_SIZE),
    ROUNDED_POSITION_Y = full_seq(c(min(ROUNDED_POSITION_Y), max(ROUNDED_POSITION_Y)), period = BIN_SIZE)
  ) %>%
  select(ROUNDED_POSITION_X, ROUNDED_POSITION_Y, DELTA_X, DELTA_Y) %>%
  distinct() %>% 
  as.data.table()

# Add derivatives
StreamTempTable$DELTA_X[is.na(StreamTempTable$DELTA_X)] = 0
StreamTempTable$DELTA_Y[is.na(StreamTempTable$DELTA_Y)] = 0




ggplot(
) +
  geom_streamline(
    data = StreamTempTable,
    aes(
      x = ROUNDED_POSITION_X,
      y = ROUNDED_POSITION_Y,
      dx = DELTA_X,
      dy = DELTA_Y,
      color = sqrt(..dx..^2 + ..dy..^2),
      size = ..step..
    ),
    arrow = NULL,
    # n = 25,
    L = 3, res = 20, lineend = "round"
  ) +
  geom_arrow(
    data = TempTable,
    aes(
      x = ROUNDED_POSITION_X,
      y = ROUNDED_POSITION_Y,
      mag = 0.5,
      angle = atan2(DELTA_Y, DELTA_X)*180/pi,
      # color = sqrt(DELTA_Y^2 + DELTA_X^2)
      color = COLOR
    ),
    size = .01, # Small arrow head
    arrow.length = .5, # Small arrow head
    # size = 1, # Big arrow head
    # arrow.length = 1,# Big arrow head
    lineend = "square"
  ) +
  scale_size(range = c(.01, 1), guide = "none") +
  scale_alpha(guide = "none") +
  scale_mag(
    # max = 3, # Small arrow head
    max = 1, # Big arrow head
    guide = 'none'
  ) +
  scale_color_viridis(
    option = "viridis"
  ) +
  scale_y_reverse()+
  labs(
    x = "um",
    y = "um",
    color = "Speed (nm/s)"
  ) +
  dark_theme_classic(
    base_size = 20
  ) +
  theme(
    legend.key.height= unit(1.25, 'cm'),
    legend.position = "right"
  ) +
  coord_fixed()




ggplot(
) +
  # geom_streamline(
  #   data = StreamTempTable,
  #   aes(
  #     x = ROUNDED_POSITION_X,
  #     y = ROUNDED_POSITION_Y,
  #     dx = DELTA_X,
  #     dy = DELTA_Y,
  #     color = sqrt(..dx..^2 + ..dy..^2),
  #     size = ..step..
  #   ),
  #   arrow = NULL,
#   # n = 25,
#   L = 3, res = 20, lineend = "round"
# ) +
geom_tile(
  data = TempTable,
  aes(
    x = ROUNDED_POSITION_X,
    y = ROUNDED_POSITION_Y,
    fill = N
  )
) +
  geom_arrow(
    data = TempTable,
    aes(
      x = ROUNDED_POSITION_X,
      y = ROUNDED_POSITION_Y,
      mag = 0.5,
      angle = atan2(DELTA_Y, DELTA_X)*180/pi,
      color = sqrt(DELTA_Y^2 + DELTA_X^2)
      # color = N
    ),
    size = .01, # Small arrow head
    arrow.length = .5, # Small arrow head
    # size = 1, # Big arrow head
    # arrow.length = 1,# Big arrow head
    lineend = "square"
  ) +
  scale_fill_distiller(
    palette = "Greys"
  ) + 
  scale_size(range = c(.01, 1), guide = "none") +
  scale_alpha(guide = "none") +
  scale_mag(
    # max = 3, # Small arrow head
    max = 1, # Big arrow head
    guide = 'none'
  ) +
  scale_color_viridis(
    option = "viridis"
  ) +
  scale_y_reverse()+
  labs(
    x = "um",
    y = "um",
    color = "Speed (nm/s)"
  ) +
  dark_theme_classic(
    base_size = 20
  ) +
  theme(
    legend.key.height= unit(1.25, 'cm'),
    legend.position = "right"
  ) +
  coord_fixed()


