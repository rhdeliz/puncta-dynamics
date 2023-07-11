library(zoo)
library(trajr)

paths <- c(
  "/Users/u_deliz/Desktop/ArchivePhasePortraitAnalysis/December 2021/Archive/Four Samples Essential.csv.gz"
  # "/Users/u_deliz/Desktop/Essential.csv.gz"
)
Table <- mclapply(paths, fread)
Table <- rbindlist(Table, fill = T)

SplitTable <-
  Table %>%
  filter(
    PROTEIN == "MyD88"
  ) %>%
  arrange(
    COHORT,
    LIGAND_DENSITY_CAT,
    PROTEIN,
    UNIVERSAL_TRACK_ID,
    FRAMES_ADJUSTED
  ) %>% 
  group_by(
    IMAGE,
    CELL
  ) %>% 
  mutate(
    FRAMES_SINCE_LANDING = FRAME - min(FRAME)
  ) %>% 
  filter(
    FRAMES_SINCE_LANDING <= 200,
    FRAMES_ADJUSTED <= 100
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    N = n()
  ) %>% 
  # filter(
  #   N >= 50
  # ) %>% 
  arrange(
    N
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    new_id = cur_group_id()
  ) %>% 
  ungroup() %>% 
  mutate(
    batch_block = round(detectCores()*5/new_id)/detectCores()*5
  ) %>% 
  as_tibble() %>% 
  group_split(
    batch_block
  )

SmoothFx <- function(table){
  
  Tracks <-
    table %>% 
    group_by(
      UNIVERSAL_TRACK_ID
    ) %>% 
    # mutate(
    #   ABSOLUTE_POSITION_X = signal::sgolayfilt(ABSOLUTE_POSITION_X, p = 1, n = 7),
    #   ABSOLUTE_POSITION_Y = signal::sgolayfilt(ABSOLUTE_POSITION_Y, p = 1, n = 7)
    # ) %>%
    # mutate(
    #   ABSOLUTE_POSITION_X = signal::sgolayfilt(ABSOLUTE_POSITION_X, p = 3, n = 5),
    #   ABSOLUTE_POSITION_Y = signal::sgolayfilt(ABSOLUTE_POSITION_Y, p = 3, n = 5)
    # ) %>%
    drop_na(
      ABSOLUTE_POSITION_X,
      ABSOLUTE_POSITION_Y
    ) %>%
    group_by(
      UNIVERSAL_TRACK_ID
    ) %>% 
    mutate(
      ABSOLUTE_POSITION_X = ABSOLUTE_POSITION_X - ABSOLUTE_POSITION_X[n()],
      ABSOLUTE_POSITION_Y = ABSOLUTE_POSITION_Y - ABSOLUTE_POSITION_Y[n()],
    ) %>% 
    mutate(
      # Check if it needs flipping
      AVG_ABSOLUTE_POSITION_X = mean(ABSOLUTE_POSITION_X),
      AVG_ABSOLUTE_POSITION_Y = mean(ABSOLUTE_POSITION_Y)
    ) %>% 
    mutate(
      # Flip if necessary
      ABSOLUTE_POSITION_X = ifelse(AVG_ABSOLUTE_POSITION_X < 0, ABSOLUTE_POSITION_X*-1, ABSOLUTE_POSITION_X),
      ABSOLUTE_POSITION_Y = ifelse(AVG_ABSOLUTE_POSITION_Y < 0, ABSOLUTE_POSITION_Y*-1, ABSOLUTE_POSITION_Y)
    ) %>% 
    mutate(
      # START_X = ABSOLUTE_POSITION_X[1],
      # START_Y = ABSOLUTE_POSITION_Y[1],
      # END_X = ABSOLUTE_POSITION_X[n()],
      # END_Y = ABSOLUTE_POSITION_Y[n()],
      
      START_X = ABSOLUTE_POSITION_X[n()-10],
      START_Y = ABSOLUTE_POSITION_Y[n()-10],
      END_X = ABSOLUTE_POSITION_X[n()],
      END_Y = ABSOLUTE_POSITION_Y[n()]
    ) %>% 
    mutate(
      DELTA_X = END_X - START_X, 
      DELTA_Y = END_Y - START_Y
    ) %>% 
    mutate(
      RAD_ANGLE = atan2(-DELTA_Y, DELTA_X) + pi
    ) %>% 
    as_tibble() %>% 
    group_split(
      UNIVERSAL_TRACK_ID
    )
  
  AlignFx <- function(track){  
    df <- NULL
    df$x <- track$ABSOLUTE_POSITION_X
    df$y <- track$ABSOLUTE_POSITION_Y
    df$id <- track$UNIVERSAL_TRACK_ID
    df <- as.data.table(df)
    
    theta <- as.numeric(track$RAD_ANGLE[1])
    
    in_traj <- TrajFromCoords(df)
    out_traj <- TrajRotate(in_traj, angle = theta, relative = F)
    out_traj <- out_traj %>% select(x, y) %>% as.data.table()
    
    track$ABSOLUTE_POSITION_X <- out_traj$x
    track$ABSOLUTE_POSITION_Y <- out_traj$y

    track <-
      track %>% 
      as.data.table() %>% 
      mutate(
        t = 1:n()
      ) %>% 
      mutate(
        t_end = n() - t
      ) %>% 
      as.data.table()
    
    return(track)
  }
  AlignedTrack <- lapply(Tracks, AlignFx)
  AlignedTrack <- rbindlist(AlignedTrack)
  
  return(AlignedTrack)
}
SmoothTrack <- mclapply(SplitTable, SmoothFx, mc.cores = detectCores(logical = F))
SmoothTrack <- rbindlist(SmoothTrack)

SmoothTrack <-
  SmoothTrack %>% 
  group_by(
    UNIVERSAL_TRACK_ID 
  ) %>% 
  mutate(
    ADJUSTED_POSITION_X = ABSOLUTE_POSITION_X/max(abs(ABSOLUTE_POSITION_X)),
    ADJUSTED_POSITION_Y = ABSOLUTE_POSITION_Y/max(abs(ABSOLUTE_POSITION_Y)),
    
    MAX_QUERY = max(COMPLEMENTARY_NORMALIZED_INTENSITY_1),
    adjusted_t = t/max(t)
  ) %>% 
  mutate(
    MAX_QUERY = round(MAX_QUERY/2)*2,
    adjusted_t = round(adjusted_t*100)/100
  ) %>% 
  group_by(
    COHORT,
    MAX_QUERY
  ) %>% 
  mutate(
    N = ifelse(FRAMES_ADJUSTED == 1, 1, 0)
  ) %>% 
  mutate(
    N = sum(N)
  ) %>% 
  filter(
    N >= 50
  ) %>% 
  as.data.table()

SmoothTrack$COHORT <- factor(SmoothTrack$COHORT, levels = paste("MyD88", PROTEIN_ORDER))

ggplot() +
  geom_path(
    data = SmoothTrack,
    aes(
      # ABSOLUTE_POSITION_X, ABSOLUTE_POSITION_Y,
      ADJUSTED_POSITION_X, ADJUSTED_POSITION_Y,
      group = UNIVERSAL_TRACK_ID
      # color = new_id
    ),
    color = "black",
    alpha = 0.1
  ) +
  # scale_color_distiller(
  #   palette = "Spectral"
  # ) +
  # scale_x_continuous(
  #   limits = c(0, 1)
  # ) +
  # facet_grid(
  #   COHORT~MAX_QUERY
  # ) +
  theme_void() +
  coord_fixed()


TimeRandom <-
  SmoothTrack %>%
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    DELTA_X = ABSOLUTE_POSITION_X - lag(ABSOLUTE_POSITION_X),
    DELTA_Y = ABSOLUTE_POSITION_Y - lag(ABSOLUTE_POSITION_Y)
  ) %>% 
  mutate(
    DISPLACEMENT = sqrt(ABSOLUTE_POSITION_X^2 + ABSOLUTE_POSITION_Y^2)
  ) %>% 
  group_by(
    COHORT,
    MAX_QUERY,
    # adjusted_t
    # t,
    t_end
  ) %>%
  summarize(
    DISPLACEMENT = median(DISPLACEMENT, na.rm = T)
  ) %>%
  as.data.table()


ggplot(
  TimeRandom,
  aes(
    x = t,
    # x = t_end,
    y = DISPLACEMENT,
    color = MAX_QUERY,
    group = MAX_QUERY
  )
) +
  geom_path() +
  scale_color_viridis(
    option = "plasma"
  ) +
  facet_wrap(
    ~COHORT
  ) +
  theme_classic()


