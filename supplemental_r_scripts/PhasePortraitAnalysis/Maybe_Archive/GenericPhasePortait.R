# Indicate parameters
X_VARIABLE = "TOTAL_INTENSITY"
X_TITLE = "Total Intensity of Cluster"
X_BIN_SIZE = 1

# X_VARIABLE = "TIME_ADJUSTED"
# X_TITLE = "Cluster Time (s)"
# X_BIN_SIZE = 8
# 
Y_VARIABLE = "RATIO"
Y_TITLE = "Stoichiometry"
Y_BIN_SIZE = 0.05

# Y_VARIABLE = "TOTAL_INTENSITY"
# Y_TITLE = "Total Intensity of Cluster"
# Y_BIN_SIZE = 2

# Frames since landing
FRAMES_SINCE_LANDING_BIN = 200
# Add flow lines or not
STREAMLINE = F

# Table Name
SaveName <-
  paste0(
    "Normalized StatTable - ",
    "LeadLag ", LEAD_LAG, " - ",
    "StepSize ", round(STEP_SIZE, 2),
    ".csv.gz"
  )

SaveName <- gsub(" - \\.", "\\.", SaveName)

# Import table
# StatTable <- fread(file = file.path(OUTPUT_DIRECTORY, SaveName))

# Pre-processing
TempTable <-
  StatTable %>% 
  filter(
    FRAMES_SINCE_LANDING_CAT <= FRAMES_SINCE_LANDING_BIN,
    # COHORT == "MyD88 TRAF6",
    LIGAND_DENSITY_CAT == 32,
    # To mirror phase portrait
    # Filter out unusually large clusters
    REFERENCE_TOTAL_INTENSITY <= 2,
    REFERENCE_TOTAL_INTENSITY >= 0,
    QUERY_TOTAL_INTENSITY <= 2,
    QUERY_TOTAL_INTENSITY >= 0,
    FRAMES_ADJUSTED > min(FRAMES_ADJUSTED),
    FRAMES_ADJUSTED < max(FRAMES_ADJUSTED)
  ) %>% 
  arrange(
    COHORT,
    LIGAND_DENSITY_CAT,
    FPS,
    IMAGE,
    UNIVERSAL_TRACK_ID,
    FRAMES_ADJUSTED
  ) %>% 
  group_by(
    COHORT,
    LIGAND_DENSITY_CAT,
    FPS
  ) %>% 
  mutate(
    REFERENCE_TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY*max(MAX_REFERENCE),
    QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY*max(MAX_QUERY)
  ) %>%
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>%
  mutate(
    TIME_ADJUSTED = TIME_ADJUSTED - min(TIME_ADJUSTED),
    # Get track parameters
    LIFETIME = max(FRAMES_ADJUSTED) - min(FRAMES_ADJUSTED) + 1,
    TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY + QUERY_TOTAL_INTENSITY,
  ) %>%
  mutate(
    RATIO = REFERENCE_TOTAL_INTENSITY/(TOTAL_INTENSITY)
  ) %>% 
  filter(
    # To mirror phase portrait
    LIFETIME >= (LEAD_LAG*2+1)
  ) %>%
  as.data.frame()

TempTable <- 
  TempTable %>% 
  mutate(
    X = !!as.name(X_VARIABLE),
    Y = !!as.name(Y_VARIABLE),
    SAMPLE_ID = UNIVERSAL_TRACK_ID,
    # R Facet
    GROUP = paste(COHORT, LIGAND_DENSITY_CAT, FPS)
  ) %>% 
  arrange(
    GROUP,
    SAMPLE_ID,
    X,
    Y
  ) %>% 
  group_by(
    SAMPLE_ID
  ) %>% 
  mutate(
    LEAD_X = lead(X, LEAD_LAG),
    LAG_X = lag(X, LEAD_LAG),
    
    LEAD_Y = lead(Y, LEAD_LAG),
    LAG_Y = lag(Y, LEAD_LAG),
    
    ROUNDED_X = round(X/X_BIN_SIZE)*X_BIN_SIZE,
    ROUNDED_Y = round(Y/Y_BIN_SIZE)*Y_BIN_SIZE,
  ) %>% 
  mutate(
    DELTA_X = LEAD_X - LAG_X,
    DELTA_Y = LEAD_Y - LAG_Y
  ) %>% 
  drop_na(
    DELTA_X,
    DELTA_Y
  ) %>% 
  group_by(
    GROUP,
    ROUNDED_X,
    ROUNDED_Y
  ) %>% 
  summarize(
    DELTA_X = median(DELTA_X),
    DELTA_Y = median(DELTA_Y),
    N = n()
  ) %>% 
  filter(
    N >= 50
  ) %>% 
  group_by(
    GROUP
  ) %>% 
  mutate(
    DELTA_X = DELTA_X/max(abs(DELTA_X)),
    DELTA_Y = DELTA_Y/max(abs(DELTA_Y))
  ) %>% 
  distinct() %>% 
  as.data.table()

# Make grid (needed for streamlines)
StreamTempTable <-
  TempTable %>%
  as_tibble() %>%
  arrange(
    GROUP,
    ROUNDED_X,
    ROUNDED_Y
  ) %>%
  group_by(
    GROUP
  ) %>%
  mutate(
    ROUNDED_X = round(ROUNDED_X/X_BIN_SIZE),
    ROUNDED_Y = round(ROUNDED_Y/Y_BIN_SIZE)
  ) %>%
  complete(
    ROUNDED_X = full_seq(c(min(ROUNDED_X), max(ROUNDED_X)), period = 1),
    ROUNDED_Y = full_seq(c(min(ROUNDED_Y), max(ROUNDED_Y)), period = 1)
  ) %>%
  mutate(
    ROUNDED_X = ROUNDED_X*X_BIN_SIZE,
    ROUNDED_Y = ROUNDED_Y*Y_BIN_SIZE
  ) %>% 
  as.data.table()

# Add missing derivatives
StreamTempTable$DELTA_X[is.na(StreamTempTable$DELTA_X)] = 0
StreamTempTable$DELTA_Y[is.na(StreamTempTable$DELTA_Y)] = 0
StreamTempTable$N[is.na(StreamTempTable$N)] = 0

if(STREAMLINE == T){
  
  # Plot
  ggplot(
  ) +
    geom_streamline(
      data = StreamTempTable,
      aes(
        x = ROUNDED_X,
        y = ROUNDED_Y,
        dx = DELTA_X,
        dy = DELTA_Y,
        color = sqrt(..dx..^2 + ..dy..^2),
        size = ..step..
        # alpha = ..step..
      ),
      arrow = NULL,
      # n = 10,
      # arrow.length = 0.3,
      # jitter = 4,
      L = 10, res = 10, lineend = "round"
    ) +
    geom_arrow(
      data = TempTable,
      aes(
        x = ROUNDED_X,
        y = ROUNDED_Y,
        mag = 0.5,
        angle = atan2(DELTA_Y, DELTA_X)*180/pi,
        color = sqrt(DELTA_Y^2 + DELTA_X^2),
        # color = FRAMES_ADJUSTED
      ),
      # size = .01, # Small arrow head
      # arrow.length = .5, # Small arrow head
      
      size = 1, # Big arrow head
      arrow.length = 1,# Big arrow head
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
      option = "plasma",
      guide = "none"
      # direction = -1,
      # trans = "log10",
      # breaks = trans_breaks("log10", function(x) 10^x),
      # labels = trans_format("log10", math_format(10^.x))
    )  +
    labs(
      x = X_TITLE,
      y = Y_TITLE
    ) +
    facet_wrap(
      ~GROUP
    ) +
    theme_classic() +
    coord_fixed(
      ratio = X_BIN_SIZE/Y_BIN_SIZE
    )
  
} else{
  
  # Plot
  ggplot(
  ) +
    geom_arrow(
      data = TempTable,
      aes(
        x = ROUNDED_X,
        y = ROUNDED_Y,
        mag = 0.5,
        angle = atan2(DELTA_Y, DELTA_X)*180/pi,
        color = sqrt(DELTA_Y^2 + DELTA_X^2),
        # color = FRAMES_ADJUSTED
      ),
      # size = .01, # Small arrow head
      # arrow.length = .5, # Small arrow head
      
      size = 1, # Big arrow head
      arrow.length = 1,# Big arrow head
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
      option = "plasma",
      guide = "none"
      # direction = -1,
      # trans = "log10",
      # breaks = trans_breaks("log10", function(x) 10^x),
      # labels = trans_format("log10", math_format(10^.x))
    )  +
    labs(
      x = X_TITLE,
      y = Y_TITLE
    ) +
    facet_wrap(
      ~GROUP
    ) +
    theme_classic() +
    coord_fixed(
      ratio = X_BIN_SIZE/Y_BIN_SIZE
    )
  
}