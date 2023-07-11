# Select data
Trajectories <-
  NormStatTable %>% 
  filter(
    COHORT == "MyD88 TRAF6",
    LIGAND_DENSITY_CAT == 32
  ) %>%
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>%
  mutate(
    LIFETIME = max(FRAMES_ADJUSTED)
  ) %>%
  filter(
    LIFETIME >= 19
  ) %>%
  mutate(
    REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 1, n = 7),
    QUERY_TOTAL_INTENSITY = signal::sgolayfilt(QUERY_TOTAL_INTENSITY, p = 1, n = 7),
  ) %>%
  mutate(
    REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 3, n = 15),
    QUERY_TOTAL_INTENSITY = signal::sgolayfilt(QUERY_TOTAL_INTENSITY, p = 3, n = 15)
  ) %>%
  mutate(
    ROUNDED_REFERENCE_TOTAL_INTENSITY = round(REFERENCE_TOTAL_INTENSITY/STEP_SIZE)*STEP_SIZE,
    ROUNDED_QUERY_TOTAL_INTENSITY = round(QUERY_TOTAL_INTENSITY/STEP_SIZE)*STEP_SIZE
  ) %>% 
  group_by(
    ROUNDED_REFERENCE_TOTAL_INTENSITY,
    ROUNDED_QUERY_TOTAL_INTENSITY
  ) %>% 
  mutate(
    N = n()
  ) %>% 
  filter(
    N >= 50
  ) %>% 
  as.data.table()

# Identify which tracks pass through particular box
TracjectoryPasses <-
  Trajectories %>% 
  select(
    ROUNDED_REFERENCE_TOTAL_INTENSITY,
    ROUNDED_QUERY_TOTAL_INTENSITY,
    UNIVERSAL_TRACK_ID
  ) %>% 
  distinct() %>% 
  as.data.table()

ExpandTable <- function(TrackX){
  TempTrajectories <-
    Trajectories %>% 
    filter(
      UNIVERSAL_TRACK_ID == TracjectoryPasses$UNIVERSAL_TRACK_ID[TrackX]
    ) %>% 
    mutate(
      REFERENCE_LABEL = TracjectoryPasses$ROUNDED_REFERENCE_TOTAL_INTENSITY[TrackX],
      QUERY_LABEL = TracjectoryPasses$ROUNDED_QUERY_TOTAL_INTENSITY[TrackX]
    ) %>% 
    as.data.table()
  
  return(TempTrajectories)
}
ExpandedTrajectories <- mclapply(1:NROW(TracjectoryPasses), ExpandTable)
ExpandedTrajectories <- rbindlist(ExpandedTrajectories)

# Plot prep (shade, grid box id)
ExpandedTrajectories <-
  ExpandedTrajectories %>% 
  group_by(
    REFERENCE_LABEL,
    QUERY_LABEL
  ) %>% 
  mutate(
    SHADE = n()
  ) %>% 
  ungroup() %>% 
  mutate(
    SHADE = max(SHADE)/SHADE,
    REFERENCE_LABEL = paste(REFERENCE_PROTEIN, REFERENCE_LABEL),
    QUERY_LABEL = paste(QUERY_PROTEIN, QUERY_LABEL)
  ) %>% 
  mutate(
    SHADE = SHADE/max(SHADE)
  ) %>% 
  as.data.table()

# Label grid properly
REFERENCE_STEPS = max(Trajectories$ROUNDED_REFERENCE_TOTAL_INTENSITY)/STEP_SIZE
REFERENCE_STEPS_SIZE = (0:REFERENCE_STEPS)*STEP_SIZE
REFERENCE_STEPS = paste(Trajectories$REFERENCE_PROTEIN[1], REFERENCE_STEPS_SIZE)

QUERY_STEPS = max(Trajectories$ROUNDED_QUERY_TOTAL_INTENSITY)/STEP_SIZE
QUERY_STEPS_SIZE = (0:QUERY_STEPS)*STEP_SIZE
QUERY_STEPS = paste(Trajectories$QUERY_PROTEIN[1], QUERY_STEPS_SIZE)

ExpandedTrajectories$REFERENCE_LABEL <- factor(ExpandedTrajectories$REFERENCE_LABEL, levels = REFERENCE_STEPS)
ExpandedTrajectories$QUERY_LABEL <- factor(ExpandedTrajectories$QUERY_LABEL, levels = rev(QUERY_STEPS))

# Phase portrait arrows
TrajectoryArrow <-
  ExpandedTrajectories %>% 
  group_by(
    ROUNDED_REFERENCE_TOTAL_INTENSITY,
    ROUNDED_QUERY_TOTAL_INTENSITY,
    REFERENCE_LABEL,
    QUERY_LABEL
  ) %>% 
  summarize(
    DELTA_REFERENCE_TOTAL_INTENSITY = median(DELTA_REFERENCE_TOTAL_INTENSITY, na.rm = T),
    DELTA_QUERY_TOTAL_INTENSITY = median(DELTA_QUERY_TOTAL_INTENSITY, na.rm = T)
  ) %>%
  group_by(
    REFERENCE_LABEL,
    QUERY_LABEL
  ) %>% 
  mutate(
    DELTA_REFERENCE_TOTAL_INTENSITY = DELTA_REFERENCE_TOTAL_INTENSITY/max(DELTA_REFERENCE_TOTAL_INTENSITY),
    DELTA_QUERY_TOTAL_INTENSITY = DELTA_QUERY_TOTAL_INTENSITY/max(DELTA_QUERY_TOTAL_INTENSITY)
  ) %>% 
  as.data.table()


# To make box of where filter was applied
TrajectoryPoints <- expand.grid(REFERENCE_STEPS, QUERY_STEPS)
names(TrajectoryPoints) <- c("REFERENCE_LABEL", "QUERY_LABEL")

TrajectoryPointsSize <- expand.grid(REFERENCE_STEPS_SIZE, QUERY_STEPS_SIZE)
names(TrajectoryPointsSize) <- c("ROUNDED_REFERENCE_TOTAL_INTENSITY", "ROUNDED_QUERY_TOTAL_INTENSITY")

TrajectoryPoints <- cbind(TrajectoryPoints, TrajectoryPointsSize)
remove(TrajectoryPointsSize)

TrajectoryPoints$REFERENCE_LABEL <- factor(TrajectoryPoints$REFERENCE_LABEL, levels = REFERENCE_STEPS)
TrajectoryPoints$QUERY_LABEL <- factor(TrajectoryPoints$QUERY_LABEL, levels = rev(QUERY_STEPS))

ARROW_LENGTH = STEP_SIZE/max(
  c(TrajectoryArrow$ROUNDED_REFERENCE_TOTAL_INTENSITY,
    TrajectoryArrow$ROUNDED_QUERY_TOTAL_INTENSITY
  ))*2


ggplot() +
  geom_rect(
    data = TrajectoryPoints,
    aes(
      xmin = ROUNDED_REFERENCE_TOTAL_INTENSITY-STEP_SIZE/2,
      ymin = ROUNDED_QUERY_TOTAL_INTENSITY-STEP_SIZE/2,
      xmax = ROUNDED_REFERENCE_TOTAL_INTENSITY+STEP_SIZE/2,
      ymax = ROUNDED_QUERY_TOTAL_INTENSITY+STEP_SIZE/2
    ),
    color = "red",
    fill = NA
  ) +
  geom_arrow(
    data = TrajectoryArrow,
    aes(
      x = ROUNDED_REFERENCE_TOTAL_INTENSITY,
      y = ROUNDED_QUERY_TOTAL_INTENSITY,
      mag = .5,
      angle = atan2(DELTA_QUERY_TOTAL_INTENSITY, DELTA_REFERENCE_TOTAL_INTENSITY)*180/pi,
      color = sqrt(DELTA_REFERENCE_TOTAL_INTENSITY^2 + DELTA_QUERY_TOTAL_INTENSITY^2)
    ),
    arrow.length = ARROW_LENGTH,
    lineend = "square"
  ) +
  scale_size(range = c(ARROW_LENGTH/10, ARROW_LENGTH), guide = "none") +
  scale_color_viridis(
    # option = "magma",
    guide = "none"
  ) +
  ggnewscale::new_scale_color() +
  geom_path(
    data = ExpandedTrajectories,
    aes(
      x = REFERENCE_TOTAL_INTENSITY,
      y = QUERY_TOTAL_INTENSITY,
      group = UNIVERSAL_TRACK_ID,
      alpha = SHADE
    ),
    size = 0.1
  ) +
  scale_color_gradient(
    low = "black",
    high = "grey"
  ) +
  labs(
    x = paste(ExpandedTrajectories$REFERENCE_PROTEIN[1], "Size"),
    y = paste(ExpandedTrajectories$QUERY_PROTEIN[1], "Size")
  ) + 
  facet_grid(
    QUERY_LABEL~REFERENCE_LABEL
  ) +
  theme_classic() +
  theme(
    
  ) +
  coord_fixed()


