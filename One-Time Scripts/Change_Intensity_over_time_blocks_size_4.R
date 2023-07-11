
bl <- colorRampPalette(c("navy","royalblue","lightskyblue"))(200)                      
re <- colorRampPalette(c("mistyrose", "red2","darkred"))(200)

DeltaRefTime <-
  NormStatTable %>% 
  filter(
    # ROUNDED_REFERENCE_TOTAL_INTENSITY == 6
  ) %>% 
  mutate(
    # ROUNDED_REFERENCE_TOTAL_INTENSITY = ifelse(ROUNDED_REFERENCE_TOTAL_INTENSITY > 0, STEP_SIZE, 0),
    # ROUNDED_QUERY_TOTAL_INTENSITY = ifelse(ROUNDED_QUERY_TOTAL_INTENSITY > 0, STEP_SIZE, 0),
    # FRAMES_SINCE_LANDING_CAT = ifelse(FRAMES_SINCE_LANDING>= 150, 150, 0)
  ) %>% 
  group_by(
    COHORT,
    REFERENCE_PROTEIN,
    QUERY_PROTEIN,
    FRAMES_SINCE_LANDING_CAT,
    ROUNDED_QUERY_TOTAL_INTENSITY,
    ROUNDED_REFERENCE_TOTAL_INTENSITY
  ) %>% 
  filter(
    # ROUNDED_REFERENCE_TOTAL_INTENSITY == STEP_SIZE,
    # ROUNDED_QUERY_TOTAL_INTENSITY == STEP_SIZE

    ROUNDED_REFERENCE_TOTAL_INTENSITY <= 15,
    ROUNDED_QUERY_TOTAL_INTENSITY <= 15
  ) %>% 
  summarize(
    DELTA_REFERENCE_TOTAL_INTENSITY = median(DELTA_REFERENCE_TOTAL_INTENSITY),
    N = n()
  ) %>% 
  filter(
    N >= 100
  ) %>% 
  group_by(
    COHORT,
    QUERY_PROTEIN
  ) %>% 
  mutate(
    N = N/sum(N)*100
  ) %>% 
  as.data.table()

ggplot(
  DeltaRefTime,
  aes(
    FRAMES_SINCE_LANDING_CAT,
    ROUNDED_QUERY_TOTAL_INTENSITY,
    # color = DELTA_REFERENCE_TOTAL_INTENSITY,
    fill = DELTA_REFERENCE_TOTAL_INTENSITY
  )
) +
  geom_tile() +
  geom_hline(
    yintercept = STEP_SIZE
  ) +
  geom_vline(
    xintercept = 150
  ) +
  scale_fill_gradientn(
    colours=c(bl,"white", re), na.value = "grey98",
    limits = c(-max(DeltaRefTime$DELTA_REFERENCE_TOTAL_INTENSITY),
               max(DeltaRefTime$DELTA_REFERENCE_TOTAL_INTENSITY))) +
  labs(
    x = "Cell Time (s)",
    y = "Query N. Int.",
    fill = "d-MyD88"
  ) +
  # scale_x_continuous(
  #   limits = c(0, 250)
  # ) +
  facet_grid(
    ROUNDED_REFERENCE_TOTAL_INTENSITY~QUERY_PROTEIN
  ) +
  theme_classic()



DeltaRefTime <-
  NormStatTable %>% 
  filter(
    # ROUNDED_REFERENCE_TOTAL_INTENSITY == 6
  ) %>% 
  mutate(
    ROUNDED_REFERENCE_TOTAL_INTENSITY = ifelse(ROUNDED_REFERENCE_TOTAL_INTENSITY > 0, STEP_SIZE, 0),
    ROUNDED_QUERY_TOTAL_INTENSITY = ifelse(ROUNDED_QUERY_TOTAL_INTENSITY > 0, STEP_SIZE, 0),
    # FRAMES_SINCE_LANDING_CAT = ifelse(FRAMES_SINCE_LANDING>= 150, 150, 0)
  ) %>% 
  group_by(
    COHORT,
    REFERENCE_PROTEIN,
    QUERY_PROTEIN,
    FRAMES_SINCE_LANDING_CAT,
    ROUNDED_QUERY_TOTAL_INTENSITY,
    ROUNDED_REFERENCE_TOTAL_INTENSITY
  ) %>% 
  filter(
    ROUNDED_REFERENCE_TOTAL_INTENSITY == STEP_SIZE,
    ROUNDED_QUERY_TOTAL_INTENSITY == STEP_SIZE
    
    # ROUNDED_REFERENCE_TOTAL_INTENSITY <= 15,
    # ROUNDED_QUERY_TOTAL_INTENSITY <= 15
  ) %>% 
  summarize(
    DELTA_REFERENCE_TOTAL_INTENSITY = median(DELTA_REFERENCE_TOTAL_INTENSITY),
    N = n()
  ) %>% 
  filter(
    N >= 100
  ) %>% 
  group_by(
    COHORT,
    QUERY_PROTEIN
  ) %>% 
  mutate(
    N = N/sum(N)*100
  ) %>% 
  as.data.table()

ggplot(
  DeltaRefTime,
  aes(
    FRAMES_SINCE_LANDING_CAT,
    DELTA_REFERENCE_TOTAL_INTENSITY,
    color = QUERY_PROTEIN,
    group = QUERY_PROTEIN
  )
) +
  geom_point(
    aes(
      size = N
    )
  ) +
  geom_path() +
  geom_hline(
    yintercept = 0
  ) +
  # geom_vline(
  #   xintercept = 150
  # ) +
  labs(
    x = "Cell Time (s)",
    y = "d-MyD88",
    color ="Cell Line"
  ) +
  # scale_x_continuous(
  #   limits = c(0, 250)
  # ) +
  scale_color_manual(
    values = ColorTable,
    # labels = PROTEIN_ORDER[2:NROW(PROTEIN_ORDER)],
    limits = force,
    # drop = F,
    # guide = "none"
  ) +
  # facet_grid(
  #   ~QUERY_PROTEIN
  # ) +
  theme_classic()


NormStatTable %>% 
  filter(
    # ROUNDED_REFERENCE_TOTAL_INTENSITY == 6
  ) %>% 
  mutate(
    ROUNDED_REFERENCE_TOTAL_INTENSITY = ifelse(ROUNDED_REFERENCE_TOTAL_INTENSITY > 0, STEP_SIZE, 0),
    ROUNDED_QUERY_TOTAL_INTENSITY = ifelse(ROUNDED_QUERY_TOTAL_INTENSITY > 0, STEP_SIZE, 0),
    # FRAMES_SINCE_LANDING_CAT = ceiling(FRAMES_SINCE_LANDING/25)
  ) %>% 
  filter(
    ROUNDED_REFERENCE_TOTAL_INTENSITY == STEP_SIZE,
    ROUNDED_QUERY_TOTAL_INTENSITY == STEP_SIZE
  ) %>% 
  group_by(
    COHORT,
    REFERENCE_PROTEIN,
    QUERY_PROTEIN,
    # FRAMES_SINCE_LANDING_CAT,
  ) %>% 
  summarize(
    DELTA_REFERENCE_TOTAL_INTENSITY = median(DELTA_REFERENCE_TOTAL_INTENSITY),
    N = n()
  ) %>% 
  filter(
    N >= 100
  ) %>% 
  as.data.table() %>% 
ggplot(
  aes(
    QUERY_PROTEIN,
    DELTA_REFERENCE_TOTAL_INTENSITY,
    color = QUERY_PROTEIN,
    group = QUERY_PROTEIN
  )
) +
  geom_point() +
  geom_path() +
  geom_hline(
    yintercept = 0
  ) +
  labs(
    x = "Cell Time (s)",
    y = "d-MyD88",
    color = "Cell Line"
  ) +
  scale_color_manual(
    values = ColorTable,
    # labels = PROTEIN_ORDER[2:NROW(PROTEIN_ORDER)],
    limits = force,
    # drop = F,
    # guide = "none"
  ) +
  theme_classic()

