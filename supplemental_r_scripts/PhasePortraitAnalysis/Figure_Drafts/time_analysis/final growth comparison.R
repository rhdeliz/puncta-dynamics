LEAD_LAG = 5
TIME_THRESHOLD = .5
START_END_FACTOR = 3
TIME_BINS = 5
PERCENT_THRESHOLD = 3
FRAMES_CUTOFF = 50

GrowthTable <-
  NormStatTable %>% 
  filter(
    # COHORT %in% c("MyD88 TRAF6", "MyD88 NEMO", "MyD88 HOIL1")
    # COHORT == "MyD88 NEMO",
    # LIGAND_DENSITY_CAT == 32
  ) %>%
  group_by(
    UNIVERSAL_TRACK_ID,
    FRAMES_ADJUSTED
  ) %>% 
  mutate(
    N = n()
  ) %>% 
  filter(
    N == 1
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    MAX_REFERENCE = max(REFERENCE_TOTAL_INTENSITY),
    MAX_QUERY = max(QUERY_TOTAL_INTENSITY),
    N = n(),
    STARTING_INTENSITY = median(REFERENCE_TOTAL_INTENSITY[1:3]) + median(QUERY_TOTAL_INTENSITY[1:3])
  ) %>% 
  filter(
    MAX_REFERENCE >= 2,
    MAX_QUERY >= 2,
    N > 12,
    STARTING_INTENSITY <= 6,
    FRAMES_ADJUSTED <= FRAMES_CUTOFF
  ) %>% 
  filter(
    MAX_REFERENCE >= 6 | MAX_QUERY >= 6
  ) %>% 
  arrange(
    FRAMES_ADJUSTED
  ) %>% 
  mutate(
    SCALED_REFERENCE_TOTAL_INTENSITY = (REFERENCE_TOTAL_INTENSITY-min(REFERENCE_TOTAL_INTENSITY))/(max(REFERENCE_TOTAL_INTENSITY)-min(REFERENCE_TOTAL_INTENSITY)),
    SCALED_QUERY_TOTAL_INTENSITY = (QUERY_TOTAL_INTENSITY-min(QUERY_TOTAL_INTENSITY))/(max(QUERY_TOTAL_INTENSITY)-min(QUERY_TOTAL_INTENSITY)),
  ) %>%
  mutate(
    DELTA_REFERENCE = lead(REFERENCE_TOTAL_INTENSITY, LEAD_LAG) - lag(REFERENCE_TOTAL_INTENSITY, LEAD_LAG),
    DELTA_QUERY = lead(QUERY_TOTAL_INTENSITY, LEAD_LAG) - lag(QUERY_TOTAL_INTENSITY, LEAD_LAG)
  ) %>%
  drop_na(
    DELTA_REFERENCE,
    DELTA_QUERY
  ) %>%
  arrange(
    FRAMES_ADJUSTED
  ) %>% 
  mutate(
    SCALED_FRAMES_ADJUSTED = (FRAMES_ADJUSTED-min(FRAMES_ADJUSTED))/(max(FRAMES_ADJUSTED)-min(FRAMES_ADJUSTED)),
    QRY_FUNCTION = cumsum(DELTA_QUERY)
  ) %>% 
  arrange(
    -FRAMES_ADJUSTED
  ) %>% 
  mutate(
    REF_FUNCTION = cumsum(DELTA_REFERENCE)
  ) %>% 
  arrange(
    FRAMES_ADJUSTED
  ) %>% 
  mutate(
    REF_FUNCTION = REF_FUNCTION/max(abs(REF_FUNCTION)),
    QRY_FUNCTION = QRY_FUNCTION/max(abs(QRY_FUNCTION))
  ) %>% 
  mutate(
    START_DELTA_REF = ifelse(SCALED_FRAMES_ADJUSTED <= TIME_THRESHOLD, REF_FUNCTION, 0),
    START_DELTA_QRY = ifelse(SCALED_FRAMES_ADJUSTED <= TIME_THRESHOLD, QRY_FUNCTION, 0),
    END_DELTA_REF = ifelse(SCALED_FRAMES_ADJUSTED > TIME_THRESHOLD, REF_FUNCTION, 0),
    END_DELTA_QRY = ifelse(SCALED_FRAMES_ADJUSTED > TIME_THRESHOLD, QRY_FUNCTION, 0),
    START_FRAME_COUNT = ifelse(SCALED_FRAMES_ADJUSTED <= TIME_THRESHOLD, 1, 0),
    END_FRAME_COUNT = ifelse(SCALED_FRAMES_ADJUSTED > TIME_THRESHOLD, 1, 0)
  ) %>% 
  mutate(
    START_FRAME_COUNT = sum(START_FRAME_COUNT)/(sum(START_FRAME_COUNT) + sum(END_FRAME_COUNT)),
    END_FRAME_COUNT = sum(END_FRAME_COUNT)/(sum(START_FRAME_COUNT) + sum(END_FRAME_COUNT))
  ) %>% 
  mutate(
    # START_FRAME_COUNT = START_FRAME_COUNT/TIME_THRESHOLD,
    # END_FRAME_COUNT = END_FRAME_COUNT/TIME_THRESHOLD
    START_FRAME_COUNT = 1,
    END_FRAME_COUNT = 1
  ) %>% 
  mutate(
    START_DELTA_REF = mean(START_DELTA_REF, na.rm = T)/START_FRAME_COUNT,
    START_DELTA_QRY = mean(START_DELTA_QRY, na.rm = T)/START_FRAME_COUNT,
    END_DELTA_REF = mean(END_DELTA_REF, na.rm = T)/END_FRAME_COUNT,
    END_DELTA_QRY = mean(END_DELTA_QRY, na.rm = T)/END_FRAME_COUNT
  ) %>% 
  mutate(
    START_DELTA_REF_BIN = round(START_DELTA_REF*START_END_FACTOR)/START_END_FACTOR,
    END_DELTA_REF_BIN = round(END_DELTA_REF*START_END_FACTOR)/START_END_FACTOR,
    START_DELTA_QRY_BIN = round(START_DELTA_QRY*START_END_FACTOR)/START_END_FACTOR,
    END_DELTA_QRY_BIN = round(END_DELTA_QRY*START_END_FACTOR)/START_END_FACTOR,
    # FRAMES_ADJUSTED = round(TIME_ADJUSTED/4/25)*4*25
    FRAMES_ADJUSTED = round(SCALED_FRAMES_ADJUSTED*TIME_BINS)/TIME_BINS
  ) %>%
  # filter(
  #   START_DELTA_REF >= 0,
  #   START_DELTA_QRY >= 0
  # ) %>% 
  # filter(
  #   END_DELTA_REF >= 0 |
  #     END_DELTA_QRY >= 0
  # ) %>% 
  group_by(
    COHORT,
    LIGAND_DENSITY_CAT,
    FRAMES_ADJUSTED,
    START_DELTA_REF_BIN,
    END_DELTA_REF_BIN,
    START_DELTA_QRY_BIN,
    END_DELTA_QRY_BIN
  ) %>% 
  summarize(
    REF_FUNCTION = mean(REF_FUNCTION),
    QRY_FUNCTION = mean(QRY_FUNCTION),
    SD_REFERENCE_TOTAL_INTENSITY = sd(REFERENCE_TOTAL_INTENSITY)/sqrt(n()),
    SD_QUERY_TOTAL_INTENSITY = sd(QUERY_TOTAL_INTENSITY)/sqrt(n()),
    
    REFERENCE_TOTAL_INTENSITY = mean(REFERENCE_TOTAL_INTENSITY),
    QUERY_TOTAL_INTENSITY = mean(QUERY_TOTAL_INTENSITY),
    
    SCALED_REFERENCE_TOTAL_INTENSITY = mean(SCALED_REFERENCE_TOTAL_INTENSITY),
    SCALED_QUERY_TOTAL_INTENSITY = mean(SCALED_QUERY_TOTAL_INTENSITY),
    
    N = NROW(unique(UNIVERSAL_TRACK_ID))
  ) %>% 
  group_by(
    COHORT,
    LIGAND_DENSITY_CAT,
    START_DELTA_REF_BIN,
    END_DELTA_REF_BIN,
    START_DELTA_QRY_BIN,
    END_DELTA_QRY_BIN,
  ) %>%
  mutate(
    N = max(N),
    GROUP = paste(
      START_DELTA_REF_BIN,
      END_DELTA_REF_BIN,
      START_DELTA_QRY_BIN,
      END_DELTA_QRY_BIN
    ),
    X_LABEL = REFERENCE_TOTAL_INTENSITY[n()],
    Y_LABEL = QUERY_TOTAL_INTENSITY[n()],
    SCALED_X_LABEL = SCALED_REFERENCE_TOTAL_INTENSITY[n()],
    SCALED_Y_LABEL = SCALED_QUERY_TOTAL_INTENSITY[n()]
  ) %>% 
  group_by(
    COHORT,
    LIGAND_DENSITY_CAT,
    GROUP
  ) %>% 
  mutate(
    N = ifelse(FRAMES_ADJUSTED == min(FRAMES_ADJUSTED), N, 0)
  ) %>% 
  group_by(
    COHORT,
    LIGAND_DENSITY_CAT
  ) %>% 
  mutate(
    N_TOTAL = sum(N)
  ) %>% 
  group_by(
    COHORT,
    LIGAND_DENSITY_CAT,
    GROUP
  ) %>% 
  mutate(
    N = round(sum(N)/N_TOTAL*100, 1),
    COHORT = factor(COHORT, levels = paste(USE_REFERENCE_PROTEIN, PROTEIN_ORDER))
  ) %>% 
  filter(
    N >= PERCENT_THRESHOLD
  ) %>% 
  as.data.table() %>% 
  ggplot(
    aes(
      # SCALED_REFERENCE_TOTAL_INTENSITY,
      # SCALED_QUERY_TOTAL_INTENSITY,
      
      REFERENCE_TOTAL_INTENSITY,
      QUERY_TOTAL_INTENSITY,
      group = GROUP,
      color = GROUP,
      size = N
    )
  ) + 
  # geom_errorbar(
  #   aes(
  #     ymin = QUERY_TOTAL_INTENSITY - SD_QUERY_TOTAL_INTENSITY,
  #     ymax = QUERY_TOTAL_INTENSITY + SD_QUERY_TOTAL_INTENSITY
  #   ),
  #   # color = "black",
  #   size = 1
  # ) +
  # geom_errorbar(
  #   aes(
  #     xmin = REFERENCE_TOTAL_INTENSITY -SD_REFERENCE_TOTAL_INTENSITY,
  #     xmax = REFERENCE_TOTAL_INTENSITY +SD_REFERENCE_TOTAL_INTENSITY,
  #   ),
  #   # color = "black",
  #   size = 1
  # ) +
  geom_path(
    alpha = .75
  ) +
  # geom_point() +
geom_text(
  aes(x = X_LABEL, y = Y_LABEL, label = paste(N, "%")),
  color = "black"
) +
  # geom_text(
  #   aes(x = SCALED_X_LABEL, y = SCALED_Y_LABEL, label = paste(N, "%")),
  #   color = "black"
  # ) +
  scale_color_viridis(
    option = "turbo",
    discrete =  T,
    guide = "none"
  ) +
  # facet_wrap(~UNIVERSAL_TRACK_ID)+
  # facet_grid(
  #   START_DELTA_QRY_BIN+END_DELTA_QRY_BIN~START_DELTA_REF_BIN+END_DELTA_REF_BIN
  # ) +
  scale_size(
    guide = "none"
  ) +
  labs(
    x = "Normalized Reference Intensity",
    y = "Normalized Query Intensity"
  ) +
  facet_wrap(
    ~ COHORT,
    # ~COHORT+LIGAND_DENSITY_CAT
    scales = "free"
  ) +
  # coord_fixed() +
  theme_classic()

GrowthTable
