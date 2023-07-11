LEAD_LAG = 5
TIME_CUTOFF = .5

GrowthTable <-
  NormStatTable %>% 
  filter(
    COHORT == "MyD88 TRAF6",
    LIGAND_DENSITY_CAT == 32
    # IMAGE == "20200110 cl069 1R myd88gfp mScarlett-traf6"
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
    STARTING_INTENSITY <= 6
  ) %>% 
  filter(
    MAX_REFERENCE >= 6 |
    MAX_QUERY >= 6
  ) %>% 
  arrange(
    FRAMES_ADJUSTED
  ) %>% 
  # mutate(
  #   REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 1, n = 7),
  #   QUERY_TOTAL_INTENSITY = signal::sgolayfilt(QUERY_TOTAL_INTENSITY, p = 1, n = 7),
  # ) %>%
  # mutate(
  #   REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 3, n = 11),
  #   QUERY_TOTAL_INTENSITY = signal::sgolayfilt(QUERY_TOTAL_INTENSITY, p = 3, n = 11)
  # ) %>%
  mutate(
    REFERENCE_TOTAL_INTENSITY = (REFERENCE_TOTAL_INTENSITY-min(REFERENCE_TOTAL_INTENSITY))/(max(REFERENCE_TOTAL_INTENSITY)-min(REFERENCE_TOTAL_INTENSITY)),
    QUERY_TOTAL_INTENSITY = (QUERY_TOTAL_INTENSITY-min(QUERY_TOTAL_INTENSITY))/(max(QUERY_TOTAL_INTENSITY)-min(QUERY_TOTAL_INTENSITY)),
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
  filter(
    FRAMES_ADJUSTED <= 50
  ) %>%
  mutate(
    FRAMES_ADJUSTED = (FRAMES_ADJUSTED-min(FRAMES_ADJUSTED))/(max(FRAMES_ADJUSTED)-min(FRAMES_ADJUSTED))
  ) %>% 
  mutate(
    QRY_FUNCTION = cumsum(DELTA_QUERY),
    START_DELTA_REF = ifelse(FRAMES_ADJUSTED <= TIME_CUTOFF, DELTA_REFERENCE, NA),
    START_DELTA_QRY = ifelse(FRAMES_ADJUSTED <= TIME_CUTOFF, DELTA_QUERY, NA),
    END_DELTA_REF = ifelse(FRAMES_ADJUSTED > TIME_CUTOFF, DELTA_REFERENCE, NA),
    END_DELTA_QRY = ifelse(FRAMES_ADJUSTED > TIME_CUTOFF, DELTA_QUERY, NA)
  ) %>% 
  mutate(
    START_DELTA_REF = mean(START_DELTA_REF, na.rm = T),
    START_DELTA_QRY = mean(START_DELTA_QRY, na.rm = T),
    END_DELTA_REF = mean(END_DELTA_REF, na.rm = T),
    END_DELTA_QRY = mean(END_DELTA_QRY, na.rm = T)
  ) %>% 
  mutate(
    START_DELTA_REF = (START_DELTA_REF - min(REFERENCE_TOTAL_INTENSITY))/(max(REFERENCE_TOTAL_INTENSITY) - min(REFERENCE_TOTAL_INTENSITY)),
    START_DELTA_QRY = (START_DELTA_QRY - min(REFERENCE_TOTAL_INTENSITY))/(max(REFERENCE_TOTAL_INTENSITY) - min(REFERENCE_TOTAL_INTENSITY)),
    END_DELTA_REF = (END_DELTA_REF - min(REFERENCE_TOTAL_INTENSITY))/(max(REFERENCE_TOTAL_INTENSITY) - min(REFERENCE_TOTAL_INTENSITY)),
    END_DELTA_QRY = (END_DELTA_QRY - min(REFERENCE_TOTAL_INTENSITY))/(max(REFERENCE_TOTAL_INTENSITY) - min(REFERENCE_TOTAL_INTENSITY))
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
    REFERENCE_TOTAL_INTENSITY = (REFERENCE_TOTAL_INTENSITY - min(REFERENCE_TOTAL_INTENSITY))/(max(REFERENCE_TOTAL_INTENSITY)-min(REFERENCE_TOTAL_INTENSITY)),
    QUERY_TOTAL_INTENSITY = (QUERY_TOTAL_INTENSITY - min(QUERY_TOTAL_INTENSITY))/(max(QUERY_TOTAL_INTENSITY)-min(QUERY_TOTAL_INTENSITY)),
    
    REF_FUNCTION = REF_FUNCTION/max(abs(REF_FUNCTION)),
    QRY_FUNCTION = QRY_FUNCTION/max(abs(QRY_FUNCTION)),
    
    # FRAMES_ADJUSTED = (FRAMES_ADJUSTED-min(FRAMES_ADJUSTED))/(max(FRAMES_ADJUSTED)-min(FRAMES_ADJUSTED))
  ) %>% 
  mutate(
    # STARTING_DELTA_REF = mean(DELTA_REFERENCE[1:10]),
    # STARTING_DELTA_QRY = mean(DELTA_QUERY[1:10]),
    
    START_CORRELATION = cor(REFERENCE_TOTAL_INTENSITY[1:10], QUERY_TOTAL_INTENSITY[1:10])
    
    # STARTING_DELTA_REF = mean(REFERENCE_TOTAL_INTENSITY[1:10]),
    # STARTING_DELTA_QRY = mean(QUERY_TOTAL_INTENSITY[1:10]),
  ) %>% 
  mutate(
    CORRELATION = cor(-REF_FUNCTION, QRY_FUNCTION, method = "spearman"),
    # CORRELATION = cor(-REF_FUNCTION, QRY_FUNCTION, method = "spearman"),
    STARTING_FRACTION = median(REFERENCE_TOTAL_INTENSITY[1:5])/(median(REFERENCE_TOTAL_INTENSITY[1:5]) + median(QUERY_TOTAL_INTENSITY[1:5])),
    # DELTA_FRACTION = STARTING_DELTA_REF/STARTING_DELTA_QRY,
    # END_FRACTION = median(REFERENCE_TOTAL_INTENSITY[(n()-5):n()])/(median(REFERENCE_TOTAL_INTENSITY[(n()-5):n()]) + median(QUERY_TOTAL_INTENSITY[(n()-5):n()])),
  ) %>%  
  mutate(
    ROUND_REFERENCE_TOTAL_INTENSITY = round(REFERENCE_TOTAL_INTENSITY/2)*2,
    ROUND_QUERY_TOTAL_INTENSITY = round(QUERY_TOTAL_INTENSITY/2)*2,

    CORRELATION = round(CORRELATION*4)/4,
    STARTING_FRACTION = round(STARTING_FRACTION*5)/5,
    
    # STARTING_DELTA_REF = round(STARTING_DELTA_REF, 1),
    # STARTING_DELTA_QRY = round(STARTING_DELTA_QRY, 1),
    
    START_CORRELATION = round(START_CORRELATION*5)/5
    # END_FRACTION = round(END_FRACTION*5)/5
    # DELTA_FRACTION = round(DELTA_FRACTION/2)*2
  ) %>%
  arrange(
    CORRELATION,
    STARTING_FRACTION,
    FRAMES_ADJUSTED
  ) %>% 
  as.data.table()

cor(GrowthTable$REF_FUNCTION, GrowthTable$QRY_FUNCTION, method = "spearman")


GrowthTable %>%
  filter(
    # CORRELATION >= -.2,
    # CORRELATION <= .2
  ) %>%
  # mutate(
  #   UNIVERSAL_TRACK_ID = factor(UNIVERSAL_TRACK_ID, levels = unique(GrowthTable$UNIVERSAL_TRACK_ID))
  # ) %>%
  # mutate(
  #   FRAMES_ADJUSTED = round(FRAMES_ADJUSTED*25)
  # ) %>%
  # group_by(
  #   CORRELATION,
  #   STARTING_FRACTION,
  #   FRAMES_ADJUSTED
  # ) %>%
  # summarise(
  #   UNIVERSAL_TRACK_ID = 1,
  #   REF_FUNCTION = median(REF_FUNCTION),
  #   QRY_FUNCTION = median(QRY_FUNCTION)
  # ) %>%
  as.data.table() %>%
ggplot(
  aes(group = UNIVERSAL_TRACK_ID)
) +
  geom_hline(
    yintercept = 0, color = "black"
  ) +
  geom_path(
    aes(
      FRAMES_ADJUSTED,
      # REFERENCE_TOTAL_INTENSITY
      REF_FUNCTION
    ),
    color = "green",
    # alpha = .5
  ) +
  geom_path(
    aes(
      FRAMES_ADJUSTED,
      QRY_FUNCTION
      # QUERY_TOTAL_INTENSITY
    ),
    color = "magenta",
    # alpha = .5
  ) +
  # facet_wrap(
  #   ~UNIVERSAL_TRACK_ID,
  # ) +
  facet_grid(
    STARTING_FRACTION~CORRELATION, 
    scales = "free"
  ) +
  theme_classic()









GrowthTable %>%
as.data.table() %>%
  ggplot(
    aes(group = UNIVERSAL_TRACK_ID)
  ) +
  geom_hline(
    yintercept = 0, color = "black"
  ) +
  geom_path(
    aes(
      REFERENCE_TOTAL_INTENSITY,
      QUERY_TOTAL_INTENSITY
      # REF_FUNCTION
    ),
    # color = "green",
    size = .1
  ) +
  # geom_path(
  #   aes(
  #     FRAMES_ADJUSTED,
  #     # QRY_FUNCTION
  #     QUERY_TOTAL_INTENSITY
  #   ),
  #   color = "magenta",
  #   size = .1
  # ) +
  # facet_wrap(
  #   ~UNIVERSAL_TRACK_ID,
  # ) +
  facet_grid(
    ~CORRELATION,
    # STARTING_FRACTION~CORRELATION,
    # STARTING_DELTA_QRY~STARTING_DELTA_REF,
    scales = "free"
  ) +
  theme_classic()
