START_END_FACTOR = 3

GrowthTable %>%
  filter(
    # UNIVERSAL_TRACK_ID %in% unique(GrowthTable$UNIVERSAL_TRACK_ID)[1:20]
    # between(START_DELTA_REF, -1, 1),
    # between(END_DELTA_REF, -1, 1),
    # between(START_DELTA_QRY, -1, 1),
    # between(END_DELTA_QRY, -1, 1),
    # START_DELTA_REF >0,
    # END_DELTA_QRY >0
  ) %>%
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>%
  mutate(
    START_DELTA_REF = ifelse(FRAMES_ADJUSTED <= TIME_CUTOFF, REF_FUNCTION, 0),
    START_DELTA_QRY = ifelse(FRAMES_ADJUSTED <= TIME_CUTOFF, QRY_FUNCTION, 0),
    END_DELTA_REF = ifelse(FRAMES_ADJUSTED > TIME_CUTOFF, REF_FUNCTION, 0),
    END_DELTA_QRY = ifelse(FRAMES_ADJUSTED > TIME_CUTOFF, QRY_FUNCTION, 0)
  ) %>% 
  mutate(
    START_DELTA_REF = mean(START_DELTA_REF, na.rm = T),
    START_DELTA_QRY = mean(START_DELTA_QRY, na.rm = T),
    END_DELTA_REF = mean(END_DELTA_REF, na.rm = T),
    END_DELTA_QRY = mean(END_DELTA_QRY, na.rm = T)
  ) %>% 
  filter(
    # START_DELTA_REF>= 0,
    # START_DELTA_QRY >=0
  ) %>% 
  # mutate(
  #   START_DELTA_REF = mean(REF_FUNCTION),
  #   START_DELTA_QRY = mean(QRY_FUNCTION)
  # ) %>% 
  # mutate(
  #   START_DELTA_REF = START_DELTA_REF/max(abs(START_DELTA_REF)),
  #   END_DELTA_REF = END_DELTA_REF/max(abs(END_DELTA_REF)),
  #   START_DELTA_QRY = START_DELTA_QRY/max(abs(START_DELTA_QRY)),
  #   END_DELTA_QRY = START_DELTA_REF/max(abs(END_DELTA_QRY))
  # ) %>%
  # mutate(
  #   START_DELTA_REF = 1/START_DELTA_REF,
  #   END_DELTA_REF = 1/END_DELTA_REF,
  #   START_DELTA_QRY = 1/START_DELTA_QRY,
  #   END_DELTA_QRY = 1/END_DELTA_QRY
  # ) %>%
  mutate(
    START_DELTA_REF_BIN = round(START_DELTA_REF*START_END_FACTOR)/START_END_FACTOR,
    END_DELTA_REF_BIN = round(END_DELTA_REF*START_END_FACTOR)/START_END_FACTOR,
    START_DELTA_QRY_BIN = round(START_DELTA_QRY*START_END_FACTOR)/START_END_FACTOR,
    END_DELTA_QRY_BIN = round(END_DELTA_QRY*START_END_FACTOR)/START_END_FACTOR,
    
    # START_DELTA_REF_BIN = ifelse(START_DELTA_REF >= 0, "+R <t", "-R <t"),
    # END_DELTA_REF_BIN = ifelse(END_DELTA_REF >= 0, "+R >t", "-R >t"),
    # START_DELTA_QRY_BIN = ifelse(START_DELTA_QRY >= 0, "+Q <t", "-Q <t"),
    # END_DELTA_QRY_BIN = ifelse(END_DELTA_QRY >= 0, "+Q >t", "-Q >t"),
    
    FRAMES_ADJUSTED = round(FRAMES_ADJUSTED*5)/5
    # FRAMES_ADJUSTED = ifelse(FRAMES_ADJUSTED <= TIME_CUTOFF, 0, 1)
  ) %>%
  group_by(
    FRAMES_ADJUSTED,
    START_DELTA_REF_BIN,
    END_DELTA_REF_BIN,
    START_DELTA_QRY_BIN,
    END_DELTA_QRY_BIN,
    # UNIVERSAL_TRACK_ID
  ) %>% 
  summarize(
    REF_FUNCTION = mean(REF_FUNCTION),
    QRY_FUNCTION = mean(QRY_FUNCTION),
    
    SD_REFERENCE_TOTAL_INTENSITY = sd(REFERENCE_TOTAL_INTENSITY)/sqrt(n()),
    SD_QUERY_TOTAL_INTENSITY = sd(QUERY_TOTAL_INTENSITY)/sqrt(n()),
    
    REFERENCE_TOTAL_INTENSITY = mean(REFERENCE_TOTAL_INTENSITY),
    QUERY_TOTAL_INTENSITY = mean(QUERY_TOTAL_INTENSITY),
    
    N = NROW(unique(UNIVERSAL_TRACK_ID))
  ) %>% 
  group_by(
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
    )
  ) %>% 
  filter(
    N >= 5
  ) %>% 
  # mutate(
  #   START_DELTA_BIN = ifelse(START_DELTA_REF >= START_DELTA_QRY, "+R", "+Q"),
  #   END_DELTA_BIN = ifelse(END_DELTA_REF >= END_DELTA_QRY, "+R", "+Q")
  # ) %>% 
  mutate(
    # START_DELTA_REF = round(START_DELTA_REF*START_END_FACTOR)/START_END_FACTOR,
    # END_DELTA_REF = round(END_DELTA_REF*START_END_FACTOR)/START_END_FACTOR,
    # START_DELTA_QRY = round(START_DELTA_QRY*START_END_FACTOR)/START_END_FACTOR,
    # END_DELTA_QRY = round(END_DELTA_QRY*START_END_FACTOR)/START_END_FACTOR
  ) %>% 
  as.data.table() %>% 
ggplot(
  aes(
    # REF_FUNCTION,
    # QRY_FUNCTION,
    
    REFERENCE_TOTAL_INTENSITY,
    QUERY_TOTAL_INTENSITY,
    group = GROUP,
    # group = UNIVERSAL_TRACK_ID,
    # color = FRAMES_ADJUSTED
    color = GROUP,
    size = N
  )
  # aes(
  #   START_DELTA_REF,
  #   # END_DELTA_REF,
  #   START_DELTA_QRY,
  #   # END_DELTA_QRY
  # )
) + 
  # geom_ribbon(
  #   aes(
      # xmin = REFERENCE_TOTAL_INTENSITY -SD_REFERENCE_TOTAL_INTENSITY,
      # xmax = REFERENCE_TOTAL_INTENSITY +SD_REFERENCE_TOTAL_INTENSITY,
      # ymin = QUERY_TOTAL_INTENSITY - SD_QUERY_TOTAL_INTENSITY,
      # ymax = QUERY_TOTAL_INTENSITY + SD_QUERY_TOTAL_INTENSITY,
  #   )
  # ) +
  geom_path(
    # size = .1
  ) +
  # geom_point(
  #   aes(size = N)
  #   # size = .1
  # ) +
  # geom_errorbar(
  #   aes(
  #     ymin = QUERY_TOTAL_INTENSITY - SD_QUERY_TOTAL_INTENSITY,
  #     ymax = QUERY_TOTAL_INTENSITY + SD_QUERY_TOTAL_INTENSITY
  #   ),
  #   color = "black"
  # ) +
  # geom_errorbar(
  #   aes(
  #     xmin = REFERENCE_TOTAL_INTENSITY -SD_REFERENCE_TOTAL_INTENSITY,
  #     xmax = REFERENCE_TOTAL_INTENSITY +SD_REFERENCE_TOTAL_INTENSITY,
  #   ),
  #   color = "black"
  # ) +
  # geom_text(
  #   aes(x = 1, y = 1, label = N)
  # ) +
  # geom_hline(
  #   yintercept = 0
  # ) +
  # geom_vline(
  #   xintercept = 0
  # ) +
  # geom_jitter() +
  scale_color_viridis(
    option = "turbo",
    discrete =  T,
    guide = "none"
  ) +
  # facet_wrap(~UNIVERSAL_TRACK_ID)+
  # facet_grid(
  #   # START_DELTA_QRY_BIN~START_DELTA_REF_BIN
  #   # START_DELTA_BIN~END_DELTA_BIN
  #   START_DELTA_QRY_BIN+END_DELTA_QRY_BIN~START_DELTA_REF_BIN+END_DELTA_REF_BIN
  # ) +
  theme_classic() + 
  coord_fixed()

