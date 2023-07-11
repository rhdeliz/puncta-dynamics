GrowthTable %>% 
  filter(
    UNIVERSAL_TRACK_ID %in% unique(GrowthTable$UNIVERSAL_TRACK_ID)[1:20],
    # FRAMES_ADJUSTED > TIME_CUTOFF
  ) %>% 
  arrange(
    UNIVERSAL_TRACK_ID,
    FRAMES_ADJUSTED
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    RATIO = DELTA_REFERENCE - DELTA_QUERY
  ) %>% 
  mutate(
    START_RATIO = ifelse(FRAMES_ADJUSTED <= TIME_CUTOFF, RATIO, 0),
    END_RATIO = ifelse(FRAMES_ADJUSTED > TIME_CUTOFF, RATIO, 0),
    # N = ifelse(FRAMES_ADJUSTED <= TIME_CUTOFF, 1, 0)
  ) %>% 
  mutate(
    START_RATIO = sum(START_RATIO),
    END_RATIO = sum(END_RATIO)
    # START_RATIO = sum(START_RATIO)/sum(N),
    # END_RATIO = sum(END_RATIO)/(n()-sum(N))
  ) %>% 
  # mutate(
  #   RATIO = RATIO/max(abs(RATIO)),
  #   START_RATIO = START_RATIO/max(abs(START_RATIO)),
  #   END_RATIO = END_RATIO/max(abs(END_RATIO))
  # ) %>%
  mutate(
    START_RATIO_BIN = ifelse(START_RATIO >= 0, "+R <t", "+Q <t"),
    END_RATIO_BIN = ifelse(END_RATIO >= 0, "+R >t", "+Q >t"),
    FRAMES_ADJUSTED = round(FRAMES_ADJUSTED*5)/5
  ) %>%
  group_by(
    START_RATIO_BIN,
    END_RATIO_BIN,
    FRAMES_ADJUSTED
  ) %>% 
  summarize(
    REFERENCE_TOTAL_INTENSITY = mean(REFERENCE_TOTAL_INTENSITY),
    QUERY_TOTAL_INTENSITY = mean(QUERY_TOTAL_INTENSITY)
  ) %>% 
  as.data.table() %>% 
  ggplot(
    aes(
      REFERENCE_TOTAL_INTENSITY,
      QUERY_TOTAL_INTENSITY,
      # group = UNIVERSAL_TRACK_ID,
      color = FRAMES_ADJUSTED
    )
  ) + 
  geom_path(
    # size = .1
  ) +
  # geom_point(
  #   # aes(size = N)
  #   # size = .1
  # ) +
  geom_hline(
    yintercept = 0
  ) +
  # geom_vline(
  #   xintercept = 0
  # ) +
  geom_point() +
  # scale_color_distiller(
  #   palette = "PiYG",
  #   direction = 1
  # ) +
  scale_color_viridis() +
  # facet_wrap(
  #   ~UNIVERSAL_TRACK_ID
  # ) +
  facet_grid(
    START_RATIO_BIN~END_RATIO_BIN
    # START_DELTA_QRY_BIN+END_DELTA_QRY_BIN~START_DELTA_REF_BIN+END_DELTA_REF_BIN
  ) +
  theme_classic() +
  coord_fixed()
