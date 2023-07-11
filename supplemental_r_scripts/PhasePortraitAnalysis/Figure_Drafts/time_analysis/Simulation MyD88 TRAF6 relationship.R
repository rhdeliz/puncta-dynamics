



SampleTrack <- NULL
SampleTrack$REFERENCE_TOTAL_INTENSITY = c(0:10/10, 10:0/10, c(0:10*0), 0:10/20, 10:0/20+.25)
SampleTrack$QUERY_TOTAL_INTENSITY = c(c(0:10*0), 0:10/10, rep(1, 11), 10:0/10, 0:10/10)

SampleTrack <-
  as_tibble(SampleTrack) %>% 
  ungroup() %>% 
  mutate(
    FRAMES_ADJUSTED = 1:n(),
    REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 1, n = 7),
    QUERY_TOTAL_INTENSITY = signal::sgolayfilt(QUERY_TOTAL_INTENSITY, p = 1, n = 7)
  ) %>% 
  mutate(
    REFERENCE_TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY - min(REFERENCE_TOTAL_INTENSITY),
    QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY - min(QUERY_TOTAL_INTENSITY)
  ) %>%
  mutate(
    REFERENCE_TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY/max(REFERENCE_TOTAL_INTENSITY),
    QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY/max(QUERY_TOTAL_INTENSITY)
  ) %>% 
  mutate(
    DELTA_REFERENCE = lead(REFERENCE_TOTAL_INTENSITY, 1) - lag(REFERENCE_TOTAL_INTENSITY,1),
    DELTA_QUERY = lead(QUERY_TOTAL_INTENSITY,1) - lag(QUERY_TOTAL_INTENSITY,1)
  ) %>% 
  drop_na(
    DELTA_REFERENCE,
    DELTA_QUERY
  ) %>% 
  arrange(
    FRAMES_ADJUSTED
  ) %>% 
  mutate(
    CUMULATIVE_DELTA_QUERY = cumsum(DELTA_QUERY),
    MEAN_CUMULATIVE_DELTA_QUERY = cummean(DELTA_QUERY),
  ) %>% 
  arrange(
    FRAMES_ADJUSTED
  ) %>% 
  mutate(
    QRY_FUNCTION = cummean(QUERY_TOTAL_INTENSITY) - QUERY_TOTAL_INTENSITY[1]
  ) %>% 
  arrange(
    -FRAMES_ADJUSTED
  ) %>% 
  mutate(
    REF_FUNCTION = cummean(REFERENCE_TOTAL_INTENSITY - REFERENCE_TOTAL_INTENSITY[n()]),
  ) %>% 
  arrange(
    FRAMES_ADJUSTED
  ) %>% 
  as.data.table()





SampleTrack %>% 
  ggplot()  +
  geom_hline(
    yintercept = 0
  ) +
  geom_path(
    aes(
      FRAMES_ADJUSTED,
      REFERENCE_TOTAL_INTENSITY
    ),
    color = "darkgreen"
  ) +
  geom_path(
    aes(
      FRAMES_ADJUSTED,
      QUERY_TOTAL_INTENSITY
    ),
    color = "darkmagenta"
  ) +
  geom_path(
    aes(
      FRAMES_ADJUSTED,
      REF_FUNCTION
    ),
    color = "green"
  ) +
  geom_path(
    aes(
      FRAMES_ADJUSTED,
      DELTA_QUERY
      # QRY_FUNCTION
    ),
    color = "magenta"
  ) +
  # geom_path(
  #   aes(
  #     FRAMES_ADJUSTED,
  #     DELTA_QUERY,
  #     QUERY_TOTAL_INTENSITY,
  #     # CUMULATIVE_DELTA_QUERY,
  # MEAN_CUMULATIVE_DELTA_QUERY,
  #   ),
  #   color = "black"
  # ) +
  # scale_fill_viridis() +
theme_classic()

