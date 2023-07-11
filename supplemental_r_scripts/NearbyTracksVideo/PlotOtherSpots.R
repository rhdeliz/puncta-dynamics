setwd(RETURN_TO_DIRECTORY)

plot_other_spots <-
  plot_spot %>%
  mutate(
    FRAME = t,
    REF_POSITION_X = POSITION_X,
    REF_POSITION_Y = POSITION_Y
  ) %>%
  select(
    FRAME,
    REF_POSITION_X,
    REF_POSITION_Y
  )

plot_other_spots <- merge(plot_tracks, plot_other_spots, by = "FRAME")

plot_other_spots <-
  plot_other_spots %>%
  filter(
    FRAME == FOCUS_TIME
  ) %>%
  mutate(
    POSITION_X = POSITION_X - REF_POSITION_X,
    POSITION_Y = POSITION_Y - REF_POSITION_Y,
    x = POSITION_X/PIXEL_SIZE,
    y = POSITION_Y/PIXEL_SIZE,
    t = FRAME
  ) %>%
  select(
    t,
    x,
    y
  ) %>%
  filter(
    x >= -20,
    y >= -20,
    x <= 20,
    y <= 20
  )
