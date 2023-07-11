RefFractionalComposition <-
  # Input table
  FilteredTable %>% 
  filter(
    # Filtering parameters, just like the phase portrait
    LIFETIME >= LIFETIME_MIN,
    FRAMES_ADJUSTED <= MAX_FRAMES,
    STARTING_NORMALIZED_INTENSITY <= MAX_STARTING_INTENSITY,
    # Limit intensities to avoid biasing fractional composition
    REFERENCE_TOTAL_INTENSITY <= 7.5,
    QUERY_TOTAL_INTENSITY <= 2.5
  ) %>% 
  group_by(
    LIGAND_DENSITY_CAT,
    REFERENCE_PROTEIN,
    QUERY_PROTEIN,
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    # Fractional composition
    TIME_CHANGE = lead(TIME, n=3L) - lag(TIME, n=3L),
    BIN = REFERENCE_TOTAL_INTENSITY/(REFERENCE_TOTAL_INTENSITY + QUERY_TOTAL_INTENSITY),
    # Growth rate of MyD88
    REF_GROWTH = lead(REFERENCE_TOTAL_INTENSITY, n=3L) - lag(REFERENCE_TOTAL_INTENSITY, n=3L),
    # Time adjustment so that it's per second
    REF_GROWTH = REF_GROWTH/TIME_CHANGE,
  )

QryFractionalComposition <-
  # Input table
  RefFractionalComposition %>% 
  group_by(
    LIGAND_DENSITY_CAT,
    REFERENCE_PROTEIN,
    QUERY_PROTEIN,
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    # Fractional composition
    BIN = QUERY_TOTAL_INTENSITY/(REFERENCE_TOTAL_INTENSITY + QUERY_TOTAL_INTENSITY),
    # Repear for IRAK (query protein_
    QRY_GROWTH = lead(QUERY_TOTAL_INTENSITY, n=3L) - lag(QUERY_TOTAL_INTENSITY, n=3L),
    QRY_GROWTH = QRY_GROWTH/TIME_CHANGE
  )

FractionalComposition <- bind_rows(RefFractionalComposition, QryFractionalComposition)

FractionalComposition <-
  FractionalComposition %>% 
  mutate(
    # Binned fractional composition to every 0.5x protein
    BIN = round(BIN*20)/20,
    # Just for plotting IRAK4 first, then IRAK1 (protein assembly order)
    QUERY_PROTEIN = factor(QUERY_PROTEIN, levels = c("IRAK4", "IRAK1"))
  ) %>% 
  group_by(
    # Fractional composition bins are grouped
    LIGAND_DENSITY_CAT,
    REFERENCE_PROTEIN,
    QUERY_PROTEIN,
    BIN
  ) %>% 
  summarize(
    # SD_REF_GROWTH = sd(REF_GROWTH, na.rm = T),
    # SD_QRY_GROWTH = sd(QRY_GROWTH, na.rm = T),
    # Mean growth of MyD88
    REF_GROWTH = mean(REF_GROWTH, na.rm = T),
    # Mean growth of IRAKs
    QRY_GROWTH = mean(QRY_GROWTH, na.rm = T),
    # Size N of bin
    N = n()
  ) %>% 
  group_by(
    LIGAND_DENSITY_CAT,
    REFERENCE_PROTEIN,
    QUERY_PROTEIN
  ) %>% 
  tidyr::drop_na() %>%
  group_by(
    LIGAND_DENSITY_CAT,
    REFERENCE_PROTEIN,
    QUERY_PROTEIN
  ) %>%
  filter(
    # Constraint to 0-15x proteins
    # A bit of flexibility since the calibration signal variance increases with protein size
    BIN >= 0,
    BIN <= 1.5
  ) %>%
  mutate(
    # SD_REF_GROWTH = SD_REF_GROWTH,##/max(REF_GROWTH),
    # SD_QRY_GROWTH = SD_QRY_GROWTH,#/max(QRY_GROWTH),
    # Scaled signals so that all peaks have a height of 1 a.u./s
    REF_GROWTH = REF_GROWTH/ max(REF_GROWTH),
    QRY_GROWTH = QRY_GROWTH/max(QRY_GROWTH),
    # # Smoothed using Savitzky-Golay filter with window n
    SMOOTH_REF_GROWTH = signal::sgolayfilt(REF_GROWTH, p = 2, n = 3),
    SMOOTH_REF_GROWTH = signal::sgolayfilt(SMOOTH_REF_GROWTH, p = 3, n = 7),
    SMOOTH_REF_GROWTH = signal::sgolayfilt(SMOOTH_REF_GROWTH, p = 2, n = 3),
    SMOOTH_QRY_GROWTH = signal::sgolayfilt(QRY_GROWTH, p = 2, n = 3),
    SMOOTH_QRY_GROWTH = signal::sgolayfilt(SMOOTH_QRY_GROWTH, p = 3, n = 7),
    SMOOTH_QRY_GROWTH = signal::sgolayfilt(SMOOTH_QRY_GROWTH, p = 2, n = 3)
  ) %>% 
  filter(
    # Constraint to 0-15x proteins
    # A bit of flexibility since the calibration signal variance increases with protein size
    BIN >= .25,
    BIN <= 1.33
  )

# Plot
ggplot() + 
  # Intercept at 0 to indicate growth/shrinking
  geom_hline(
    yintercept = 0
  ) +
  # Plots smoothed data as curve
  geom_path(
    data = FractionalComposition,
    aes(
      x = BIN,
      y = SMOOTH_REF_GROWTH,
      color = REFERENCE_PROTEIN
    )
  ) +
  geom_path(
    data = FractionalComposition,
    aes(
      x = BIN,
      y = SMOOTH_QRY_GROWTH,
      color = QUERY_PROTEIN
    )
  ) +
# geom_ribbon(
#   data = FractionalComposition,
#   aes(
#     x = BIN,
#     ymin = SD_REF_GROWTH - REF_GROWTH,
#     ymax = SD_REF_GROWTH + REF_GROWTH,
#     fill = REFERENCE_PROTEIN
#   ),
#   # size = .5,
#   alpha = .5
# )+
  # geom_ribbon(
  #   data = FractionalComposition,
  #   aes(
  #     x = BIN,
  #     ymin = SD_QRY_GROWTH - QRY_GROWTH,
  #     ymax = SD_QRY_GROWTH + QRY_GROWTH,
  #     fill = QUERY_PROTEIN
  #   ),
  #   # size = .5,
  #   alpha = .5
  # )+
  # Plots actual data as point
  geom_point(
    data = FractionalComposition,
    aes(
      x = BIN,
      y = REF_GROWTH,
      color = REFERENCE_PROTEIN
    ),
    # size = .5,
    alpha = .5
  ) +
  geom_point(
    data = FractionalComposition,
    aes(
      x = BIN,
      y = QRY_GROWTH,
      color = QUERY_PROTEIN
    ),
    # size = .5,
    alpha = .5
  ) +
  scale_color_brewer(
    palette = "Set1"
  ) +
  # Labels
  labs(
    x = "Fractional Composition",
    y = "Scaled Growth Rate (a.u./s)",
    color = "Protein"
  ) +
  # Make plots per cell line
  facet_grid(
    ~QUERY_PROTEIN
  ) +
  dark_theme_classic() +
  theme(
    legend.position = "bottom"
  )
  ggsave(
    "fraction.pdf",
    height = 9/2,
    width = 16/2
  )

