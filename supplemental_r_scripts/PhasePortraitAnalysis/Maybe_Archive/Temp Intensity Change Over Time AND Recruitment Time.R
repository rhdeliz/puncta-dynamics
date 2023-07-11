OldEssential <- fread("/Users/u_deliz/Desktop/OldEssential.csv.gz")

Test <-
  NormStatTable %>% 
  filter(
    QUERY_PROTEIN == "NEMO",
    LIGAND_DENSITY_CAT==32,
    FRAMES_SINCE_LANDING <= 100
  ) %>% 
  as.data.table()

IntSummary <-
  NormStatTable %>% 
  filter(
    QUERY_PROTEIN == "NEMO",
    LIGAND_DENSITY_CAT==32,
    FRAMES_SINCE_LANDING <= 100
  ) %>% 
  group_by(
    REFERENCE_PROTEIN,
    QUERY_PROTEIN,
    FRAMES_ADJUSTED
  ) %>% 
  summarize(
    TIME_ADJUSTED = median(TIME_ADJUSTED),
    REFERENCE_TOTAL_INTENSITY = mean(REFERENCE_TOTAL_INTENSITY, na.rm = T),
    QUERY_TOTAL_INTENSITY = mean(QUERY_TOTAL_INTENSITY, na.rm = T)
  ) %>% 
  group_by(
    REFERENCE_PROTEIN,
    QUERY_PROTEIN
  ) %>% 
  mutate(
    SCALED_REFERENCE_TOTAL_INTENSITY = (REFERENCE_TOTAL_INTENSITY-min(REFERENCE_TOTAL_INTENSITY))/(max(REFERENCE_TOTAL_INTENSITY)-min(REFERENCE_TOTAL_INTENSITY)),
    SCALED_QUERY_TOTAL_INTENSITY = (QUERY_TOTAL_INTENSITY-min(QUERY_TOTAL_INTENSITY))/(max(QUERY_TOTAL_INTENSITY)-min(QUERY_TOTAL_INTENSITY))
  ) %>% 
  as.data.table()

ggplot(
  IntSummary,
  aes(
    REFERENCE_TOTAL_INTENSITY,
    QUERY_TOTAL_INTENSITY,
    color = TIME_ADJUSTED
  )
) +
  geom_path() +
  scale_color_viridis() +
  labs(
    x = "MyD88 Size",
    y = "NEMO Size",
    color = "Time (s)"
  ) +
  theme_classic()


ggplot() +
  geom_path(
    data = IntSummary,
    aes(
      y = REFERENCE_TOTAL_INTENSITY,
      x = TIME_ADJUSTED,
      color = REFERENCE_PROTEIN
    )
  ) +
  geom_path(
    data = IntSummary,
    aes(
      y = QUERY_TOTAL_INTENSITY,
      x = TIME_ADJUSTED,
      color = QUERY_PROTEIN
    )
  ) +
  scale_color_manual(
    values = c("green", "magenta")
  ) +
  labs(
    x = "MyD88 Size",
    y = "NEMO Size",
    color = "Time (s)"
  ) +
  theme_classic()

RecruitmentTime <-
  NormStatTable %>% 
  filter(
    QUERY_PROTEIN == "NEMO",
    LIGAND_DENSITY_CAT==32,
    FRAMES_SINCE_LANDING <= 100
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    REFERENCE_TOTAL_INTENSITY = rollmean(REFERENCE_TOTAL_INTENSITY, 9, align = c("right"), fill = NA),
    QUERY_TOTAL_INTENSITY = rollmean(QUERY_TOTAL_INTENSITY, 9, align = c("right"), fill = NA)
  ) %>%
  drop_na(
    REFERENCE_TOTAL_INTENSITY,
    QUERY_TOTAL_INTENSITY,
  ) %>% 
  mutate(
    SCALED_REFERENCE_TOTAL_INTENSITY = (REFERENCE_TOTAL_INTENSITY-min(REFERENCE_TOTAL_INTENSITY))/(max(REFERENCE_TOTAL_INTENSITY)-min(REFERENCE_TOTAL_INTENSITY)),
    SCALED_QUERY_TOTAL_INTENSITY = (QUERY_TOTAL_INTENSITY-min(QUERY_TOTAL_INTENSITY))/(max(QUERY_TOTAL_INTENSITY)-min(QUERY_TOTAL_INTENSITY))
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    REF_50 = which(SCALED_REFERENCE_TOTAL_INTENSITY >= 0.5)[1],
    REF_10 = which(SCALED_REFERENCE_TOTAL_INTENSITY >= 0.1)[1],
    QRY_50 = which(SCALED_QUERY_TOTAL_INTENSITY >= 0.5)[1]
  ) %>% 
  mutate(
    REF_50_TIME = TIME_ADJUSTED[REF_50],
    REF_10_TIME = TIME_ADJUSTED[REF_10],
    QRY_50_TIME = TIME_ADJUSTED[QRY_50],
    ID = 1:n()
  ) %>% 
  filter(
    ID == 1
  ) %>% 
  mutate(
    DELTA_TIME = QRY_50_TIME - REF_50_TIME,
    DELTA_TIME_10 = QRY_50_TIME - REF_10_TIME
  ) %>% 
  distinct() %>% 
  as.data.table()

ggplot(
  RecruitmentTime,
  aes(
    x = DELTA_TIME
  )
) +
  geom_density() +
  theme_classic()


RecruitmentTime


ggplot() +
  geom_path(
    data = IntSummary,
    aes(
      y = SCALED_REFERENCE_TOTAL_INTENSITY,
      x = TIME_ADJUSTED,
      color = REFERENCE_PROTEIN
    )
  ) +
  geom_path(
    data = IntSummary,
    aes(
      y = SCALED_QUERY_TOTAL_INTENSITY,
      x = TIME_ADJUSTED,
      color = QUERY_PROTEIN
    )
  ) +
  scale_color_manual(
    values = c("green", "magenta")
  ) +
  labs(
    x = "Scaled MyD88",
    y = "Scaled NEMO",
    color = "Time (s)"
  ) +
  theme_classic()


StoichiometrySummary <-
  NormStatTable %>% 
  filter(
    QUERY_PROTEIN == "NEMO",
    LIGAND_DENSITY_CAT==32,
    FRAMES_SINCE_LANDING <= 100
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    STOICHIOMETRY = REFERENCE_TOTAL_INTENSITY/(REFERENCE_TOTAL_INTENSITY+QUERY_TOTAL_INTENSITY)
  ) %>% 
  mutate(
    STARTING_STOICHIOMETRY = mean(STOICHIOMETRY[1:4], na.rm = T)
  ) %>% 
  mutate(
    ROUNDED_STARTING_STOICHIOMETRY = round(STARTING_STOICHIOMETRY*10)/10
  ) %>% 
  group_by(
    FRAMES_ADJUSTED,
    ROUNDED_STARTING_STOICHIOMETRY
  ) %>% 
  summarize(
    TIME_ADJUSTED = median(TIME_ADJUSTED),
    REFERENCE_TOTAL_INTENSITY = mean(REFERENCE_TOTAL_INTENSITY, na.rm = T),
    QUERY_TOTAL_INTENSITY = mean(QUERY_TOTAL_INTENSITY, na.rm = T)
  ) %>% 
  group_by(
    ROUNDED_STARTING_STOICHIOMETRY
  ) %>% 
  mutate(
    REFERENCE_TOTAL_INTENSITY = rollmean(REFERENCE_TOTAL_INTENSITY, 9, align = c("center"), fill = NA),
    QUERY_TOTAL_INTENSITY = rollmean(QUERY_TOTAL_INTENSITY, 9, align = c("center"), fill = NA),
    N = n()
  ) %>% 
  filter(
    N > 5
  ) %>% 
  drop_na() %>%
  as.data.table()


# 
# AddColumn <- function(x){
#   TempIntSummary <- IntSummary
#   TempIntSummary$ROUNDED_STARTING_STOICHIOMETRY <- x
#   return(TempIntSummary)
# }
# AvailStoich <- unique(StoichiometrySummary$ROUNDED_STARTING_STOICHIOMETRY)
# StoichiometrySummaryAll <- lapply(AvailStoich, AddColumn)
# StoichiometrySummaryAll <- rbindlist(StoichiometrySummaryAll)

ggplot(
) +
  # geom_path(
  #   data = StoichiometrySummaryAll,
  #   aes(
  #     REFERENCE_TOTAL_INTENSITY,
  #     QUERY_TOTAL_INTENSITY
  #   ),
  #   color = "black"
  # ) +
  geom_path(
    data = StoichiometrySummary,
    aes(
      REFERENCE_TOTAL_INTENSITY,
      QUERY_TOTAL_INTENSITY,
      color = TIME_ADJUSTED
    )
  ) +
  scale_color_viridis() +
  labs(
    x = "MyD88 Size",
    y = "NEMO Size",
    color = "Time (s)"
  ) +
  facet_wrap(
    ~ROUNDED_STARTING_STOICHIOMETRY
  ) +
  theme_classic()




LifetimeTable <-
  Test %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    REFERENCE_TOTAL_INTENSITY = rollmedian(REFERENCE_TOTAL_INTENSITY, 11, align = c("right"), fill = NA),
    QUERY_TOTAL_INTENSITY = rollmedian(QUERY_TOTAL_INTENSITY, 11, align = c("right"), fill = NA)
  ) %>%
  # tidyr::drop_na(
  #   REFERENCE_TOTAL_INTENSITY,
  #   QUERY_TOTAL_INTENSITY
  # ) %>%
  mutate(
    STOICHIOMETRY = REFERENCE_TOTAL_INTENSITY/(REFERENCE_TOTAL_INTENSITY+QUERY_TOTAL_INTENSITY),
    MAX_REFERENCE_TOTAL_INTENSITY = max(REFERENCE_TOTAL_INTENSITY, na.rm = T),
    MAX_QUERY_TOTAL_INTENSITY = max(QUERY_TOTAL_INTENSITY, na.rm = T),
    DELTA_QUERY_TOTAL_INTENSITY = median(DELTA_QUERY_TOTAL_INTENSITY),
    DELTA_REFERENCE_TOTAL_INTENSITY = median(DELTA_REFERENCE_TOTAL_INTENSITY)
  ) %>% 
  mutate(
    STARTING_REFERENCE_TOTAL_INTENSITY = median(REFERENCE_TOTAL_INTENSITY[1:4], na.rm = T),
    STARTING_STOICHIOMETRY = median(STOICHIOMETRY[1:4], na.rm = T)
  ) %>% 
  filter(
    REFERENCE_TOTAL_INTENSITY == MAX_REFERENCE_TOTAL_INTENSITY
  ) %>% 
  mutate(
    # ID = 1:n(),
    ROUNDED_TIME_ADJUSTED = round(TIME_ADJUSTED/25)*25,
    ROUNDED_DELTA = round(DELTA_QUERY_TOTAL_INTENSITY*10^5)/10^5
  ) %>%
  group_by(
    # PLOT_FACETS,
    ROUNDED_TIME_ADJUSTED,
    # ROUNDED_DELTA
  ) %>%
  summarize(
    N = n(),
    DELTA_QUERY_TOTAL_INTENSITY = median(DELTA_QUERY_TOTAL_INTENSITY, na.rm = T),
    DELTA_REFERENCE_TOTAL_INTENSITY = median(DELTA_REFERENCE_TOTAL_INTENSITY, na.rm = T),
    MAX_QUERY_TOTAL_INTENSITY = median(MAX_QUERY_TOTAL_INTENSITY, na.rm = T),
    MAX_REFERENCE_TOTAL_INTENSITY=  median(MAX_REFERENCE_TOTAL_INTENSITY, na.rm = T)
  ) %>%
  # group_by(
  #   ROUNDED_TIME_ADJUSTED
  # ) %>%
  # mutate(
  #   N = (N-min(N)+1)/(max(N) - min(N))
  # ) %>% 
  # filter(
  #   N >= 50
  # ) %>% 
  as.data.table()



ggplot(
  LifetimeTable,
  aes(
    x = ROUNDED_TIME_ADJUSTED,
    y = DELTA_QUERY_TOTAL_INTENSITY,
    size = N,
    color = MAX_REFERENCE_TOTAL_INTENSITY
  )
) +
  geom_point() +
  # scale_fill_viridis() +
  # scale_x_continuous(
  #   limits = c(-10^-3, 10^-3)
  # ) +
  theme_classic()


