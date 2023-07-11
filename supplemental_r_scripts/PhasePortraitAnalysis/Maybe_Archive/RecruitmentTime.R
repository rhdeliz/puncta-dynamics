# StatTableBackup <- fread("/Users/u_deliz/Desktop/PhasePortraitAnalysis/Normalized StatTable - LeadLag 5 - StepSize 0.05.csv.gz")

StatTable <-
  StatTableBackup %>% 
  filter(
    !grepl("rid", COHORT),
    LIGAND_DENSITY_CAT == 32
    # COHORT == "MyD88 TRAF6"
  ) %>% 
  as.data.table()

# Recruitment Time
Normalization <-
  StatTable %>%
  select(
    COHORT, LIGAND_DENSITY_CAT,
    IMAGE,
    MIN_REFERENCE,
    MAX_REFERENCE,
    MIN_QUERY,
    MAX_QUERY
  ) %>%
  distinct() %>%
  group_by(
    COHORT, LIGAND_DENSITY_CAT
  ) %>% 
  mutate(
    RANK_MAX_QUERY = rank(MAX_QUERY),
    RANK_MAX_REFERENCE = rank(MAX_REFERENCE),
    SUM = RANK_MAX_QUERY + RANK_MAX_REFERENCE,
    EQUALIZER = median(SUM),
    ID = 1:n(),
    TEST = which.min(abs(SUM - EQUALIZER))
  ) %>%
  filter(
    ID == TEST
  ) %>%
  as.data.table()

NormalizedTable <-
  StatTable %>%
  mutate(
    MIN_REFERENCE = Normalization$MIN_REFERENCE[1],
    MAX_REFERENCE = Normalization$MAX_REFERENCE[1],
    MIN_QUERY = Normalization$MIN_QUERY[1],
    MAX_QUERY = Normalization$MAX_QUERY[1],
  ) %>%
  group_by(
    COHORT, LIGAND_DENSITY_CAT
  ) %>% 
  mutate(
    REFERENCE_TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY*MAX_REFERENCE+MIN_REFERENCE,
    DELTA_REFERENCE_TOTAL_INTENSITY = DELTA_REFERENCE_TOTAL_INTENSITY*MAX_REFERENCE,
    ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY = ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY*MAX_REFERENCE,
    
    QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY*MAX_REFERENCE+MIN_REFERENCE,
    DELTA_QUERY_TOTAL_INTENSITY = DELTA_QUERY_TOTAL_INTENSITY*MAX_REFERENCE,
    ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY = ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY*MAX_REFERENCE
  ) %>%
  as.data.table()

# Recruitment time
TimeTable <-
  NormalizedTable %>% 
  drop_na(REFERENCE_TOTAL_INTENSITY, QUERY_TOTAL_INTENSITY) %>%
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    LIFETIME = max(FRAMES_ADJUSTED)
  ) %>% 
  filter(
    LIFETIME >= 41
  ) %>% 
  arrange(
    UNIVERSAL_TRACK_ID,
    TIME_ADJUSTED 
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    TIME_ADJUSTED = TIME_ADJUSTED - min(TIME_ADJUSTED),
    TIME_ADJUSTED = round(TIME_ADJUSTED*FPS)/FPS
  ) %>% 
  drop_na(REFERENCE_TOTAL_INTENSITY, QUERY_TOTAL_INTENSITY) %>%
  group_by(
    COHORT,
    LIGAND_DENSITY_CAT,
    REFERENCE_PROTEIN,
    QUERY_PROTEIN,
    TIME_ADJUSTED
  ) %>%
  summarize(
    N = n(),
    REFERENCE_TOTAL_INTENSITY = median(REFERENCE_TOTAL_INTENSITY, na.rm = T),
    QUERY_TOTAL_INTENSITY = median(QUERY_TOTAL_INTENSITY, na.rm = T)
  ) %>%
  filter(
    N >=5
  ) %>%
  group_by(
    COHORT,
    LIGAND_DENSITY_CAT,
    REFERENCE_PROTEIN,
    QUERY_PROTEIN
  ) %>%
  mutate(
    REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 1, n = 11),
    QUERY_TOTAL_INTENSITY =  signal::sgolayfilt(QUERY_TOTAL_INTENSITY, p = 1, n = 11),
    # REFERENCE_TOTAL_INTENSITY = roll_mean(REFERENCE_TOTAL_INTENSITY, 7, fill = NA, align = "left"),
    # REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 3, n = 5),
    # QUERY_TOTAL_INTENSITY = roll_mean(QUERY_TOTAL_INTENSITY, 7, fill = NA, align = "left")
  ) %>% 
  drop_na(REFERENCE_TOTAL_INTENSITY, QUERY_TOTAL_INTENSITY) %>%
  group_by(
    COHORT,
    LIGAND_DENSITY_CAT,
    REFERENCE_PROTEIN,
    QUERY_PROTEIN
  ) %>%
  mutate(
    REFERENCE_TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY-min(REFERENCE_TOTAL_INTENSITY),
    QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY-min(QUERY_TOTAL_INTENSITY),
    TIME_ADJUSTED = TIME_ADJUSTED - min(TIME_ADJUSTED)
  ) %>% 
  mutate(
    REFERENCE_TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY/max(REFERENCE_TOTAL_INTENSITY),
    QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY/max(QUERY_TOTAL_INTENSITY)
  ) %>%
  as.data.table()
  
Order <-
  NormalizedTable %>% 
  drop_na(REFERENCE_TOTAL_INTENSITY, QUERY_TOTAL_INTENSITY) %>% 
  # filter(
  #   REFERENCE_TOTAL_INTENSITY >= 1.5 || QUERY_TOTAL_INTENSITY >= 1.5
  # ) %>%
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>%
  mutate(
    ALIGNMENT_TIME = TIME_ADJUSTED[which.min(REFERENCE_TOTAL_INTENSITY >= 1)],
    RECRUITMENT_TIME = TIME_ADJUSTED[which.min(QUERY_TOTAL_INTENSITY >= 2)]
  ) %>%
  group_by(
    COHORT, LIGAND_DENSITY_CAT
  ) %>% 
  summarize(
    RECRUITMENT_TIME = median(RECRUITMENT_TIME, na.rm = T)
  ) %>% 
  arrange(
    RECRUITMENT_TIME
  ) %>% 
  as.data.table()

StackedOrder <-
  TimeTable %>% 
  arrange(
    COHORT, LIGAND_DENSITY_CAT,
    TIME_ADJUSTED
  ) %>% 
  group_by(
    COHORT, LIGAND_DENSITY_CAT, QUERY_PROTEIN
  ) %>% 
  summarize(
    RECRUITMENT_TIME = TIME_ADJUSTED[which.min(QUERY_TOTAL_INTENSITY<=.49)],
    QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY[which.min(QUERY_TOTAL_INTENSITY<=.49)]
  ) %>% 
  # summarize(
  #   RECRUITMENT_TIME = TIME_ADJUSTED[which.min(QUERY_TOTAL_INTENSITY == .5)]
  # ) %>%
  arrange(
    RECRUITMENT_TIME
  ) %>% 
  as.data.table()

TimeTable$COHORT <- factor(TimeTable$COHORT, levels = unique(StackedOrder$COHORT))
StackedOrder$COHORT <- factor(StackedOrder$COHORT, levels = unique(StackedOrder$COHORT))

ggplot(
  TimeTable
) +
  geom_point(
    data = StackedOrder,
    aes(
      x = RECRUITMENT_TIME,
      y = QUERY_TOTAL_INTENSITY
    ),
    size = 3
  ) +
  geom_path(
    aes(
      x = TIME_ADJUSTED,
      y = REFERENCE_TOTAL_INTENSITY,
      color = REFERENCE_PROTEIN,
      group = COHORT
    )
  ) +
  geom_path(
    aes(
      x = TIME_ADJUSTED,
      y = QUERY_TOTAL_INTENSITY,
      color = QUERY_PROTEIN,
      group = COHORT
    )
  ) +
  scale_color_manual(
    labels = c("MyD88", unique(TimeTable$QUERY_PROTEIN)),
    breaks = c("MyD88", unique(TimeTable$QUERY_PROTEIN)),
    values = c("darkgreen", rep("purple", NROW(unique(TimeTable$QUERY_PROTEIN))))
  ) +
  labs(
    x = "Time (s)",
    # y = "Normalized Intensity",
    y = "Scaled Intensity",
    color = "Protein"
  ) +
  facet_wrap(
    ~COHORT,
    scales = "free_y",
    nrow = 2
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank()
  )
  
ggsave(
  # Save vector image
  file.path(OUTPUT_DIRECTORY, "ScaledRecruitment.pdf"),
  # file.path(OUTPUT_DIRECTORY, "Recruitment.pdf"),
  # height = 4.76*1.25,
  # width = 11.5*1.25
  height = 4.76,
  width = 11.5
)

TimeTable$QUERY_PROTEIN <- factor(TimeTable$QUERY_PROTEIN, levels = unique(StackedOrder$QUERY_PROTEIN))


ggplot(
  TimeTable
) +
  # geom_hline(
  #   yintercept = .45
  # ) +
  # geom_vline(
  #   xintercept = StackedOrder$RECRUITMENT_TIME
  # ) +
  geom_path(
    aes(
      x = TIME_ADJUSTED,
      y = QUERY_TOTAL_INTENSITY,
      color = QUERY_PROTEIN,
      group = COHORT
    )
  ) +
  labs(
    x = "Time (s)",
    y = "Scaled Intensity",
    color = "Protein"
  ) +
  scale_color_viridis(
    discrete = T,
    option = "turbo"
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank()
  )


Test <-
  NormalizedTable %>% 
  mutate(
    REFERENCE_TOTAL_INTENSITY = round(REFERENCE_TOTAL_INTENSITY*2)/2,
    QUERY_TOTAL_INTENSITY = round(QUERY_TOTAL_INTENSITY*2)/2
  ) %>% 
  group_by(
    REFERENCE_TOTAL_INTENSITY,
    QUERY_TOTAL_INTENSITY
  ) %>% 
  summarize(
    TIME_ADJUSTED = mean(TIME_ADJUSTED),
    N = n()
  ) %>% 
  filter(
    N >= 5
  ) %>% 
  arrange(
    TIME_ADJUSTED
  ) %>% 
  as.data.table()

ggplot(
  Test,
  aes(
    x = TIME_ADJUSTED,
    y = REFERENCE_TOTAL_INTENSITY
  )
) +
  geom_path()

