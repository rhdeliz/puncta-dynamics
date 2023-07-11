ClusterTime <- 
  NormStatTable %>% 
  group_by(
    COHORT,
    LIGAND_DENSITY_CAT
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>%
  mutate(
    N = n()
  ) %>%
  filter(
    N >= 16
  ) %>%
  as_tibble() %>% 
  group_split(
    COHORT
  )

THRESHOLD = 0.5

# Parallelize for speed
Analyze <- function(df){
  
  TempClusterRecruitmentTime <-
    df %>% 
    group_by(
      UNIVERSAL_TRACK_ID
    ) %>%
    mutate(
      REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 1, n = 7),
      QUERY_TOTAL_INTENSITY = signal::sgolayfilt(QUERY_TOTAL_INTENSITY, p = 1, n = 7),
    ) %>%
    mutate(
      REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 3, n = 15),
      QUERY_TOTAL_INTENSITY = signal::sgolayfilt(QUERY_TOTAL_INTENSITY, p = 3, n = 15)
    ) %>%
    group_by(
      UNIVERSAL_TRACK_ID
    ) %>%
    mutate(
      REFERENCE_TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY - min(REFERENCE_TOTAL_INTENSITY),
      QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY - min(QUERY_TOTAL_INTENSITY)
    ) %>%
    mutate(
      MAX_REFERENCE_TOTAL_INTENSITY = max(REFERENCE_TOTAL_INTENSITY),
      MAX_QUERY_TOTAL_INTENSITY = max(QUERY_TOTAL_INTENSITY)
    ) %>% 
    filter(
      MAX_REFERENCE_TOTAL_INTENSITY >= 1,
      MAX_QUERY_TOTAL_INTENSITY >= 1
    ) %>%
    # filter(
    #   MAX_REFERENCE_TOTAL_INTENSITY >= 6 |
    #   MAX_QUERY_TOTAL_INTENSITY >= 6
    # ) %>%
    mutate(
      REFERENCE_TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY/MAX_REFERENCE_TOTAL_INTENSITY,
      QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY/MAX_QUERY_TOTAL_INTENSITY
    ) %>% 
    mutate(
      REFERENCE_TOTAL_INTENSITY_50 = ifelse(REFERENCE_TOTAL_INTENSITY >= THRESHOLD, 1, 0),
      QUERY_TOTAL_INTENSITY_50 = ifelse(QUERY_TOTAL_INTENSITY >= THRESHOLD, 1, 0)
      
      # REFERENCE_TOTAL_INTENSITY_50 = abs(REFERENCE_TOTAL_INTENSITY-THRESHOLD),
      # QUERY_TOTAL_INTENSITY_50 = abs(QUERY_TOTAL_INTENSITY - THRESHOLD)
    ) %>% 
    group_by(
      COHORT,
      QUERY_PROTEIN,
      REFERENCE_PROTEIN,
      UNIVERSAL_TRACK_ID
    ) %>%
    summarize(
      REFERENCE_TIME_ADJUSTED = min(case_when(REFERENCE_TOTAL_INTENSITY_50 == 1 ~ TIME_ADJUSTED), na.rm = T),
      QUERY_TIME_ADJUSTED = min(case_when(QUERY_TOTAL_INTENSITY_50 == 1 ~ TIME_ADJUSTED), na.rm = T)
      # REFERENCE_TIME_ADJUSTED = min(case_when(REFERENCE_TOTAL_INTENSITY_50 == min(REFERENCE_TOTAL_INTENSITY_50) ~ TIME_ADJUSTED), na.rm = T),
      # QUERY_TIME_ADJUSTED = min(case_when(QUERY_TOTAL_INTENSITY_50 == min(QUERY_TOTAL_INTENSITY_50) ~ TIME_ADJUSTED), na.rm = T),
    ) %>%
    mutate(
      DELTA_TIME_ADJUSTED = QUERY_TIME_ADJUSTED - REFERENCE_TIME_ADJUSTED,
      COHORT = factor(COHORT, levels = paste(USE_REFERENCE_PROTEIN, PROTEIN_ORDER)),
      QUERY_PROTEIN = factor(QUERY_PROTEIN, levels = PROTEIN_ORDER),
      REFERENCE_PROTEIN = factor(REFERENCE_PROTEIN, levels = PROTEIN_ORDER)
    ) %>% 
    as.data.table()
  
  return(TempClusterRecruitmentTime)
}
ClusterRecruitmentTime <- mclapply(ClusterTime, Analyze)
ClusterRecruitmentTime <- rbindlist(ClusterRecruitmentTime)


ClusterRecruitmentTimeBar <-
  ClusterRecruitmentTime %>% 
  group_by(
    COHORT, QUERY_PROTEIN
  ) %>% 
  summarize(
    N = n(),
    MAD_DELTA_TIME_ADJUSTED = mad(DELTA_TIME_ADJUSTED),
    DELTA_TIME_ADJUSTED = median(DELTA_TIME_ADJUSTED),
    REFERENCE_TIME_ADJUSTED = median(REFERENCE_TIME_ADJUSTED),
    QUERY_TIME_ADJUSTED = median(QUERY_TIME_ADJUSTED)
  ) %>% 
  arrange(
    DELTA_TIME_ADJUSTED
  ) %>% 
  mutate(
    SEM = MAD_DELTA_TIME_ADJUSTED/sqrt(N)
  ) %>% 
  as_tibble()

# By protein
ggplot(ClusterRecruitmentTime) +
  geom_violin(
    aes(
      x = QUERY_PROTEIN,
      y = QUERY_TIME_ADJUSTED,
      fill = QUERY_PROTEIN
    ),
    scale = "width"
  ) +
  geom_boxplot(
    aes(
      x = QUERY_PROTEIN,
      y = QUERY_TIME_ADJUSTED,
    ),
    fill = "white",
    # alpha = 0,
    width = .1
  ) +
  scale_fill_manual(
    values = ColorTable,
    limits = force,
    guide = "none"
  ) +
  labs(
    y = "Cluster Time\nto Mid-Intensity (s)"
  ) +
  theme_classic() +
  theme(
    axis.title.x=element_blank(),
    legend.position = "none"
  )


ggsave(
  file.path(OUTPUT_DIRECTORY, "cluster_recruitment_time_all.pdf"),
  height = 2,
  width = 9.5/2.25
)


# Time difference
ggplot() +
  geom_density(
    data = ClusterRecruitmentTime,
    aes(
      x = DELTA_TIME_ADJUSTED,
      y = ..ndensity..,
      group = COHORT,
      fill  = QUERY_PROTEIN
    ),
    adjust=2
  ) +
  geom_vline(
    data = ClusterRecruitmentTimeBar,
    aes(
      xintercept = DELTA_TIME_ADJUSTED,
    ),
    color = "black"
  ) +
  scale_color_manual(
    values = ColorTable,
    limits = force,
    guide = "none"
  ) +
  scale_fill_manual(
    values = ColorTable,
    limits = force,
    guide = "none"
  ) +
  facet_wrap(
    ~QUERY_PROTEIN,
    nrow = 1
  ) +
  labs(
    x = "Difference in Cluster Time to Mid-Intensity (s)",
    y = "Scaled Density"
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  )

ggsave(
  file.path(OUTPUT_DIRECTORY, "cluster_recruitment_time.pdf"),
  height = 2,
  width = 9.5
)


# By protein
ggplot(ClusterRecruitmentTime) +
  geom_violin(
    aes(
      x = QUERY_PROTEIN,
      y = DELTA_TIME_ADJUSTED,
      fill  = QUERY_PROTEIN
    ),
    scale = "width"
  ) +
  geom_boxplot(
    aes(
      x = QUERY_PROTEIN,
      y = DELTA_TIME_ADJUSTED
    ),
    fill = "white",
    # alpha = 0,
    width = .1
  ) +
  geom_hline(
    yintercept = 0
  ) +
  scale_fill_manual(
    values = ColorTable,
    limits = force,
    guide = "none"
  ) +
  labs(
    y = "Difference in Cluster Time\nto Mid-Intensity (s)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x=element_blank()
  )

ggsave(
  file.path(OUTPUT_DIRECTORY, "cluster_recruitment_time_violin.pdf"),
  height = 2,
  width = 9.5/2.25
)


ClusterRecruitmentTimeBar %>% 
  ungroup() %>% 
  arrange(
    DELTA_TIME_ADJUSTED
  ) %>% 
  mutate(
    QUERY_PROTEIN = factor(QUERY_PROTEIN, levels = rev(unique(QUERY_PROTEIN)))
  ) %>% 
  as.data.table() %>% 
ggplot(
  aes(
    x = QUERY_PROTEIN,
    y = DELTA_TIME_ADJUSTED,
    # color = QUERY_PROTEIN,
    fill = QUERY_PROTEIN
  )
) +
  geom_col(
    color = "black"
  ) +
  geom_errorbar(
    aes(
      x = QUERY_PROTEIN,
      ymin = (DELTA_TIME_ADJUSTED - SEM),
      ymax = (DELTA_TIME_ADJUSTED + SEM),
    ),
    color = "black",
    width = 0.5
  ) +
  scale_color_manual(
    values = ColorTable,
    limits = force,
    guide = "none"
  ) +
  scale_fill_manual(
    values = ColorTable,
    limits = force,
    guide = "none"
  ) +
  labs(
    y = "Difference in Cluster Time to Mid-Intensity (s)",
    x = "Cell Line"
  ) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.title.y = element_blank()
  )


ggsave(
  file.path(OUTPUT_DIRECTORY, "cluster_recruitment_time_bar.pdf"),
  height = 2,
  width = 9.5/2.25
)


ClusterRecruitmentTimeBar %>% 
  ungroup() %>% 
  arrange(
    MAD_DELTA_TIME_ADJUSTED
  ) %>% 
  mutate(
    QUERY_PROTEIN = factor(QUERY_PROTEIN, levels = rev(unique(QUERY_PROTEIN)))
  ) %>% 
  as.data.table() %>% 
  ggplot(
    aes(
      x = QUERY_PROTEIN,
      y = MAD_DELTA_TIME_ADJUSTED,
      fill = QUERY_PROTEIN
    )
  ) +
  geom_col(
    color = "black"
  ) +
  scale_color_manual(
    values = ColorTable,
    limits = force,
    guide = "none"
  ) +
  scale_fill_manual(
    values = ColorTable,
    limits = force,
    guide = "none"
  ) +
  labs(
    y = "Variability of Cluster Time to Mid-Intensity (s)",
    x = "Cell Line"
  ) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.title.y = element_blank()
  )


ggsave(
  file.path(OUTPUT_DIRECTORY, "cluster_recruitment_time_bar_variance.pdf"),
  height = 2,
  width = 9.5/2.25
)
