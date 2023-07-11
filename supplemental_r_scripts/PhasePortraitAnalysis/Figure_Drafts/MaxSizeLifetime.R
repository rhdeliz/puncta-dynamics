setwd(OUTPUT_DIRECTORY)

# Table Name
SaveName <-
  paste0(
    "NormStatTable - ",
    "LeadLag ", LEAD_LAG,
    ".gz.parquet"
  )
SaveName <- gsub(" - \\.", "\\.", SaveName)
NormStatTable <- read_parquet(file.path(OUTPUT_DIRECTORY, SaveName))

NormStatTable$QUERY_PROTEIN = factor(NormStatTable$QUERY_PROTEIN, levels = PROTEIN_ORDER)

MaxInt_Life <-
  NormStatTable %>% 
  filter(
    FRAMES_SINCE_LANDING_CAT == 200,
    QUERY_TOTAL_INTENSITY >= 0,
    REFERENCE_TOTAL_INTENSITY >=0
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    REFERENCE_TOTAL_INTENSITY = max(REFERENCE_TOTAL_INTENSITY),
    QUERY_TOTAL_INTENSITY = max(QUERY_TOTAL_INTENSITY),
    LIFETIME = max(FRAMES_ADJUSTED)/FPS
  ) %>% 
  mutate(
    QUERY_TOTAL_INTENSITY = ifelse(QUERY_TOTAL_INTENSITY >= 2, 1, 0)
  ) %>% 
  filter(
    FRAMES_ADJUSTED == 5
  ) %>% 
  as_tibble() %>% 
  mutate(
    COLOCALIZED = ifelse(QUERY_TOTAL_INTENSITY == 1, "+ Coloc.", "- Coloc.")
  )

MaxInt_Life2 <- MaxInt_Life
MaxInt_Life2$COLOCALIZED <- "All"
MaxInt_Life2$QUERY_TOTAL_INTENSITY = 0
MaxInt_Life <- rbind(MaxInt_Life, MaxInt_Life2)
MaxInt_Life$COLOCALIZED <- factor(MaxInt_Life$COLOCALIZED, levels = c("All", "- Coloc.", "+ Coloc."))

Summary_MaxInt_Life <-
  MaxInt_Life %>% 
  group_by(
    QUERY_PROTEIN,
    QUERY_TOTAL_INTENSITY,
    COLOCALIZED,
    IMAGE
  ) %>% 
  summarize(
    REFERENCE_TOTAL_INTENSITY = median(REFERENCE_TOTAL_INTENSITY),
    LIFETIME = median(LIFETIME),
    N = n()
  ) %>% 
  group_by(
    QUERY_PROTEIN,
    QUERY_TOTAL_INTENSITY,
    COLOCALIZED
  ) %>% 
  summarize(
    MAD_N = mad(N),
    N = median(N),
    MAD_REFERENCE_TOTAL_INTENSITY = mad(REFERENCE_TOTAL_INTENSITY),
    MAD_LIFETIME = mad(LIFETIME),
    REFERENCE_TOTAL_INTENSITY = median(REFERENCE_TOTAL_INTENSITY),
    LIFETIME = median(LIFETIME)
  ) %>% 
  mutate(
    LIMIT_N = MAD_N + N,
    LIMIT_MAX = REFERENCE_TOTAL_INTENSITY + MAD_REFERENCE_TOTAL_INTENSITY,
    LIMIT_LIFE = MAD_LIFETIME + LIFETIME
  ) %>% 
  as.data.table()
  
ggplot(
  MaxInt_Life,
  aes(
    LIFETIME,
    REFERENCE_TOTAL_INTENSITY
  )
) +
  stat_density2d(geom="tile", aes(fill = ..ndensity..), adjust = 1.5, contour = FALSE) +
  geom_point(
    data = Summary_MaxInt_Life,
    aes(
      LIFETIME,
      REFERENCE_TOTAL_INTENSITY
    ),
    shape = 7,
    color = "cyan"
  ) +
  geom_hline(
    yintercept = 6
  ) +
  geom_vline(
    xintercept = 100
  ) +
  scale_fill_viridis(
    option = "plasma",
    direction = -1
  ) +
  scale_x_continuous(
    limits = c(0, quantile(MaxInt_Life$LIFETIME, .95))
  ) +
  scale_y_continuous(
    limits = c(0, quantile(MaxInt_Life$REFERENCE_TOTAL_INTENSITY, .95))
  ) +
  facet_grid(
    QUERY_PROTEIN~COLOCALIZED
  ) +
  labs(
    x = "Lifetime (s)",
    y = "Max MyD88 Size (x)",
    fill = "Scaled\nDensity"
  ) +
  theme_classic()+
  theme(
    strip.background = element_blank(),
    legend.position = "none"
  )

ggsave(
  "maxsize_life.pdf",
  width = 4,
  height = 9
)


PROTEINS = c("MyD88", "IRAK4", "IRAK1", "TRAF6",
             "HOIL1", "NEMO","TAB2",  "A20", "RelA")

ColorTable <- c("#000000", turbo(NROW(PROTEINS), direction = -1))[1:NROW(PROTEINS)]
names(ColorTable) <- PROTEINS


ggplot(
  Summary_MaxInt_Life,
  aes(
    QUERY_PROTEIN,
    REFERENCE_TOTAL_INTENSITY,
    fill = QUERY_PROTEIN
  )
) +
  geom_col(
    color = "black"
  ) +
  geom_errorbar(
    aes(
      ymin = REFERENCE_TOTAL_INTENSITY - MAD_REFERENCE_TOTAL_INTENSITY,
      ymax = REFERENCE_TOTAL_INTENSITY + MAD_REFERENCE_TOTAL_INTENSITY,
    ),
    width = 0.2,
    color = "black"
  ) +
  scale_x_discrete(
    limits=rev
  ) +
  scale_y_continuous(
    limits = c(0, max(Summary_MaxInt_Life$LIMIT_MAX))
  ) +
  scale_color_manual(
    values = ColorTable,
    # labels = PROTEIN_ORDER[2:NROW(PROTEIN_ORDER)],
    limits = force,
    # drop = F,
    guide = "none"
  ) +
  scale_fill_manual(
    values = ColorTable,
    # labels = PROTEIN_ORDER[2:NROW(PROTEIN_ORDER)],
    limits = force,
    # drop = F,
    guide = "none"
  ) +
  facet_wrap(
    ~COLOCALIZED,
    ncol = 3,
    scales = "free"
  ) +
  labs(
    x = "",
    y = "Max MyD88 Size",
  ) +
  theme_classic(
    base_size = 15
  ) +
  theme(
    strip.background = element_blank(),
    # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  ) +
  coord_flip()


ggsave(
  "max_coloc.pdf",
  width = 2.5*3.5,
  height = 2.5
)


ggplot(
  Summary_MaxInt_Life,
  aes(
    QUERY_PROTEIN,
    LIFETIME,
    fill = QUERY_PROTEIN
  )
) +
  geom_col(
    color = "black"
  ) +
  geom_errorbar(
    aes(
      ymin = LIFETIME - MAD_LIFETIME,
      ymax = LIFETIME + MAD_LIFETIME,
    ),
    width = 0.2,
    color = "black"
  ) +
  scale_x_discrete(
    limits=rev
  ) +
  scale_y_continuous(
    limits = c(0, max(Summary_MaxInt_Life$LIMIT_LIFE))
  ) +
  scale_color_manual(
    values = ColorTable,
    # labels = PROTEIN_ORDER[2:NROW(PROTEIN_ORDER)],
    limits = force,
    # drop = F,
    guide = "none"
  ) +
  scale_fill_manual(
    values = ColorTable,
    # labels = PROTEIN_ORDER[2:NROW(PROTEIN_ORDER)],
    limits = force,
    # drop = F,
    guide = "none"
  ) +
  facet_wrap(
    ~COLOCALIZED,
    ncol = 3,
    scales = "free"
  ) +
  labs(
    x = "",
    y = "Lifetime (s)",
  ) +
  theme_classic(
    base_size = 18
  ) +
  theme(
    strip.background = element_blank(),
    # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  ) +
  coord_flip()


ggsave(
  "life_coloc.pdf",
  width = 3*3.5,
  height = 3
)



FacetOrder <- paste(
  rep( c("Short-Lived", "Long-Lived"), 2),
  rep(c("Dim MyD88", "Bright MyD88"), each = 2),
  sep = ", ")


All_Summaries <-
  NormStatTable %>% 
  filter(
    FRAMES_SINCE_LANDING_CAT == 200,
    QUERY_TOTAL_INTENSITY >= 0,
    REFERENCE_TOTAL_INTENSITY >=0
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    REFERENCE_TOTAL_INTENSITY = max(REFERENCE_TOTAL_INTENSITY),
    QUERY_TOTAL_INTENSITY = max(QUERY_TOTAL_INTENSITY),
    LIFETIME = max(FRAMES_ADJUSTED)/FPS
  ) %>% 
  mutate(
    REFERENCE_TOTAL_INTENSITY = ifelse(REFERENCE_TOTAL_INTENSITY >= 6, 1, 0),
    QUERY_TOTAL_INTENSITY = ifelse(QUERY_TOTAL_INTENSITY >= 2, 1, 0),
    LIFETIME = ifelse(LIFETIME >= 100, 1, 0)
  ) %>% 
  filter(
    FRAMES_ADJUSTED == 5
  ) %>% 
  group_by(
    QUERY_PROTEIN,
    REFERENCE_TOTAL_INTENSITY,
    LIFETIME,
    IMAGE
  ) %>% 
  summarize(
    QUERY_TOTAL_INTENSITY = mean(QUERY_TOTAL_INTENSITY),
    N = n()
  ) %>% 
  as_tibble() %>% 
  group_by(
    QUERY_PROTEIN,
    IMAGE
  ) %>% 
  mutate(
    N = N/sum(N)
  ) %>% 
  as_tibble() %>% 
  group_by(
    QUERY_PROTEIN,
    REFERENCE_TOTAL_INTENSITY,
    LIFETIME
  ) %>% 
  summarize(
    MAD_N = mad(N),
    MAD_QUERY_TOTAL_INTENSITY = mad(QUERY_TOTAL_INTENSITY),
    QUERY_TOTAL_INTENSITY = mean(QUERY_TOTAL_INTENSITY),
    N = median(N)
  ) %>% 
  mutate(
    LIMIT = QUERY_TOTAL_INTENSITY + MAD_QUERY_TOTAL_INTENSITY
  ) %>% 
  group_by(
    QUERY_PROTEIN
  ) %>% 
  mutate(
    LIMIT_N = MAD_N + N,
    # N = N/sum(N),
    REFERENCE_TOTAL_INTENSITY = ifelse(REFERENCE_TOTAL_INTENSITY == 1, "Bright MyD88", "Dim MyD88"),
    LIFETIME = ifelse(LIFETIME==1, "Long-Lived", "Short-Lived")
  ) %>% 
  mutate(
    PLOT_FACET = paste(LIFETIME, REFERENCE_TOTAL_INTENSITY, sep = ", "),
    # N = round(N)
  ) %>%
  mutate(
    PLOT_FACET = factor(PLOT_FACET, levels = FacetOrder),
    # N = paste0( N, "% N")
  )


ggplot(
  All_Summaries,
  aes(
    x = QUERY_PROTEIN,
    y = QUERY_TOTAL_INTENSITY,
    fill = QUERY_PROTEIN
  )
) +
  geom_col(
    color = "black"
  ) +
  geom_errorbar(
    aes(
      ymin = QUERY_TOTAL_INTENSITY - MAD_QUERY_TOTAL_INTENSITY,
      ymax = QUERY_TOTAL_INTENSITY + MAD_QUERY_TOTAL_INTENSITY,
    ),
    width = 0.2,
    color = "black"
  ) +
  scale_x_discrete(
    limits=rev
  ) +
  scale_y_continuous(
    limits = c(0, max(All_Summaries$LIMIT)),
    labels = percent
  ) +
  scale_color_manual(
    values = ColorTable,
    # labels = PROTEIN_ORDER[2:NROW(PROTEIN_ORDER)],
    limits = force,
    # drop = F,
    guide = "none"
  ) +
  scale_fill_manual(
    values = ColorTable,
    # labels = PROTEIN_ORDER[2:NROW(PROTEIN_ORDER)],
    limits = force,
    # drop = F,
    guide = "none"
  ) +
  facet_wrap(
    ~PLOT_FACET,
    ncol = 1,
    scales = "free"
  ) +
  labs(
    x = "",
    y = "Colocalized",
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  ) +
  coord_flip()

ggsave(
  "colocalized_cats.pdf",
  width = 3,
  height = 7
)


ggplot(
  All_Summaries,
  aes(
    x = QUERY_PROTEIN,
    y = N,
    fill = QUERY_PROTEIN
  )
) +
  geom_col(
    color = "black"
  ) +
  geom_errorbar(
    aes(
      ymin = N - MAD_N,
      ymax = N + MAD_N,
    ),
    width = 0.2,
    color = "black"
  ) +
  scale_x_discrete(
    limits=rev
  ) +
  scale_y_continuous(
    limits = c(0, max(All_Summaries$LIMIT_N)),
    labels = percent
  ) +
  scale_color_manual(
    values = ColorTable,
    # labels = PROTEIN_ORDER[2:NROW(PROTEIN_ORDER)],
    limits = force,
    # drop = F,
    guide = "none"
  ) +
  scale_fill_manual(
    values = ColorTable,
    # labels = PROTEIN_ORDER[2:NROW(PROTEIN_ORDER)],
    limits = force,
    # drop = F,
    guide = "none"
  ) +
  facet_wrap(
    ~PLOT_FACET,
    ncol = 1,
    scales = "free"
  ) +
  labs(
    x = "",
    y = "N",
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  ) +
  coord_flip()

ggsave(
  "colocalized_n.pdf",
  width = 3,
  height = 7
)

