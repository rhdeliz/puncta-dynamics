R2Loop <- function(x){
  ExpandedNormStatTable %>% 
    group_by(
      QUERY_PROTEIN
    ) %>% 
    filter(
      REFERENCE_TOTAL_INTENSITY >= STEP_SIZE,
      QUERY_TOTAL_INTENSITY <= x,
      FRAMES_SINCE_LANDING_CAT == 100
    ) %>% 
    group_by(
      QUERY_PROTEIN,
      ROUNDED_QUERY_TOTAL_INTENSITY
    ) %>% 
    mutate(
      N = n()
    ) %>% 
    filter(
      N >= 175
    ) %>% 
    mutate(
      ROUNDED_QUERY_TOTAL_INTENSITY = x
    ) %>% 
    group_by(
      # FRAMES_SINCE_LANDING_CAT,
      QUERY_PROTEIN,
      ROUNDED_QUERY_TOTAL_INTENSITY
    ) %>% 
    summarize(
      N = n(),
      dR_Q = cor(DELTA_REFERENCE_TOTAL_INTENSITY, QUERY_TOTAL_INTENSITY)#, method = "spearman")
    ) %>% 
    filter(
      N >= 175
    ) %>% 
    as.data.table()
}
R2Table <- mclapply(0:33*3, R2Loop)
R2Table <- rbindlist(R2Table)

ggplot(
  R2Table,
  aes(
    ROUNDED_QUERY_TOTAL_INTENSITY,
    dR_Q,
    group = QUERY_PROTEIN,
    color = QUERY_PROTEIN
    # group = FRAMES_SINCE_LANDING_CAT*4,
    # color = FRAMES_SINCE_LANDING_CAT*4
  )
) +
  geom_hline(
    yintercept = 0,
    size = 2
  ) +
  geom_path() +
  labs(
    x = "Query Size",
    y = "R^2 dMyD88 vs Query Size",
    color = "Cell t (s)"
  ) +
  scale_color_manual(
    values = ColorTable,
    # labels = PROTEIN_ORDER[2:NROW(PROTEIN_ORDER)],
    limits = force,
    # drop = F,
    # guide = "none"
  ) +
  # scale_color_viridis() +
  # facet_wrap(
  #   ~QUERY_PROTEIN
  #   # ~FRAMES_SINCE_LANDING_CAT
  # ) +
  theme_classic()

