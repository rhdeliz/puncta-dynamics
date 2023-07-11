library(ggbeeswarm)

ExpandFx <- function(x){
  TempTable <-
    NormStatTable %>% 
    mutate(
      # FRAMES_SINCE_LANDING_CAT = x,
      FRAMES_SINCE_LANDING_CAT = ceiling(FRAMES_SINCE_LANDING/50)*50,
    ) %>% 
    filter(
      FRAMES_SINCE_LANDING <= x
      # FRAMES_SINCE_LANDING_CAT == x
    ) %>% 
    mutate(
      FRAMES_SINCE_LANDING_CAT = x
    ) %>% 
    as.data.table()
  
  return(TempTable)
}
ExpandedNormStatTable <- lapply(unique(NormStatTable$FRAMES_SINCE_LANDING_CAT), ExpandFx)
ExpandedNormStatTable <- rbindlist(ExpandedNormStatTable)


CellTable <- 
  ExpandedNormStatTable %>% 
  group_by(
    QUERY_PROTEIN,
    # ROUNDED_QUERY_TOTAL_INTENSITY
    # FRAMES_SINCE_LANDING_CAT
  ) %>% 
  mutate(
    # QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY/quantile(QUERY_TOTAL_INTENSITY, .75)
  ) %>% 
  filter(
    REFERENCE_TOTAL_INTENSITY >= STEP_SIZE,
    QUERY_TOTAL_INTENSITY >= STEP_SIZE,
    # FRAMES_SINCE_LANDING_CAT >= 100,
    # FRAMES_SINCE_LANDING_CAT <= 200
  ) %>% 
  group_by(
    QUERY_PROTEIN,
    IMAGE,
    CELL
  ) %>%
  summarize(
    N = n(),
    dR_dQ = cor(DELTA_REFERENCE_TOTAL_INTENSITY, DELTA_QUERY_TOTAL_INTENSITY),#, method = "spearman"),
    dR_Q = cor(DELTA_REFERENCE_TOTAL_INTENSITY, QUERY_TOTAL_INTENSITY)#, method = "spearman")
  ) %>% 
  filter(
    N >= 100
  ) %>% 
  as.data.table()

ImageTable <-
  CellTable %>% 
  group_by(
    QUERY_PROTEIN,
    IMAGE
  ) %>%
  summarize(
    dR_Q = median(dR_Q)
  ) %>% 
  as.data.table()


ProteinTable <-
  CellTable %>% 
  group_by(
    QUERY_PROTEIN
  ) %>%
  summarize(
    dR_Q = median(dR_Q)
  ) %>% 
  arrange(
    -dR_Q
  ) %>% 
  as.data.table()

CellTable$QUERY_PROTEIN = factor(CellTable$QUERY_PROTEIN, levels = ProteinTable$QUERY_PROTEIN)
ImageTable$QUERY_PROTEIN = factor(ImageTable$QUERY_PROTEIN, levels = ProteinTable$QUERY_PROTEIN)
ProteinTable$QUERY_PROTEIN = factor(ProteinTable$QUERY_PROTEIN, levels = ProteinTable$QUERY_PROTEIN)

ggplot() +
  geom_violin(
    data = CellTable,
    aes(
      QUERY_PROTEIN,
      dR_Q,
    ),
    fill = "white",
    scale = "width",
    width = .67
  ) +
  geom_hline(
    yintercept = 0,
    color = "red",
    size = 2
  ) +
  geom_violin(
    data = CellTable,
    aes(
      QUERY_PROTEIN,
      dR_Q,
      fill = QUERY_PROTEIN
    ),
    alpha = .15,
    scale = "width",
    width = .67
  ) +
  geom_beeswarm(
    data = ImageTable,
    aes(
      QUERY_PROTEIN,
      dR_Q,
      fill = QUERY_PROTEIN
    ),
    shape = 21,
    size = 4,
    cex = 2
  ) +
  geom_crossbar(
    data = ProteinTable,
    aes(
      QUERY_PROTEIN,
      dR_Q,
      ymin = dR_Q,
      ymax = dR_Q
    ),
    width = .5
  ) +
  scale_fill_manual(
    values = ColorTable,
    limits = force,
    guide = "none"
  ) +
  labs(
    x = "Cell Line",
    y = "R (∆, mol.)"
  ) +
  theme_classic(
    base_size = 20
  ) +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank()
  )


SaveName <- "correlation.pdf"

ggsave(
  file.path(OUTPUT_DIRECTORY, SaveName),
  height = 3*1,
  width = 3*sqrt(2)*2,
  device = cairo_pdf
)
