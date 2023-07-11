OUTPUT_DIRECTORY = "/Users/u_deliz/Desktop/figure_drafts"
setwd(OUTPUT_DIRECTORY)

REF_STEP_SIZE = 6
QRY_STEP_SIZE = 4
TIME_BLOCK = 3
MAX_SIZE = 1

# Transform max size to limits
MAX_REF_SIZE = MAX_SIZE*REF_STEP_SIZE
MAX_QRY_SIZE = MAX_SIZE*QRY_STEP_SIZE

Alluvium_Data <-
  NormStatTable %>% 
  filter(
    COHORT == "MyD88 IRAK4",
    FRAMES_ADJUSTED <= round(TIME_BLOCK/FPS)
  ) %>% 
  arrange(
    UNIVERSAL_TRACK_ID,
    FRAMES_ADJUSTED
  ) %>% 
  mutate(
    FRAMES_ADJUSTED = floor(FRAMES_ADJUSTED/TIME_BLOCK)*TIME_BLOCK
    # FRAMES_ADJUSTED = floor(FRAMES_ADJUSTED/max(FRAMES_ADJUSTED)*TIME_BLOCK)/TIME_BLOCK/FPS*max(FRAMES_ADJUSTED)
  ) %>% 
  mutate(
    FRAMES_ADJUSTED = round(FRAMES_ADJUSTED/FPS)
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    MAX_REFERENCE_TOTAL_INTENSITY = max(REFERENCE_TOTAL_INTENSITY),
    MAX_QUERY_TOTAL_INTENSITY = max(QUERY_TOTAL_INTENSITY)
  ) %>% 
  filter(
    MAX_REFERENCE_TOTAL_INTENSITY <= MAX_REF_SIZE,
    MAX_QUERY_TOTAL_INTENSITY <= MAX_QRY_SIZE
  ) %>% 
  as.data.table()

Img_Period <- round(1/Alluvium_Data$FPS[1])
RefProt <- Alluvium_Data$REFERENCE_PROTEIN[1]
QryProt <- Alluvium_Data$QUERY_PROTEIN[1]

Alluvium_Data <-
  Alluvium_Data %>% 
  group_by(
    REFERENCE_PROTEIN,
    QUERY_PROTEIN,
    UNIVERSAL_TRACK_ID,
    FRAMES_ADJUSTED
  ) %>% 
  summarize(
    REFERENCE_TOTAL_INTENSITY = median(REFERENCE_TOTAL_INTENSITY),
    QUERY_TOTAL_INTENSITY = median(QUERY_TOTAL_INTENSITY),

    MAX_REFERENCE_TOTAL_INTENSITY = mean(MAX_REFERENCE_TOTAL_INTENSITY),
    MAX_QUERY_TOTAL_INTENSITY = mean(MAX_QUERY_TOTAL_INTENSITY)
  ) %>% 
  mutate(
    REFERENCE_TOTAL_INTENSITY = round(REFERENCE_TOTAL_INTENSITY/REF_STEP_SIZE)*REF_STEP_SIZE,
    QUERY_TOTAL_INTENSITY = round(QUERY_TOTAL_INTENSITY/QRY_STEP_SIZE)*QRY_STEP_SIZE,
    
    MAX_REFERENCE_TOTAL_INTENSITY = round(MAX_REFERENCE_TOTAL_INTENSITY/REF_STEP_SIZE)*REF_STEP_SIZE,
    MAX_QUERY_TOTAL_INTENSITY = round(MAX_QUERY_TOTAL_INTENSITY/QRY_STEP_SIZE)*QRY_STEP_SIZE
  ) %>% 
  mutate(
    START_REFERENCE = REFERENCE_TOTAL_INTENSITY[1],
    START_QUERY = QUERY_TOTAL_INTENSITY[1],
    END_REFERENCE = REFERENCE_TOTAL_INTENSITY[n()],
    END_QUERY = QUERY_TOTAL_INTENSITY[n()]
  ) %>% 
  filter(
    START_REFERENCE <= 0,
    START_QUERY <= 0
  ) %>%
  mutate(
    INTENSITY = paste0(REFERENCE_TOTAL_INTENSITY, "x ", RefProt, " + ", QUERY_TOTAL_INTENSITY, "x ", QryProt),
    END_INTENSITY = paste0(END_REFERENCE, "x ", RefProt, " + ", END_QUERY, "x ", QryProt)
  ) %>% 
  as.data.table()
  
MaxRef <- max(Alluvium_Data$END_REFERENCE)
MaxQry <- max(Alluvium_Data$END_QUERY)

Alluvium_Data <-
  Alluvium_Data %>% 
  select(-c(
    REFERENCE_TOTAL_INTENSITY,
    QUERY_TOTAL_INTENSITY,
    MAX_REFERENCE_TOTAL_INTENSITY,
    MAX_QUERY_TOTAL_INTENSITY,
    END_REFERENCE,
    END_QUERY
  )) %>% 
  pivot_wider(
    names_from = FRAMES_ADJUSTED,
    id_cols = c(
      "UNIVERSAL_TRACK_ID",
      "END_INTENSITY"
    ),
    values_from = INTENSITY
  ) %>% 
  ungroup() %>% 
  select(-c(
    UNIVERSAL_TRACK_ID
  )) %>% 
  group_by_all(
  ) %>% 
  summarize(
    N = n()
  ) %>% 
  as.data.table()

Alluvium_Data[is.na(Alluvium_Data)] <- "Aborted"

MaxRefScale <- 0:(MaxRef/REF_STEP_SIZE)*REF_STEP_SIZE
MaxQryScale <- 0:(MaxQry/QRY_STEP_SIZE)*QRY_STEP_SIZE
MaxRefScale <- rep(MaxRefScale, each = NROW(MaxQryScale))
MaxQryScale <- rep(MaxQryScale, NROW(MaxRefScale))

ColorTable <- NULL
ColorTable$END_REFERENCE <- MaxRefScale
ColorTable$END_QUERY <- MaxQryScale 
ColorTable <- as.data.table(ColorTable)

ColorTable <-
  ColorTable %>% 
  mutate(
    SCALED_END_REFERENCE = END_REFERENCE/(max(END_REFERENCE)*100/85)+.15,
    SCALED_END_QUERY = END_QUERY/(max(END_QUERY)*100/85)+.15
  ) %>% 
  mutate(
    FILL_COLOR = rgb(1-SCALED_END_REFERENCE, 1-SCALED_END_QUERY, 1-SCALED_END_REFERENCE),
    END_INTENSITY = paste0(END_REFERENCE, "x ", RefProt, " + ", END_QUERY, "x ", QryProt)
  ) %>% 
  arrange(
    END_QUERY,
    END_REFERENCE
  ) %>% 
  select(
    END_INTENSITY,
    FILL_COLOR
  ) %>% 
  distinct() %>% 
  as.data.table()

ColorTable$END_INTENSITY <- factor(ColorTable$END_INTENSITY, levels = unique(ColorTable$END_INTENSITY))

order_int <- function(x){factor(x, levels = c(rev(as.character(ColorTable$END_INTENSITY)), "Aborted"))}

Alluvium_Data <-
  Alluvium_Data %>%
  mutate(
    across(
      names(Alluvium_Data)[names(Alluvium_Data)!="N"],
      # names(Alluvium_Data)[grepl(" s", names(Alluvium_Data))],
      order_int
    )
  ) %>%
  select(
    END_INTENSITY, 
    N,
    everything()
  ) %>% 
  as.data.table()

AxisTimes <- colnames(Alluvium_Data)[3:NCOL(Alluvium_Data)]
AxisTimes <- as.numeric(AxisTimes)
AxisTimes <- round(AxisTimes)
AxisTimes <- paste(AxisTimes, "s")

library(ggrepel)
library(ggalluvial)

# Fix column order
ggplot(
  data = Alluvium_Data,
  aes(
    y = N,
    axis1 = !!as.name(names(Alluvium_Data)[3]),
    axis2 = !!as.name(names(Alluvium_Data)[4]),
    axis3 = !!as.name(names(Alluvium_Data)[5]),
    axis4 = !!as.name(names(Alluvium_Data)[6]),
    # axis5 = !!as.name(names(Alluvium_Data)[7]),
    # axis6 = !!as.name(names(Alluvium_Data)[8]),
    # axis7 = !!as.name(names(Alluvium_Data)[9]),
    # axis8 = !!as.name(names(Alluvium_Data)[10]),
    # axis9 = !!as.name(names(Alluvium_Data)[11]),
    # axis10 = !!as.name(names(Alluvium_Data)[12])
  )
) +
  geom_alluvium(
    aes(fill = END_INTENSITY),
    alpha = 1, width = 1/6
  ) +
  geom_alluvium(
    aes(fill = END_INTENSITY),
    alpha = .5, width = 1/6
  ) +
  geom_stratum(width = 1/12, fill = "black", color = "grey", alpha = .5) +
  geom_text_repel(stat = "stratum", aes(label = after_stat(stratum)),
                  reverse = T, color = "black",
                  direction = "y", nudge_x = .5) +
  scale_fill_manual(
    values = ColorTable$FILL_COLOR,
    breaks = ColorTable$END_INTENSITY
  ) +
  scale_x_continuous(
    breaks = 1:NROW(AxisTimes),
    labels = AxisTimes
  ) +
  # scale_y_log10() +
  labs(
    x = "Cluster Time (s)",
    y = "N",
    # y = "log(N)",
    fill = "Ending\nStoichiometry"
  ) +
  theme_classic(
    base_size = 20
  ) +
  # guides(fill=guide_legend(nrow=2)) +
  theme(
    legend.position = "right",
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank()
  )

ggsave(
  "alluvium.pdf",
  height = 2.5*2,
  width = 8.66*1.5
)



# Fix column order
ggplot(
  data = Alluvium_Data,
  aes(
    y = N,
    axis1 = !!as.name(names(Alluvium_Data)[3]),
    axis2 = !!as.name(names(Alluvium_Data)[4]),
    axis3 = !!as.name(names(Alluvium_Data)[5]),
    axis4 = !!as.name(names(Alluvium_Data)[6]),
    # axis5 = !!as.name(names(Alluvium_Data)[7]),
    # axis6 = !!as.name(names(Alluvium_Data)[8]),
    # axis7 = !!as.name(names(Alluvium_Data)[9]),
    # axis8 = !!as.name(names(Alluvium_Data)[10]),
    # axis9 = !!as.name(names(Alluvium_Data)[11]),
    # axis10 = !!as.name(names(Alluvium_Data)[12])
  )
) +
  geom_alluvium(
    aes(fill = END_INTENSITY),
    alpha = 1, width = 1/6
  ) +
  geom_alluvium(
    aes(fill = END_INTENSITY),
    alpha = .5, width = 1/6
  ) +
  geom_stratum(width = 1/12, fill = "black", color = "grey", alpha = .5) +
  geom_text_repel(stat = "stratum", aes(label = after_stat(stratum)),
                  reverse = T, color = "black",
                  direction = "y", nudge_x = .5) +
  scale_fill_manual(
    values = ColorTable$FILL_COLOR,
    breaks = ColorTable$END_INTENSITY
  ) +
  scale_x_continuous(
    breaks = 1:NROW(AxisTimes),
    labels = AxisTimes
  ) +
  scale_y_log10() +
  labs(
    x = "Cluster Time (s)",
    y = "log(N)",
    fill = "Ending\nStoichiometry"
  ) +
  theme_classic(
    base_size = 20
  ) +
  # guides(fill=guide_legend(nrow=2)) +
  theme(
    legend.position = "right",
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank()
  )

ggsave(
  "alluvium_log.pdf",
  height = 2.5*2,
  width = 8.66*1.5
)