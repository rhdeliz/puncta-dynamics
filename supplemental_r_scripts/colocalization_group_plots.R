Table <- 
  ColocalizationNeeded %>% 
  filter(
    cell %in% Cells[3]
  )
Tables <- mclapply(Table$table, fread, mc.cores = detectCores(logical = F))
Tables <- rbindlist(Tables)

Groups <- mclapply(CostTablePaths, fread)
Groups <- rbindlist(Groups)

Test <- merge(Tables, Groups, by = "UNIVERSAL_TRACK_ID")


Test2 <-
  Test %>%
  ungroup() %>% 
  mutate(
    TRACK_ID = group_indices(., IMAGE, CELL, PROTEIN, COLOCALIZATION_GROUP, UNIVERSAL_TRACK_ID)
  ) %>% 
  group_by(
    IMAGE, CELL, PROTEIN, COLOCALIZATION_GROUP
  ) %>% 
  mutate(
    TRACK_ID = TRACK_ID - min(TRACK_ID)
  )


library(ggnewscale)
library(ggdark)
library(ggforce)

ggplot() +
  geom_path(
    data = Test2,
    aes(
      x = POSITION_X/CALIBRATION_UM,
      y = POSITION_Y/CALIBRATION_UM,
      color = PROTEIN,
      group = UNIVERSAL_TRACK_ID
    )
  ) +
  scale_color_brewer(
    palette = "Set1"
  ) +
  new_scale_color()+
  geom_path(
    data = Test2,
    aes(
      x = POSITION_X/CALIBRATION_UM,
      y = POSITION_Y/CALIBRATION_UM,
      color = COLOCALIZATION_GROUP,
      group = UNIVERSAL_TRACK_ID
    ),
    size = 0.5,
    alpha = .5
  ) +
  scale_color_distiller(
    palette = "Spectral"
  )+
  dark_theme_classic() +
  theme(
    legend.position = "none"
  ) +
  coord_fixed()


Pages <- unique(Test2$COLOCALIZATION_GROUP)
PagePlot <- function(PageX){
  Plot <-
    ggplot(
      data = Test2 %>% filter(COLOCALIZATION_GROUP==PageX),
      aes(
        x = POSITION_X/CALIBRATION_UM,
        y = POSITION_Y/CALIBRATION_UM,
        color = PROTEIN,
        alpha = TRACK_ID,
        group = UNIVERSAL_TRACK_ID,
        label = FRAME
      )
    ) +
    geom_path() +
    geom_text() +
    scale_color_brewer(
      palette = "Set1"
    )+
    scale_alpha(
      range = c(1, .25)
    ) +
    labs(
      title = paste("Colocalization Group", PageX),
      x = "pixels",
      y = "pixels",
      alpha = "Track ID",
      color = "Protein"
    ) +
    dark_theme_classic() +
    theme(
      legend.position = "right"
    ) +
    coord_fixed()
}
Plots <- mclapply(Pages, PagePlot)

pdf("plots.pdf", height = 9, width = 16)
for(page in Pages){
  print(Plots[page])
}
dev.off()
