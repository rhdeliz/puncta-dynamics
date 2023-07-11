BoxColor <- NULL
BoxColor$t <- unique(plot_img_spot$t)
BoxColor$t <- rep(BoxColor$t, 3)
BoxColor$protein <- rep(c(REFERENCE_PROTEIN, QUERY_PROTEIN, "Merge"), each = NROW(unique(plot_img_spot$t)))
BoxColor$protein <- factor(BoxColor$protein, levels = c(REFERENCE_PROTEIN, QUERY_PROTEIN, "Merge"))
BoxColor$color <- rep(c("green", "magenta", "white"), each = NROW(unique(plot_img_spot$t)))
BoxColor$x <- 1
BoxColor$y <- 1
BoxColor <- as.data.table(BoxColor)

ggplot(
  plot_img_spot
) +
  with_blur(
    geom_raster(
      aes(  
        x = x,
        y = y,
        fill = color
      ),
      interpolate = F
    ),
    sigma = 0
  ) +
  geom_point(
    data = BoxColor,
    aes(
      x = x,
      y = y,
      color = color
    ),
    shape = 22,
    size = BOX_SIZE
  ) +
  facet_grid(
    protein~real_t,
    switch = "y"
  ) +
  scale_color_identity( 
    guide = "none"
  ) +
  scale_fill_identity(
    guide = "none"
  ) +
  labs(
    title = paste(round((BOX_SIZE*2+1)*PIXEL_SIZE, 2), "µm ×", round((BOX_SIZE*2+1)*PIXEL_SIZE, 2), "µm")
  ) +
  dark_theme_classic(
    base_size = 18
  ) +
  coord_equal() +
  theme(
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank()
  )

setwd(OUTPUT_PATH)

ggsave(
  "Montage.pdf",
  height = 3*2,
  width = 6
)
