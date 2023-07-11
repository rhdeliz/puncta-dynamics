# 1D Phase Portrait

library(data.table)
library(dplyr)
library(ggplot2)
library(ggdark)
library(scales)

PhasePortrait <- fread("/Users/u_deliz/Desktop/lin_statTable.csv.gz")

PhasePortrait <-
  PhasePortrait %>% 
  ungroup() %>% 
  mutate(
    QUERY_PROTEIN = factor(QUERY_PROTEIN, levels = c("IRAK4", "IRAK1")),
    IMAGENUMBER = group_indices(., LIGAND_DENSITY_CAT, REFERENCE_PROTEIN, QUERY_PROTEIN, IMAGE)
  ) %>% 
  group_by(
    LIGAND_DENSITY_CAT, REFERENCE_PROTEIN, QUERY_PROTEIN
  ) %>% 
  mutate(
    IMAGENUMBER = IMAGENUMBER - min(IMAGENUMBER) + 1
  )

MyD88 <-
  PhasePortrait %>% 
  mutate(
    ASSEMBLED = REFERENCE_TOTAL_INTENSITY/6
  ) %>% 
  filter(
    ASSEMBLED <= 1
  ) %>% 
  group_by(
    LIGAND_DENSITY_CAT,
    QUERY_PROTEIN,
    IMAGENUMBER,
    ASSEMBLED
  ) %>% 
  summarize(
    DELTA_REFERENCE_TOTAL_INTENSITY = median(DELTA_REFERENCE_TOTAL_INTENSITY, na.rm = T),
    DELTA_QUERY_TOTAL_INTENSITY = median(DELTA_QUERY_TOTAL_INTENSITY, na.rm = T)
  ) %>% 
  group_by(
    LIGAND_DENSITY_CAT,
    QUERY_PROTEIN,
    IMAGENUMBER
  ) %>% 
  mutate(
    DELTA_REFERENCE_TOTAL_INTENSITY = DELTA_REFERENCE_TOTAL_INTENSITY/max(abs(DELTA_REFERENCE_TOTAL_INTENSITY)),
    DELTA_QUERY_TOTAL_INTENSITY = DELTA_QUERY_TOTAL_INTENSITY/max(abs(DELTA_QUERY_TOTAL_INTENSITY))
  ) %>% 
  group_by(
    LIGAND_DENSITY_CAT,
    QUERY_PROTEIN
  ) %>% 
  mutate(
    IMAGENUMBER = -IMAGENUMBER/max(IMAGENUMBER)
  )

Myddosome <-
  PhasePortrait %>% 
  filter(
    REFERENCE_TOTAL_INTENSITY <= 6,
    QUERY_TOTAL_INTENSITY <= 4
  ) %>% 
  mutate(
    ASSEMBLED = REFERENCE_TOTAL_INTENSITY*.1 + QUERY_TOTAL_INTENSITY*.1,
    ASSEMBLED = round(ASSEMBLED*20)/20
  ) %>% 
  group_by(
    LIGAND_DENSITY_CAT,
    QUERY_PROTEIN,
    IMAGENUMBER,
    ASSEMBLED
  ) %>% 
  summarize(
    N = n(),
    DELTA_REFERENCE_TOTAL_INTENSITY = median(DELTA_REFERENCE_TOTAL_INTENSITY, na.rm = T),
    DELTA_QUERY_TOTAL_INTENSITY = median(DELTA_QUERY_TOTAL_INTENSITY, na.rm = T)
  ) %>% 
  group_by(
    LIGAND_DENSITY_CAT,
    QUERY_PROTEIN,
    IMAGENUMBER
  ) %>% 
  mutate(
    DELTA_REFERENCE_TOTAL_INTENSITY = DELTA_REFERENCE_TOTAL_INTENSITY/max(abs(DELTA_REFERENCE_TOTAL_INTENSITY)),
    DELTA_QUERY_TOTAL_INTENSITY = DELTA_QUERY_TOTAL_INTENSITY/max(abs(DELTA_QUERY_TOTAL_INTENSITY))
  ) %>% 
  group_by(
    LIGAND_DENSITY_CAT,
    QUERY_PROTEIN
  ) %>% 
  mutate(
    IMAGENUMBER = -IMAGENUMBER/max(IMAGENUMBER)
  )

ggplot(
  MyD88
) + 
  geom_hline(
    yintercept = 0
  ) +
  geom_line(
    aes(
      x = ASSEMBLED,
      y = DELTA_REFERENCE_TOTAL_INTENSITY,
      color = "MyD88",
      alpha = -IMAGENUMBER,
      group = -IMAGENUMBER
    ),
  ) +
  geom_line(
    aes(
      x = ASSEMBLED,
      y = DELTA_QUERY_TOTAL_INTENSITY,
      color = QUERY_PROTEIN,
      alpha = -IMAGENUMBER,
      group = -IMAGENUMBER
    )
  ) +
  labs(
    x = "Assembled MyD88",
    y = "Relative Growth",
    color = "Protein",
    alpha = "Replicate"
  ) +
  scale_x_continuous(
    labels = percent
  ) +
  scale_color_manual(
    values = c("magenta", "pink", "green")
  ) +
  facet_grid(
    QUERY_PROTEIN~LIGAND_DENSITY_CAT,
    scales = "free"
  ) +
  dark_theme_classic() +
  theme(
    legend.position = "bottom"
  )


ggsave(
  "~/Desktop/Myd88Replicates.pdf",
  height = 4.76,
  width = 11.5
)






# Replicates combined

MyD88 <-
  PhasePortrait %>% 
  mutate(
    ASSEMBLED = REFERENCE_TOTAL_INTENSITY/6
  ) %>% 
  filter(
    # ASSEMBLED <= 1,
    # DELTA_REFERENCE_TOTAL_INTENSITY >= 0,
    # DELTA_QUERY_TOTAL_INTENSITY >= 0
  ) %>% 
  group_by(
    LIGAND_DENSITY_CAT,
    QUERY_PROTEIN,
    ASSEMBLED
  ) %>% 
  summarize(
    DELTA_REFERENCE_TOTAL_INTENSITY = median(DELTA_REFERENCE_TOTAL_INTENSITY, na.rm = T),
    DELTA_QUERY_TOTAL_INTENSITY = median(DELTA_QUERY_TOTAL_INTENSITY, na.rm = T)
  ) %>% 
  group_by(
    LIGAND_DENSITY_CAT,
    QUERY_PROTEIN
  ) %>% 
  mutate(
    DELTA_REFERENCE_TOTAL_INTENSITY = DELTA_REFERENCE_TOTAL_INTENSITY/max(abs(DELTA_REFERENCE_TOTAL_INTENSITY)),
    DELTA_QUERY_TOTAL_INTENSITY = DELTA_QUERY_TOTAL_INTENSITY/max(abs(DELTA_QUERY_TOTAL_INTENSITY))
  )

Myddosome <-
  PhasePortrait %>% 
  filter(
    REFERENCE_TOTAL_INTENSITY <= 6,
    QUERY_TOTAL_INTENSITY <= 4,
    DELTA_REFERENCE_TOTAL_INTENSITY >= 0,
    DELTA_QUERY_TOTAL_INTENSITY >= 0
  ) %>% 
  mutate(
    ASSEMBLED = ORIG_REFERENCE_TOTAL_INTENSITY*.1 + ORIG_QUERY_TOTAL_INTENSITY*.1,
    ASSEMBLED = round(ASSEMBLED*40)/40
  ) %>% 
  group_by(
    LIGAND_DENSITY_CAT,
    QUERY_PROTEIN,
    ASSEMBLED
  ) %>% 
  summarize(
    N = n(),
    DELTA_REFERENCE_TOTAL_INTENSITY = median(DELTA_REFERENCE_TOTAL_INTENSITY, na.rm = T),
    DELTA_QUERY_TOTAL_INTENSITY = median(DELTA_QUERY_TOTAL_INTENSITY, na.rm = T)
  ) %>% 
  group_by(
    LIGAND_DENSITY_CAT,
    QUERY_PROTEIN
  ) %>% 
  mutate(
    DELTA_REFERENCE_TOTAL_INTENSITY = DELTA_REFERENCE_TOTAL_INTENSITY/max(abs(DELTA_REFERENCE_TOTAL_INTENSITY)),
    DELTA_QUERY_TOTAL_INTENSITY = DELTA_QUERY_TOTAL_INTENSITY/max(abs(DELTA_QUERY_TOTAL_INTENSITY))
  )

ggplot(
  Myddosome
) + 
  geom_hline(
    yintercept = 0
  ) +
  # geom_point(
  #   data = MyddosomePoints,
  #   aes(
  #     x = ASSEMBLED,
  #     y = DELTA_REFERENCE_TOTAL_INTENSITY,
  #     color = "MyD88"
  #   ),
  #   alpha = 0.25
  # ) +
  # geom_point(
  #   data = MyddosomePoints,
  #   aes(
  #     x = ASSEMBLED,
  #     y = DELTA_QUERY_TOTAL_INTENSITY,
  #     color = QUERY_PROTEIN
  #   ),
  #   alpha = 0.25
  # ) +
  geom_line(
    aes(
      x = ASSEMBLED,
      y = DELTA_REFERENCE_TOTAL_INTENSITY,
      color = "MyD88"
    ),
  ) +
  geom_line(
    aes(
      x = ASSEMBLED,
      y = DELTA_QUERY_TOTAL_INTENSITY,
      color = QUERY_PROTEIN
    )
  ) +
  labs(
    x = "Assembled Myddosome",
    y = "Relative Growth",
    color = "Protein"
  ) +
  scale_x_continuous(
    labels = percent
  ) +
  scale_color_manual(
    values = c("magenta", "pink", "green")
  ) +
  facet_grid(
    QUERY_PROTEIN~LIGAND_DENSITY_CAT,
    scales = "free"
  ) +
  dark_theme_classic() +
  theme(
    legend.position = "bottom"
  )

ggsave(
  "~/Desktop/MyddosomeQuarterAssembly.pdf",
  height = 4.76,
  width = 11.5
)

