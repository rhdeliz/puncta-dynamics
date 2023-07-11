# Make query table
QryTrackTable <-
  Table %>% 
  mutate(
    NORMALIZED_INTENSITY = COMPLEMENTARY_NORMALIZED_INTENSITY_1,
    PROTEIN = COMPLEMENTARY_PROTEIN_1,
    NORMALIZED_CHANGEPOINT_INTENSITY = COMPLEMENTARY_CHANGEPOINT_NORMALIZED_INTENSITY_1
  )

# Add new rows
TrackTable <- rbindlist(list(QryTrackTable, Table))

TrackTable <- 
  TrackTable %>% 
  mutate(
    # Rescale time
    TIME = TIME - min(TIME),
    # Define protein order
    PROTEIN = factor(PROTEIN, levels = c(REFERENCE_PROTEIN, QUERY_PROTEIN))
  )
