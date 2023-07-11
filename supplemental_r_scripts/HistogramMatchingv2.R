pacman::p_load(ijtiff, parallel, dtplyr, data.table, dplyr)

IMAGE_PATH = "/Users/u_deliz/Desktop/WellB02_ChannelX RFP Cy5,X GFP,X DAPI_Seq0000_X DAPI.tif"
IMAGE_BIN_SIZE = 2
BITS = 16
# INTENSITY_BIN_SIZE = 25


tictoc::tic()
# Binning matrix function
mat <- function(n, r) {
  suppressWarnings(matrix(c(rep(1, r), rep(0, n)), n, n/r))
}

# Import image function
ReadImageFx <- function(FrameX){
  # Start table for collating data
  Image <- NULL
  # Read image
  GFP_Image <- read_tif(IMAGE_PATH, FrameX, msg = F)
  GFP_Image <- matrix(GFP_Image, nrow = NROW(GFP_Image))
  
  # Bin, if applicable
  if(IMAGE_BIN_SIZE > 1){
    x <- mat(NROW(GFP_Image), IMAGE_BIN_SIZE)
    GFP_Image <- t(x) %*% GFP_Image %*% x
    GFP_Image <- GFP_Image/IMAGE_BIN_SIZE^2
  }
  
  # Make table from pixels
  Image$GFP <- as.vector(GFP_Image)
  Image$FRAME <- FrameX
  
  # # Scale images
  Image$GFP_Norm <- Image$GFP - min(Image$GFP)
  Image$GFP_Norm <- Image$GFP_Norm/max(Image$GFP_Norm)*(2^BITS-1)
  
  # Get ratio of bright to dim
  # To be used later for filtering dark and dusty images
  IMG_HIGH = quantile(Image$GFP_Norm, .9)
  IMG_LOW = quantile(Image$GFP_Norm, .1)
  Image$RATIO = IMG_HIGH/IMG_LOW
  
  # Export
  Image <- as.data.table(Image)
  return(Image)
}

# Parallelize
FRAME_COUNT = count_frames(IMAGE_PATH)[1]
ALlImages <- mclapply(1:FRAME_COUNT, ReadImageFx, mc.cores = detectCores(logical = T))
ALlImages <- rbindlist(ALlImages)

# Eliminate dark and dusty images
ALlImages <-
  ALlImages %>%
  filter(
    RATIO >= 2
  ) %>% 
  as.data.table()

# Find reference frame
FRAME_COUNT = NROW(unique(ALlImages$FRAME))
MATRIX_SIZE <- NROW(ALlImages)/FRAME_COUNT
SUBSET_INDEX <- (1:FRAME_COUNT)*MATRIX_SIZE
RefFrame <- quantile(ALlImages$RATIO[SUBSET_INDEX], .9)
RefFrame <- abs(ALlImages$RATIO[SUBSET_INDEX] - RefFrame)
RefFrame <- which(RefFrame == min(RefFrame))
RefFrame <- min(RefFrame)

# Extract reference frame data
RefImage <-  
  ALlImages %>% 
  filter(
    FRAME == RefFrame
  ) %>% 
  mutate(
    ID = 1:n()
  ) %>% 
  arrange(
    GFP
  ) %>% 
  group_by(
    GFP
  ) %>% 
  mutate(
    Count = n(),
  ) %>% 
  ungroup() %>% 
  mutate(
    Frequency = Count/sum(Count),
    Cummulative = cumsum(Frequency),
    SK = round(Cummulative*(2^BITS-1))
  ) %>%
  as.data.table()

# Reference histogram
RefTable <-
  RefImage %>% 
  group_by(
    SK
  ) %>% 
  mutate(
    Transformed = min(GFP),
    Transformed = round(Transformed)
  ) %>% 
  select(
    Transformed,
    SK
  ) %>%
  arrange(
    SK
  ) %>% 
  distinct() %>% 
  as.data.table()

# Split frames for speed
SplitTable <- ALlImages %>% as_tibble() %>% group_by(FRAME) %>% group_split()
tictoc::toc()

tictoc::tic()
NormalizeFx <- function(QryImage){
  # Get table
  QryImage <-
    QryImage %>% 
    mutate(
      ID = 1:n()
    ) %>% 
    arrange(
      GFP
    ) %>% 
    group_by(
      GFP
    ) %>% 
    mutate(
      Count = n(),
    ) %>% 
    ungroup() %>% 
    mutate(
      Frequency = Count/sum(Count),
      Cummulative = cumsum(Frequency),
      SK = round(Cummulative*(2^BITS-1))
    ) %>%
    as.data.table()
  
  # Mapping
  # Match Each Query image sk value with closest SK value of referenc image and then
  # find the corresponding grey value of the reference image
  x = RefTable$SK
  dt = data.table(x, val = RefTable$Transformed)
  setattr(dt, "sorted", "x")
  setkey(dt, x)
  Index <- dt[J(QryImage$SK), roll = "nearest"]
  QryImage$Transformed <- Index$val
  QryImage <- QryImage %>% arrange(ID) %>% as.data.table()
  
  return(QryImage)
}
NormalizedImg <- mclapply(SplitTable, NormalizeFx, mc.cores = detectCores(logical = F))
NormalizedImg <- rbindlist(NormalizedImg)
tictoc::toc()

FrameOrder <- NULL
FrameOrder$ID <- 1:FRAME_COUNT
FrameOrder <- as_tibble(FrameOrder)

FrameOrder <-
  FrameOrder %>% 
  mutate(
    N_Row = round(sqrt(max(ID))),
    Row = ceiling(ID/N_Row),
    Column = ID-((Row-1)*(N_Row)),
    Column = ifelse(
      (Row %% 2) == 0,
      Column,
      N_Row - Column + 1
    ),
    FRAME = ID
  ) %>% 
  arrange(Row, Column) %>% 
  mutate(NEW_FRAME = 1:n()) %>% 
  select(FRAME, NEW_FRAME)

# Write multipage tiff
NormalizedImg <- merge(NormalizedImg, FrameOrder, by = "FRAME")
NormalizedImg <- NormalizedImg %>% arrange(NEW_FRAME) %>% as.data.table()
# Resort frames
ROWS = sqrt(MATRIX_SIZE)
COLUMNS = sqrt(MATRIX_SIZE)
write_tif(array(NormalizedImg$Transformed, dim = c(ROWS, COLUMNS,FRAME_COUNT)), "test.tif", bits_per_sample = BITS, overwrite = TRUE)


ijtiff::display(matrix(QryImage$GFP, nrow = 1200))
# 
# 
# 
# 
# ReadTiff <- function(FrameX){
#   # Start table for collating data
#   Image <- NULL
#   
#   # Read image
#   GFP_Image <- read_tif(IMAGE_PATH, FrameX, msg = F)
#   GFP_Image <- matrix(GFP_Image, nrow = NROW(GFP_Image))
#   
#   # Bin, if applicable
#   if(IMAGE_BIN_SIZE > 1){
#     x <- mat(NROW(GFP_Image), IMAGE_BIN_SIZE)
#     GFP_Image <- t(x) %*% GFP_Image %*% x
#     GFP_Image <- GFP_Image/IMAGE_BIN_SIZE^2
#   }
#   
#   # Make table from pixels
#   Image$GFP <- as.vector(GFP_Image)
#   Image$FRAME <- FrameX
#   return(Image)
# }
# Image <- mclapply(1:169, ReadTiff)
# Image <- rbindlist(Image)
# Image <- merge(Image, FrameOrder, by = "FRAME")
# Image <- Image %>% arrange(NEW_FRAME) %>% as.data.table()
# write_tif(array(Image$GFP, dim = c(ROWS, COLUMNS,FRAME_COUNT)), "quantile.tif", bits_per_sample = BITS, overwrite = TRUE)
# 
# 
# IMAGE_PATH="/Users/u_deliz/Desktop/WellB02_ChannelX RFP Cy5,X GFP,X DAPI_Seq0000_X DAPI-normalized.tif"
