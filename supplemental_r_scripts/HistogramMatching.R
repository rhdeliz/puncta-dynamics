library(ijtiff)
library(parallel)
library(data.table)
library(dplyr)

# IMAGE_PATH = "/Users/u_deliz/Desktop/StainTest/EL4-WT/p-p38 Stimulated 20201011 WellC02_EL4-WT_p-p38_X DAPI.tif"
IMAGE_PATH = "/Users/u_deliz/Desktop/WellB02_ChannelX RFP Cy5,X GFP,X DAPI_Seq0000_X DAPI.tif"
IMAGE_BIN_SIZE = 2
BITS = 8
# INTENSITY_BIN_SIZE = 25


tictoc::tic()
# For binning matrix
mat <- function(n, r) {
  suppressWarnings(matrix(c(rep(1, r), rep(0, n)), n, n/r))
}

# Import image
ReadImageFx <- function(ImageFrameX){
  # Pull image
  ImageX <- ceiling(ImageFrameX/FRAME_COUNT)
  # Pull frame
  FrameX <- ImageFrameX - (ImageX*FRAME_COUNT) + FRAME_COUNT
  
  # Start table for collating data
  Image <- NULL
  
  # Read image
  GFP_Image <- read_tif(IMAGE_PATH, FrameX, msg = F)
  GFP_Image <- matrix(GFP_Image, nrow = NROW(GFP_Image))
  
  # Bin, if applicable
  if(IMAGE_BIN_SIZE > 1){
    x <- mat(NROW(GFP_Image), IMAGE_BIN_SIZE)
    GFP_Image <- t(x) %*% GFP_Image %*% x
  }
  
  # Make table from pixels
  Image$GFP <- as.vector(GFP_Image)
  Image$FRAME <- ImageFrameX
  
  # Scale images
  # Image$GFP <- Image$GFP - min(Image$GFP)
  # Image$GFP <- Image$GFP/max(Image$GFP)*(2^BITS-1)
  
  # Get ratio of bright to dim
  # To be used later for filtering dark and dusty images
  IMG_HIGH = quantile(Image$GFP, .9)
  IMG_LOW = quantile(Image$GFP, .1)
  Image$RATIO = IMG_HIGH/IMG_LOW
  
  # Export
  Image <- as.data.table(Image)
  return(Image)
}
# Parallelize
FRAME_COUNT = count_frames(IMAGE_PATH)[1]
ALlImages <- mclapply(1:FRAME_COUNT, ReadImageFx)
ALlImages <- rbindlist(ALlImages)

# Eliminate dark and dusty images
# ALlImages <-
#   ALlImages %>% 
#   filter(
#     RATIO >= 20,
#     RATIO <= 100
#   )

# Get reference frame
# A bright image

RefFrame <- quantile(ALlImages$RATIO, .9)
RefFrame <- abs(ALlImages$RATIO - RefFrame)
RefFrame <- which(RefFrame == min(RefFrame))
RefFrame <- ALlImages$FRAME[RefFrame]
RefFrame <- min(RefFrame)
# Extract reference frame
RefImage <- 
  ALlImages %>% 
  filter(
    FRAME == RefFrame
  )

# Histogram match images
NormalizeFx <- function(FrameX){
  
  # Pull table
  QueryImage <- 
    ALlImages %>% 
    filter(
      FRAME == FrameX
    )
  # Normalize
  NormImg <- landsat::histmatch(RefImage$GFP, QueryImage$GFP, minval = 0, maxval = 2^BITS-1, by = 1)
  NormImg <- NormImg$newimage
  
  return(NormImg)
}
# Parallelize
NormImages <- mclapply(unique(ALlImages$FRAME), NormalizeFx, mc.cores = detectCores(logical = F))
NormImages <- unlist(NormImages)

# Write multipage tiff
NormImages[is.na(NormImages)] = 0
FRAME_COUNT = NROW(unique(ALlImages$FRAME))
# Resort frames
write_tif(array(NormImages, dim = c(600,600,FRAME_COUNT)), "test.tif", bits_per_sample = 16, overwrite = TRUE)
tictoc::toc()
