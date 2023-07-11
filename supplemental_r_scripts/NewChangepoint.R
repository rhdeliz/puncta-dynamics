library(changepoint)

ChangepointFx <- function(x){
  tryCatch({
    
    df <- NULL
    df$x <- 1:NROW(x)
    # df$y <- log(x+1)
    df$y <- x
    
    # Run changepoint
    m.pelt <- cpt.mean(df$y, method = "PELT", penalty = "AIC", minseglen = 3)
    # m.pelt <- cpt.mean(df$y, method = "PELT", penalty = "Manual", pen.value = .7, minseglen = 3)
    # Get y
    y <- param.est(m.pelt)$mean
    # Transform y back to linear
    # y <- exp(y)-1
    # Expand y to match time
    len <- seg.len(m.pelt)
    y <- rep(y, times = len)
    
    return(y)
  }, error = function(e) {print("Error with ChangepointFx")})}

Test <- rpois(20, lambda = 5)+1
Test <- zoo::rollmean(Test, 4)

ExpandFx <- function(x){
  Length <- 0.5/rexp(1, x)+5
  Length <- round(Length)
  Result <- rnorm(Length, x, .1*x)
  Ground <- rep(x, Length)
  df <- NULL
  df$Result = Result
  df$Ground = Ground
  df$CP = ChangepointFx(Result)
  return(df)
}
NoiseTable <- lapply(Test, ExpandFx)
NoiseTable <- rbindlist(NoiseTable)
NoiseTable$t <- 1:NROW(NoiseTable)


ggplot(
  NoiseTable
) +
  geom_path(
    aes(
      x = t,
      y = Result,
      color = "Noise"
    )
  ) +
  geom_path(
    aes(
      x = t,
      y = Ground,
      color = "Ground Truth"
    )
  ) +
  geom_path(
    aes(
      x = t,
      y = CP,
      color = "Change Point"
    )
  ) +
  scale_color_manual(
    values = c("blue", "red", "#404040")
  ) +
  dark_theme_classic()

NoiseTable <-
  NoiseTable %>% 
  mutate(
    Difference = Ground - CP,
    Difference_pct = Difference/Ground*100
  ) %>% 
  as.data.table()

plot(NoiseTable$Difference_pct)
