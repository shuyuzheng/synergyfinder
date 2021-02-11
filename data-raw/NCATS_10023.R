meta <- read.csv("data-raw/NCATS_10023/metadata.csv", stringsAsFactors = FALSE)
response <- read.csv("data-raw/NCATS_10023/responses.csv", stringsAsFactors = FALSE)
meta$Conc3 <- c(
  0.75, 0.375, 0.1875, 0.938, 0.469, 0.0234, 0.0117, 0.0059,
  0.0029, 0.0015, 0.0007, 0
)
meta$Drug3 <- rep("Piperaquine", 12)

NCAST_10023 <- NULL

for (i in 1:12) {
  Conc1 <- as.numeric(unlist(strsplit(meta$RowConcs[i], ",")))
  Conc2 <- as.numeric(unlist(strsplit(meta$ColConcs[i], ",")))
  tmp <- response[which(response$BlockId == i), ]
  df <- data.frame(BlockId = rep(1, nrow(df)), stringsAsFactors = FALSE)
  df$Drug1 <- rep(meta$RowName[i], nrow(df))
  df$Drug2 <- rep(meta$ColName[i], nrow(df))
  df$Drug3 <- rep(meta$Drug3[i], nrow(df))
  df$Conc1 <- rep(Conc1, each = 10)
  df$Conc2 <- rep(Conc2, times = 10)
  df$Conc3 <- rep(meta$Conc3[i], nrow(df))
  df$Response <- tmp$Value
  df$ConcUnit1 <- rep(meta$RowConcUnit[i], nrow(df))
  df$ConcUnit2 <- rep(meta$ColConcUnit[i], nrow(df))
  df$ConcUnit3 <- rep("uM", nrow(df))
  NCAST_10023 <- rbind.data.frame(NCAST_10023, df)
}

usethis::use_data(NCAST_10023, overwrite = TRUE, internal = TRUE)
