NCATS_screening_data<- read.csv("data-raw/NCATS_screening_data.csv", stringsAsFactors = FALSE)
ONEIL_screening_data <- read.csv("data-raw/ONEIL_screening_data.csv")

usethis::use_data(NCATS_screening_data, overwrite = TRUE)
usethis::use_data(ONEIL_screening_data, overwrite = TRUE)