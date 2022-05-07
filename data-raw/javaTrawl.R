## code to prepare `DATASET` dataset goes here

df.javaTrawl <- data.frame(year=c(1969:1977),
                        catch=c(50,49,47.5,45,51,56,66,58,52),
                        effort=c(623,628,520,513,661,919,1158,1970,1317))

usethis::use_data(df.javaTrawl, overwrite = TRUE)
