## code to prepare `DATASET` dataset goes here

r <- 0.2
K <- 1000
q <- 0.00025
Bo <- K
nYears <- 20
effort <- c(seq(1,500, length.out = nYears/2), rev(seq(1,500, length.out = nYears / 2)))

B <- CPUE <- C <- rep(NA, nYears)
procError <- 0.05
catchError <- 0.05
for (i in 1:nYears) {
  if (i == 1) B[i] <- Bo
  if (i>1) B[i] <- B[i-1] + B[i-1] * r * (1 - B[i-1]/K) - C[i-1]
  C[i] <- q * effort[i] * B[i]
  CPUE[i] <- C[i] / effort[i]
}

df.goodcontrast0 <- data.frame(year=1981:2000,
                            catch = C,
                            cpue = CPUE,
                            biomass = B)

usethis::use_data(df.goodcontrast0, overwrite = TRUE)

