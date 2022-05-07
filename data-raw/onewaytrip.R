## code to prepare `DATASET` dataset goes here

r <- 0.2
K <- 1000
q <- 0.00025
Bo <- K
nYears <- 20
effort <- seq(1,500, length.out = nYears)

B <- CPUE <- C <- rep(NA, nYears)
procError <- 0.05
catchError <- 0.05
for (i in 1:nYears) {
  if (i == 1) B[i] <- Bo
  if (i>1) B[i] <- B[i-1] + rlnorm(1, meanlog = -procError^2/2, sdlog = procError) * B[i-1] * r * (1 - B[i-1]/K) - C[i-1]
  C[i] <- q * effort[i] * B[i] * rlnorm(1, meanlog = -catchError^2/2, sdlog = catchError)
  CPUE[i] <- C[i] / effort[i]
}

df.onewaytrip <- data.frame(year=1981:2000,
                         catch = C,
                         effort = effort)

usethis::use_data(df.onewaytrip, overwrite = TRUE)
