library(mipfp)

# 2020 True population

true_2020 <- array(c(
  400, 150,   # Young Male
  380, 170,   # Young Female
  700, 200,   # Adult Male
  720, 180    # Adult Female
),
dim = c(2,2,2),
dimnames = list(
  Age = c("Young","Adult"),
  Sex = c("Male","Female"),
  Employment = c("Employed","NotEmployed")
))

true_2020

# Outdated seed 2010

seed_2010 <- array(c(
  500, 100,
  500, 100,
  600, 100,
  600, 100
),
dim = c(2,2,2),
dimnames = list(
  Age = c("Young","Adult"),
  Sex = c("Male","Female"),
  Employment = c("Employed","NotEmployed")
))

seed_2010

#Marginals (updated data - constraints)

margin_age_sex <- apply(true_2020, c(1,2), sum)
margin_age_sex

margin_employment <- apply(true_2020, 3, sum)
margin_employment

# IPF

target.list <- list(
  c(1,2),
  c(3)
)

target.data <- list(
  margin_age_sex,
  margin_employment
)

result <- Ipfp(
  seed = seed_2010,
  target.list = target.list,
  target.data = target.data,
  iter = 1000,
  tol = 1e-8
)

ipf_result <- result$x.hat
ipf_result

# Comparison

comparison <- data.frame(
  True = as.vector(true_2020),
  IPF  = as.vector(ipf_result)
)

comparison

error <- abs(ipf_result - true_2020)
error

mean(abs(error))
