library(mipfp)

# Dimensions:
# Age: Young, Adult
# Sex: Male, Female
# Employment: Employed, NotEmployed

#2020

seed <- array(c(
  500, 200,   # Young Male
  480, 220,   # Young Female
  600, 150,   # Adult Male
  620, 180    # Adult Female
),
dim = c(2, 2, 2),
dimnames = list(
  Age = c("Young", "Adult"),
  Sex = c("Male", "Female"),
  Employment = c("Employed", "NotEmployed")
))

View(as.data.frame(seed))

# 2025

margin_age_sex <- matrix(c(
  900, 950,   # Young: Male, Female
  1100, 1150  # Adult: Male, Female
),
nrow = 2,
byrow = TRUE)

dimnames(margin_age_sex) <- list(
  Age = c("Young", "Adult"),
  Sex = c("Male", "Female")
)

View(as.data.frame(margin_age_sex))

margin_employment <- c(
  Employed = 3200,
  NotEmployed = 900
)

View(as.data.frame(margin_employment))


# IPF targets

target.list <- list(
  c(1,2),  # Age × Sex
  c(3)     # Employment
)

target.data <- list(
  margin_age_sex,
  margin_employment
)

result <- Ipfp(
  seed = seed,
  target.list = target.list,
  target.data = target.data,
  iter = 1000,
  tol = 1e-8
)

updated_array <- result$x.hat
View(as.data.frame(updated_array))




apply(updated_array, c(1,2), sum)

apply(updated_array, 3, sum)

