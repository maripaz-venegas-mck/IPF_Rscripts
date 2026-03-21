# Costa Rica - INEC

library(tidyr)
library(dplyr)
library(mipfp)

# Variable table
# This is a table in which the reference values are stored. 
# Each row represents a grup of people or households
seed <- expand.grid(
  size = c("1-2","3-4","5+"),
  income = c("low","mid","high")
)

#Frequency table

# Initial uniform frequency and population
seed$freq <- 1
pop <- seed[rep(1:nrow(seed), round(seed$freq*10000)), ]


# Constraints
# This is the information as aggregates of the variable that will help simulate a real population
marg_size <- c("1-2" = 0.35,
               "3-4" = 0.40,
               "5+"  = 0.25)
marg_income <- c("low"  = 0.30,
                 "mid"  = 0.45,
                 "high" = 0.25)

# N - total population to be built
N <- 10000
# Constraints are applied
marg_size_n <- c("1-2"=3500, "3-4"=4000, "5+"=2500)
marg_income_n <- c("low"=3000, "mid"=4500, "high"=2500)

#Data as matrix
seed_matrix <- xtabs(freq ~ size + income, data = seed)

#Apply IPF
result_matrix <- Estimate(
  seed = seed_matrix,
  target.data = c(marg_size_n, marg_income_n),
  target.list = list(1, 2)
)


#Microdata
result_df <- as.data.frame(as.table(result_matrix$x.hat))
colnames(result_df) <- c("size","income","freq")

pop <- result_df |>
  slice(rep(1:n(), round(freq))) |>
  select(-freq)

sum(result_df$freq)
rowSums(result_matrix$x.hat)
colSums(result_matrix$x.hat)

#Demand example
pop$can_buy <- pop$income == "high" |
  (pop$income == "mid" & pop$size %in% c("1-2","3-4"))

mean(pop$can_buy)


# Validation
prop.table(table(pop$size))
prop.table(table(pop$income))


