#Packages

library(mipfp)
library(ggplot2)
library(reshape2)

#Parameters

n_sim      <- 200        # Monte Carlo repetitions
drift_low  <- 0.05       # 5% structural drift
drift_med  <- 0.15       # 15% drift
drift_high <- 0.30       # 30% drift

scenarios <- c("Low","Medium","High")


# Base Census Structure

generate_base <- function() {
  base <- array(c(
    300, 200,   # Young LowEdu
    250, 150,   # Young HighEdu
    400, 100,   # Adult LowEdu
    350, 150    # Adult HighEdu
  ),
  dim = c(2,2,2),
  dimnames = list(
    Age = c("Young","Adult"),
    Education = c("Low","High"),
    Employment = c("Employed","NotEmployed")
  ))
  return(base)
}

# Structural Drift

apply_drift <- function(base, drift_intensity) {
  
  true_pop <- base
  
  # Modify employment probabilities conditional on education
  for(a in 1:2){
    for(e in 1:2){
      
      employed   <- true_pop[a,e,1]
      unemployed <- true_pop[a,e,2]
      total      <- employed + unemployed
      
      # Shift employment probability
      shift <- drift_intensity * runif(1,-1,1)
      
      new_prob <- (employed/total) + shift
      new_prob <- max(min(new_prob,0.95),0.05)
      
      true_pop[a,e,1] <- total * new_prob
      true_pop[a,e,2] <- total * (1-new_prob)
    }
  }
  
  return(true_pop)
}

# Margins (observed surveys, constraints)

get_margins <- function(true_pop) {
  
  margin_age_edu <- apply(true_pop, c(1,2), sum)
  margin_emp     <- apply(true_pop, 3, sum)
  
  return(list(margin_age_edu, margin_emp))
}

# IPF

run_ipf <- function(seed, margins){
  
  target.list <- list(c(1,2), c(3))
  target.data <- margins
  
  result <- Ipfp(
    seed = seed,
    target.list = target.list,
    target.data = target.data,
    iter = 1000,
    tol = 1e-8
  )
  
  return(result$x.hat)
}

# Error metrics

compute_metrics <- function(true_pop, ipf_pop){
  
  mae <- mean(abs(true_pop - ipf_pop))
  
  p_true <- as.vector(true_pop/sum(true_pop))
  p_ipf  <- as.vector(ipf_pop/sum(ipf_pop))
  
  kl <- sum(p_true * log(p_true/p_ipf))
  
  relative_mae <- mean(abs(true_pop - ipf_pop)) / mean(true_pop)
  
  return(c(MAE = mae, KL = kl, RelativeMAE=relative_mae))
}


# Monte carlo loop

run_simulation <- function(drift_intensity, label){
  
  results <- matrix(NA, n_sim, 3)
  colnames(results) <- c("MAE","KL", "RelativeMAE")   # <- important fix
  
  for(i in 1:n_sim){
    
    base      <- generate_base()
    true_pop  <- apply_drift(base, drift_intensity)
    margins   <- get_margins(true_pop)
    ipf_pop   <- run_ipf(base, margins)
    metrics   <- compute_metrics(true_pop, ipf_pop)
    
    results[i,] <- metrics
  }
  
  df <- as.data.frame(results)
  df$Scenario <- label
  
  return(df)
}
  

# Run scenarios

set.seed(12231233)

res_low  <- run_simulation(drift_low,  "Low")
res_med  <- run_simulation(drift_med,  "Medium")
res_high <- run_simulation(drift_high, "High")

results_all <- rbind(res_low, res_med, res_high)


# Plots

ggplot(results_all, aes(x=Scenario, y=MAE)) +
  geom_boxplot() +
  ggtitle("IPF Reliability Across Drift Scenarios") +
  theme_minimal()

ggplot(results_all, aes(x=Scenario, y=KL)) +
  geom_boxplot() +
  ggtitle("KL Divergence Across Drift Scenarios") +
  theme_minimal()

  # | KL Value  | Interpretation             |
  # | --------- | -------------------------- |
  # | < 0.01    | Almost identical           |
  # | 0.01–0.05 | Small distortion           |
  # | 0.05–0.15 | Moderate divergence        |
  # | 0.15–0.30 | Significant distortion     |
  # | > 0.30    | Strong structural mismatch |



ggplot(results_all, aes(x=Scenario, y=RelativeMAE)) +
  geom_boxplot() +
  ggtitle("KL Divergence Across Drift Scenarios") +
  theme_minimal()

