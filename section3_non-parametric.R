# ============================================================
# MA416 Project 1 - Section 3: Non-Parametric Simulation of Pitching Volatility
# Subject: Tarik Skubal Contact Volatility under Pressure
# Custom Metric: Trajectory-Weighted Exit Velocity (TEV)
# ============================================================

ab$launch_speed <- as.numeric(ab$launch_speed)
ab$launch_angle <- as.numeric(ab$launch_angle)

# Filter to Tarik Skubal's contact at-bats
ts <- ab[!is.na(ab$launch_speed) & !is.na(ab$launch_angle) & 
           ab$player_name == "Skubal, Tarik", ]

cat("\nTotal contact at-bats found for Tarik Skubal:", nrow(ts), "\n")

# Apply the Custom Angle Weight (TEV)
# Centers the optimal angle at 20 degrees using cosine
ts$tev <- ts$launch_speed * cos((ts$launch_angle - 20) * pi / 180)

empty_1b <- is.na(ts$on_1b) | ts$on_1b == "" | ts$on_1b == "null"
empty_2b <- is.na(ts$on_2b) | ts$on_2b == "" | ts$on_2b == "null"
empty_3b <- is.na(ts$on_3b) | ts$on_3b == "" | ts$on_3b == "null"

ts$leverage <- ifelse(empty_1b & empty_2b & empty_3b, "Bases Empty", "Runners On")

cat("\n--- Tarik Skubal Contact At-Bats by Leverage ---\n")
print(table(ts$leverage))

# Calculate the Observed Custom Test Statistic
tev_on    <- ts$tev[ts$leverage == "Runners On"]
tev_empty <- ts$tev[ts$leverage == "Bases Empty"]

obs_tau <- var(tev_on, na.rm=TRUE) - var(tev_empty, na.rm=TRUE)

cat("\n--- Observed Volatility Statistic ---\n")
cat("Variance (Runners On): ", round(var(tev_on, na.rm=TRUE), 4), "\n")
cat("Variance (Bases Empty):", round(var(tev_empty, na.rm=TRUE), 4), "\n")
cat("Observed Tau:          ", round(obs_tau, 4), "\n")


# ============================================================
# PART A: MC Permutation Test (For the p-value)
# H0: Leverage does not affect contact volatility.
# ============================================================
set.seed(416)
nmc <- 10000

# Pre-allocate numeric vector for speed
mc_tau <- numeric(nmc)

for(k in 1:nmc) {
  shuffled_leverage <- sample(ts$leverage)
  sim_on    <- ts$tev[shuffled_leverage == "Runners On"]
  sim_empty <- ts$tev[shuffled_leverage == "Bases Empty"]
  
  mc_tau[k] <- var(sim_on, na.rm=TRUE) - var(sim_empty, na.rm=TRUE)
}

emp_pval <- mean(abs(mc_tau) >= abs(obs_tau), na.rm=TRUE)

cat("\n--- MC Permutation Test (H0) ---\n")
cat("Empirical p-value: ", round(emp_pval, 6), "\n")


# ============================================================
# PART B: Bootstrapping (For the 95% Confidence Interval)
# ============================================================
set.seed(416)
nboot <- 10000

# Pre-allocate numeric vector for speed
boot_tau <- numeric(nboot)

for(k in 1:nboot) {
  boot_on    <- sample(tev_on, length(tev_on), replace = TRUE)
  boot_empty <- sample(tev_empty, length(tev_empty), replace = TRUE)
  
  boot_tau[k] <- var(boot_on, na.rm=TRUE) - var(boot_empty, na.rm=TRUE)
}

ci_tau <- quantile(boot_tau, c(0.025, 0.975), na.rm=TRUE)

cat("\n--- Bootstrap 95% CI (L09) ---\n")
cat("Lower Bound:", round(ci_tau[1], 4), "\n")
cat("Upper Bound:", round(ci_tau[2], 4), "\n")


# ============================================================
# GRAPHICAL SUMMARY
# ============================================================

# Save the plot directly to your current working directory
png("volatility_null.png", width=800, height=600, res=100)
plot(density(mc_tau, na.rm=TRUE), col='black', lwd=2, 
     main='Permutation Null Distribution of Volatility Statistic (Skubal)',
     xlab='Difference in TEV Variance', ylab='Density')
abline(v=obs_tau, col='red', lwd=2, lty=1)
legend('topright', legend=c('Null Distribution', 'Observed Tau'), 
       col=c('black', 'red'), lwd=2, lty=c(1,1))
dev.off()

cat("\nSuccess! The graph has been saved as 'volatility_null.png' in:\n")
cat(getwd(), "\n")
