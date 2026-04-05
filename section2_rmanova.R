# MA416 Project 1 - Section 2 and Section 3 (RM-ANOVA)
# Subject: Pitcher, Within-Factor: TTO (1, 2, 3)
# Response: mean xwOBA per (pitcher, TTO)

dir.create('C:/Users/aasmi/ma 416/Project/Plots', showWarnings=FALSE)

ab <- read.csv('C:/Users/aasmi/ma 416/Project/atbat_2024.csv', stringsAsFactors=FALSE)

tto_per_pitcher <- tapply(ab$TTO, ab$pitcher, function(x) sort(unique(x)))
eligible        <- names(tto_per_pitcher)[sapply(tto_per_pitcher, function(x) all(c(1,2,3) %in% x))]
ab_elig         <- ab[ab$pitcher %in% eligible & ab$TTO %in% c(1,2,3) & !is.na(ab$estimated_woba_using_speedangle), ]

pitcher_tto_mean <- tapply(ab_elig$estimated_woba_using_speedangle,
                           list(ab_elig$pitcher, ab_elig$TTO), mean)
pitcher_tto_mean <- pitcher_tto_mean[complete.cases(pitcher_tto_mean), ]

T1 <- pitcher_tto_mean[, 1]
T2 <- pitcher_tto_mean[, 2]
T3 <- pitcher_tto_mean[, 3]

X      <- cbind(T1, T2, T3)
n      <- nrow(X)
g      <- ncol(X)
alpha  <- 0.05

cat("Eligible pitchers (all 3 TTO levels):", n, "\n")
cat("TTO column means:\n"); print(round(colMeans(X), 6))

# Assumption 1 - Normality

mysf = function(X){
  n      = length(X)
  X      = sort(X)
  mk     = qnorm((1:n - 3/8) / (n + 1/4), lower.tail=T)
  ak     = mk / sqrt(sum(mk^2))
  SST    = sum((X - mean(X))^2)
  tau    = (sum(ak * X))^2
  SFstat = tau / SST
  return(SFstat)
}

sf_obs <- c(mysf(T1), mysf(T2), mysf(T3))
names(sf_obs) <- c('TTO1', 'TTO2', 'TTO3')
cat("\nSF statistics per TTO level:\n"); print(round(sf_obs, 6))

set.seed(416)
nmc <- 10000

sf_pvals <- sapply(list(T1, T2, T3), function(X){
  obs   <- mysf(X)
  n     <- length(X)
  mc_sf <- c()
  for(k in 1:nmc){
    mcsam <- rnorm(n, 0, 1)
    mc_sf <- c(mc_sf, mysf(mcsam))
  }
  mean(mc_sf <= obs)
})
names(sf_pvals) <- c('TTO1', 'TTO2', 'TTO3')
cat("\nSF MC p-values per TTO level (nmc=10000, left-tail):\n")
print(round(sf_pvals, 6))

png('C:/Users/aasmi/ma 416/Project/Plots/qq_plots_rm_normality.png', width=1200, height=500)
par(mfrow=c(1, 3))
for(j in 1:g){
  Xj       <- sort(X[, j])
  nj       <- length(Xj)
  Fx       <- seq(1, nj, 1) / (nj + 1)
  S_ideal  <- qnorm(Fx, 0, 1, lower.tail=T)
  S_actual <- (Xj - mean(Xj)) / sd(Xj)
  plot(S_ideal, S_actual, pch=16, cex=0.5, col='black',
       xlab='Theoretical Quantiles', ylab='Sample Quantiles',
       main=paste('QQ Plot - TTO', j))
  abline(0, 1, lty=2, col='red')
}
par(mfrow=c(1, 1))
dev.off()

# Assumption 2 - Sphericity (Mauchly's Test)

sighatx  <- var(X)
ones     <- matrix(1, nrow=g, ncol=g)
C        <- diag(g) - (1/g)*ones
sighatc  <- C %*% sighatx %*% C
sighat2  <- (1/g) * sum(diag(sighatx))
W        <- det(sighatx) / (sighat2^g)

M_stat_o <- -(n-1) * log(W)
M_stat_b <- -((n-1) - (2*g+1)/6) * log(W)
mpval_o  <- pchisq(M_stat_o, g*(g-1)/2,     lower.tail=F)
mpval_b  <- pchisq(M_stat_b, g*(g-1)/2 - 1, lower.tail=F)

cat("\n--- Mauchly's Sphericity Test ---\n")
cat("W:         ", round(W, 6), "\n")
cat("M_stat_o:  ", round(M_stat_o, 4), "  p-val:", round(mpval_o, 6), "\n")
cat("M_stat_b:  ", round(M_stat_b, 4), "  p-val:", round(mpval_b, 6), "\n")

eps_gg <- (sum(diag(sighatc)))^2 / ((g-1) * sum(diag(sighatc %*% sighatc)))
eps_hf <- (n*(g-1)*eps_gg - 2) / ((g-1)*(n-1) - (g-1)*eps_gg)

cat("eps_GG:", round(eps_gg, 5), "\n")
cat("eps_HF:", round(eps_hf, 5), "\n")

# Normality: All three SF p-values equal zero at n = 311, mirroring the Section 1 result.
# The SF statistic for TTO3 (0.745) is substantially lower than for TTO1 (0.928) and
# TTO2 (0.942), indicating a more pronounced departure from normality at the third time
# through the order — likely because the distribution of mean xwOBA at TTO3 is noisier
# and more right-skewed than at earlier TTO levels. With n = 311, the SF test is
# hypersentive and will reject for even minor structural departures. The QQ plots provide
# the more honest qualitative picture. This normality violation motivates the non-parametric
# follow-up in Section 3.
#
# Sphericity: Both Mauchly statistics are extreme (M_o = 416.72, M_b = 415.16, both
# p-values = 0), providing overwhelming evidence against sphericity. The variances of the
# pairwise TTO differences are not equal. Since eps_GG = 0.687 < 0.75, the
# Greenhouse-Geisser correction is applied, adjusting the numerator df to v_msm = 1.37 and
# denominator df to v_mse = 426.00 for a more conservative F-test.
#
# Pairwise dependence: Each subject (pitcher) contributes one observation per TTO level —
# the within-subject mean xwOBA across all their contact at-bats at that TTO. The same
# pitcher is measured at TTO 1, 2, and 3, so the three groups are pairwisely dependent by
# design. RM-ANOVA accounts for this structure by partitioning out SSS.
#
# Independence within groups: Across pitchers, observations are independent. One pitcher's
# mean xwOBA at a given TTO level has no direct influence on another's. The 311 subjects
# were selected on the basis of having contact outcomes at all three TTO levels in 2024.
#
# Independence (qualitative): Each pitcher contributes exactly one mean per TTO level.
# The within-subject correlation is captured by SSS, not treated as noise. Across pitchers,
# the observations are independent — eligibility and performance of one pitcher do not
# influence another's measurements.
#
# Normality (qualitative + quantitative): SF rejects normality at all three TTO levels
# (p = 0, n = 311). The most severe departure is at TTO3 (SF = 0.745), where the
# distribution of pitcher means is right-skewed due to the smaller and less stable sample
# of at-bats. Per the professor's guidance, at n = 311 the SF test is hypersensitive and
# even minor deviations produce p = 0; the QQ plots confirm non-normality but not a
# catastrophic departure. We proceed with the corrected F-test while addressing normality
# formally through the non-parametric Section 3.
#
# Sphericity: Mauchly's test is decisively rejected (M_o = 416.72, p ≈ 0). The pairwise
# TTO differences have unequal variances. eps_GG = 0.687 < 0.75, so the Greenhouse-Geisser
# correction is applied with v_msm = 1.37 and v_mse = 426.00.

# Section 2 - RM-ANOVA

xbarj   <- colMeans(X)
xbari   <- rowMeans(X)
xbarbar <- mean(X)
s2all   <- var(c(T1, T2, T3))

SST <- s2all * (n*g - 1)
SSM <- n * sum((xbarj - xbarbar)^2)
SSS <- g * sum((xbari - xbarbar)^2)
E   <- X - xbari - matrix(xbarj, n, g, byrow=T) + xbarbar
SSE <- sum(E^2)

cat("\nSST - (SSM + SSS + SSE):", SST - (SSM + SSS + SSE), " (should be ~0)\n")

MSM   <- SSM / (g-1)
MSE   <- SSE / ((n-1)*(g-1))
Fstat <- MSM / MSE
Fcrit <- qf(alpha, g-1, (n-1)*(g-1), lower.tail=F)
pval  <- pf(Fstat, g-1, (n-1)*(g-1), lower.tail=F)

cat("\n--- RM-ANOVA (Uncorrected) ---\n")
cat("SSM:", round(SSM,6), "  SSE:", round(SSE,6), "  SST:", round(SST,6), "\n")
cat("MSM:", round(MSM,6), "  MSE:", round(MSE,6), "\n")
cat("Fstat:", round(Fstat,4), "  Fcrit:", round(Fcrit,4), "\n")
cat("p-val:", round(pval,6), "\n")

if(eps_gg >= 0.75){
  eps_use <- eps_hf
  cat("\neps_GG >= 0.75: applying HF correction\n")
} else {
  eps_use <- eps_gg
  cat("\neps_GG < 0.75: applying GG correction\n")
}

v_msm         <- (g-1) * eps_use
v_mse         <- (n-1) * (g-1) * eps_use
Fcrit_corr    <- qf(alpha, v_msm, v_mse, lower.tail=F)
pval_corrected <- pf(Fstat, v_msm, v_mse, lower.tail=F)

cat("\n--- RM-ANOVA (Corrected) ---\n")
cat("epsilon used:", round(eps_use, 5), "\n")
cat("v_msm:", round(v_msm, 4), "  v_mse:", round(v_mse, 4), "\n")
cat("Fstat:", round(Fstat, 4), "  Fcrit:", round(Fcrit_corr, 4), "\n")
cat("p-val (corrected):", round(pval_corrected, 6), "\n")

# Section 3 - Non-Parametric Follow-Up (RM)
#
# H0: mu_TTO1 = mu_TTO2 = mu_TTO3  (no time-through-order effect on xwOBA)
# Ha: At least one TTO level mean is significantly different
#
# Custom test statistic: T = max(xbarj) - min(xbarj), the range of TTO column means.
# Under H0 the TTO label carries no information, so shuffling within each pitcher row
# preserves the subject structure while destroying any TTO effect. A large observed T
# relative to the permutation null distribution gives evidence against H0.

myrange_rm = function(X){
  xbarj = colMeans(X)
  return(max(xbarj) - min(xbarj))
}

obs_stat <- myrange_rm(X)
cat("\nObserved test statistic (range of TTO means):", round(obs_stat, 6), "\n")

set.seed(416)
nmc      <- 10000
mc_stats <- c()
for(k in 1:nmc){
  Xperm    <- t(apply(X, 1, sample))
  mc_stats <- c(mc_stats, myrange_rm(Xperm))
}

mycrit   <- quantile(mc_stats, 1-alpha)
emp_pval <- mean(mc_stats >= obs_stat)

cat("\n--- MC Within-Subject Permutation Test ---\n")
cat("Critical value (alpha=0.05):", round(mycrit, 6), "\n")
cat("Observed statistic:         ", round(obs_stat, 6), "\n")
cat("Empirical p-value:          ", round(emp_pval, 6), "\n")
cat("Decision: Reject H0?", obs_stat > mycrit, "\n")

set.seed(416)
nboot   <- 10000
ci_tto  <- list()
tto_names <- c('TTO1', 'TTO2', 'TTO3')
for(j in 1:g){
  Xcol   <- X[, j]
  mymean <- c()
  for(k in 1:nboot){
    sboot  <- sample(Xcol, length(Xcol), replace=T)
    mymean <- c(mymean, mean(sboot))
  }
  ci_tto[[j]] <- quantile(mymean, c(alpha/2, 1-alpha/2))
}

cat("\n--- Bootstrap 95% CIs for TTO Group Means ---\n")
for(j in 1:g){
  cat(tto_names[j], ": (", round(ci_tto[[j]][1], 6),
      ",", round(ci_tto[[j]][2], 6), ")\n")
}

set.seed(416)
nmc_pow   <- 1000
mystat_ha <- c()

for(k in 1:nmc_pow){
  Xsim <- matrix(0, nrow=n, ncol=g)
  for(j in 1:g){
    Xsim[, j] <- rnorm(n, mean=xbarj[j], sd=sd(X[, j]))
  }
  mystat_ha <- c(mystat_ha, myrange_rm(Xsim))
}

type2_hat    <- mean(mystat_ha <= mycrit)
power_hat    <- mean(mystat_ha > mycrit)
se_powerhat  <- sqrt(power_hat * (1-power_hat) / nmc_pow)
zcrit_ci     <- qnorm(alpha/2, lower.tail=F)
eps_powerhat <- zcrit_ci * se_powerhat
mcci_power   <- c(power_hat - eps_powerhat, power_hat + eps_powerhat)

cat("\n--- Power Estimation (MC under Ha*) ---\n")
cat("power_hat:", round(power_hat, 4), "\n")
cat("type2_hat:", round(type2_hat, 4), "\n")
cat("Wald 95% CI for power: (",
    round(mcci_power[1], 4), ",", round(mcci_power[2], 4), ")\n")

png('C:/Users/aasmi/ma 416/Project/Plots/rm_permutation_null_dist.png', width=900, height=600)
plot(density(mc_stats), col='black', lwd=2,
     main='Within-Subject Permutation Null Distribution (TTO)',
     xlab='Range of TTO Means (T)', ylab='Density')
abline(v=obs_stat, col='red',  lwd=2, lty=1)
abline(v=mycrit,   col='blue', lwd=2, lty=2)
legend('topright',
       legend=c('Null Distribution', 'Observed T', 'Critical Value'),
       col=c('black','red','blue'), lwd=2, lty=c(1,1,2))
dev.off()

# Both the GG-corrected RM-ANOVA (Fstat = 7.618, Fcrit = 3.448, p = 0.0025) and the
# within-subject permutation test (T_obs = 0.0276, critical value = 0.0166, p = 0.0001)
# reject H0 at alpha = 0.05, providing consistent evidence that mean xwOBA allowed
# changes significantly as pitchers face the batting order additional times within a game.
#
# The TTO column means increase monotonically: TTO1 = 0.317, TTO2 = 0.330, TTO3 = 0.345.
# The bootstrap 95% CIs reveal the structure of this effect: TTO1 (0.312, 0.323) and
# TTO3 (0.329, 0.361) do not overlap, confirming a practically meaningful increase in
# xwOBA from the first to the third time through the order. TTO2 (0.324, 0.336) sits
# between them with partial overlap on both sides, suggesting a gradual rather than
# sudden degradation. This pattern is consistent with the well-documented time-through-order
# penalty in baseball: batters gain a measurable advantage on repeated looks at a starting
# pitcher as they accumulate information on pitch selection and velocity.
#
# The power of the permutation test is estimated at 89.2% (Wald 95% CI: 0.873–0.911),
# indicating that if the true TTO means are as observed in our 2024 sample, the test would
# correctly reject H0 in approximately 89 out of 100 repeated experiments. This high power
# is consistent with n = 311 and the clear monotonic trend in TTO means. The agreement
# between the parametric and non-parametric results strengthens the conclusion: pitchers
# allow meaningfully more hard contact the further they advance through the batting order.
