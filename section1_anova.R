# MA416 Project 1 - Section 1 and Section 3
# Groups: AL Central Teams (KC, DET, MIN, CWS, CLE)
# Response: estimated_woba_using_speedangle (xwOBA on contact)

dir.create('C:/Users/aasmi/ma 416/Project/Plots', showWarnings=FALSE)

ab <- read.csv('C:/Users/aasmi/ma 416/Project/atbat_2024.csv', stringsAsFactors=FALSE)

al_central <- c("KC", "DET", "MIN", "CWS", "CLE")

ab$pitcher_team <- ifelse(ab$inning_topbot == "Top", ab$home_team, ab$away_team)

ab_alc <- ab[(ab$home_team %in% al_central | ab$away_team %in% al_central) &
               ab$pitcher_team %in% al_central, ]

pitcher_teams <- tapply(ab_alc$pitcher_team, ab_alc$pitcher, function(x) unique(x))
multi_team    <- pitcher_teams[sapply(pitcher_teams, length) > 1]
ab_alc        <- ab_alc[!ab_alc$pitcher %in% names(multi_team), ]
ab_alc        <- ab_alc[!is.na(ab_alc$estimated_woba_using_speedangle), ]

Y     <- ab_alc$estimated_woba_using_speedangle
group <- ab_alc$pitcher_team

cat("At-bats per team (xwOBA non-NA):\n")
print(table(group))

# ============================================================
# ASSUMPTION 1: NORMALITY
# Qualitative: QQ plots per group
# Quantitative: Shapiro-Francia test + MC p-value per group
# ============================================================

mysf = function(X){
  n  = length(X)
  X  = sort(X)
  mk = qnorm((1:n - 3/8) / (n + 1/4), lower.tail=T)
  ak = mk / sqrt(sum(mk^2))
  SST  = sum((X - mean(X))^2)
  tau  = (sum(ak * X))^2
  SFstat = tau / SST
  return(SFstat)
}

dsplit <- split(Y, group)
nj     <- sapply(dsplit, length)

sf_obs <- sapply(dsplit, mysf)
cat("\nSF statistics per team:\n")
print(round(sf_obs, 6))

set.seed(416)
nmc <- 1000

sf_pvals <- sapply(names(dsplit), function(team){
  X   <- dsplit[[team]]
  n   <- length(X)
  obs <- mysf(X)
  mc_sf <- c()
  for(k in 1:nmc){
    mcsam <- rnorm(n, 0, 1)
    mc_sf <- c(mc_sf, mysf(mcsam))
  }
  mean(mc_sf <= obs)
})

cat("\nSF MC p-values per team (nmc=1000, left-tail):\n")
print(round(sf_pvals, 6))

png('C:/Users/aasmi/ma 416/Project/Plots/qq_plots_normality.png', width=1200, height=800)
par(mfrow=c(2, 3))
for(team in names(dsplit)){
  X <- sort(dsplit[[team]])
  n <- length(X)
  Fx       <- seq(1, n, 1) / (n + 1)
  S_ideal  <- qnorm(Fx, 0, 1, lower.tail=T)
  S_actual <- (X - mean(X)) / sd(X)
  plot(S_ideal, S_actual, pch=16, cex=0.3, col='black',
       xlab='Theoretical Quantiles', ylab='Sample Quantiles',
       main=paste('QQ Plot -', team))
  abline(0, 1, lty=2, col='red')
}
par(mfrow=c(1, 1))
dev.off()

# ============================================================
# ASSUMPTION 2: HOMOSCEDASTICITY (EQUAL VARIANCES)
# Quantitative: Bartlett Test (L06)
# ============================================================

s2j  <- sapply(dsplit, var)
vj   <- nj - 1
n    <- sum(nj)
g    <- length(dsplit)

sp2    <- sum(vj * s2j) / (n - g)
Batop  <- (n - g) * log(sp2) - sum(vj * log(s2j))
Babot  <- 1 + 1/(3*(g-1)) * (sum(1/vj) - 1/(n-g))
Bastat <- Batop / Babot

alpha  <- 0.05
Bacrit <- qchisq(alpha, g-1, lower.tail=F)
Bapval <- pchisq(Bastat, g-1, lower.tail=F)

cat("\n--- Bartlett Test (H0: equal variances) ---\n")
cat("Bastat:", round(Bastat, 4), "\n")
cat("Bacrit:", round(Bacrit, 4), "\n")
cat("p-val: ", round(Bapval, 6), "\n")

cat("\nGroup variances:\n");  print(round(s2j, 6))
cat("Group means:\n");        print(round(sapply(dsplit, mean), 6))

# ============================================================
# INTERPRETATION AND PATH FORWARD
#
# Normality: All five SF p-values = 0, meaning normality is
# heavily violated across every AL Central group. However, with
# n ~ 6,200 per group, the SF test (and all normality tests) are
# extremely sensitive — even minor departures from perfect
# normality will produce p = 0 at this sample size. The QQ plots
# provide the more honest qualitative picture: systematic
# deviation from the diagonal indicates the xwOBA distribution
# is non-normal (likely right-skewed with a mass at 0 from weak
# contact). Per the professor's guidance, heavy normality
# violation motivates a non-parametric follow-up in Section 3.
#
# Homoscedasticity: Bartlett rejects equal variances
# (p = 7e-06). Again, large n inflates Bartlett's sensitivity.
# Group variances are practically similar (0.127 to 0.144),
# suggesting the violation is more statistical than practically
# meaningful. Regardless, we present Welch's heteroscedastic
# ANOVA as the appropriate alternative, which relaxes the
# equal-variance assumption and uses a Satterthwaite-corrected
# denominator df.
#
# Path forward:
# - Section 1 ANOVA: homoscedasticity is violated, so we
#   proceed directly with Welch's heteroscedastic ANOVA as
#   the appropriate test.
# - Section 3 Non-Parametric: motivated by both the SF test
#   (all p-values = 0) and the QQ plots (systematic deviation
#   from the diagonal confirms non-normality qualitatively);
#   will apply an MC-based permutation/rank-style test on the
#   same groups as a distribution-free alternative.
# ============================================================

# ============================================================
# SECTION 1: HETEROSCEDASTIC (WELCH'S) ONE-FACTOR ANOVA (L06)
#
# Independence (qualitative): Each at-bat is assigned to
# exactly one pitcher's team. No pitcher appeared for more
# than one AL Central team in 2024 (verified: 0 multi-team
# pitchers). Observations across teams come from different
# pitchers in separate games. AL Central teams play a roughly
# balanced schedule against one another, so no team faces a
# systematically different batter pool. At-bats within a team
# share pitchers (some clustering), but with 6,200+ at-bats
# spread across a full pitching roster, this effect is minor.
# Independence is a reasonable qualitative assumption.
#
# Normality (qualitative + quantitative): SF test rejects for
# all five teams (p = 0). However, with n ~ 6,200 per group,
# normality tests are hyper-sensitive — the Glivenko-Cantelli
# theorem guarantees power → 1 as n → ∞, so even trivial
# departures produce p = 0 at this scale. Per the professor's
# guidance, in the super-large-n regime we slightly inflate
# our effective alpha threshold and proceed. The QQ plots
# confirm some deviation from normality but no catastrophic
# departure.
#
# Homoscedasticity: Rejected by Bartlett (p = 7e-06). We
# therefore use Welch's heteroscedastic ANOVA, which relaxes
# the equal-variance assumption and uses a corrected
# denominator df (v*) per the lecture formula.
# ============================================================

wj      <- nj / s2j
xbarH   <- sum(wj * xbarj) / sum(wj)
SSMstar <- sum(wj * (xbarj - xbarH)^2)
MSMstar <- SSMstar / (g - 1)

vstar   <- sum(wj)^2 / sum(wj^2 / vj)
Fstar   <- MSMstar / (1 + 2*(g-2) / vstar)
Fcrit_h <- qf(alpha, g-1, vstar, lower.tail=F)
pval_h  <- pf(Fstar, g-1, vstar, lower.tail=F)

cat("\n--- Heteroscedastic ANOVA (Welch's, L06) ---\n")
cat("MSMstar:", round(MSMstar, 6), "\n")
cat("vstar:  ", round(vstar, 4), "\n")
cat("Fstar:  ", round(Fstar, 4), "\n")
cat("Fcrit:  ", round(Fcrit_h, 4), "\n")
cat("p-val:  ", round(pval_h, 6), "\n")

# ============================================================
# SECTION 3: NON-PARAMETRIC FOLLOW-UP
# Motivation: Normality heavily violated in assumption testing
# (all SF p-values = 0, QQ plots confirm non-normality).
# A distribution-free permutation test is used as follow-up.
#
# H0: mu_CLE = mu_CWS = mu_DET = mu_KC = mu_MIN
# Ha: At least one pair of group means is significantly different
#
# Custom Test Statistic (designed by group):
#   T = max(xbarj) - min(xbarj)   [range of group means]
# Under H0, all group means are equal so T ≈ 0.
# Large T gives evidence against H0 (right-tail test).
#
# Method: MC permutation test — permute group labels across all
# observations (keeping group sizes fixed), recompute T each
# iteration to build the null distribution. No distributional
# assumption required. Bootstrap sampling test per L09/L11.
# ============================================================

myrange = function(Y, group){
  ds    = split(Y, group)
  xbarj = sapply(ds, mean)
  return(max(xbarj) - min(xbarj))
}

obs_stat <- myrange(Y, group)
cat("\nObserved test statistic (range of group means):", round(obs_stat, 6), "\n")

# ============================================================
# STEP 1: MC UNDER H0 — Permutation Null Distribution (L11)
# ============================================================

set.seed(416)
nmc      <- 10000
mc_stats <- c()
for(k in 1:nmc){
  group_perm <- sample(group, length(group), replace=F)
  mc_stats   <- c(mc_stats, myrange(Y, group_perm))
}

mycrit   <- quantile(mc_stats, 1-alpha)
emp_pval <- mean(mc_stats >= obs_stat)

cat("\n--- MC Permutation Test (H0) ---\n")
cat("Critical value (alpha=0.05):", round(mycrit, 6), "\n")
cat("Observed statistic:         ", round(obs_stat, 6), "\n")
cat("Empirical p-value:          ", round(emp_pval, 6), "\n")
cat("Decision: Reject H0?", obs_stat > mycrit, "\n")

# ============================================================
# STEP 2: BOOTSTRAP CIs FOR EACH GROUP MEAN (L09)
# ============================================================

set.seed(416)
nboot    <- 10000
ci_means <- list()

for(team in names(dsplit)){
  X      <- dsplit[[team]]
  mymean <- c()
  for(k in 1:nboot){
    sboot  <- sample(X, length(X), replace=T)
    mymean <- c(mymean, mean(sboot))
  }
  ci_means[[team]] <- quantile(mymean, c(alpha/2, 1-alpha/2))
}

cat("\n--- Bootstrap 95% CIs for Group Means (L09) ---\n")
for(team in names(ci_means)){
  cat(team, ": (", round(ci_means[[team]][1], 6),
      ",", round(ci_means[[team]][2], 6), ")\n")
}

# ============================================================
# STEP 3: MC UNDER Ha* — POWER ESTIMATION (L11)
# Ha*: true group means and sds = observed sample values.
# Generate data from each group under Ha*, apply permutation
# test, estimate power = P(reject H0 | Ha* true).
# Note: using nmc_pow=1000 due to large group sizes.
# ============================================================

set.seed(416)
nmc_pow   <- 1000
mystat_ha <- c()

for(k in 1:nmc_pow){
  Ysim <- c()
  gsim <- c()
  for(team in names(dsplit)){
    ni   <- nj[team]
    mui  <- xbarj[team]
    si   <- sqrt(s2j[team])
    Ysim <- c(Ysim, rnorm(ni, mui, si))
    gsim <- c(gsim, rep(team, ni))
  }
  mystat_ha <- c(mystat_ha, myrange(Ysim, gsim))
}

type2_hat <- mean(mystat_ha <= mycrit)
power_hat <- mean(mystat_ha > mycrit)

cat("\n--- Power Estimation (MC under Ha*, L11) ---\n")
cat("power_hat:", round(power_hat, 4), "\n")
cat("type2_hat:", round(type2_hat, 4), "\n")

se_powerhat  <- sqrt(power_hat * (1-power_hat) / nmc_pow)
zcrit_ci     <- qnorm(alpha/2, lower.tail=F)
eps_powerhat <- zcrit_ci * se_powerhat
mcci_power   <- c(power_hat - eps_powerhat, power_hat + eps_powerhat)

cat("Wald 95% CI for power: (",
    round(mcci_power[1], 4), ",", round(mcci_power[2], 4), ")\n")

# ============================================================
# STEP 4: GRAPHICAL SUMMARY (L09/L11 base R)
# ============================================================

png('C:/Users/aasmi/ma 416/Project/Plots/permutation_null_dist.png', width=900, height=600)
plot(density(mc_stats), col='black', lwd=2,
     main='Permutation Null Distribution of Range Statistic',
     xlab='Range of Group Means (T)', ylab='Density')
abline(v=obs_stat, col='red',  lwd=2, lty=1)
abline(v=mycrit,   col='blue', lwd=2, lty=2)
legend('topright',
       legend=c('Null Distribution', 'Observed T', 'Critical Value'),
       col=c('black','red','blue'), lwd=2, lty=c(1,1,2))
dev.off()

# ============================================================
# INTERPRETATION
#
# Both the Welch's heteroscedastic ANOVA (Fstar = 4.49,
# p = 0.0013) and the non-parametric permutation test
# (T_obs = 0.0251, p = 0.0021) reject H0 at alpha = 0.05,
# providing consistent evidence that mean xwOBA allowed is
# not equal across all five AL Central teams in 2024.
#
# The bootstrap 95% CIs reveal where the differences lie:
# DET pitchers allowed the lowest mean xwOBA (~0.297,
# CI: 0.288–0.306) while CWS pitchers allowed the highest
# (~0.322, CI: 0.313–0.332). These two intervals do not
# overlap, indicating a practically meaningful gap. CLE, KC,
# and MIN fall in between, with partially overlapping CIs.
#
# The power of the permutation test is estimated at 92.7%
# (Wald 95% CI: 0.911–0.943), meaning that if the true group
# means are as observed in our sample, the test would
# correctly reject H0 in approximately 93 out of 100
# repeated experiments. This high power is expected given
# the very large sample sizes (~6,200 per team).
#
# The agreement between the parametric (Welch's) and
# non-parametric (permutation) approaches strengthens the
# conclusion: AL Central pitching staffs differ meaningfully
# in the quality of contact they allow, with Detroit's
# rotation and bullpen suppressing hard contact most
# effectively and Chicago's allowing the most.
# ============================================================

