# ============================================================
# MA416 Project - At-Bat Level Aggregation (2024 Season)
# ============================================================

filepath <- "C:/Users/aasmi/ma 416/Project/2020_24Training.csv"
outpath  <- "C:/Users/aasmi/ma 416/Project/atbat_2024.csv"

# ---- 1. Load & Filter to 2024 ----
cat("Reading data...\n")
raw <- read.csv(filepath, stringsAsFactors = FALSE)

cat("Filtering to 2024...\n")
raw <- raw[raw$game_year == 2024, ]
cat("Rows in 2024:", nrow(raw), "\n")

# ---- 2. Sort by game > at-bat > pitch ----
raw <- raw[order(raw$game_pk, raw$at_bat_number, raw$pitch_number), ]

# ---- 3. Get last pitch of each at-bat (carries outcome fields) ----
cat("Getting last pitch per at-bat...\n")
ab_key <- paste(raw$game_pk, raw$at_bat_number, sep = "_")
last_idx <- tapply(seq_len(nrow(raw)), ab_key, function(i) i[length(i)])
last_idx <- as.integer(unlist(last_idx))
ab <- raw[last_idx, ]

# n_pitches = pitch_number on last pitch (pitch_number counts from 1 within each AB)
ab$n_pitches <- ab$pitch_number

# ---- 4. Compute TTO ----
# Sort by game > pitcher > at_bat_number so cumsum runs in chronological order
cat("Computing TTO...\n")
ab <- ab[order(ab$game_pk, ab$pitcher, ab$at_bat_number), ]

tto_key    <- paste(ab$game_pk, ab$pitcher, ab$batter, sep = "_")
ab$TTO     <- as.integer(ave(rep(1L, nrow(ab)), tto_key, FUN = cumsum))

# ---- 5. Binary Outcome Indicators ----
ev <- ab$events

ab$is_single    <- as.integer(ev == "single")
ab$is_double    <- as.integer(ev == "double")
ab$is_triple    <- as.integer(ev == "triple")
ab$is_hr        <- as.integer(ev == "home_run")
ab$is_walk      <- as.integer(ev %in% c("walk", "intent_walk"))
ab$is_strikeout <- as.integer(ev %in% c("strikeout", "strikeout_double_play"))
ab$is_hbp       <- as.integer(ev == "hit_by_pitch")
ab$is_hit       <- as.integer(ev %in% c("single", "double", "triple", "home_run"))
ab$is_error     <- as.integer(ev == "field_error")
ab$is_sac       <- as.integer(ev %in% c("sac_fly", "sac_bunt",
                                         "sac_fly_double_play", "sac_bunt_double_play"))

# ---- 6. Select Final Columns ----
keep <- c(
  # Identifiers / metadata
  "game_pk", "game_date",
  "pitcher", "player_name",          # player_name = pitcher name
  "batter",
  "stand", "p_throws",               # batter/pitcher handedness (L/R)
  "home_team", "away_team",
  "inning", "inning_topbot",
  "outs_when_up",
  "on_1b", "on_2b", "on_3b",         # runner IDs (NA = base empty)
  "at_bat_number", "n_pitches",
  "TTO",

  # Raw outcome label
  "events",

  # Binary outcomes
  "is_hit", "is_single", "is_double", "is_triple", "is_hr",
  "is_walk", "is_strikeout", "is_hbp", "is_error", "is_sac",

  # Batted-ball / contact (NA when no contact)
  "launch_speed", "launch_angle", "hit_distance_sc",

  # Expected stats (NA when no contact)
  "estimated_ba_using_speedangle",
  "estimated_woba_using_speedangle",
  "estimated_slg_using_speedangle"
)

keep     <- keep[keep %in% names(ab)]
ab_final <- ab[, keep]

# Remove rows with no event (e.g. truncated_pa mid-game)
ab_final <- ab_final[!is.na(ab_final$events) & ab_final$events != "", ]

# ---- 7. Save ----
cat("Saving to", outpath, "\n")
write.csv(ab_final, outpath, row.names = FALSE)
cat("Done! At-bats written:", nrow(ab_final), "\n")

# ---- 8. Count pitchers with all 3 TTO levels ----
cat("\n--- TTO Coverage ---\n")
tto_per_pitcher <- tapply(ab_final$TTO, ab_final$pitcher, function(x) sort(unique(x)))
has_all_three   <- sapply(tto_per_pitcher, function(x) all(c(1, 2, 3) %in% x))

cat("Total unique pitchers:         ", length(tto_per_pitcher), "\n")
cat("Pitchers with TTO 1, 2, and 3:", sum(has_all_three), "\n")
cat("Pitchers missing TTO 3:        ", sum(sapply(tto_per_pitcher, function(x) !(3 %in% x))), "\n")

# At-bat counts per TTO group for the 316 eligible pitchers
eligible_pitchers <- names(has_all_three)[has_all_three]
ab_eligible <- ab_final[ab_final$pitcher %in% eligible_pitchers, ]
cat("\nAt-bats per TTO group (eligible pitchers only):\n")
print(table(ab_eligible$TTO))

# ---- 9. Handedness matchup group counts ----
cat("\n--- Handedness Matchup Groups ---\n")
ab_final$matchup <- paste(ab_final$p_throws, ab_final$stand, sep="-")
cat("Total at-bats:", nrow(ab_final), "\n")
print(table(ab_final$matchup))

# ---- 10. AL Central one-factor ANOVA setup ----
al_central <- c("KC", "DET", "MIN", "CWS", "CLE")
ab_alc <- ab_final[ab_final$home_team %in% al_central | ab_final$away_team %in% al_central, ]

# Assign at-bat to pitcher's team
ab_alc$pitcher_team <- ifelse(ab_alc$inning_topbot == "Top", ab_alc$home_team, ab_alc$away_team)
ab_alc <- ab_alc[ab_alc$pitcher_team %in% al_central, ]

# Check for pitchers who threw for more than one AL Central team
pitcher_teams <- tapply(ab_alc$pitcher_team, ab_alc$pitcher, function(x) unique(x))
multi_team    <- pitcher_teams[sapply(pitcher_teams, length) > 1]

cat("\n--- AL Central ANOVA ---\n")
cat("Pitchers who appeared for more than one AL Central team:", length(multi_team), "\n")
if(length(multi_team) > 0) print(multi_team)

# Exclude those pitchers
ab_alc_clean <- ab_alc[!ab_alc$pitcher %in% names(multi_team), ]

cat("\nTotal at-bats (after exclusions):", nrow(ab_alc_clean), "\n")
cat("At-bats per team:\n")
print(table(ab_alc_clean$pitcher_team))
