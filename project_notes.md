# MA416 Project Notes

## Dataset
- File: `2020_24Training.csv` — MLB Statcast pitch-by-pitch data (2020–2024), ~3.3M rows
- Filtered to: **2024 season only**
- Aggregated to: **at-bat level** via `build_atbat_2024.R` → output: `atbat_2024.csv`
- At-bat = last pitch of each (game_pk, at_bat_number) carries the outcome

## TTO (Time Through the Order)
- Computed from scratch (ignore original `n_thruorder_pitcher` column)
- Logic: per (game_pk, pitcher), sort by at_bat_number, count cumulative appearances of each batter
- TTO = 1st time facing that batter = 1, 2nd = 2, 3rd = 3, etc.
- Uses `ave(..., FUN = cumsum)` on key = paste(game_pk, pitcher, batter)

---

## Section 1 — One-Factor ANOVA

### Design
- **Groups (5):** AL Central teams — KC, DET, MIN, CWS, CLE
- **At-bat assignment:** pitcher's team (if inning_topbot == "Top" → home team is pitching, else away team)
- **Response:** `estimated_woba_using_speedangle` (xwOBA on contact)
- **Observations:** 33,037 total at-bats, ~6,500–6,700 per team (well balanced)
- **Multi-team pitchers:** 0 — no exclusions needed

### Independence Note
- Each at-bat belongs to exactly one pitcher's team → clean independence
- No pitcher threw for more than one AL Central team in 2024

---

## Section 2 — Repeated-Measures ANOVA

### Design
- **Subject:** pitcher
- **Within-subject factor:** TTO (1, 2, 3)
- **Response:** mean `estimated_woba_using_speedangle` per (pitcher, TTO)
- **Eligible pitchers:** 316 (those with at-bats in all 3 TTO levels)
- **At-bat counts (eligible pitchers only):**
  - TTO 1: 54,986
  - TTO 2: 45,902
  - TTO 3: 22,857
- **Hypothesis:** Does pitcher performance (xwOBA allowed) change as TTO increases?

### Key Assumptions to Test
- Normality of differences between TTO levels (SF or SW test)
- Sphericity (Mauchly's test — both O and B statistics)
- If sphericity violated: GG correction (eps_GG < 0.75) or HF correction (eps_GG > 0.75)

---

## Section 3 — Non-Parametric Test
- TBD

---

## atbat_2024.csv Column Reference

| Column | Description |
|---|---|
| game_pk | Game ID |
| game_date | Date of game |
| pitcher | Pitcher ID |
| player_name | Pitcher name |
| batter | Batter ID |
| stand | Batter handedness (L/R) |
| p_throws | Pitcher handedness (L/R) |
| home_team / away_team | Team abbreviations |
| inning / inning_topbot | Inning info |
| outs_when_up | Outs at start of AB |
| on_1b / on_2b / on_3b | Runner IDs (NA = empty) |
| at_bat_number | Sequential AB number within game |
| n_pitches | Pitches in the AB |
| TTO | Time through the order (1/2/3) |
| events | Raw outcome label |
| is_hit / is_single / is_double / is_triple / is_hr | Binary hit outcomes |
| is_walk / is_strikeout / is_hbp / is_error / is_sac | Binary other outcomes |
| launch_speed | Exit velocity (NA if no contact) |
| launch_angle | Launch angle (NA if no contact) |
| hit_distance_sc | Hit distance (NA if no contact) |
| estimated_ba_using_speedangle | xBA on contact (NA if no contact) |
| estimated_woba_using_speedangle | xwOBA on contact (NA if no contact) |
| estimated_slg_using_speedangle | xSLG on contact (NA if no contact) |
