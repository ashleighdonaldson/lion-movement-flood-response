# Run from a folder containing the original data files listed below.
# Required: Combined collar data with weather astronomy and movement(in).csv,
# elevation.tif, rivers.shp (+ sidecar files), Finalweatherstation.csv,
# gps_with_road_metrics.gpkg, DynamicWorld_WoodyBinary_2025_Satara.tif,
# DynamicWorld_WoodyProbability_2025_Satara.tif,
# DynamicWorld_TreeProbability_2025_Satara.tif,
# DynamicWorld_ShrubProbability_2025_Satara.tif.

required_files <- c(
  'Combined collar data with weather astronomy and movement(in).csv',
  'elevation.tif', 'rivers.shp', 'Finalweatherstation.csv',
  'gps_with_road_metrics.gpkg',
  'DynamicWorld_WoodyBinary_2025_Satara.tif',
  'DynamicWorld_WoodyProbability_2025_Satara.tif',
  'DynamicWorld_TreeProbability_2025_Satara.tif',
  'DynamicWorld_ShrubProbability_2025_Satara.tif'
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) stop('Missing input files: ', paste(missing_files, collapse=', '))

# Set optional analyses to FALSE if you only need a subset.
RUN_WOODY_SELECTION <- TRUE
RUN_ROAD_SELECTION <- TRUE
RUN_VEGETATION_MAPS <- TRUE
RUN_RAINFALL_FIGURE <- TRUE

library(dplyr)
library(lubridate)
library(sf)
library(terra)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)
library(nlme)
library(amt)
library(survival)

# GPS data and flood periods
lions <- c('8387', '8388', '8389', '8391', '8393', '8395')
flood_start <- as.Date('2026-01-10')
flood_end   <- as.Date('2026-01-19')
before_start <- flood_start - 10
after_end <- flood_end + 10
period_levels <- c('Before', 'During', 'After')

gps <- read.csv('Combined collar data with weather astronomy and movement(in).csv') %>%
  mutate(
    TagId = as.character(TagId),
    Timestamp = ymd_hms(sub(' SAST$', '', Timestamp_SAST), tz = 'Africa/Johannesburg'),
    Date = as.Date(Timestamp, tz = 'Africa/Johannesburg'),
    Period = case_when(
      Date >= before_start & Date < flood_start ~ 'Before',
      Date >= flood_start & Date <= flood_end ~ 'During',
      Date > flood_end & Date <= after_end ~ 'After',
      TRUE ~ NA_character_
    ),
    Period = factor(Period, levels = period_levels)
  ) %>%
  filter(TagId %in% lions, !is.na(Period), !is.na(Timestamp),
         !is.na(longitude), !is.na(latitude)) %>%
  arrange(TagId, Timestamp)

print(gps %>% count(TagId, Period))

# Project coordinates once, for all movement and spatial analyses.
gps_sf <- st_as_sf(gps, coords = c('longitude', 'latitude'), crs = 4326,
                   remove = FALSE) %>% st_transform(32736)
xy <- st_coordinates(gps_sf)
gps_xy <- gps_sf %>% st_drop_geometry() %>%
  mutate(X = xy[, 1], Y = xy[, 2])

# Hourly movement
hourly_movement <- gps_xy %>%
  group_by(TagId, Period) %>%
  arrange(Timestamp, .by_group = TRUE) %>%
  mutate(
    TimeDifference = as.numeric(difftime(Timestamp, lag(Timestamp), units = 'hours')),
    Distance_km = sqrt((X - lag(X))^2 + (Y - lag(Y))^2) / 1000,
    Hour = hour(Timestamp),
    SameDayStep = Date == lag(Date)
  ) %>%
  ungroup() %>%
  filter(TimeDifference >= 0.75, TimeDifference <= 1.25,
         is.finite(Distance_km))
print(hourly_movement %>% count(TagId, Period))

# Daily movement
gps_steps <- gps_xy %>%
  group_by(TagId) %>%
  arrange(Timestamp, .by_group = TRUE) %>%
  mutate(
    Previous_Time = lag(Timestamp),
    Previous_Date = lag(Date),
    TimeDifference = as.numeric(difftime(Timestamp, Previous_Time,
                                         units = 'hours')),
    Distance_km = sqrt((X - lag(X))^2 + (Y - lag(Y))^2) / 1000,
    SameDayStep = Date == Previous_Date
  ) %>%
  ungroup() %>%
  filter(!is.na(TimeDifference), TimeDifference > 0,
         is.finite(Distance_km))

daily_movement <- gps_steps %>%
  filter(SameDayStep, TimeDifference <= 6) %>%
  group_by(TagId, Date, Period) %>%
  summarise(
    DailyDistance_km = sum(Distance_km),
    N_Steps = n(),
    Hours_Covered = sum(TimeDifference),
    Max_Gap = max(TimeDifference),
    .groups = 'drop'
  ) %>%
  mutate(
    TagId = factor(TagId),
    Period = factor(Period, levels = period_levels),
    DayNumber = as.numeric(Date)
  ) %>%
  arrange(TagId, Date)

print(daily_movement %>% count(TagId, Period))
print(daily_movement %>%
        group_by(Period) %>%
        summarise(n_days = n(), mean_distance = mean(DailyDistance_km),
                  mean_hours = mean(Hours_Covered),
                  min_hours = min(Hours_Covered), .groups = 'drop'))

# Gaussian LME with period-specific variance
movement_model_var <- lme(
  DailyDistance_km ~ Period,
  random = ~1 | TagId,
  correlation = corCAR1(form = ~DayNumber | TagId),
  weights = varIdent(form = ~1 | Period),
  data = daily_movement,
  method = "REML"
)

summary(movement_model_var)

anova(movement_model_var)

pairs(emmeans(movement_model_var, ~Period), adjust = "none")

movement_plot <- ggplot(daily_movement, aes(Period, DailyDistance_km)) +
  geom_boxplot(fill = 'grey85', outlier.shape = NA) +
  geom_jitter(width = 0.12, alpha = 0.35, size = 1.5) +
  labs(x = 'Flood period', y = 'Daily distance travelled (km)') +
  theme_classic(base_size = 12)
print(movement_plot)

daily_movement <- daily_movement %>%
  mutate(
    RowID = row_number(),
    PointColour = "Other"
  )

top4_ids <- daily_movement %>%
  filter(Period == "During") %>%
  slice_max(DailyDistance_km, n = 4, with_ties = FALSE) %>%
  pull(RowID)

daily_movement <- daily_movement %>%
  mutate(
    PointColour = if_else(RowID %in% top4_ids, "Top 4", "Other")
  )

movement_plot <- ggplot(
  daily_movement,
  aes(Period, DailyDistance_km)
) +
  geom_boxplot(
    fill = "grey85",
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(colour = PointColour),
    width = 0.12,
    alpha = 0.6,
    size = 1.5
  ) +
  scale_colour_manual(
    values = c(
      "Other" = "grey50",
      "Top 4" = "red"
    )
  ) +
  labs(
    x = "Flood period",
    y = "Daily distance travelled (km)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 17),
    axis.title.y = element_text(size = 17),
    axis.text.x = element_text(size = 15, colour = "black"),
    axis.text.y = element_text(size = 15, colour = "black")
  )

print(movement_plot)


# Elevation and river use

elevation <- terra::rast('elevation.tif')

rivers <- sf::st_read('rivers.shp', quiet = TRUE) %>%
  sf::st_transform(32736)

river_network <- sf::st_union(sf::st_geometry(rivers))

gps_sf <- sf::st_transform(gps_sf, 32736)

gps_sf <- sf::st_as_sf(
  gps_xy,
  coords = c("X", "Y"),
  crs = 32736,
  remove = FALSE
)

gps_xy$Elevation_m <- terra::extract(
  elevation,
  terra::vect(sf::st_transform(gps_sf, terra::crs(elevation)))
)[[2]]

# Distance to nearest river
gps_xy$Distance_to_river_m <- as.numeric(
  sf::st_distance(
    gps_sf,
    river_network
  )
)

# Descriptive use models
# Calculate daily mean elevation per lion
daily_elevation <- gps_xy %>%
  group_by(TagId, Date, Period) %>%
  summarise(
    Elevation_m = mean(Elevation_m, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    TagId = factor(TagId),
    Period = factor(Period, levels = c("Before", "During", "After")),
    DayNumber = as.numeric(as.Date(Date))
  )

elevation_daily_model_var <- lme(
  Elevation_m ~ Period,
  random = ~1 | TagId,
  correlation = corCAR1(form = ~DayNumber | TagId),
  weights = varIdent(form = ~1 | Period),
  data = daily_elevation,
  method = "REML"
)

anova(elevation_daily_model_var)
summary(elevation_daily_model_var)
summary(emmeans(elevation_daily_model_var, ~Period))

pairs(
  emmeans(elevation_daily_model_var, ~Period),
  adjust = "tukey"
)

daily_river <- gps_xy %>%
  group_by(TagId, Date, Period) %>%
  summarise(
    Distance_to_river_m = mean(Distance_to_river_m, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    TagId = factor(TagId),
    Period = factor(Period, levels = c("Before", "During", "After")),
    DayNumber = as.numeric(as.Date(Date))
  )

river_daily_model <- lme(
  log1p(Distance_to_river_m) ~ Period,
  random = ~1 | TagId,
  correlation = corCAR1(form = ~DayNumber | TagId),
  data = daily_river,
  method = "REML"
)

anova(river_daily_model)
summary(river_daily_model)
summary(emmeans(river_daily_model, ~Period))

pairs(
  emmeans(river_daily_model, ~Period),
  adjust = "tukey"
)

daily_river %>%
  group_by(Period) %>%
  summarise(
    Mean = mean(Distance_to_river_m, na.rm = TRUE),
    SD = sd(Distance_to_river_m, na.rm = TRUE),
    SE = SD / sqrt(sum(!is.na(Distance_to_river_m))),
    .groups = "drop"
  )

# Step selection: elevation and rivers

library(dplyr)
library(amt)
library(sf)
library(terra)
library(survival)

make_ssf_steps <- function(df) {
  
  if (nrow(df) < 5) return(NULL)
  
  trk <- amt::make_track(
    df %>% arrange(Timestamp),
    X, Y, Timestamp,
    crs = 32736,
    all_cols = TRUE
  )
  
  trk <- amt::track_resample(
    trk,
    rate = lubridate::hours(1),
    tolerance = lubridate::minutes(15)
  )
  
  steps <- amt::steps_by_burst(trk)
  
  if (nrow(steps) < 3) return(NULL)
  
  amt::random_steps(steps, n_control = 15)
}

set.seed(123)

# Generate observed and available steps
ssf_groups <- gps_xy %>%
  group_by(TagId, Period) %>%
  group_split()

ssf_list <- lapply(seq_along(ssf_groups), function(i) {
  
  df <- ssf_groups[[i]]
  steps <- make_ssf_steps(df)
  
  if (is.null(steps)) return(NULL)
  
  steps %>%
    mutate(
      TagId = df$TagId[1],
      Period = as.character(df$Period[1]),
      Step_ID_unique = paste(i, step_id_, sep = "_")
    )
})

ssf_list <- Filter(Negate(is.null), ssf_list)

if (!length(ssf_list)) {
  stop("No valid iSSF steps generated.")
}

ssf_steps <- bind_rows(ssf_list)

step_points <- sf::st_as_sf(
  ssf_steps,
  coords = c("x2_", "y2_"),
  crs = 32736,
  remove = FALSE
)

step_elevation <- terra::extract(
  elevation,
  terra::vect(
    sf::st_transform(step_points, terra::crs(elevation))
  )
)[[2]]

river_projected <- sf::st_transform(
  river_network,
  sf::st_crs(step_points)
)

nearest_river <- sf::st_nearest_feature(
  step_points,
  river_projected
)

step_river_distance <- as.numeric(
  sf::st_distance(
    step_points,
    river_projected[nearest_river, ],
    by_element = TRUE
  )
)

ssf_steps <- ssf_steps %>%
  mutate(
    Elevation_m = step_elevation,
    Distance_to_river_m = step_river_distance,
    Period = factor(
      Period,
      levels = c("Before", "During", "After")
    ),
    Elevation_scaled = as.numeric(scale(Elevation_m)),
    River_scaled = as.numeric(
      scale(log1p(Distance_to_river_m))
    ),
    LogStepLength = log(pmax(sl_, 0.001)),
    CosTurnAngle = cos(ta_)
  ) %>%
  group_by(Step_ID_unique) %>%
  filter(
    all(
      is.finite(Elevation_scaled) &
        is.finite(River_scaled) &
        is.finite(sl_) &
        is.finite(LogStepLength) &
        is.finite(CosTurnAngle)
    )
  ) %>%
  ungroup()

issf_model <- survival::clogit(
  case_ ~ Elevation_scaled * Period + River_scaled * Period +
    sl_ + LogStepLength + CosTurnAngle +
    strata(Step_ID_unique) + cluster(TagId),
  data = ssf_steps, method = 'efron'
)
print(summary(issf_model))
saveRDS(issf_model, 'Lion_iSSF_model.rds')


# Elevation and river distance figures 
publication_theme <- ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    legend.position = "none",
    axis.title = ggplot2::element_text(size = 16, colour = "black"),
    axis.text = ggplot2::element_text(size = 14, colour = "black")
  )

elevation_plot <- ggplot2::ggplot(gps_xy,
    ggplot2::aes(x = Period, y = Elevation_m, fill = Period)) +
  ggplot2::geom_boxplot(width = 0.55, alpha = 0.6,
                        outlier.shape = NA, colour = "black") +
  ggplot2::geom_jitter(width = 0.10, size = 0.7, alpha = 0.15,
                       colour = "black") +
  ggplot2::scale_fill_manual(values = c(Before = "grey", During = "grey", After = "grey")) +
  ggplot2::scale_y_continuous(breaks = seq(200, 350, 50)) +
  ggplot2::coord_cartesian(ylim = c(200, 350)) +
  ggplot2::labs(x = "Flood period", y = "Elevation (m a.s.l.)") +
  publication_theme

river_plot <- ggplot2::ggplot(gps_xy,
    ggplot2::aes(x = Period, y = Distance_to_river_m, fill = Period)) +
  ggplot2::geom_boxplot(width = 0.55, alpha = 0.6,
                        outlier.shape = NA, colour = "black") +
  ggplot2::geom_jitter(width = 0.10, size = 0.7, alpha = 0.15,
                       colour = "black") +
  ggplot2::scale_fill_manual(values = c(Before = "grey", During = "grey", After = "grey")) +
  ggplot2::scale_y_continuous(breaks = scales::breaks_width(500)) +
  ggplot2::labs(x = "Flood period", y = "Distance to nearest river (m)") +
  publication_theme

print(elevation_plot)
print(river_plot)
ggplot2::ggsave("Lion_Elevation_Use.svg", elevation_plot, width = 8, height = 6)
ggplot2::ggsave("Lion_River_Distance_Use.svg", river_plot, width = 8, height = 6)

# Diel movement
for (pkg in c('dplyr','ggplot2','mgcv','nlme','MASS','tidyr')) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop('Please install: ', pkg)
}

set.seed(123)
period_levels <- c('Before','During','After')

required <- c('TagId','Period','Timestamp','Distance_km')
if (!exists('hourly_movement')) stop('Run the movement-preparation section first: hourly_movement is missing.')
if (!all(required %in% names(hourly_movement))) {
  stop('Missing hourly_movement columns: ', paste(setdiff(required,names(hourly_movement)),collapse=', '))
}

diel_data_all <- hourly_movement |>
  dplyr::mutate(
    TagId = factor(TagId),
    Period = factor(as.character(Period), levels = period_levels),
    HourOfDay = as.numeric(format(Timestamp, '%H', tz = 'Africa/Johannesburg')) +
      as.numeric(format(Timestamp, '%M', tz = 'Africa/Johannesburg')) / 60 +
      as.numeric(format(Timestamp, '%S', tz = 'Africa/Johannesburg')) / 3600,
    TimeHours = as.numeric(Timestamp) / 3600
  ) |>
  dplyr::filter(!is.na(Period), is.finite(Distance_km),
                is.finite(HourOfDay), is.finite(TimeHours)) |>
  dplyr::arrange(TagId, Timestamp)

cat('Total hourly steps:', nrow(diel_data_all), '\n')
cat('Zero-distance steps excluded from Gamma GAMM:', sum(diel_data_all$Distance_km == 0), '\n')
cat('Negative-distance steps excluded:', sum(diel_data_all$Distance_km < 0), '\n')

diel_data <- dplyr::filter(diel_data_all, Distance_km > 0)
stopifnot(nrow(diel_data) > 0)
print(dplyr::count(diel_data, TagId, Period))

diel_gamm <- mgcv::gamm(
  Distance_km ~ Period + s(HourOfDay, by = Period, bs = 'cc', k = 18),
  random = list(TagId = ~ 1),
  correlation = nlme::corCAR1(form = ~ TimeHours | TagId),
  family = Gamma(link = 'log'),
  data = diel_data,
  knots = list(HourOfDay = c(0, 24)),
  method = 'REML',
  control = nlme::lmeControl(opt = 'optim', maxIter = 200, msMaxIter = 200,
                             niterEM = 0, returnObject = TRUE)
)

print(summary(diel_gamm$gam))
print(summary(diel_gamm$lme))

prediction_hours <- seq(0, 24, by = 0.05)
diel_prediction <- expand.grid(
  HourOfDay = prediction_hours,
  Period = factor(period_levels, levels = period_levels)
) |>
  dplyr::arrange(Period, HourOfDay)

Xpred <- predict(diel_gamm$gam, newdata = diel_prediction, type = 'lpmatrix')
beta <- stats::coef(diel_gamm$gam)
Vbeta <- stats::vcov(diel_gamm$gam)
eta <- as.vector(Xpred %*% beta)
se_eta <- sqrt(pmax(0, rowSums((Xpred %*% Vbeta) * Xpred)))

diel_prediction <- diel_prediction |>
  dplyr::mutate(
    Predicted_km = exp(eta),
    Lower95_km = exp(eta - 1.96 * se_eta),
    Upper95_km = exp(eta + 1.96 * se_eta)
  )

diel_plot <- ggplot2::ggplot(
  diel_prediction,
  ggplot2::aes(x = HourOfDay, y = Predicted_km,
               linetype = Period, group = Period)
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = Lower95_km, ymax = Upper95_km, fill = Period),
    linetype = 0, alpha = 0.18
  ) +
  ggplot2::geom_line(colour = 'black', linewidth = 0.9) +
  ggplot2::scale_linetype_manual(
    values = c(Before = 'solid', During = 'longdash', After = 'dotted')
  ) +
  ggplot2::scale_fill_manual(
    values = c(Before = 'grey50', During = 'grey65', After = 'grey80')
  ) +
  ggplot2::scale_x_continuous(
    breaks = c(0, 6, 12, 18, 24),
    labels = c('00:00','06:00','12:00','18:00','00:00'),
    limits = c(0,24), expand = c(0,0)
  ) +
  ggplot2::labs(x = 'Time of day', y = 'Predicted hourly movement (km)') +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    legend.position = 'bottom',
    legend.title = ggplot2::element_blank(),
    axis.title = ggplot2::element_text(size = 16, colour = 'black'),
    axis.text = ggplot2::element_text(size = 14, colour = 'black'),
    legend.text = ggplot2::element_text(size = 12)
  )
print(diel_plot)
ggplot2::ggsave('Lion_Diel_GAMM_95CI.svg', diel_plot,
                width = 9, height = 6, units = 'in')

peak_windows <- list(Morning = c(0,8), Evening = c(17,24))
extract_peaks <- function(values, period_name, window_name) {
  rows <- which(as.character(diel_prediction$Period) == period_name &
                  diel_prediction$HourOfDay >= peak_windows[[window_name]][1] &
                  diel_prediction$HourOfDay <= peak_windows[[window_name]][2])
  winner <- rows[which.max(values[rows])]
  c(PeakHour = diel_prediction$HourOfDay[winner],
    PeakDistance_km = values[winner])
}

peak_estimates <- do.call(rbind, lapply(period_levels, function(p) {
  do.call(rbind, lapply(names(peak_windows), function(w) {
    out <- extract_peaks(diel_prediction$Predicted_km, p, w)
    data.frame(Period = p, Window = w,
               PeakHour = unname(out['PeakHour']),
               PeakDistance_km = unname(out['PeakDistance_km']))
  }))
}))
print(peak_estimates, row.names = FALSE)

n_sim <- 5000L
sim_beta <- MASS::mvrnorm(n = n_sim, mu = beta, Sigma = Vbeta)
eta_sim <- Xpred %*% t(sim_beta)  # prediction rows x simulations

peak_draws <- do.call(rbind, lapply(period_levels, function(p) {
  do.call(rbind, lapply(names(peak_windows), function(w) {
    rows <- which(as.character(diel_prediction$Period) == p &
                    diel_prediction$HourOfDay >= peak_windows[[w]][1] &
                    diel_prediction$HourOfDay <= peak_windows[[w]][2])
    best <- max.col(t(eta_sim[rows, , drop = FALSE]), ties.method = 'first')
    data.frame(Simulation = seq_len(n_sim), Period = p, Window = w,
               PeakHour = diel_prediction$HourOfDay[rows][best])
  }))
}))

peak_uncertainty <- peak_draws |>
  dplyr::group_by(Period, Window) |>
  dplyr::summarise(
    MedianHour = stats::median(PeakHour),
    Lower95Hour = unname(stats::quantile(PeakHour, 0.025)),
    Upper95Hour = unname(stats::quantile(PeakHour, 0.975)),
    .groups = 'drop'
  )
print(peak_uncertainty)

peak_comparisons <- peak_draws |>
  tidyr::pivot_wider(names_from = Period, values_from = PeakHour) |>
  dplyr::mutate(
    `During - Before` = During - Before,
    `After - Before` = After - Before,
    `After - During` = After - During
  ) |>
  tidyr::pivot_longer(cols = c('During - Before','After - Before','After - During'),
                      names_to = 'Comparison', values_to = 'DifferenceHours') |>
  dplyr::group_by(Window, Comparison) |>
  dplyr::summarise(
    MedianDifferenceHours = stats::median(DifferenceHours),
    Lower95 = unname(stats::quantile(DifferenceHours, 0.025)),
    Upper95 = unname(stats::quantile(DifferenceHours, 0.975)),
    ExcludesZero = Lower95 > 0 | Upper95 < 0,
    .groups = 'drop'
  )
print(peak_comparisons, n = Inf)

utils::write.csv(diel_prediction, 'Lion_Diel_GAMM_Predictions.csv', row.names = FALSE)
utils::write.csv(peak_estimates, 'Lion_Diel_Peak_Estimates.csv', row.names = FALSE)
utils::write.csv(peak_uncertainty, 'Lion_Diel_Peak_Uncertainty.csv', row.names = FALSE)
utils::write.csv(peak_comparisons, 'Lion_Diel_Peak_Comparisons.csv', row.names = FALSE)
saveRDS(diel_gamm, 'Lion_Diel_GAMM_Model.rds')


# Rainfall and hourly movement
rain_movement_car <- hourly_movement %>%
  filter(is.finite(rain), rain >= 0) %>%
  mutate(
    TagId = factor(TagId),
    LionPeriod = interaction(TagId, Period, drop = TRUE),
    TimeHours = as.numeric(Timestamp) / 3600,
    SinHour = sin(2 * pi * Hour / 24),
    CosHour = cos(2 * pi * Hour / 24),
    LogDistance = log1p(Distance_km)
  ) %>% arrange(LionPeriod, Timestamp)

rain_model_car <- nlme::lme(
  fixed = LogDistance ~ rain * Period + SinHour + CosHour,
  random = ~ 1 | TagId,
  correlation = nlme::corCAR1(form = ~ TimeHours | TagId/LionPeriod),
  data = rain_movement_car,
  method = 'ML',
  control = nlme::lmeControl(opt = 'optim')
)
summary(rain_model_car)
coef(rain_model_car$modelStruct$corStruct, unconstrained = FALSE)
rain_trends_car <- emtrends(rain_model_car, ~ Period, var = 'rain')
summary(rain_trends_car, infer = TRUE)
pairs(rain_trends_car)

prediction_data <- rain_movement_car %>% group_by(Period) %>%
  summarise(min_rain = min(rain), max_rain = max(rain), .groups = 'drop') %>%
  group_by(Period) %>%
  reframe(rain = seq(min_rain, max_rain, length.out = 100)) %>%
  mutate(SinHour = 0, CosHour = 0)
prediction_data$PredictedLog <- predict(rain_model_car,
                                        newdata = prediction_data, level = 0)
X <- model.matrix(~ rain * Period + SinHour + CosHour, prediction_data)
V <- vcov(rain_model_car)
prediction_data <- prediction_data %>%
  mutate(
    SE = sqrt(rowSums((X %*% V) * X)),
    Predicted_km = expm1(PredictedLog),
    Lower_km = expm1(PredictedLog - 1.96 * SE),
    Upper_km = expm1(PredictedLog + 1.96 * SE)
  )


rain_plot_during <- ggplot2::ggplot() +
  ggplot2::geom_point(
    data = dplyr::filter(rain_movement_car, Period == "During"),
    ggplot2::aes(x = rain, y = Distance_km),
    colour = "grey45", alpha = 0.12, size = 0.7
  ) +
  ggplot2::geom_ribbon(
    data = dplyr::filter(prediction_data, Period == "During"),
    ggplot2::aes(x = rain, ymin = Lower_km, ymax = Upper_km),
    fill = "grey75", alpha = 0.5
  ) +
  ggplot2::geom_line(
    data = dplyr::filter(prediction_data, Period == "During"),
    ggplot2::aes(x = rain, y = Predicted_km),
    colour = "black", linewidth = 0.9
  ) +
  ggplot2::scale_x_continuous(breaks = scales::breaks_width(2)) +
  ggplot2::labs(x = "Hourly rainfall (mm)", y = "Hourly movement distance (km)") +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 16),
    axis.text = ggplot2::element_text(size = 14, colour = "black")
  )
print(rain_plot_during)
ggplot2::ggsave("Lion_Hourly_Rainfall_Movement_During.svg",
                rain_plot_during, width = 9, height = 6)

# Road distance
if (RUN_ROAD_SELECTION) {
for (p in c('sf','terra','dplyr','ggplot2','nlme','emmeans','survival')) {
  if (!requireNamespace(p, quietly=TRUE)) stop('Install package: ',p)
}
stopifnot(exists('gps_xy'), exists('ssf_steps'))
stopifnot(all(c('X','Y','TagId','Date','Period') %in% names(gps_xy)))
stopifnot(all(c('x2_','y2_','Step_ID_unique','case_','Period','TagId',
                'Elevation_scaled','River_scaled','sl_','LogStepLength','CosTurnAngle') %in% names(ssf_steps)))

# Read road network
road_file <- 'gps_with_road_metrics.gpkg'
if (!file.exists(road_file)) stop('File not found: ', road_file,
                                   '. Place it in your working directory or edit road_file.')
if (!'roads_used' %in% sf::st_layers(road_file)$name)
  stop('GeoPackage does not contain the roads_used layer.')
roads_sf <- sf::st_read(road_file, layer='roads_used', quiet=TRUE)
roads_sf <- sf::st_transform(roads_sf,32736)
roads_sf <- sf::st_collection_extract(sf::st_make_valid(roads_sf),'LINESTRING',warn=FALSE)
roads_sf <- roads_sf[!sf::st_is_empty(roads_sf),]
if(nrow(roads_sf)==0) stop('No valid road lines in ',road_file)

if ('analysis_area' %in% sf::st_layers(road_file)$name) {
  area_sf <- sf::st_read(road_file, layer='analysis_area', quiet=TRUE) |>
    sf::st_transform(32736)
  xy_check <- rbind(as.matrix(gps_xy[,c('X','Y')]),
                    as.matrix(ssf_steps[,c('x2_','y2_')]))
  xy_check <- xy_check[stats::complete.cases(xy_check),,drop=FALSE]
  check_pts <- sf::st_as_sf(data.frame(x=xy_check[,1],y=xy_check[,2]),
                             coords=c('x','y'),crs=32736)
  inside <- lengths(sf::st_intersects(check_pts,area_sf))>0
  message(sprintf('GPS and available endpoints inside GeoPackage analysis area: %.1f%%',
                  100*mean(inside)))
  if(any(!inside)) warning('Some endpoints fall outside the supplied road coverage. ',
                           'Check road network completeness before interpreting distances.')
}

# Road raster
road_raster_file <- 'Satara_Roads_10m_GPKG.tif'
if (!file.exists(road_raster_file)) {
  r_template <- terra::rast(terra::ext(terra::vect(roads_sf)),
                            resolution=10,crs='EPSG:32736')
  r_roads <- terra::rasterize(terra::vect(roads_sf),r_template,field=1,
                             background=NA,touches=TRUE)
  terra::writeRaster(r_roads,road_raster_file,overwrite=TRUE,
                     gdal=c('COMPRESS=LZW'))
}

road_union <- sf::st_union(sf::st_geometry(roads_sf))
road_distance <- function(x,y) {
  pts <- sf::st_as_sf(data.frame(x=x,y=y),coords=c('x','y'),crs=32736)
  as.numeric(sf::st_distance(sf::st_geometry(pts),road_union))
}

# GPS distances to roads
road_gps <- gps_xy |>
  dplyr::mutate(Distance_to_road_m=road_distance(X,Y),
                DayNumber=as.numeric(as.Date(Date)),
                TagId=factor(TagId),
                Period=factor(as.character(Period),levels=c('Before','During','After')))

road_daily <- road_gps |>
  dplyr::group_by(TagId,Date,Period,DayNumber) |>
  dplyr::summarise(Distance_to_road_m=mean(Distance_to_road_m,na.rm=TRUE),
                   N_GPS=dplyr::n(),.groups='drop') |>
  dplyr::filter(is.finite(Distance_to_road_m))

road_daily_model <- nlme::lme(
  fixed=log1p(Distance_to_road_m) ~ Period,
  random=~1|TagId,
  correlation=nlme::corCAR1(form=~DayNumber|TagId),
  data=road_daily, method='REML',na.action=na.omit)
cat('\n--- ROAD DISTANCE: period effect ---\n')
print(anova(road_daily_model))
print(summary(road_daily_model))
print(emmeans::emmeans(road_daily_model,~Period))
print(pairs(emmeans::emmeans(road_daily_model,~Period),adjust='tukey'))

road_distance_summary <- road_daily |>
  dplyr::group_by(Period) |>
  dplyr::summarise(N_LionDays=dplyr::n(),Mean_m=mean(Distance_to_road_m),
                   SD_m=sd(Distance_to_road_m),Median_m=median(Distance_to_road_m),
                   .groups='drop')
print(road_distance_summary)
write.csv(road_distance_summary,'Lion_Road_Distance_Daily_Summary.csv',row.names=FALSE)

road_plot <- ggplot2::ggplot(
  road_gps,
  ggplot2::aes(x = Period, y = Distance_to_road_m, fill = Period)
) +
  ggplot2::geom_boxplot(
    width = 0.55,
    alpha = 0.6,
    outlier.shape = NA,
    colour = "black"
  ) +
  ggplot2::geom_jitter(
    width = 0.10,
    size = 0.7,
    alpha = 0.15,
    colour = "black"
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Before" = "grey",
      "During" = "grey",
      "After" = "grey"
    )
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq(0, 5000, by = 1000)
  ) +
  ggplot2::coord_cartesian(ylim = c(0, 5000)) +
  ggplot2::labs(
    x = "Flood period",
    y = "Distance to nearest road (m)"
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    legend.position = "none",
    axis.title.x = ggplot2::element_text(size = 16, colour = "black"),
    axis.title.y = ggplot2::element_text(size = 16, colour = "black"),
    axis.text.x = ggplot2::element_text(size = 14, colour = "black"),
    axis.text.y = ggplot2::element_text(size = 14, colour = "black")
  )

road_plot
print(road_plot)
ggplot2::ggsave('Lion_Road_Distance_Use.svg',road_plot,width=8,height=5)


# Road step selection 
road_ssf <- ssf_steps
road_ssf$Period <- factor(as.character(road_ssf$Period), levels=c('Before','During','After'))
road_ssf$Distance_to_road_m <- road_distance(road_ssf$x2_,road_ssf$y2_)
road_ssf$Road_scaled <- as.numeric(scale(log1p(road_ssf$Distance_to_road_m)))
road_ssf <- road_ssf |>
  dplyr::group_by(Step_ID_unique) |>
  dplyr::filter(all(is.finite(Road_scaled))) |>
  dplyr::ungroup()
stopifnot(nrow(road_ssf)>0)
road_sets <- road_ssf |>
  dplyr::group_by(Step_ID_unique) |>
  dplyr::summarise(n=dplyr::n(),used=sum(as.integer(case_)),.groups='drop')
if (!all(road_sets$used == 1)) stop('Some step sets do not contain exactly one used step.')
if (length(unique(road_sets$n)) != 1) warning('Step sets have varying numbers of steps; inspect road_sets.')

road_issf_model <- survival::clogit(
  case_ ~ Road_scaled*Period + sl_ + LogStepLength + CosTurnAngle +
    strata(Step_ID_unique) + cluster(TagId),
  data=road_ssf, method='efron')
cat('\n--- ROAD iSSF ---\n')
print(summary(road_issf_model))

road_coefs <- stats::coef(road_issf_model)
print(road_coefs[grepl('Road_scaled',names(road_coefs))])
road_int <- grep('Road_scaled:Period|Period.*:Road_scaled',names(road_coefs),value=TRUE)
if(length(road_int)==2) {
  b <- road_coefs[road_int]
  V <- stats::vcov(road_issf_model)[road_int,road_int,drop=FALSE]
  W <- as.numeric(t(b)%*%solve(V,b))
  cat('\nJoint test of road-selection changes: chi-square=',W,
      ', df=2, p=',stats::pchisq(W,2,lower.tail=FALSE),'\n')
} else warning('Could not identify both road-by-period interaction terms.')

all_b <- stats::coef(road_issf_model)
all_V <- stats::vcov(road_issf_model)
road_effects <- do.call(rbind,lapply(c('Before','During','After'),function(p) {
  L <- setNames(rep(0,length(all_b)),names(all_b))
  L['Road_scaled'] <- 1
  if(p!='Before') {
    nm <- intersect(c(paste0('Road_scaled:Period',p),
                      paste0('Period',p,':Road_scaled')),names(all_b))
    if(length(nm)!=1) stop('Missing road interaction for ',p)
    L[nm] <- 1
  }
  est <- sum(L*all_b)
  se <- sqrt(as.numeric(t(L)%*%all_V%*%L))
  data.frame(Period=p,Beta=est,SE=se,Lower95=est-1.96*se,
             Upper95=est+1.96*se,SelectionRatio=exp(est))
}))
print(road_effects,row.names=FALSE)
write.csv(road_effects,'Lion_Road_iSSF_Period_Effects.csv',row.names=FALSE)
saveRDS(road_issf_model,'Lion_Road_iSSF_Model.rds')

}

# Woody vegetation selection
if (RUN_WOODY_SELECTION) {

packages <- c(
  "dplyr",
  "tidyr",
  "lubridate",
  "readr",
  "readxl",
  "sf",
  "terra",
  "amt",
  "survival",
  "ggplot2",
  "purrr"
)

missing_packages <- setdiff(
  packages,
  rownames(installed.packages())
)

if (length(missing_packages) > 0) {
  stop("Install missing packages: ", paste(missing_packages, collapse=", "))
}

invisible(lapply(packages, library, character.only = TRUE))

set.seed(123)

n_available <- 20
target_crs <- "EPSG:32736"

# Flood periods: 10 days each

before_start <- as.Date("2025-12-31")
before_end   <- as.Date("2026-01-09")

during_start <- as.Date("2026-01-10")
during_end   <- as.Date("2026-01-19")

after_start  <- as.Date("2026-01-20")
after_end    <- as.Date("2026-01-29")

period_levels <- c("Before", "During", "After")


gps_files <- list.files(
  pattern = paste0(
    "^Combined collar data with weather astronomy ",
    "and movement\\(in\\)\\.(csv|xlsx|xls|rds)$"
  ),
  ignore.case = TRUE
)

print(gps_files)

if (length(gps_files) != 1) {
  stop(
    "Expected one original GPS file. ",
    "Check its filename using list.files()."
  )
}

gps_file <- gps_files[1]

extension <- tolower(tools::file_ext(gps_file))

if (extension == "csv") {
  
  gps_raw <- readr::read_csv(
    gps_file,
    show_col_types = FALSE
  )
  
} else if (extension %in% c("xlsx", "xls")) {
  
  gps_raw <- readxl::read_excel(gps_file)
  
} else if (extension == "rds") {
  
  gps_raw <- readRDS(gps_file)
  
}

cat("\nOriginal GPS observations:", nrow(gps_raw), "\n")


required <- c(
  "TagId",
  "Timestamp",
  "latitude",
  "longitude"
)

if (!all(required %in% names(gps_raw))) {
  stop("One or more required GPS columns are missing.")
}

gps <- gps_raw %>%
  dplyr::transmute(
    TagId = as.character(TagId),
    TimestampRaw = Timestamp,
    longitude = as.numeric(longitude),
    latitude = as.numeric(latitude)
  )


if (inherits(gps$TimestampRaw, "POSIXt")) {
  
  gps$Timestamp <- as.POSIXct(
    gps$TimestampRaw,
    tz = "UTC"
  )
  
} else if (is.numeric(gps$TimestampRaw)) {
  
  gps$Timestamp <- as.POSIXct(
    gps$TimestampRaw * 86400,
    origin = "1899-12-30",
    tz = "UTC"
  )
  
} else {
  
  gps$Timestamp <- lubridate::parse_date_time(
    as.character(gps$TimestampRaw),
    orders = c(
      "ymd HMS",
      "ymd HM",
      "dmy HMS",
      "dmy HM",
      "mdy HMS",
      "mdy HM"
    ),
    tz = "UTC"
  )
  
}

gps <- gps %>%
  dplyr::filter(
    !is.na(TagId),
    !is.na(Timestamp),
    !is.na(longitude),
    !is.na(latitude),
    dplyr::between(longitude, -180, 180),
    dplyr::between(latitude, -90, 90),
    
    # Exclude lion 8394
    TagId != "8394"
  ) %>%
  dplyr::arrange(TagId, Timestamp) %>%
  dplyr::distinct(
    TagId,
    Timestamp,
    .keep_all = TRUE
  )

cat("\nLions included:\n")
print(sort(unique(gps$TagId)))

stopifnot(!("8394" %in% gps$TagId))


gps <- gps %>%
  dplyr::mutate(
    LocalDate = as.Date(
      lubridate::with_tz(
        Timestamp,
        "Africa/Johannesburg"
      )
    ),
    
    Period = dplyr::case_when(
      LocalDate >= before_start &
        LocalDate <= before_end ~ "Before",
      
      LocalDate >= during_start &
        LocalDate <= during_end ~ "During",
      
      LocalDate >= after_start &
        LocalDate <= after_end ~ "After",
      
      TRUE ~ NA_character_
    ),
    
    Period = factor(
      Period,
      levels = period_levels
    )
  ) %>%
  dplyr::filter(!is.na(Period))

cat("\nGPS observations per lion and period:\n")

print(
  gps %>%
    dplyr::count(TagId, Period, .drop = FALSE),
  n = Inf
)

if (nrow(gps) == 0) {
  stop("No GPS observations fall within the flood periods.")
}


gps_sf <- sf::st_as_sf(
  gps,
  coords = c("longitude", "latitude"),
  crs = 4326
)

gps_sf <- sf::st_transform(
  gps_sf,
  32736
)

xy <- sf::st_coordinates(gps_sf)

gps$x <- xy[, 1]
gps$y <- xy[, 2]


raster_file <- "DynamicWorld_WoodyBinary_2025_Satara.tif"

if (!file.exists(raster_file)) {
  stop("Dynamic World woody raster not found.")
}

landcover <- terra::rast(raster_file)

if (terra::nlyr(landcover) != 1) {
  stop("Expected one woody vegetation raster layer.")
}

print(landcover)

raster_freq <- terra::freq(landcover)

print(raster_freq)

if (!all(raster_freq$value %in% c(0, 1))) {
  stop(
    "Woody vegetation raster must contain only 0 and 1."
  )
}

if (!terra::same.crs(landcover, target_crs)) {
  
  landcover <- terra::project(
    landcover,
    target_crs,
    method = "near"
  )
  
}


buffer_radius <- 50

res_x <- terra::res(landcover)[1]
res_y <- terra::res(landcover)[2]

nx <- ceiling(buffer_radius / res_x)
ny <- ceiling(buffer_radius / res_y)

dx <- seq(-nx, nx) * res_x
dy <- seq(-ny, ny) * res_y

kernel <- outer(
  dy,
  dx,
  function(y, x) {
    as.numeric(
      sqrt(x^2 + y^2) <= buffer_radius
    )
  }
)

kernel[kernel == 0] <- NA_real_

woody_cover_raster <- terra::focal(
  landcover,
  w = kernel,
  fun = "mean",
  na.rm = TRUE
) * 100

names(woody_cover_raster) <- "WoodyCover"

terra::writeRaster(
  woody_cover_raster,
  "WoodyCover_50m.tif",
  overwrite = TRUE
)


gps_split <- gps %>%
  dplyr::group_by(TagId, Period) %>%
  dplyr::group_split(.keep = TRUE)

names(gps_split) <- vapply(
  gps_split,
  function(x) {
    paste(
      unique(x$TagId),
      unique(as.character(x$Period)),
      sep = "_"
    )
  },
  character(1)
)


make_steps <- function(dat) {
  
  if (nrow(dat) < 5) {
    return(NULL)
  }
  
  track <- amt::make_track(
    dat,
    x,
    y,
    Timestamp,
    crs = 32736
  )
  
  # Split movement into approximately hourly bursts
  
  track <- amt::track_resample(
    track,
    rate = lubridate::hours(1),
    tolerance = lubridate::minutes(10)
  )
  
  steps <- amt::steps_by_burst(track)
  
  if (nrow(steps) < 5) {
    return(NULL)
  }
  
  
  steps <- steps %>%
    dplyr::filter(
      is.finite(sl_),
      sl_ > 0,
      is.finite(ta_)
    )
  
  if (nrow(steps) < 5) {
    return(NULL)
  }
  
  
  steps <- amt::random_steps(
    steps,
    n_control = n_available
  )
  
  steps$TagId <- unique(dat$TagId)
  
  steps$Period <- as.character(
    unique(dat$Period)
  )
  
  steps
}

steps_list <- lapply(
  gps_split,
  make_steps
)

steps_list <- Filter(
  Negate(is.null),
  steps_list
)

if (length(steps_list) == 0) {
  stop("No valid hourly movement steps were created.")
}

steps_all <- dplyr::bind_rows(
  steps_list,
  .id = "TrackID"
)

steps_all <- steps_all %>%
  dplyr::mutate(
    Period = factor(
      Period,
      levels = period_levels
    ),
    
    StepGroup = interaction(
      TrackID,
      step_id_,
      drop = TRUE
    )
  )

stopifnot(!("8394" %in% steps_all$TagId))


endpoints <- terra::vect(
  steps_all,
  geom = c("x2_", "y2_"),
  crs = target_crs
)

extracted <- terra::extract(
  woody_cover_raster,
  endpoints
)

steps_all$WoodyCover <- extracted$WoodyCover

steps_all <- steps_all %>%
  dplyr::mutate(
    WoodyProp = WoodyCover / 100,
    log_sl = log(sl_),
    cos_ta = cos(ta_)
  )


steps_all <- steps_all %>%
  dplyr::group_by(StepGroup) %>%
  dplyr::filter(
    dplyr::n() == n_available + 1,
    sum(case_ == TRUE) == 1,
    
    all(
      is.finite(WoodyProp),
      is.finite(sl_),
      is.finite(log_sl),
      is.finite(cos_ta)
    )
  ) %>%
  dplyr::ungroup()

if (nrow(steps_all) == 0) {
  stop(
    "No complete matched step sets remain. ",
    "Check raster coverage and generated endpoints."
  )
}

step_summary <- steps_all %>%
  dplyr::filter(case_ == TRUE) %>%
  dplyr::count(
    TagId,
    Period,
    name = "UsedSteps"
  )

cat("\nValid used steps per lion and period:\n")
print(step_summary, n = Inf)

write.csv(
  step_summary,
  "Lion_iSSF_step_summary.csv",
  row.names = FALSE
)


fit_lion <- function(dat) {
  
  survival::clogit(
    case_ ~
      0 +
      WoodyProp:Period +
      sl_ +
      log_sl +
      cos_ta +
      strata(StepGroup),
    
    data = dat,
    method = "exact"
  )
  
}

lion_data <- steps_all %>%
  dplyr::group_by(TagId) %>%
  dplyr::group_split()

names(lion_data) <- vapply(
  lion_data,
  function(x) unique(x$TagId),
  character(1)
)

# Only fit lions with sufficient matched steps
# and woody-cover variation.

eligible <- vapply(
  lion_data,
  function(x) {
    
    used <- x[x$case_ == TRUE, ]
    
    nrow(used) >= 10 &&
      length(unique(x$WoodyProp)) > 1
    
  },
  logical(1)
)

cat("\nLion model eligibility:\n")
print(eligible)

lion_data <- lion_data[eligible]

if (length(lion_data) == 0) {
  stop("No lions have sufficient data for model fitting.")
}

lion_models <- lapply(
  lion_data,
  fit_lion
)

cat("\nIndividual-lion iSSF models fitted.\n")


get_period_effects <- function(model, lion_id) {
  
  b <- stats::coef(model)
  V <- stats::vcov(model)
  
  available_terms <- grep(
    "^WoodyProp:Period",
    names(b),
    value = TRUE
  )
  
  if (length(available_terms) == 0) {
    stop("No woody selection coefficients for ", lion_id)
  }
  
  data.frame(
    Lion = lion_id,
    
    Period = factor(
      sub(
        "^WoodyProp:Period",
        "",
        available_terms
      ),
      levels = period_levels
    ),
    
    Estimate = unname(b[available_terms]),
    
    SE = sqrt(
      diag(V)[available_terms]
    )
  ) %>%
    dplyr::mutate(
      Lower95 = Estimate - 1.96 * SE,
      Upper95 = Estimate + 1.96 * SE
    )
}

lion_effects <- purrr::imap_dfr(
  lion_models,
  get_period_effects
) %>%
  dplyr::filter(
    !is.na(Period),
    is.finite(Estimate),
    is.finite(SE)
  )

cat("\nIndividual-lion selection coefficients:\n")
print(lion_effects)


population_summary <- lion_effects %>%
  dplyr::group_by(Period) %>%
  dplyr::summarise(
    N = dplyr::n(),
    
    MeanEstimate = mean(Estimate),
    
    SD = if (N > 1) {
      sd(Estimate)
    } else {
      NA_real_
    },
    
    SE = if (N > 1) {
      SD / sqrt(N)
    } else {
      NA_real_
    },
    
    Lower95 = if (N > 1) {
      MeanEstimate -
        qt(0.975, df = N - 1) * SE
    } else {
      NA_real_
    },
    
    Upper95 = if (N > 1) {
      MeanEstimate +
        qt(0.975, df = N - 1) * SE
    } else {
      NA_real_
    },
    
    .groups = "drop"
  )

cat("\nPopulation-level selection:\n")
print(population_summary)


effects_wide <- lion_effects %>%
  dplyr::select(
    Lion,
    Period,
    Estimate
  ) %>%
  tidyr::pivot_wider(
    names_from = Period,
    values_from = Estimate
  )

paired_contrast <- function(data, p1, p2) {
  
  valid <- data %>%
    dplyr::filter(
      !is.na(.data[[p1]]),
      !is.na(.data[[p2]])
    )
  
  differences <- valid[[p1]] - valid[[p2]]
  
  if (length(differences) < 2) {
    
    return(
      data.frame(
        Comparison = paste(p1, "-", p2),
        N = length(differences),
        MeanDifference = if (length(differences) == 1) {
          differences
        } else {
          NA_real_
        },
        Lower95 = NA_real_,
        Upper95 = NA_real_,
        P = NA_real_
      )
    )
    
  }
  
  test <- stats::t.test(
    differences,
    mu = 0
  )
  
  data.frame(
    Comparison = paste(p1, "-", p2),
    N = length(differences),
    MeanDifference = mean(differences),
    Lower95 = test$conf.int[1],
    Upper95 = test$conf.int[2],
    P = test$p.value
  )
}

period_contrasts <- dplyr::bind_rows(
  
  paired_contrast(
    effects_wide,
    "During",
    "Before"
  ),
  
  paired_contrast(
    effects_wide,
    "After",
    "Before"
  ),
  
  paired_contrast(
    effects_wide,
    "After",
    "During"
  )
  
)

period_contrasts$AdjustedP <- p.adjust(
  period_contrasts$P,
  method = "holm"
)

cat("\nPeriod comparisons:\n")
print(period_contrasts)

# Individual lion selection

selection_plot <- ggplot(
  lion_effects,
  aes(
    x = Period,
    y = Estimate,
    group = Lion,
    colour = Lion
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey40"
  ) +
  geom_line(
    linewidth = 0.7
  ) +
  geom_point(
    size = 2.5
  ) +
  theme_classic(
    base_size = 12
  ) +
  labs(
    x = "Flood period",
    y = "Selection coefficient for woody vegetation",
    colour = "Lion"
  )

print(selection_plot)


# Population selection

population_plot <- ggplot(
  population_summary,
  aes(
    x = Period,
    y = MeanEstimate
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey40"
  ) +
  geom_errorbar(
    aes(
      ymin = Lower95,
      ymax = Upper95
    ),
    width = 0.15
  ) +
  geom_point(
    size = 3
  ) +
  theme_classic(
    base_size = 12
  ) +
  labs(
    x = "Flood period",
    y = "Mean woody vegetation selection coefficient"
  )

print(population_plot)


saveRDS(
  lion_models,
  "Lion_iSSF_models.rds"
)

saveRDS(
  steps_all,
  "Lion_iSSF_steps.rds"
)

write.csv(
  lion_effects,
  "Lion_iSSF_individual_effects.csv",
  row.names = FALSE
)

write.csv(
  population_summary,
  "Lion_iSSF_population_summary.csv",
  row.names = FALSE
)

write.csv(
  period_contrasts,
  "Lion_iSSF_period_contrasts.csv",
  row.names = FALSE
)

}


# Woody vegetation selection plot
# Each point is one lion's period-specific selection coefficient.
if (RUN_WOODY_SELECTION) {
  woody_boxplot <- ggplot2::ggplot(
    lion_effects, ggplot2::aes(x = Period, y = Estimate, fill = Period)
  ) +
    ggplot2::geom_boxplot(width = 0.55, alpha = 0.6,
                          outlier.shape = NA, colour = "black") +
    ggplot2::geom_jitter(width = 0.10, size = 2.5,
                         alpha = 0.85, colour = "black") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    ggplot2::scale_fill_manual(values = c(Before = "grey", During = "grey", After = "grey")) +
    ggplot2::labs(x = "Flood period", y = "Woody vegetation selection coefficient") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      legend.position = "none",
      axis.title = ggplot2::element_text(size = 16, colour = "black"),
      axis.text = ggplot2::element_text(size = 14, colour = "black")
    )
  print(woody_boxplot)
  ggplot2::ggsave("Lion_Woody_Selection_Boxplot.svg", woody_boxplot,
                  width = 9, height = 7)
}

# Daily rainfall
if (RUN_RAINFALL_FIGURE) {
library(dplyr)
library(ggplot2)
library(readr)

weather <- read_csv("Finalweatherstation.csv")

start_date <- as.Date("2025-12-31")
end_date   <- as.Date("2026-01-29")

rainfall_daily <- weather %>%
  mutate(
    Date = as.Date(Date, format = "%Y/%m/%d"),
    rain = as.numeric(rain)
  ) %>%
  filter(
    Date >= start_date,
    Date <= end_date
  ) %>%
  group_by(Date) %>%
  summarise(
    Rainfall_mm = sum(rain, na.rm = TRUE),
    .groups = "drop"
  )

print(rainfall_daily, n = 30)

max(rainfall_daily$Rainfall_mm, na.rm = TRUE)


rainfall_plot <- ggplot(
  rainfall_daily,
  aes(x = Date, y = Rainfall_mm)
) +
  
  geom_col(
    fill = "black",
    width = 0.9
  ) +
  
  # Lines separating flood periods
  geom_vline(
    xintercept = as.Date(
      c("2026-01-10", "2026-01-20")
    ) - 0.5,
    linetype = "dashed",
    colour = "grey40",
    linewidth = 0.6
  ) +
  
  annotate(
    "text",
    x = as.Date(
      c("2026-01-05", "2026-01-15", "2026-01-25")
    ),
    y = max(rainfall_daily$Rainfall_mm, na.rm = TRUE) * 1.15,
    label = c("Before", "During", "After"),
    size = 4
  ) +
  
  scale_x_date(
    breaks = seq(
      start_date,
      end_date,
      by = "1 day"
    ),
    date_labels = "%d %b",
    expand = expansion(add = 0.5)
  ) +
  
  scale_y_continuous(
    breaks = scales::breaks_width(25),
    expand = expansion(mult = c(0, 0.12))
  ) +
  
  labs(
    x = "Date",
    y = "Total daily rainfall (mm)"
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    axis.title = element_text(size = 14)
  )

print(rainfall_plot)

ggsave(
  filename = "Daily_Rainfall_Flood_Periods.svg",
  plot = rainfall_plot,
  width = 11,
  height = 7,
  units = "in",
  device = "svg"
)

}

# Vegetation probability maps
if (RUN_VEGETATION_MAPS) {

packages <- c(
  "readr",
  "dplyr",
  "sf",
  "terra",
  "ggplot2",
  "patchwork"
)

missing_packages <- setdiff(
  packages,
  rownames(installed.packages())
)

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

invisible(lapply(packages, library, character.only = TRUE))

gps_data <- readr::read_csv(
  "Combined collar data with weather astronomy and movement(in).csv",
  show_col_types = FALSE
)

gps <- gps_data %>%
  mutate(
    longitude = as.numeric(longitude),
    latitude = as.numeric(latitude)
  ) %>%
  filter(
    is.finite(longitude),
    is.finite(latitude),
    between(longitude, -180, 180),
    between(latitude, -90, 90)
  )

if (nrow(gps) == 0) {
  stop("No valid GPS coordinates found.")
}

cat("Total GPS locations:", nrow(gps), "\n")


woody <- terra::rast(
  "DynamicWorld_WoodyProbability_2025_Satara.tif"
)

tree <- terra::rast(
  "DynamicWorld_TreeProbability_2025_Satara.tif"
)

shrub <- terra::rast(
  "DynamicWorld_ShrubProbability_2025_Satara.tif"
)

# Confirm rasters have matching geometry
if (!terra::compareGeom(
  woody,
  tree,
  shrub,
  stopOnError = FALSE
)) {
  stop("Raster geometries do not match.")
}


satara <- data.frame(
  Name = "Satara Rest Camp",
  longitude = 31.78,
  latitude = -24.39
)

satara_sf <- sf::st_as_sf(
  satara,
  coords = c("longitude", "latitude"),
  crs = 4326
)

# Transform to raster CRS
satara_sf <- sf::st_transform(
  satara_sf,
  crs = terra::crs(woody)
)

satara_xy <- sf::st_coordinates(satara_sf)

satara_plot <- data.frame(
  x = satara_xy[, 1],
  y = satara_xy[, 2],
  Name = "Satara Rest Camp"
)

gps_sf <- sf::st_as_sf(
  gps,
  coords = c("longitude", "latitude"),
  crs = 4326
)

gps_projected <- sf::st_transform(
  gps_sf,
  crs = terra::crs(woody)
)

gps_bbox <- sf::st_bbox(gps_projected)

# Include both the GPS extent and Satara Rest Camp
all_x <- c(
  as.numeric(gps_bbox[c("xmin", "xmax")]),
  satara_plot$x
)

all_y <- c(
  as.numeric(gps_bbox[c("ymin", "ymax")]),
  satara_plot$y
)

# Buffer in metres
buffer_distance <- 500

crop_extent <- terra::ext(
  min(all_x) - buffer_distance,
  max(all_x) + buffer_distance,
  min(all_y) - buffer_distance,
  max(all_y) + buffer_distance
)


woody_crop <- terra::crop(
  woody,
  crop_extent
)

tree_crop <- terra::crop(
  tree,
  crop_extent
)

shrub_crop <- terra::crop(
  shrub,
  crop_extent
)

cat(
  "Cropped raster:",
  terra::ncol(woody_crop), "columns x",
  terra::nrow(woody_crop), "rows\n"
)

to_percent <- function(r) {
  
  maximum <- terra::global(
    r,
    "max",
    na.rm = TRUE
  )[1, 1]
  
  if (!is.finite(maximum)) {
    stop("Raster has no valid probability values.")
  }
  
  if (maximum <= 1) {
    r <- r * 100
  }
  
  return(r)
}

woody_crop <- to_percent(woody_crop)
tree_crop <- to_percent(tree_crop)
shrub_crop <- to_percent(shrub_crop)


# Combined tree and shrub probability

combined_crop <- terra::clamp(
  tree_crop + shrub_crop,
  lower = 0,
  upper = 100
)

names(woody_crop) <- "Woody"
names(tree_crop) <- "Tree"
names(shrub_crop) <- "Shrub"
names(combined_crop) <- "Combined"


# Save cropped rasters

terra::writeRaster(
  woody_crop,
  "Woody_Probability_Cropped.tif",
  overwrite = TRUE
)

terra::writeRaster(
  tree_crop,
  "Tree_Probability_Cropped.tif",
  overwrite = TRUE
)

terra::writeRaster(
  shrub_crop,
  "Shrub_Probability_Cropped.tif",
  overwrite = TRUE
)

terra::writeRaster(
  combined_crop,
  "TreeShrub_Combined_Cropped.tif",
  overwrite = TRUE
)

prepare_plot_raster <- function(
    r,
    max_cells = 500000,
    fill_iterations = 10
) {
  
  # Treat non-finite values as missing
  r <- terra::ifel(is.finite(r), r, NA)
  
  # Reduce raster size for plotting
  aggregation_factor <- max(
    1,
    ceiling(
      sqrt(terra::ncell(r) / max_cells)
    )
  )
  
  if (aggregation_factor > 1) {
    r <- terra::aggregate(
      r,
      fact = aggregation_factor,
      fun = "mean",
      na.rm = TRUE
    )
  }
  
  # Fill small missing gaps using neighbouring cells
  for (i in seq_len(fill_iterations)) {
    
    neighbourhood <- terra::focal(
      r,
      w = matrix(1, 3, 3),
      fun = "mean",
      na.rm = TRUE,
      na.policy = "only"
    )
    
    r <- terra::cover(r, neighbourhood)
  }
  
  return(r)
}

woody_plot_r <- prepare_plot_raster(woody_crop)
tree_plot_r <- prepare_plot_raster(tree_crop)
shrub_plot_r <- prepare_plot_raster(shrub_crop)
combined_plot_r <- prepare_plot_raster(combined_crop)

raster_to_df <- function(r) {
  
  df <- as.data.frame(
    r,
    xy = TRUE,
    na.rm = FALSE
  )
  
  names(df)[3] <- "Probability"
  
  return(df)
}

woody_df <- raster_to_df(woody_plot_r)
tree_df <- raster_to_df(tree_plot_r)
shrub_df <- raster_to_df(shrub_plot_r)
combined_df <- raster_to_df(combined_plot_r)

map_extent <- terra::ext(woody_plot_r)

make_heatmap <- function(data) {
  
  ggplot(
    data,
    aes(
      x = x,
      y = y,
      fill = Probability
    )
  ) +
    
    geom_raster(
      interpolate = FALSE
    ) +
    
    scale_fill_gradientn(
      colours = c(
        "#F7FCF5",
        "#C7E9C0",
        "#74C476",
        "#238B45",
        "#00441B"
      ),
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100),
      name = "Probability (%)",
      na.value = "white"
    ) +
    
    geom_point(
      data = satara_plot,
      aes(x = x, y = y),
      inherit.aes = FALSE,
      shape = 21,
      fill = "white",
      colour = "black",
      size = 3,
      stroke = 0.8
    ) +
    
    geom_text(
      data = satara_plot,
      aes(
        x = x,
        y = y,
        label = Name
      ),
      inherit.aes = FALSE,
      hjust = -0.05,
      vjust = -0.7,
      colour = "black",
      size = 4,
      fontface = "bold"
    ) +
    
    coord_fixed(
      ratio = 1,
      xlim = c(map_extent$xmin, map_extent$xmax),
      ylim = c(map_extent$ymin, map_extent$ymax),
      expand = FALSE
    ) +
    
    labs(
      x = "Easting (m)",
      y = "Northing (m)"
    ) +
    
    theme_classic(
      base_size = 12
    ) +
    
    theme(
      plot.title = element_blank(),
      panel.grid = element_blank(),
      axis.text = element_text(colour = "black"),
      legend.position = "right"
    )
}


woody_heatmap <- make_heatmap(woody_df)

tree_heatmap <- make_heatmap(tree_df)

shrub_heatmap <- make_heatmap(shrub_df)

combined_heatmap <- make_heatmap(combined_df)

print(woody_heatmap)
print(tree_heatmap)
print(shrub_heatmap)
print(combined_heatmap)

vegetation_combined <- (
  woody_heatmap | tree_heatmap
) / (
  shrub_heatmap | combined_heatmap
) &
  theme(
    plot.title = element_blank(),
    legend.position = "bottom"
  )

print(vegetation_combined)

ggsave(
  "Woody_Probability_Heatmap.png",
  woody_heatmap,
  width = 9,
  height = 8,
  dpi = 300
)

ggsave(
  "Tree_Probability_Heatmap.png",
  tree_heatmap,
  width = 9,
  height = 8,
  dpi = 300
)

ggsave(
  "Shrub_Probability_Heatmap.png",
  shrub_heatmap,
  width = 9,
  height = 8,
  dpi = 300
)

ggsave(
  "Combined_TreeShrub_Heatmap.png",
  combined_heatmap,
  width = 9,
  height = 8,
  dpi = 300
)

ggsave(
  "Vegetation_Probability_Combined.png",
  vegetation_combined,
  width = 14,
  height = 12,
  dpi = 300
)

ggsave(
  "Vegetation_Probability_Combined.svg",
  vegetation_combined,
  width = 14,
  height = 12
)

cat("\nAll four vegetation heat maps saved successfully.\n")
}

# Minimum convex polygon area
calculate_mcp <- function(df) {
  if (nrow(df) < 4) return(NA_real_)
  centre_x <- mean(df$X)
  centre_y <- mean(df$Y)
  df <- df %>% mutate(dist_centre = sqrt((X - centre_x)^2 +
                                           (Y - centre_y)^2)) %>%
    arrange(dist_centre) %>% slice_head(n = floor(nrow(df) * 0.95))
  pts <- st_as_sf(df, coords = c('X', 'Y'), crs = 32736)
  as.numeric(st_area(st_convex_hull(st_union(pts)))) / 1e6
}
mcp_results <- gps_xy %>% group_by(TagId, Period) %>%
  group_modify(~ tibble(Area_km2 = calculate_mcp(.x))) %>%
  ungroup() %>% filter(is.finite(Area_km2), Area_km2 > 0) %>%
  mutate(LogArea = log(Area_km2))
mcp_model <- lmer(LogArea ~ Period + (1 | TagId), data = mcp_results)
summary(mcp_model)
anova(mcp_model)

mcp_results %>%
  dplyr::select(TagId, Period, Area_km2) %>%
  arrange(TagId, Period) %>%
  print(n = Inf)

library(dplyr)

mcp_means <- mcp_results %>%
  group_by(Period) %>%
  summarise(
    Mean_km2 = mean(Area_km2, na.rm = TRUE),
    SD_km2 = sd(Area_km2, na.rm = TRUE),
    SE_km2 = SD_km2 / sqrt(n()),
    .groups = "drop"
  )

print(mcp_means)
