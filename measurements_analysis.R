## Load Required Libraries
library(readxl)      # Read Excel files
library(dplyr)       # Data manipulation
library(ggplot2)     # Data visualization
library(mgcv)        # Generalized additive models
library(reshape2)    # Data reshaping (for correlation heatmaps)
library(FSA)         # Dunn's test for multiple comparisons
library(lme4)        # Linear mixed-effects models
library(lmerTest)    # p-values for lme4
library(performance) # Model performance metrics (ICC, R2)
library(MuMIn)       # R^2 for mixed models
library(ggeffects)   # Predictions from mixed models

#### Load Data ####
# Load measurements data
measurements <- read_excel("Database.xlsx", sheet = "Measurements") %>%
    filter(!is.na(EW))  # Remove missing values in EW column

# Load temperature and oxygen model reconstructions
Temp <- read_excel("temp_o2.xlsx", sheet = "Temp")
O2 <- read_excel("temp_o2.xlsx", sheet = "O2")

#### Data Preparation and Log-Transformation ####
## Log-transform external width
measurements <- measurements %>%
    mutate(
        logEW = log(EW),                        # Natural log transformation
        Abs_Paleolat = abs(Paleolat),           # Distance from equator
        Setting = factor(Setting),          # Convert to factor for modelling
        Zone = factor(Zone),
        Map = factor(Map)
    )

# Filter data by time slice for temporal analyses
measurements_280 <- measurements %>% filter(Time == 280)  # Late Cisuralian
measurements_300 <- measurements %>% filter(Time == 300)  # LCa-EPe
measurements_320 <- measurements %>% filter(Time == 320)  # MCa
measurements_340 <- measurements %>% filter(Time == 340)  # ECa
measurements_360 <- measurements %>% filter(Time == 360)  # EMi

#### Summary Statistics Functions ####
## Calculate summary statistics by grouping variables
#' @param data Data frame containing EW measurements
#' @param time_value Numeric age for the time slice
#' @return Summary table with mean, SD, N, min, max, median per group
calculate_summary <- function(data, time_value) {
    data %>%
        group_by(Unit, Setting, Category, Basin, Paleolat, Zone) %>%
        summarise(
            Time = first(time_value),
            Mean = mean(EW, na.rm = TRUE),
            SD = sd(EW, na.rm = TRUE),
            N = n(),
            Max = max(EW, na.rm = TRUE),
            Min = min(EW, na.rm = TRUE),
            Median = median(EW, na.rm = TRUE),
            .groups = "drop"
        )
}

## Generic summary statistics function for any grouping variable
#' @param data Data frame containing EW measurements
#' @param group_var Character vector of column names to group by
#' @return Summary table with mean, SD, N, min, max, median per group
summary_stats <- function(data, group_var) {
    data %>%
        group_by(across(all_of(group_var))) %>%
        summarise(
            mean = mean(EW, na.rm = TRUE),
            sd = sd(EW, na.rm = TRUE),
            n = n(),
            max = max(EW, na.rm = TRUE),
            min = min(EW, na.rm = TRUE),
            median = median(EW, na.rm = TRUE),
            .groups = "drop"
        )
}

# Generate summary tables for different groupings
summary_ichnogenus <- summary_stats(measurements, "Ichnogenus")
summary_map <- summary_stats(measurements, "Map")
summary_zone <- summary_stats(measurements, "Zone")
summary_setting <- summary_stats(measurements, "Setting")
summary_category <- summary_stats(measurements, "Category")
summary_category_setting <- summary_stats(measurements, c("Setting", "Category"))
summary_setting_zone <- summary_stats(measurements, c("Zone", "Setting"))

# Display summaries
print(summary_zone)
print(summary_setting)
print(summary_category)

#### Descriptive Visualizations: Boxplots by Category ####
# External width by climate zone
ggplot(measurements, aes(x = Zone, y = EW, fill = Zone)) +
    geom_boxplot() +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +
    labs(y = "External Width (mm)") +
    scale_y_log10() +  # Log scale for better visualization of range
    theme_bw()

# External width by depositional setting
ggplot(measurements, aes(x = Setting, y = EW, fill = Setting)) +
    geom_boxplot() +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 5) +
    labs(y = "External width (mm)") +
    scale_y_log10() +
    coord_cartesian(ylim = c(1, 800)) +
    theme_bw()

# External width by category
ggplot(measurements, aes(x = Category, y = EW, fill = Category)) +
    geom_boxplot() +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 5) +
    labs(y = "External width (mm)") +
    scale_y_log10() +
    coord_cartesian(ylim = c(1, 800)) +
    theme_bw()

#### Statistical Tests: Kruskal-Wallis and Dunn's Post-Hoc ####
# Climate zones
kruskal.test(EW ~ Zone, data = measurements)
dunnTest(EW ~ Zone, data = measurements, method = "bonferroni")

# Depositional settings
kruskal.test(EW ~ Setting, data = measurements)
dunnTest(EW ~ Setting, data = measurements, method = "bonferroni")

# Categories within settings
kruskal.test(EW ~ Category, data = measurements)
dunnTest(EW ~ Category, data = measurements, method = "bonferroni")

#### Latitudinal Patterns: Size vs. Absolute Paleolatitude ####
# Global dataset: size vs. absolute paleolatitude
ggplot(measurements, aes(x = Abs_Paleolat, y = logEW, color = Setting)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", color = "red") +
    geom_smooth(method = "gam", color = "black") +
    labs(x = "Absolute Paleolatitude", y = "log(External Width mm)") +
    scale_x_reverse(breaks = c(0, 20, 40, 60, 80, 90)) +
    theme_bw()

## Correlation tests
# Spearman rank correlation
cor.test(measurements$Abs_Paleolat, measurements$logEW,
         use = "complete.obs", method = "spearman")

# Linear model (LM)
lm_global <- lm(logEW ~ Abs_Paleolat, data = measurements)
summary(lm_global)

# Generalized Additive Model (GAM) for non-linear relationships
gam_global <- gam(logEW ~ s(Abs_Paleolat), data = measurements)
summary(gam_global)

#### Temporal Patterns: Size vs. Absolute Paleolatitude by Time Slice ####
# Middle Carboniferous (320 Ma)
ggplot(measurements_320, aes(y = logEW, x = Abs_Paleolat, color = Setting)) +
    geom_point() +
    geom_smooth(method = "lm", color = "red") +
    geom_smooth(method = "gam", color = "black") +
    labs(x = "Absolute Paleolatitude", y = "log(External Width mm)") +
    scale_x_reverse(breaks = c(0, 20, 40, 60, 80, 90)) +
    theme_bw()

lm_320 <- lm(logEW ~ Abs_Paleolat, data = measurements_320)
summary(lm_320)

gam_320 <- gam(logEW ~ s(Abs_Paleolat), data = measurements_320)
summary(gam_320)

cor.test(measurements_320$Abs_Paleolat, measurements_320$logEW,
         use = "complete.obs", method = "spearman")

# Late Carboniferous–Early Permian (300 Ma)
ggplot(measurements_300, aes(y = logEW, x = Abs_Paleolat, color = Setting)) +
    geom_point() +
    geom_smooth(method = "lm", color = "red") +
    geom_smooth(method = "gam", color = "black") +
    labs(x = "Absolute Paleolatitude", y = "log(External Width mm)") +
    scale_x_reverse(breaks = c(0, 20, 40, 60, 80, 90)) +
    theme_bw()

lm_300 <- lm(logEW ~ Abs_Paleolat, data = measurements_300)
summary(lm_300)

gam_300 <- gam(logEW ~ s(Abs_Paleolat), data = measurements_300)
summary(gam_300)

cor.test(measurements_300$Abs_Paleolat, measurements_300$logEW,
         use = "complete.obs", method = "spearman")

# Late Cisuralian (280 Ma)
ggplot(measurements_280, aes(y = logEW, x = Abs_Paleolat, color = Setting)) +
    geom_point() +
    geom_smooth(method = "lm", color = "red") +
    geom_smooth(method = "gam", color = "black") +
    labs(x = "Absolute Paleolatitude", y = "log(External Width mm)") +
    scale_x_reverse(breaks = c(0, 20, 40, 60, 80, 90)) +
    theme_bw()

lm_280 <- lm(logEW ~ Abs_Paleolat, data = measurements_280)
summary(lm_280)

gam_280 <- gam(logEW ~ s(Abs_Paleolat), data = measurements_280)
summary(gam_280)

cor.test(measurements_280$Abs_Paleolat, measurements_280$logEW,
         use = "complete.obs", method = "spearman")

#### Temporal Trends: Size Through Geological Time ####
# External width over time
ggplot(measurements, aes(x = Mean_time, y = EW, color = Zone)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "gam", color = "black") +
    geom_smooth(method = "lm", color = "red") +
    scale_x_reverse(breaks = c(360, 340, 320, 300, 280)) +
    coord_cartesian(ylim = c(1, 1000), xlim = c(354, 277)) +
    labs(x = "Age (Ma)", y = "External Width (mm)") +
    theme_bw()

# log-transformed version
ggplot(measurements, aes(x = Mean_time, y = logEW, color = Zone)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "gam", color = "black") +
    geom_smooth(method = "lm", color = "red") +
    scale_x_reverse(breaks = c(360, 340, 320, 300, 280)) +
    labs(x = "Age (Ma)", y = "log(External Width mm)") +
    theme_bw()

# Correlation and models
cor.test(measurements$Mean_time, measurements$logEW,
         use = "complete.obs", method = "spearman")

lm_time <- lm(logEW ~ Mean_time, data = measurements)
summary(lm_time)

gam_time <- gam(logEW ~ s(Mean_time), data = measurements)
summary(gam_time)

#### Temperature and Oxygen Reconstructions Over Time ####
# Plot temperature and oxygen models against geological time
ggplot() +
    # Temperature models (red shades)
    geom_line(aes(x = Age_1, y = Scotese), color = "red", data = Temp) +
    geom_line(aes(x = Age_2, y = Song), color = "brown1", data = Temp) +
    geom_line(aes(x = Age_3, y = Geocarbsulf), color = "brown", data = Temp) +
    # Oxygen models (blue shades)
    geom_line(aes(x = Age_2, y = Royer), color = "blue", data = O2) +
    geom_line(aes(x = Age_3, y = Mils), color = "cornflowerblue", data = O2) +
    geom_smooth(aes(x = Age_1, y = Song), color = "darkblue", 
                method = "gam", se = FALSE, data = O2) +
    scale_x_reverse(breaks = c(360, 340, 320, 300, 280)) +
    coord_cartesian(ylim = c(10, 40), xlim = c(354, 277)) +
    labs(x = "Age (Ma)", y = "°C and pO₂ (%)") +
    theme_bw()

#### Correlation Between Trace Fossil Size and Environmental Variables ####
# This section interpolates all data to a common time scale and calculates
# correlations between the GAM-smoothed size curve and temperature/oxygen models.

# Extract GAM-smoothed size curve
p_size <- ggplot(measurements, aes(x = Mean_time, y = logEW)) +
    scale_x_reverse() +
    geom_smooth(method = "gam", color = "black") +
    labs(title = "Smoothed size through time",
         x = "Age (Ma)", y = "log(External Width)") +
    theme_bw()

# Extract the smoothed values from the ggplot object
size_smooth <- ggplot_build(p_size)$data[[1]]
size_smooth$x <- abs(size_smooth$x)
names(size_smooth)[names(size_smooth) == "x"] <- "Mean_time"
names(size_smooth)[names(size_smooth) == "y"] <- "EW_smoothed"

# Fit GAM to Song's oxygen data for interpolation
gam_song <- gam(Song ~ s(Age_1), data = O2)
age_values <- seq(354, 277, length.out = 100)
song_predictions <- predict(gam_song, newdata = data.frame(Age_1 = age_values))
song_smooth <- data.frame(Age_1 = age_values, Song_smoothed = song_predictions)

# Interpolate all variables to common time scale
common_age <- seq(277, 354, by = 1)

# Interpolate size curve
size_interp <- approx(size_smooth$Mean_time, size_smooth$EW_smoothed, 
                      xout = common_age)

# Interpolate temperature models
temp_Scotese_interp <- approx(Temp$Age_1, Temp$Scotese, xout = common_age)
temp_Song_interp <- approx(Temp$Age_2, Temp$Song, xout = common_age)
temp_Geoc_interp <- approx(Temp$Age_3, Temp$Geocarbsulf, xout = common_age)

# Interpolate oxygen models
o2_Royer_interp <- approx(O2$Age_2, O2$Royer, xout = common_age)
o2_Mils_interp <- approx(O2$Age_3, O2$Mils, xout = common_age)
o2_Song_interp <- approx(song_smooth$Age_1, song_smooth$Song_smoothed, 
                         xout = common_age)

# Combine all interpolated data
combined_data <- data.frame(
    Age = common_age,
    Size_smoothed = size_interp$y,
    Temp_Scotese = temp_Scotese_interp$y,
    Temp_Song = temp_Song_interp$y,
    Temp_Geoc = temp_Geoc_interp$y,
    O2_Mils = o2_Mils_interp$y,
    O2_Song = o2_Song_interp$y,
    O2_Royer = o2_Royer_interp$y
)

head(combined_data)

# Correlation matrix (Spearman)
cor_matrix <- cor(combined_data[, -1], use = "complete.obs", method = "spearman")
print(cor_matrix)

# Correlation heatmap
cor_melt <- reshape2::melt(cor_matrix)

ggplot(cor_melt, aes(Var1, Var2, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", midpoint = 0, 
                         limit = c(-1, 1)) +
    labs(title = "Correlation Heatmap: Size vs. Environmental Variables") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Individual correlation tests
# Temperature models
cor.test(combined_data$Size_smoothed, combined_data$Temp_Scotese,
         method = "spearman")
cor.test(combined_data$Size_smoothed, combined_data$Temp_Song,
         method = "spearman")
cor.test(combined_data$Size_smoothed, combined_data$Temp_Geoc,
         method = "spearman")

# Oxygen models
cor.test(combined_data$Size_smoothed, combined_data$O2_Mils,
         method = "spearman")
cor.test(combined_data$Size_smoothed, combined_data$O2_Song,
         method = "spearman")
cor.test(combined_data$Size_smoothed, combined_data$O2_Royer,
         method = "spearman")

# Scatter plots for visual inspection
ggplot(combined_data, aes(x = Temp_Scotese, y = Size_smoothed)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", color = "blue") +
    labs(x = "Temperature (°C)", y = "Smoothed log(EW)") +
    theme_minimal()

ggplot(combined_data, aes(x = O2_Mils, y = Size_smoothed)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", color = "blue") +
    labs(x = "pO₂ (%)", y = "Smoothed log(EW)") +
    theme_minimal()

#### Environmental Partitioning: Subaqueous vs. Non-Subaqueous
# Hypothesis h3 predicts that terrestrial and aquatic arthropods respond
# differently to temperature and oxygen. We test this by partitioning the data.

# Coastal subaqueous traces only
subaqueous <- measurements %>% 
    filter(Setting == "Coastal" & Category == "Subaqueous")

# Repeat the GAM smoothing and correlation for this subset
p_sub <- ggplot(subaqueous, aes(x = Mean_time, y = logEW)) +
    scale_x_reverse() +
    geom_smooth(method = "loess", color = "black", span = 0.8) +
    labs(x = "Age (Ma)", y = "log(External Width)") +
    theme_bw()

size_sub_smooth <- ggplot_build(p_sub)$data[[1]]
size_sub_smooth$x <- abs(size_sub_smooth$x)
names(size_sub_smooth)[names(size_sub_smooth) == "x"] <- "Mean_time"
names(size_sub_smooth)[names(size_sub_smooth) == "y"] <- "EW_smoothed"

# Interpolate and correlate for subaqueous data
size_sub_interp <- approx(size_sub_smooth$Mean_time, size_sub_smooth$EW_smoothed, 
                          xout = common_age)

combined_sub <- data.frame(
    Age = common_age,
    Size_smoothed = size_sub_interp$y,
    Temp_Song = temp_Song_interp$y,
    O2_Song = o2_Song_interp$y
)

# Correlations for subaqueous data
cor.test(combined_sub$Size_smoothed, combined_sub$Temp_Song, method = "spearman")
cor.test(combined_sub$Size_smoothed, combined_sub$O2_Song, method = "spearman")

# Non-subaqueous traces (terrestrial and transitional)
subaerial <- measurements %>% filter(!(Category == "Subaqueous"))

# Repeat the GAM smoothing and correlation for this subset
p_sub_terr <- ggplot(subaerial, aes(x = Mean_time, y = logEW)) +
    scale_x_reverse() +
    geom_smooth(method = "loess", color = "black", span = 1) +
    labs(x = "Age (Ma)", y = "log(External Width)") +
    theme_bw()

size_terr_smooth <- ggplot_build(p_sub_terr)$data[[1]]
size_terr_smooth$x <- abs(size_terr_smooth$x)
names(size_terr_smooth)[names(size_terr_smooth) == "x"] <- "Mean_time"
names(size_terr_smooth)[names(size_terr_smooth) == "y"] <- "EW_smoothed"

# Interpolate and correlate for non-subaqueous data
size_terr_interp <- approx(size_terr_smooth$Mean_time, size_terr_smooth$EW_smoothed,
                           xout = common_age)

combined_terr <- data.frame(
    Age = common_age,
    Size_smoothed = size_terr_interp$y,
    Temp_Scotese = temp_Scotese_interp$y,
    Temp_Geoc = temp_Geoc_interp$y,
    O2_Mils = o2_Mils_interp$y,
    O2_Royer = o2_Royer_interp$y
)

# Correlations for non-subaqueous data
cor.test(combined_terr$Size_smoothed, combined_terr$Temp_Scotese, method = "spearman")
cor.test(combined_terr$Size_smoothed, combined_terr$Temp_Geoc, method = "pearson")
cor.test(combined_terr$Size_smoothed, combined_terr$O2_Mils, method = "spearman")
cor.test(combined_terr$Size_smoothed, combined_terr$O2_Royer, method = "spearman")

#### Linear Mixed-Effects Models (LMMs) ####
# Full model with all random effects
lmm_full <- lmer(
    logEW ~ Abs_Paleolat +
        (1 | Setting) +
        (1 | Zone) +
        (1 | Map),
    data = measurements
)

summary(lmm_full)

# Variance components
VarCorr(lmm_full)
performance::icc(lmm_full)  # Intraclass correlation coefficient

## Marginal and conditional R²
MuMIn::r.squaredGLMM(lmm_full)

# Model comparison via AIC
m0 <- lm(logEW ~ Abs_Paleolat, data = measurements)                    # Null model
m1 <- lmer(logEW ~ Abs_Paleolat + (1 | Setting), data = measurements)  # Setting only
m2 <- lmer(logEW ~ Abs_Paleolat + (1 | Zone), data = measurements)     # Zone only
m3 <- lmer(logEW ~ Abs_Paleolat + (1 | Map), data = measurements)      # Time only
m4 <- lmer(logEW ~ Abs_Paleolat + (1 | Setting) + (1 | Zone) +         # Full model
               (1 | Map), data = measurements)

AIC(m0, m1, m2, m3, m4)

# Predictions from the full model
pred <- ggpredict(lmm_full, terms = "Abs_Paleolat")

ggplot() +
    geom_point(data = measurements, aes(Abs_Paleolat, logEW), alpha = 0.5) +
    geom_line(data = pred, aes(x, predicted), colour = "red", linewidth = 1.2) +
    geom_ribbon(data = pred, aes(x = x, ymin = conf.low, ymax = conf.high),
                alpha = 0.25) +
    scale_x_reverse() +
    theme_bw() +
    labs(x = "Absolute paleolatitude (°)", y = "log(External Width mm)")