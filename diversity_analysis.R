## Load Required Libraries
library(readxl)      # Read Excel files
library(dplyr)       # Data manipulation
library(tidyr)       # Data reshaping
library(purrr)       # Functional programming
library(ggplot2)     # Data visualization
library(vegan)       # Community ecology (diversity, PERMANOVA)
library(iNEXT)       # Coverage-based rarefaction
library(mgcv)        # Generalized additive models
library(xgboost)     # Gradient boosting
library(caret)       # Machine learning utilities
library(lme4)        # Linear mixed-effects models
library(lmerTest)    # p-values for lme4
library(MuMIn)       # R^2 for mixed models

#### Load Data Sheets ####
# Load all data sheets from the Excel database
diversity <- read_excel("Database.xlsx", sheet = "Diversity")
lat       <- read_excel("Database.xlsx", sheet = "Lat")
ichnotaxa <- read_excel("Database.xlsx", sheet = "Ichnotaxa")
deposits  <- read_excel("Database.xlsx", sheet = "Deposits")

##### Descriptive Statistics and Basic Visualizations ####
# Boxplot of ichnogenus richness by climate zone
ggplot(diversity, aes(x = Zone, y = Ichnogenus, fill = Zone)) +
    geom_boxplot() +
    stat_summary(fun = mean, geom = "point", shape = 18, 
                 size = 3, color = "black") +
    labs(y = "Number of Ichnogenera") +
    theme_bw()

## Correlation between ichnodiversity and ichnodisparity
# Spearman rank correlation
cor.test(diversity$Ichnogenus, diversity$Designs,
         method = "spearman", use = "complete.obs")

# Linear regression: ichnodisparity as a function of ichnodiversity
lm_summary <- lm(Ichnogenus ~ Designs, data = diversity)
summary(lm_summary)

# Scatter plot with regression line
ggplot(diversity, aes(x = Ichnogenus, y = Designs)) +
    geom_point() +
    geom_smooth(method = "lm", colour = "green") +
    labs(x = "Ichnogenus Richness", y = "Architectural Designs") +
    theme_bw()

#### Coverage-Based Rarefaction ####
# Prepare metadata and occurrence matrix
metadata <- deposits %>%
    select(Basin, Unit, Time, Paleolat, Zone, Setting, Category)

# Extract occurrence matrix
occ <- deposits %>%
    select(9:ncol(deposits)) %>%
    mutate(across(everything(), as.numeric))

# Coverage-based rarefaction by climate zone
zone.list <- split(occ, metadata$Zone)

# Calculate incidence frequencies for each zone
zone.incidence <- lapply(zone.list, function(x) {
    c(nrow(x), colSums(x))  # First element = number of sites, then incidences
})

# Run iNEXT for incidence-frequency data
zone.iNEXT <- iNEXT(zone.incidence, datatype = "incidence_freq", q = 0)

# Extract richness at 95% sample coverage
zone.richness <- zone.iNEXT$iNextEst$coverage_based %>%
    group_by(Assemblage) %>%
    slice_min(abs(SC - 0.95))

print(zone.richness)

# Coverage-based rarefaction by depositional setting
setting.list <- split(occ, metadata$Setting)

setting.incidence <- lapply(setting.list, function(x) {
    c(nrow(x), colSums(x))
})

setting.iNEXT <- iNEXT(setting.incidence, datatype = "incidence_freq", q = 0)

setting.richness <- setting.iNEXT$iNextEst$coverage_based %>%
    group_by(Assemblage) %>%
    slice_min(abs(SC - 0.95))


print(setting.richness)

# Coverage-based rarefaction by category
category.list <- split(occ, metadata$Category)

category.incidence <- lapply(category.list, function(x) {
    c(nrow(x), colSums(x))
})

category.iNEXT <- iNEXT(category.incidence, datatype = "incidence_freq", q = 0)

category.richness <- category.iNEXT$iNextEst$coverage_based %>%
    group_by(Assemblage) %>%
    slice_min(abs(SC - 0.95))

print(category.richness)

# Coverage-based rarefaction by time slice
time.list <- split(occ, metadata$Time)

time.incidence <- lapply(time.list, function(x) {
    c(nrow(x), colSums(x))
})

time.iNEXT <- iNEXT(time.incidence, datatype = "incidence_freq", q = 0)

time.richness <- time.iNEXT$iNextEst$coverage_based %>%
    group_by(Assemblage) %>%
    slice_min(abs(SC - 0.95))

print(time.richness)

## Visualize standardized richness
# By climate zone
ggplot(zone.richness, aes(Assemblage, qD, fill = Assemblage)) +
    geom_col() +
    coord_cartesian(ylim = c(0, 150)) +
    theme_bw() +
    scale_fill_viridis_d() +
    labs(y = "Coverage-standardized richness", x = "Climate Zone")

# By depositional setting
ggplot(setting.richness, aes(Assemblage, qD, fill = Assemblage)) +
    geom_col() +
    coord_cartesian(ylim = c(0, 150)) +
    theme_bw() +
    scale_fill_viridis_d() +
    labs(y = "Coverage-standardized richness", x = "Depositional Setting")

# By time slice
ggplot(time.richness, aes(Assemblage, qD, fill = Assemblage)) +
    geom_col() +
    coord_cartesian(ylim = c(0, 150)) +
    theme_bw() +
    scale_fill_viridis_d() +
    labs(y = "Coverage-standardized richness", x = "Time Slice")

# By category
ggplot(category.richness, aes(Assemblage, qD, fill = Assemblage)) +
    geom_col() +
    coord_cartesian(ylim = c(0, 150)) +
    theme_bw() +
    scale_fill_viridis_d() +
    labs(y = "Coverage-standardized richness", x = "Category")

## Collector curves
# Sample coverage-based rarefaction curves
ggiNEXT(zone.iNEXT, type = 1)

# Sample-size based rarefaction curves
p_size <- ggiNEXT(zone.iNEXT, type = 1)
p_size +
    theme_bw(base_size = 14) +
    labs(x = "Number of deposits", y = "Estimated richness")

#### XGBoost ####
# Prepare training data
training <- metadata %>%
    distinct(Zone, Paleolat, Setting, Category) %>%
    left_join(zone.richness, by = c("Zone" = "Assemblage"))

# Create predictor matrix (one-hot encoding)
X <- model.matrix(
    ~ Paleolat + Setting + Category + Zone - 1,
    data = training
)
y <- training$qD  # Response variable = coverage-standardized richness

# Cross-validation
dtrain <- xgb.DMatrix(data = X, label = y)

cv <- xgb.cv(
    params = list(
        objective = "reg:squarederror",
        eta = 0.05,
        max_depth = 3,
        subsample = 0.8,
        colsample_bytree = 0.8
    ),
    data = dtrain,
    nfold = 10,
    nrounds = 500,
    early_stopping_rounds = 20,
    metrics = "rmse",
    verbose = 0
)

best_iter <- cv$best_iteration

# Train final model
model <- xgb.train(
    params = list(
        objective = "reg:squarederror",
        eta = 0.05,
        max_depth = 3,
        subsample = 0.8,
        colsample_bytree = 0.8
    ),
    data = dtrain,
    nrounds = best_iter
)

# Extract and visualize feature importance
importance <- xgb.importance(model = model)
xgb.plot.importance(importance)
print(importance)

#### Latitudinal Diversity Patterns ####
# Prepare data: observed richness per deposit
metadata <- deposits %>% select(N, Paleolat)
occ <- deposits %>% select(9:ncol(deposits)) %>%
    mutate(across(everything(), as.numeric))

raw_data <- metadata %>%
    mutate(
        Paleolat = abs(Paleolat),  # Use absolute paleolatitude
        Richness = rowSums(occ)
    )

# Bin data
metadata <- metadata %>%
    mutate(
        Paleolat = abs(Paleolat),
        Bin = floor(Paleolat / 5) * 5
    )

# Calculate pooled richness per bin
bin_occ <- occ %>%
    mutate(Bin = metadata$Bin) %>%
    group_by(Bin) %>%
    summarise(across(everything(), sum))

# Convert to presence-absence per bin
bin_occ_pa <- bin_occ %>%
    mutate(across(-Bin, ~ ifelse(. > 0, 1, 0)))

bin_richness <- bin_occ_pa %>%
    mutate(Richness = rowSums(across(-Bin)))

## Statistical tests
# Spearman correlation
cor.test(bin_richness$Bin, bin_richness$Richness, method = "spearman")

# Linear model
lm.bin <- lm(Richness ~ Bin, data = bin_richness)
summary(lm.bin)

# Generalized additive model
gam.bin <- gam(Richness ~ s(Bin), data = bin_richness)
summary(gam.bin)

# Visualize latitudinal diversity gradient
bin_n <- metadata %>%
    group_by(Bin) %>%
    summarise(nDeposits = n(), .groups = "drop")

bin_richness <- bin_richness %>%
    left_join(bin_n, by = "Bin")

ggplot(bin_richness, aes(x = Bin, y = Richness, size = nDeposits)) +
    geom_point(alpha = 0.8, colour = "black", shape = 20) +
    geom_smooth(method = "lm", colour = "blue", linewidth = 1, se = TRUE) +
    geom_smooth(method = "gam", colour = "red", linewidth = 1, se = TRUE) +
    scale_size_continuous(name = "Number of deposits", range = c(2.5, 10)) +
    scale_x_reverse() +
    scale_y_binned(limits = c(0, 50), n.breaks = 5) +
    labs(x = "Absolute paleolatitude", y = "Observed ichnogeneric richness") +
    theme_bw()

# Ichnodisparity vs. paleolatitude
cor.test(lat$Designs, lat$Paleolat, method = "spearman", use = "complete.obs")

lm.designs <- lm(Designs ~ Paleolat, data = lat)
summary(lm.designs)

gam.designs <- gam(Designs ~ s(Paleolat), data = lat)
summary(gam.designs)

#### Latitudinal Patterns by Depositional Setting ####
# Prepare data
metadata <- deposits %>% select(N, Paleolat, Setting)
occ <- deposits %>% select(9:ncol(deposits)) %>%
    mutate(across(everything(), as.numeric))

dados <- bind_cols(metadata, occ) %>%
    mutate(
        Paleolat = abs(Paleolat),
        Bin = floor(Paleolat / 5) * 5
    )

# Calculate richness per bin for each setting
bin_list <- split(dados, dados$Setting)

richness_list <- lapply(names(bin_list), function(s) {
    df <- bin_list[[s]]
    
    occ_mat <- df %>%
        select(-(N:Setting), -Paleolat, -Bin)
    
    bin_occ <- occ_mat %>%
        mutate(Bin = df$Bin) %>%
        group_by(Bin) %>%
        summarise(across(everything(), sum), .groups = "drop")
    
    bin_occ_pa <- bin_occ %>%
        mutate(across(-Bin, ~ ifelse(. > 0, 1, 0)))
    
    richness <- bin_occ_pa %>%
        mutate(Richness = rowSums(across(-Bin))) %>%
        select(Bin, Richness)
    
    ndep <- df %>%
        count(Bin, name = "nDeposits")
    
    richness %>%
        left_join(ndep, by = "Bin") %>%
        mutate(Setting = s)
})

setting_richness <- bind_rows(richness_list)

## Statistical tests per setting
# Spearman correlations
cor_results <- setting_richness %>%
    group_by(Setting) %>%
    group_modify(~ {
        test <- cor.test(.x$Bin, .x$Richness, method = "spearman")
        tibble(rho = unname(test$estimate), p.value = test$p.value)
    })
print(cor_results)

# Linear models
lm_results <- lapply(split(setting_richness, setting_richness$Setting),
                     function(x) summary(lm(Richness ~ Bin, data = x)))
print(lm_results)

# Generalized additive models
gam_results <- lapply(split(setting_richness, setting_richness$Setting),
                      function(x) summary(gam(Richness ~ s(Bin), data = x)))
print(gam_results)

# Visualize latitudinal patterns by setting
ggplot(setting_richness, aes(Bin, Richness, colour = Setting)) +
    geom_point(aes(size = nDeposits), alpha = 0.8) +
    geom_smooth(method = "lm", linewidth = 1) +
    geom_smooth(method = "gam", linewidth = 1) +
    facet_wrap(~ Setting, scales = "free_y") +
    scale_size(range = c(3, 9), name = "Deposits") +
    scale_x_reverse() +
    theme_bw() +
    labs(x = "Absolute palaeolatitude", y = "Observed ichnogeneric richness")

#### Alpha and Beta Diversity by Climate Zone ####
# Prepare data
metadata <- deposits %>% select(Unit, Zone)
occ <- deposits %>% select(9:ncol(deposits))
ichno <- bind_cols(metadata, occ)

# Alpha diversity (Shannon, Simpson, richness)
alpha_div <- ichno %>%
    rowwise() %>%
    mutate(
        Richness = sum(c_across(-(Unit:Zone))),
        Shannon = diversity(c_across(-(Unit:Zone)), index = "shannon"),
        Simpson = diversity(c_across(-(Unit:Zone)), index = "simpson")
    ) %>%
    ungroup()

# Summary by climate zone
alpha_summary <- alpha_div %>%
    group_by(Zone) %>%
    summarise(
        n = n(),
        Richness_mean = mean(Richness),
        Richness_sd   = sd(Richness),
        Shannon_mean  = mean(Shannon),
        Shannon_sd    = sd(Shannon),
        Simpson_mean  = mean(Simpson),
        Simpson_sd    = sd(Simpson)
    )
print(alpha_summary)

# Beta diversity (Jaccard dissimilarity)
jaccard_dist <- vegdist(occ, method = "jaccard", binary = TRUE)

# Multivariate dispersion (distance to centroid)
beta_disp <- betadisper(jaccard_dist, group = metadata$Zone)
anova(beta_disp)
permutest(beta_disp)

# Summary by climate zone
beta_summary <- tibble(
    Zone = metadata$Zone,
    Distance = beta_disp$distances
) %>%
    group_by(Zone) %>%
    summarise(
        Beta = mean(Distance),
        SD = sd(Distance),
        N = n(),
        .groups = "drop"
    )
print(beta_summary)

# PERMANOVA
adonis2(occ ~ Zone, data = metadata, method = "jaccard",
        binary = TRUE, permutations = 999)

## Visualize alpha and beta diversity
# Richness
ggplot(alpha_summary, aes(x = Zone, y = Richness_mean, fill = Zone)) +
    geom_col(width = 0.7, colour = "black") +
    geom_errorbar(aes(ymin = Richness_mean - Richness_sd,
                      ymax = Richness_mean + Richness_sd),
                  width = 0.15, linewidth = 0.8) +
    theme_bw(base_size = 14) +
    labs(x = "", y = "Mean richness")

# Shannon diversity
ggplot(alpha_summary, aes(x = Zone, y = Shannon_mean, fill = Zone)) +
    coord_cartesian(ylim = c(0, 2.3)) +
    geom_col(colour = "black") +
    geom_errorbar(aes(ymin = Shannon_mean - Shannon_sd,
                      ymax = Shannon_mean + Shannon_sd)) +
    theme_bw(base_size = 14) +
    labs(x = "", y = "Mean Shannon diversity")

# Simpson diversity
ggplot(alpha_summary, aes(x = Zone, y = Simpson_mean, fill = Zone)) +
    geom_col(colour = "black") +
    coord_cartesian(ylim = c(0, 1)) +
    geom_errorbar(aes(ymin = Simpson_mean - Simpson_sd,
                      ymax = Simpson_mean + Simpson_sd),
                  width = 0.50) +
    theme_bw() +
    labs(x = "", y = "Mean Simpson diversity")

# Beta diversity
beta_plot <- tibble(Zone = metadata$Zone, Distance = beta_disp$distances)

ggplot(beta_plot, aes(x = Zone, y = Distance, fill = Zone)) +
    coord_cartesian(ylim = c(0, 1)) +
    stat_summary(fun = mean, geom = "col", width = 0.65, colour = "black") +
    stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
                 geom = "errorbar", width = 0.5) +
    theme_bw() +
    labs(x = "", y = "Distance to centroid (β diversity)")

#### Temporal Trends in Diversity ####
# Diversity through time by climate zone
metadata <- deposits %>% select(Time, Zone)
occ <- deposits %>% select(9:ncol(deposits)) %>%
    mutate(across(everything(), as.numeric))

dados <- bind_cols(metadata, occ)

# Convert Time to numeric ages
dados <- dados %>%
    mutate(Time = recode(Time,
                         "Early Mississippian" = 360,
                         "Early Carb." = 340,
                         "Mid Carb." = 320,
                         "Late Carb.-Early Perm." = 300,
                         "Late Cisuralian" = 280))

occ.cols <- names(occ)

# Pooled richness per time slice and climate zone
zone_time <- dados %>%
    group_by(Zone, Time) %>%
    summarise(across(all_of(occ.cols), max), .groups = "drop") %>%
    mutate(Richness = rowSums(across(all_of(occ.cols))))

# Visualize
ggplot(zone_time, aes(x = Time, y = Richness, colour = Zone, group = Zone)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    scale_x_reverse(breaks = c(360, 340, 320, 300, 280)) +
    theme_bw() +
    labs(x = "Age (Ma)", y = "Observed ichnogeneric richness",
         colour = "Climatic zone")

# Diversity through time by depositional setting
metadata <- deposits %>% select(Time, Setting)
dados <- bind_cols(metadata, occ)

dados <- dados %>%
    mutate(Time = recode(Time,
                         "Early Mississippian" = 360,
                         "Early Carb." = 340,
                         "Mid Carb." = 320,
                         "Late Carb.-Early Perm." = 300,
                         "Late Cisuralian" = 280))

setting_time <- dados %>%
    group_by(Setting, Time) %>%
    summarise(across(all_of(occ.cols), max), .groups = "drop") %>%
    mutate(Richness = rowSums(across(all_of(occ.cols))))

# Visualize
ggplot(setting_time, aes(x = Time, y = Richness, colour = Setting, 
                         group = Setting)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    scale_x_reverse(breaks = c(360, 340, 320, 300, 280)) +
    theme_bw(base_size = 14) +
    labs(x = "Age (Ma)", y = "Observed ichnogeneric richness", 
         colour = "Depositional setting")