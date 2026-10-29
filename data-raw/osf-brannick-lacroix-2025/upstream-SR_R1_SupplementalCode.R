####ND Time Series x Music x Attention ####
# Date last modified: 19 June 2025

# Requires the following files to run:
# audio_happy <- readWave("Mendelson10.wav")
# audio_sad <- readWave("Shostavich10.wav")
#m <- music eye data <- read.csv("MusicSaccadeReport.csv") Music Database 
# d <- "8-25_Database.csv" ANT Database 

# Scipt Written by: S. Fissel Brannick (TAG Lab, Midwestern University Glendale) 
# Blinking indexes dynamic attending during and after music listening
# Fissel Brannick & LaCroix (Purdue)
# https://doi.org/10.31234/osf.io/9md3g

#### ii. Set current directory ####
setwd("./RCODE/NDts") # set working directory
#### i. prep workspace ####
rm(list=ls())
#LOAD PACKAGES 

# Acoustic Analysis Packages 
library(tuneR)
library(seewave)
library(warbleR)
library(av)
library(audio)
library(RcppRoll)
#packages <- c("tuneR", "seewave", "warbleR", "av", "audio", "RcppRoll")
#lapply(packages, citation)

# Datafgrame manipulation packages 
library(dplyr)
library(scales)
library(stats)
library(fuzzyjoin)
library(zoo)
library(tidyr)
library(pracma)  
library(purrr)
library(stringr)
library(broom)
#packages <- c("tuneR", "seewave", "warbleR", "av", "audio", "RcppRoll")
#lapply(packages, citation)

# MEMs
library(lme4)
library(lmerTest)
library(emmeans)

#Plotting Package 
library(ggplot2)
#packages <- c("lme4", "lmerTest", "emmeans", "ggplot2")
#lapply(packages, citation)

#### 1. Load and Process Audio Files ####
happy_audio <- readWave("Mendelson10.wav")
sad_audio   <- readWave("Shostavich10.wav")

# Convert to mono 
happy_mono <- mono(happy_audio, "both")
sad_mono   <- mono(sad_audio, "both")

# Normalize amplitude
normalize_wave <- function(wave) {
  wave@left <- wave@left / 2^(wave@bit - 1)
  wave@left <- wave@left / max(abs(wave@left))
  return(wave)
}

happy_norm <- normalize_wave(happy_mono)
sad_norm   <- normalize_wave(sad_mono)

# Loop & Trim to 600 sec
ensure_duration_600 <- function(wave, target_sec = 600) {
  current_duration <- length(wave@left) / wave@samp.rate
  
  while (current_duration < target_sec) {
    wave <- bind(wave, wave)
    current_duration <- length(wave@left) / wave@samp.rate
  }
  
  if (current_duration > target_sec) {
    wave <- cutw(wave, from = 0, to = target_sec, output = "Wave")
  }
  
  return(wave)
}

#Check check 
happy_norm <- ensure_duration_600(happy_norm)
sad_norm   <- ensure_duration_600(sad_norm)

#time vectors 
happy_time <- seq(0, length(happy_norm@left)-1) / happy_norm@samp.rate
sad_time   <- seq(0, length(sad_norm@left)-1) / sad_norm@samp.rate

#check 
length(happy_norm@left) / happy_norm@samp.rate  # should return 600
length(sad_norm@left) / sad_norm@samp.rate      # should return 600

#### 2. Compute Spectral Centroid ########
get_spectral_centroid_ts <- function(wave, wl = 2048, step = 882) {
  spectro_out <- seewave::spectro(wave, f = wave@samp.rate, wl = wl, ovlp = 100 * (1 - step / wl), plot = FALSE)
  
  amp_matrix <- spectro_out$amp    
  freqs <- spectro_out$freq        
  
 centroids <- apply(amp_matrix, 2, function(col) {
    if (all(is.na(col)) || sum(col, na.rm = TRUE) == 0) {
      return(NA)
    }
    sum(freqs * col, na.rm = TRUE) / sum(col, na.rm = TRUE)
  })
  
  centroid_smoothed <- RcppRoll::roll_mean(centroids, n = 5, fill = NA)
  centroid_filled <- zoo::na.approx(centroid_smoothed, na.rm = FALSE)
  
  return(centroid_filled)
}

#run Fx for happy
happy_centroid <- get_spectral_centroid_ts(happy_norm)

#downsample
happy_time <- seq(0, by = 1 / 44.1e3 * 882, length.out = length(happy_centroid))

happy_cent_df <- data.frame(Time_sec = happy_time, Spectral_Centroid = happy_centroid) %>%
  dplyr::mutate(Bin = floor(Time_sec / 0.25)) %>%
  dplyr::group_by(Bin) %>%
  dplyr::reframe(
    Time_sec = Bin[1] * 0.25,  # anchor bin time
    Spectral_Centroid = mean(Spectral_Centroid, na.rm = TRUE)
  ) %>%
  dplyr::mutate(
    Music = "HAPPY",
    Spectral_Centroid_z = as.numeric(scale(Spectral_Centroid))
  )

##### check check 
nrow(happy_cent_df)        # 2400
range(happy_cent_df$Time_sec)  #0 to 600

sum(is.na(happy_cent_df$Spectral_Centroid_z))  #0

str(happy_cent_df)
dplyr::glimpse(happy_cent_df)

#run Fx for sad
sad_centroid <- get_spectral_centroid_ts(sad_norm)

# 2. Time vector 
sad_time <- seq(0, by = 1 / 44.1e3 * 882, length.out = length(sad_centroid))

# 3. Create binned dataframe and scale
sad_cent_df <- data.frame(Time_sec = sad_time, Spectral_Centroid = sad_centroid) %>%
  dplyr::mutate(Bin = floor(Time_sec / 0.25)) %>%
  dplyr::group_by(Bin) %>%
  dplyr::reframe(
    Time_sec = Bin[1] * 0.25,  # anchor bin time
    Spectral_Centroid = mean(Spectral_Centroid, na.rm = TRUE)
  ) %>%
  dplyr::mutate(
    Music = "SAD",
    Spectral_Centroid_z = as.numeric(scale(Spectral_Centroid))
  )

##### check check 
nrow(sad_cent_df)        # 2400
range(sad_cent_df$Time_sec)  #0 to 600

sum(is.na(sad_cent_df$Spectral_Centroid_z))  #0

str(sad_cent_df)
dplyr::glimpse(sad_cent_df)

#### 3. Load and Process Blink Prob ~ Music #####
#columns: Time_ms, BlinkLogit, Music, P

m <- read.csv("MusicSaccadeReport.csv")

# Filter to HAPPY and SAD 
music_data <- dplyr::filter(m, Music %in% c("HAPPY", "SAD"))

# Bin to 250ms intervals - calculate blink probability per ###
music_data <- music_data %>%
  dplyr::mutate(Bin = floor(Time_ms / 250)) %>%
  dplyr::group_by(Music, P, Bin) %>%
  dplyr::summarise(Blink_Prob = mean(BlinkLogit, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(Time_sec = Bin * 0.25) %>%
  dplyr::group_by(Music) %>%
  dplyr::mutate(Blink_Prob_z = scale(Blink_Prob)) %>%
  dplyr::ungroup()

########### combine all ############
spec_centroid_combined <- dplyr::bind_rows(happy_cent_df, sad_cent_df) %>%
  dplyr::select(Music, Time_sec, Spectral_Centroid_z)

# Round to avoid floating point mismatch
spec_centroid_combined <- dplyr::mutate(spec_centroid_combined, Time_sec = round(Time_sec, 2))
music_data <- dplyr::mutate(music_data, Time_sec = round(Time_sec, 2))

# left join 
plot_data_scaled <- dplyr::left_join(music_data, spec_centroid_combined, by = c("Music", "Time_sec"))

####### check check check 
str(plot_data_scaled)
dplyr::glimpse(plot_data_scaled)
summary(plot_data_scaled$Spectral_Centroid_z)   
sum(is.na(plot_data_scaled$Spectral_Centroid_z))  

# change from matrix to linear numeric 
plot_data_scaled <- dplyr::mutate(plot_data_scaled, Blink_Prob_z = as.numeric(Blink_Prob_z))

#### $. Sliding Window CCF ####

# Check bins
win_check <- plot_data_scaled %>%
  dplyr::filter(!is.na(Spectral_Centroid_z), !is.na(Blink_Prob_z)) %>%
  dplyr::group_by(P, Music) %>%
  dplyr::summarise(n_bins = dplyr::n(), .groups = "drop") %>%
  dplyr::arrange(n_bins)

summary(win_check$n_bins)
dplyr::filter(win_check, n_bins < 200)  # Only one case will show

# CCF
compute_windowed_ccf <- function(data, win_size = 200, max_lag = 10) {
  results <- list()
  
  unique_combinations <- unique(data[c("P", "Music")])
  
  for (i in 1:nrow(unique_combinations)) {
    subj <- unique_combinations$P[i]
    mus  <- unique_combinations$Music[i]
    
    df <- dplyr::filter(data, P == subj, Music == mus)
    df <- df[complete.cases(df[, c("Spectral_Centroid_z", "Blink_Prob_z")]), ]
    
    if (nrow(df) < win_size) next
    
    for (start in seq(1, nrow(df) - win_size, by = 20)) {
      end <- start + win_size - 1
      window <- df[start:end, ]
      
      if (nrow(window) == win_size) {
        ccf_obj <- ccf(window$Spectral_Centroid_z, window$Blink_Prob_z,
                       plot = FALSE, lag.max = max_lag)
        max_idx <- which.max(abs(ccf_obj$acf))
        peak_lag <- ccf_obj$lag[max_idx]
        peak_val <- ccf_obj$acf[max_idx]
        
        results[[length(results) + 1]] <- data.frame(
          P = subj,
          Music = mus,
          Window_Start = window$Time_sec[1],
          Window_End = window$Time_sec[nrow(window)],
          Peak_Lag = peak_lag,
          Peak_CCF = peak_val
        )
      }
    }
  }
  return(dplyr::bind_rows(results))
}

# Run CCF
sliding_ccf_df <- compute_windowed_ccf(plot_data_scaled)
head(sliding_ccf_df)

# plot2 
ggplot(sliding_ccf_df, aes(x = Window_Start, y = Peak_CCF, group = P, color = P)) +
  geom_line(alpha = 0.6) +
  facet_wrap(~ Music, scales = "free_x") +
  labs(
    title = "Sliding-Window CCF: Blink ↔ Spectral Centroid",
    x = "Time in Music (sec)",
    y = "Peak Cross-Correlation (CCF)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

# plot 3 
# Create segment labels based on Window_Start
sliding_ccf_df$Segment <- cut(
  sliding_ccf_df$Window_Start,
  breaks = c(0, 200, 400, 600),
  labels = c("0–200s", "200–400s", "400–600s"),
  right = FALSE
)

# Rename music conditions
sliding_ccf_df$Music <- dplyr::recode(
  sliding_ccf_df$Music,
  "HAPPY" = "High Dynamic",
  "SAD" = "Low Dynamic"
)

# Plot
ggplot(sliding_ccf_df, aes(x = Window_Start, y = Peak_CCF, group = P, color = P)) +
  geom_line(alpha = 0.6) +
  facet_grid(Music ~ Segment, scales = "free_x") +
  labs(
    title = "Sliding-Window CCF: Blink ↔ Spectral Centroid",
    x = "Time in Music (sec)",
    y = "Peak Cross-Correlation (CCF)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

###### plot with smoother
ggplot(sliding_ccf_df, aes(x = Window_Start, y = Peak_CCF, group = P)) +
  geom_line(alpha = 0.4) +
  geom_smooth(aes(group = Music), method = "loess", color = "black", se = FALSE, linewidth = 1.2) +
  facet_wrap(~Music) +
  labs(
    title = "Sliding-Window CCF: Blink ↔ Spectral Centroid",
    x = "Time in Music (sec)",
    y = "Peak Cross-Correlation (CCF)"
  ) +
  theme_minimal()

##### check ######
dplyr::glimpse(sliding_ccf_df)

# Time segment
sliding_ccf_df <- sliding_ccf_df %>%
  dplyr::mutate(
    Segment = dplyr::case_when(
      Window_Start < 300 ~ "First 5 min",
      Window_Start >= 300 ~ "Last 5 min"
    )
  )

#Plot by time (poly2)
ggplot(sliding_ccf_df, aes(x = Window_Start, y = Peak_CCF, color = Music)) +
  geom_point(alpha = 0.3) +
  geom_smooth(se = FALSE, method = "loess", span = 0.3, linewidth = 1.2) +
  facet_wrap(~Segment, scales = "free_x") +
  labs(
    title = "Sliding-Window CCF Trends by Music Segment",
    x = "Time in Music (sec)",
    y = "Peak Cross-Correlation (CCF)"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

# Assign Segment labels
sliding_ccf_df <- sliding_ccf_df %>%
  dplyr::filter(Window_Start >= 0, Window_Start < 600) %>%
  dplyr::mutate(
    Segment = dplyr::case_when(
      Window_Start < 300 ~ "0–300s",
      Window_Start >= 300 ~ "300–600s"
    ),
    Segment = factor(Segment, levels = c("0–300s", "300–600s")),
    Music = factor(Music)
  )

# Fit the LMER model
model_ccf <- lmer(
  Peak_CCF ~ Segment * Music + (1 | P),
  data = sliding_ccf_df
)

model_poly2 <- lmer(Peak_CCF ~ poly(Window_Start, 2) * Music + (1 | P), data = sliding_ccf_df)

model_poly3 <- lmer(Peak_CCF ~ poly(Window_Start, 3) * Music + (1 | P), data = sliding_ccf_df)

model_poly4 <- lmer(Peak_CCF ~ poly(Window_Start, 4) * Music + (1 | P), data = sliding_ccf_df)

anova(model_ccf, model_poly2, model_poly3, model_poly4) #poly 3 is best 
summary(model_poly3)
performance::r2(model_poly3)

model_poly3lm <- lm(Peak_CCF ~ poly(Window_Start, 3) * Music, data = sliding_ccf_df)
summary(model_poly3lm)

# plot 
# prediction (0 to 600s)
newdata <- expand.grid(
  Window_Start = seq(0, 600, by = 1),
  Music = c("HAPPY", "SAD")
)

# Poly 3 
poly_basis <- poly(sliding_ccf_df$Window_Start, 3)
attr(poly_basis, "coefs") <- attr(poly(sliding_ccf_df$Window_Start, 3), "coefs")

# Poly3 -- data 
newdata <- newdata %>%
  mutate(
    poly1 = predict(poly_basis, newdata$Window_Start)[, 1],
    poly2 = predict(poly_basis, newdata$Window_Start)[, 2],
    poly3 = predict(poly_basis, newdata$Window_Start)[, 3]
  )

#Pred
coefs <- fixef(model_poly3)

newdata <- newdata %>%
  rowwise() %>%
  mutate(
    Predicted_CCF = coefs["(Intercept)"] +
      poly1 * coefs["poly(Window_Start, 3)1"] +
      poly2 * coefs["poly(Window_Start, 3)2"] +
      poly3 * coefs["poly(Window_Start, 3)3"] +
      ifelse(Music == "SAD", coefs["MusicSAD"], 0) +
      ifelse(Music == "SAD", poly1 * coefs["poly(Window_Start, 3)1:MusicSAD"], 0) +
      ifelse(Music == "SAD", poly2 * coefs["poly(Window_Start, 3)2:MusicSAD"], 0) +
      ifelse(Music == "SAD", poly3 * coefs["poly(Window_Start, 3)3:MusicSAD"], 0)
  ) %>%
  ungroup()

newdata$Music <- factor(newdata$Music,
                        levels = c("HAPPY", "SAD"),
                        labels = c("High Dynamic", "Low Dynamic"))

# Plot
ggplot(newdata, aes(x = Window_Start, y = Predicted_CCF, color = Music)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "Predicted Sliding-Window CCF Over Time by Music Dynamics",
    x = "Time in Music (sec)",
    y = "Predicted Peak CCF"
  ) +
  theme_minimal(base_size = 14)

############### 5. Process ANT Data ######################## 
# Load & Clean ANT Dataset ###
d <- read.csv("8-25_Database.csv")

d$P <- as.factor(d$P)
d$BlinkLogistic <- as.factor(d$BlinkLogistic)
d$DV <- factor(d$DV, levels = c("RT", "BlinkStart"))
d$Congruency <- factor(toupper(as.character(d$Congruency)), levels = c("NEUTRAL", "INCONGRUENT"))
d$Music_Cx <- factor(d$Music_Cx, levels = c("SILENCE", "HAPPY", "SAD"))
d$ANT_Cx <- as.factor(d$ANT_Cx)

# Identify last timestamp per trial
last_times <- d %>%
  group_by(TRIAL) %>%
  summarise(LastTime = max(TimeStamp), .groups = "drop")

# Join to original data and filter from 2903 to trial-specific end
filtered_data.flanker <- d %>%
  dplyr::inner_join(last_times, by = "TRIAL") %>%
  dplyr::filter(TimeStamp >= 2903 & TimeStamp <= LastTime)

# Filter for Incon - Neutral Executive Control trials
filtered_data.flanker.ec <- dplyr::filter(
  filtered_data.flanker, 
  EC_IncongruentNeutral == "Endogenous_EC"
)

######## 6. EC Models -- Incongruent Neutral ######### 
flankeronset.ec <- glmer(BlinkLogistic ~ (ANT_Cx * Congruency * Music_Cx) + 
                           (1 + ANT_Cx | P), 
                         data = filtered_data.flanker.ec, 
                         family = "binomial", 
                         control = glmerControl(optimizer = "bobyqa"))
summary(flankeronset.ec)
performance::check_collinearity(flankeronset.ec)
performance::r2(flankeronset.ec)

#relevel
filtered_data.flanker.ec$Music_Cxsadhsi <- factor(filtered_data.flanker.ec$Music_Cx, levels = c("SAD", "HAPPY", "SILENCE"))
flankeronset.ec.releveled <- glmer(BlinkLogistic ~ (ANT_Cx * Congruency * Music_Cxsadhsi) + 
                                     (1 | P), 
                                   data = filtered_data.flanker.ec, 
                                   family = "binomial", 
                                   control = glmerControl(optimizer = "bobyqa"))
summary(flankeronset.ec.releveled)

# ec x ccf x silence as happy #
filtered_data.flanker <- d %>%
  dplyr::group_by(P, TRIAL) %>%
  dplyr::filter(TimeStamp >= 2903 & TimeStamp <= max(TimeStamp, na.rm = TRUE)) %>%
  dplyr::ungroup()

# Filter 
filtered_data.flanker.ec <- filtered_data.flanker %>%
  dplyr::filter(
    EC_IncongruentNeutral == "Endogenous_EC",
    Congruency %in% c("NEUTRAL", "INCONGRUENT"),
    Music_Cx %in% c("SILENCE", "HAPPY", "SAD"),
    DV %in% c("RT", "BlinkStart"),
    !is.na(BlinkLogistic)
  ) %>%
  dplyr::mutate(
    P = as.character(P),
    Music_Cx = factor(Music_Cx, levels = c("SILENCE", "HAPPY", "SAD")),
    DV = factor(DV, levels = c("RT", "BlinkStart"))
  )

ccf_model_data <- sliding_ccf_df %>%
  dplyr::mutate(
    P_num = stringr::str_extract(P, "\\d{4}"),
    Music_Cx = Music  # rename for compatibility with ANT dataset
  )

ccf_summary <- ccf_model_data %>%
  dplyr::group_by(P_num, Music_Cx) %>%
  dplyr::summarise(
    mean_CCF = mean(Peak_CCF, na.rm = TRUE),
    peak_CCF = max(Peak_CCF, na.rm = TRUE),
    .groups = "drop"
  )

# Placeholder
ccf_placeholder <- ccf_summary %>% dplyr::filter(Music_Cx == "HAPPY") %>% dplyr::slice(1)

silence_real <- filtered_data.flanker.ec %>%
  dplyr::filter(Music_Cx == "SILENCE") %>%
  dplyr::distinct(P, Music_Cx) %>%
  dplyr::mutate(
    P_num = as.character(P),
    mean_CCF = ccf_placeholder$mean_CCF,
    peak_CCF = ccf_placeholder$peak_CCF
  )

ccf_summary_real <- ccf_summary %>% dplyr::rename(P = P_num)

ccf_summary_combined <- dplyr::bind_rows(ccf_summary_real, silence_real)
str(ccf_summary_combined)

ccf_summary_combined <- ccf_summary_combined %>%
  dplyr::select(P, Music_Cx, mean_CCF, peak_CCF) %>%  
  dplyr::mutate(P = as.character(P))              

blink_summary <- filtered_data.flanker.ec %>%
  dplyr::group_by(P, Music_Cx, ANT_Cx, Congruency) %>%
  dplyr::summarise(
    BlinkProb = mean(as.numeric(as.character(BlinkLogistic)), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = ANT_Cx,
    values_from = BlinkProb,
    names_prefix = "BlinkProb_"
  ) %>%
  dplyr::mutate(
    Blink_change = BlinkProb_ANT_2 - BlinkProb_ANT_1,
    P = as.character(P)
  )

str(blink_summary)

filtered_data_joined <- filtered_data.flanker.ec %>%
  dplyr::mutate(P = as.character(P)) %>%
  dplyr::left_join(ccf_summary_combined, by = c("P", "Music_Cx")) %>%
  dplyr::left_join(blink_summary, by = c("P", "Music_Cx", "Congruency")) %>%
  dplyr::filter(!is.na(mean_CCF), !is.na(Blink_change)) %>%
  dplyr::mutate(
    mean_CCF_s = scale(mean_CCF)[, 1],
    Blink_change_s = scale(Blink_change)[, 1]
  )

####### checkity check check #######
table(filtered_data_joined$Music_Cx)
table(filtered_data_joined$DV)
table(filtered_data_joined$Congruency)

summary(filtered_data_joined$TimeStamp)

m.music_ccf_blink_2way <- lmer(
  TimeStamp ~ DV * Music_Cx + ANT_Cx + Blink_change_s + mean_CCF_s + 
                 DV * mean_CCF_s +
                 DV * ANT_Cx +
                 Music_Cx:mean_CCF_s +
                 Blink_change_s:Music_Cx +
                 DV * Congruency +
                 (1 + ANT_Cx | P),
  data = filtered_data_joined
)

summary(m.music_ccf_blink_2way)

m.music_ccf_blink_2way_2 <- lmer(
  TimeStamp ~ DV * Music_Cx + ANT_Cx + Blink_change_s + 
    DV * ANT_Cx +
    Blink_change_s:Music_Cx +
    DV * Congruency +
    (1 + ANT_Cx | P),
  data = filtered_data_joined
)

anova(m.music_ccf_blink_2way, m.music_ccf_blink_2way_2) #w CCF is best

filtered_data_joined$Music_Cx <- factor(
  filtered_data_joined$Music_Cx,
  levels = c("SAD", "HAPPY", "SILENCE")
)

performance::r2(m.music_ccf_blink_2way)

############### ec x ccf x silence as sad ###########
# M/Peak CCF 
sad_ccf_summary <- ccf_summary_combined %>%
  dplyr::filter(Music_Cx == "SAD") %>%
  dplyr::summarise(
    mean_CCF = mean(mean_CCF, na.rm = TRUE),
    peak_CCF = mean(peak_CCF, na.rm = TRUE)
  )

silence_real <- filtered_data_joined %>%
  dplyr::filter(Music_Cx == "SILENCE") %>%
  dplyr::distinct(P, Music_Cx)

ccf_summary_silence_sad <- silence_real %>%
  dplyr::mutate(
    mean_CCF = sad_ccf_summary$mean_CCF,
    peak_CCF = sad_ccf_summary$peak_CCF
  )

ccf_summary_sad_happy <- ccf_summary_combined %>%
  dplyr::filter(Music_Cx %in% c("SAD", "HAPPY"))

ccf_summary_final <- dplyr::bind_rows(ccf_summary_sad_happy, ccf_summary_silence_sad)

filtered_data_joined_updated <- filtered_data_joined %>%
  dplyr::select(-mean_CCF, -peak_CCF) %>%  # remove old versions if present
  dplyr::left_join(ccf_summary_final, by = c("P", "Music_Cx")) %>%
  dplyr::filter(!is.na(mean_CCF)) %>%
  dplyr::mutate(
    mean_CCF_s = scale(mean_CCF)[, 1],
    Blink_change_s = scale(Blink_change)[, 1]
  )

######## lmers ########

m.music_ccf_blink_1 <- lmer(
  TimeStamp ~ (DV + Music_Cx + ANT_Cx + Congruency + Blink_change_s + mean_CCF_s)^2 +
    (1 + ANT_Cx | P),
  data = filtered_data_joined_updated
)

summary(m.music_ccf_blink_1)

m.music_ccf_blink_2 <- lmer(
  TimeStamp ~ (DV + Music_Cx + ANT_Cx + Congruency + Blink_change_s)^2 +
    (1 + ANT_Cx | P),
  data = filtered_data_joined_updated
)

summary(m.music_ccf_blink_2)

anova(m.music_ccf_blink_1, m.music_ccf_blink_2) #w CCF is best

#relevel 
filtered_data_joined_updated$Music_Cx <- factor(
  filtered_data_joined_updated$Music_Cx,
  levels = c("SAD", "HAPPY", "SILENCE")
)

