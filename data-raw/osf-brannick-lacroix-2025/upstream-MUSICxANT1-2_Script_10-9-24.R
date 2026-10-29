####ND Time Series x Music x Attention ####
# Date last modified: 18 Nov 2024

# Requires the following files to run:
# d <- "8-25_Database.csv" 
# m <- "MusicTrialBinnedData(100ms)_8-26-24.csv"

# This file contains Data Analysis for the Blink data separated into several parts:
# A. Music Listening Preliminaries 
# B. Visuals of Raw Music Listening Data
# C. Music Listening Models, Model Fit (ROC/AUC), Music Plots 
# D. ANT1-2 ANALYSES  
# E. EXECUTIVE CONTROL (RT & BLINKING - Time Series) - RT Window (3350-4000ms)

# Scipt Written by: S. Fissel Brannick (TAG Lab, Midwestern University Glendale) 
# Blinking indexes dynamic attending during and after music listening
# Fissel Brannick & LaCroix (Purdue)
# https://doi.org/10.31234/osf.io/9md3g

#### i. prep workspace ####
rm(list=ls())

#### ii. Set current directory ####
setwd("./RCODE/NDts") # set working directory

#### iii. LOAD LIBRARIES ####
library(rmarkdown)
library(lme4) ##mixed effects modeling
library(CLME) #Constrainted Linear Mixed Effects
library(ggplot2) ##graphic interface for MEM
library(sjPlot) #for plotting lmer and glmer mods
library(sjmisc) 
library(sjstats) #use for table generation functions
library(effects)
library(lmerTest) #p Values for lmers via 
library(emmeans) #for post hocs
library(plyr)
library(dplyr) #for filtering data
library(psych)
library(grid)
library(vcd)
library(see)
library(jtools)
library(interactions)
library(Matrix)
library(tidyr)
library(zoo) #rolling averages/corr
library(performance) #for r2 functions
library(pROC) #for model fit 
library(ggtext) #custom facet labels

#### A. MUSIC Listening Preliminaries ####
# load & format music data #
m <- read.csv("MusicTrialBinnedData(100ms)_8-26-24.csv")
m$Blink_Logistic <- as.factor(m$Blink_Logistic)
m$Blink_Frequency <- as.numeric(m$Blink_Frequency)
m$Music_NoMusic <- as.factor(m$Music_NoMusic)
m$Music <- as.factor(m$Music)
m$P <- as.factor(m$P)

m$ORDER <- as.ts(m$ORDER)
m$ORDER.s <- scale(m$ORDER)

#check
str(m)

# Filter data sets by listening condition (dplyr) #
happy <- filter(m, Music == "Happy")
sad <- filter(m, Music == "Sad")
silence <- filter(m, Music == "Silence")

########### B. Visuals of Raw Music Listening Data #########
long.music.data <- m %>%
  pivot_longer(cols = c(ORDER.s), 
               names_to = "Variable", 
               values_to = "Value")

# Define custom labels for the facets
custom_labels <- c(
  "Mendelson10.mp3" = "Happy Music (Mendelssohn)",
  "Shostavich10.mp3" = "Sad Music (Shostakovich)",
  "UNDEFINED" = "Control (No Music)"
)

# Density Plot of Blink Logistic
ggplot(long.music.data, aes(x = Value, fill = Blink_Logistic, color = Blink_Logistic)) + 
  geom_density(alpha = 0.2) +  
  facet_grid(~Music) +
  labs(title = "Blinking Probability Distributions during Listening Conditions",
       x = "Order of (scaled) Time (in 100ms bins)",
       y = "Blinking Density") +
  theme_minimal()

########### C. Music Listening Models Models ###########

####### C1a. Happy Music GLMERS #########
#happy <- filter(m, Music == "Happy")
glmer_happybase <-glmer(Blink_Logistic ~ (ORDER.s ) 
                        + (1 | P), data=happy, family="binomial") 
#converged

glmer_happymax1 <-glmer(Blink_Logistic ~ (ORDER.s ) 
                        + (1 + ORDER.s | P), data=happy, family="binomial") 
#converged

glmer_happymax2 <-glmer(Blink_Logistic ~ poly(ORDER.s, 2) 
                        + (1 + ORDER.s | P), data=happy, family="binomial") 
#converged

glmer_happymax3 <-glmer(Blink_Logistic ~ poly(ORDER.s, 3) 
                        + (1 + ORDER.s | P), data=happy, family="binomial") 
#converged

glmer_happymax4 <-glmer(Blink_Logistic ~ poly(ORDER.s, 4) 
                        + (1 + ORDER.s | P), data=happy, family="binomial") 
#converged

anova(glmer_happybase, glmer_happymax1, glmer_happymax2, glmer_happymax3, glmer_happymax4) #3 and 4 are best 
#4 wins 

summary(glmer_happymax4)
glmer_happymax4_r2 <- performance::r2(glmer_happymax4)
print(glmer_happymax4_r2)

########### C1b. Happy Music -- ROC Curve #########
# Blink_Logistic as numeric
happy$Blink_Logistic <- as.numeric(happy$Blink_Logistic)

# Predicted probabilities
happy$predicted_prob <- predict(glmer_happymax4, type = "response")

#NA Checks
happy <- happy %>%
  filter(!is.na(Blink_Logistic) & !is.na(predicted_prob))

#ROC curve
roc_obj <- roc(happy$Blink_Logistic, happy$predicted_prob)

#ROC plot
plot(roc_obj, col = "blue", lwd = 2, main = "ROC Curve ~ Blinking Probability to happy music")
abline(a = 0, b = 1, lty = 2, col = "gray")  #Diagonal line (50/50 chance)
auc_val <- auc(roc_obj) #AUC 
legend("bottomright", legend = paste("AUC =", round(auc_val, 3)), col = "blue", lwd = 2)

##area under curve = 0.804 ~Good

############## C1c. Happy Music -- Plot ###########
new_data <- data.frame(ORDER.s = seq(min(happy$ORDER.s), max(happy$ORDER.s), length.out = 100))

#polynomial terms
poly_terms <- poly(new_data$ORDER.s, 4)
new_data <- cbind(new_data, as.data.frame(poly_terms))
colnames(new_data)[2:5] <- c("Poly1", "Poly2", "Poly3", "Poly4")

#predicted probabilities
new_data$predicted_prob <- predict(glmer_happymax4, newdata = new_data, type = "response", re.form = NA)

#long format for plotting
new_data_long <- new_data %>%
  pivot_longer(cols = starts_with("Poly"), names_to = "Polynomial", values_to = "Value")

#Smooth
ggplot(new_data, aes(x = ORDER.s, y = predicted_prob)) +
  geom_point(color = "blue", alpha = 0.5) +
  geom_smooth(method = "loess", color = "red", size = 1) +
  labs(
    title = "Smoothed Oscillatory Pattern of Blinking Probability",
    x = "ORDER.s (standardized)",
    y = "Predicted Probability of Blinking"
  ) +
  theme_minimal()

######## C2a. Sad Music GLMERS  ############
#sad <- filter(m, Music == "Sad")
sad_glmer_base <-glmer(Blink_Logistic ~ (ORDER.s) 
                       + (1 | P), data=sad, family="binomial") 
#converged

sad_glmer_max1 <-glmer(Blink_Logistic ~ (ORDER.s) 
                       + (1 + ORDER.s | P), data=sad, family="binomial") 
#converged

sad_glmer_max2 <-glmer(Blink_Logistic ~ poly(ORDER.s, 2)  
                       + (1 + ORDER.s | P), data=sad, family="binomial") 
#converged 

sad_glmer_max3 <-glmer(Blink_Logistic ~ poly(ORDER.s, 3)
                       + (1 + ORDER.s | P), data=sad, family="binomial") 
#did not converge 

sad_glmer3 <-glmer(Blink_Logistic ~ poly(ORDER.s, 3)
                   + (1 | P), data=sad, family="binomial") 
#converged 

sad_glmer_max4 <-glmer(Blink_Logistic ~ poly(ORDER.s, 4) 
                       + (1 + ORDER.s | P), data=sad, family="binomial") 
#converged 

anova(sad_glmer_base, sad_glmer_max1, sad_glmer_max2, sad_glmer3, sad_glmer_max4) #2 & 4 are significant
#4 is best 
summary(sad_glmer_max4)
r2_sadglmer_max4 <- performance::r2(sad_glmer_max4) 
print(r2_sadglmer_max4)

######## C2b. Sad Music -- ROC Curve ##########
#Blink_Logistic as numeric
sad$Blink_Logistic <- as.numeric(sad$Blink_Logistic)

#predicted probabilities
sad$predicted_prob <- predict(sad_glmer_max4, type = "response")

#NA Checks 
sad <- sad %>%
  filter(!is.na(Blink_Logistic) & !is.na(predicted_prob))

# ROC curve
roc_obj <- roc(sad$Blink_Logistic, sad$predicted_prob)

plot(roc_obj, col = "blue", lwd = 2, main = "ROC Curve ~ Blinking Probability to Sad Music")
abline(a = 0, b = 1, lty = 2, col = "gray")  #Diagonal (50/50 chance)
auc_val <- auc(roc_obj) 
legend("bottomright", legend = paste("AUC =", round(auc_val, 3)), col = "blue", lwd = 2)

##area under curve = 0.729 ~Moderate

############ C2c. Sad Music Plot #############
new_data.sad <- data.frame(ORDER.s = seq(min(sad$ORDER.s), max(sad$ORDER.s), length.out = 100))

#polynomial terms
poly_terms <- poly(new_data.sad$ORDER.s, 4)
new_data.sad <- cbind(new_data.sad, as.data.frame(poly_terms))
colnames(new_data.sad)[2:5] <- c("Poly1", "Poly2", "Poly3", "Poly4")

#predicted probabilities
new_data.sad$predicted_prob.sad <- predict(sad_glmer_max4, newdata = new_data.sad, type = "response", re.form = NA)

#Smooth
ggplot(new_data.sad, aes(x = ORDER.s, y = predicted_prob.sad)) +
  geom_point(color = "blue", alpha = 0.5) +
  geom_smooth(method = "loess", color = "red", size = 1) +
  labs(
    title = "Smoothed Oscillatory Pattern of Blinking Probability",
    x = "ORDER.s (standardized)",
    y = "Predicted Probability of Blinking"
  ) +
  theme_minimal()

########## C3b. Silence GLMERS #######
modelsilence_base <-glmer(Blink_Logistic ~ (ORDER.s) 
                          + (1 | P), data=silence, family="binomial") 
#converged 

modelsilence_max1 <-glmer(Blink_Logistic ~ (ORDER.s) 
                          + (1 + ORDER.s | P), data=silence, family="binomial") 
#converged 
modelsilence_max2 <-glmer(Blink_Logistic ~ poly(ORDER.s, 2) 
                          + (1 + ORDER.s | P), data=silence, family="binomial") 
#converged 
modelsilence_max3 <-glmer(Blink_Logistic ~ poly(ORDER.s, 3) 
                          + (1 + ORDER.s | P), data=silence, family="binomial") 
#converged
modelsilence_max4 <-glmer(Blink_Logistic ~ poly(ORDER.s, 4) 
                          + (1 + ORDER.s | P), data=silence, family="binomial") 

#did not converge 

modelsilence4 <-glmer(Blink_Logistic ~ poly(ORDER.s, 4) 
                      + (1 | P), data=silence, family="binomial") 
#converged 

anova(modelsilence_base, modelsilence_max1, modelsilence_max2, modelsilence_max3, modelsilence4) #fourth is bad
#model 3

summary(modelsilence_max3)
r2_silenceglmer_max3 <- performance::r2(modelsilence_max3) 
print(r2_silenceglmer_max3)

######## C3b. Silence -- ROC Curve ##########
#Blink_Logistic as numeric
silence$Blink_Logistic <- as.numeric(silence$Blink_Logistic)

#predicted probabilities
silence$predicted_prob <- predict(modelsilence_max3, type = "response")

#NA Checks 
silence <- silence %>%
  filter(!is.na(Blink_Logistic) & !is.na(predicted_prob))

# ROC curve -- Calculate
roc_obj.si <- roc(silence$Blink_Logistic, silence$predicted_prob)

# ROC curve -- Calculate
plot(roc_obj.si, col = "blue", lwd = 2, main = "ROC Curve ~ Blinking Probability to Silence")
abline(a = 0, b = 1, lty = 2, col = "gray")  # diagonal (50/50 chance)
auc_val <- auc(roc_obj.si) 
legend("bottomright", legend = paste("AUC =", round(auc_val, 3)), col = "blue", lwd = 2)

##area under curve = 0.824 ~Good

############### C3c. Silence -- Plot ########
new_data.si <- data.frame(ORDER.s = seq(min(silence$ORDER.s), max(silence$ORDER.s), length.out = 100))

#polynomial terms
poly_terms <- poly(new_data.si$ORDER.s, 3)
new_data.si <- cbind(new_data.si, as.data.frame(poly_terms))
colnames(new_data.si)[2:4] <- c("Poly1", "Poly2", "Poly3")

#predicted probabilities
new_data.si$predicted_prob.si <- predict(modelsilence_max3, newdata = new_data.si, type = "response", re.form = NA)

#Smooth
ggplot(new_data.si, aes(x = ORDER.s, y = predicted_prob.si)) +
  geom_point(color = "blue", alpha = 0.5) +
  geom_smooth(method = "loess", color = "red", size = 1) +
  labs(
    title = "Smoothed Oscillatory Pattern of Blinking Probability",
    x = "ORDER.s (standardized)",
    y = "Predicted Probability of Blinking"
  ) +
  theme_minimal()

######### C4a. Plot -- pp by raw data -- all listening conditions##########
#Modify ORDER.s
reduced_ORDER.s <- seq(min(m$ORDER.s), max(m$ORDER.s), length.out = 110618)

# Happy music predicted probabilities
new_data_happy <- data.frame(ORDER.s = reduced_ORDER.s)
new_data_happy$predicted_prob <- predict(glmer_happymax4, newdata = new_data_happy, type = "response", re.form = NA)
new_data_happy$Music <- "Happy"

# Sad music predicted probabilities
new_data_sad <- data.frame(ORDER.s = reduced_ORDER.s)
new_data_sad$predicted_prob <- predict(sad_glmer_max4, newdata = new_data_sad, type = "response", re.form = NA)
new_data_sad$Music <- "Sad"

# Silence predicted probabilities
new_data_silence <- data.frame(ORDER.s = reduced_ORDER.s)
new_data_silence$predicted_prob <- predict(modelsilence_max3, newdata = new_data_silence, type = "response", re.form = NA)
new_data_silence$Music <- "Silence"

# Combine datasets for all music conditions
predicted_data <- rbind(new_data_happy, new_data_sad, new_data_silence)

# Create bins for ORDER.s using the reduced sequence
m_binned <- m %>%
  mutate(ORDER_bin = cut(ORDER.s, breaks = reduced_ORDER.s, include.lowest = TRUE)) %>%
  group_by(ORDER_bin, Music) %>%
  summarise(
    Blink_Logistic = list(Blink_Logistic), # Keep all raw values as a list
    Density = n()
  ) %>%
  ungroup()

# Map bins back to the reduced ORDER.s
bin_centers <- reduced_ORDER.s[-length(reduced_ORDER.s)] + diff(reduced_ORDER.s) / 2
m_binned <- m_binned %>%
  mutate(ORDER.s = bin_centers[as.numeric(ORDER_bin)])

# Expand Blink_Logistic for density calculations
m_binned_expanded <- m_binned %>%
  unnest(Blink_Logistic)

# Merge with predicted probabilities
final_data <- left_join(m_binned_expanded, predicted_data, by = c("ORDER.s", "Music"))

# Find the nearest ORDER.s in predicted_data for each value in m_binned
m_binned$ORDER.s <- sapply(m_binned$ORDER.s, function(x) {
  predicted_data$ORDER.s[which.min(abs(predicted_data$ORDER.s - x))]
})

# Merge binned data with predicted probabilities
final_data <- left_join(m_binned, predicted_data, by = c("ORDER.s", "Music"))

# Check for non-NA predicted probabilities
summary(final_data$predicted_prob)

# Confirm unique ORDER.s values
unique(final_data$ORDER.s)

facet_labels <- c(
  "Happy" = "<span style='color:red;'>Happy</span>",
  "Sad" = "<span style='color:darkgreen;'>Sad</span>",
  "Silence" = "<span style='color:blue;'>Control</span>"
)

# Updated plot code
ggplot(m, aes(x = ORDER.s)) +
  # Raw blinking density with matching outline and fill for Blink_Logistic
  geom_density(aes(fill = Blink_Logistic, color = Blink_Logistic), alpha = 0.3, position = "identity") +
  # Overlay predicted probabilities with finer granularity
  geom_line(data = final_data, aes(y = predicted_prob, color = Music), 
            linetype = "dotted", size = 0.50) +
  # Facet by Music condition with custom labels
  facet_wrap(~Music, labeller = labeller(Music = facet_labels)) +
  # Adjust titles, labels, and legend
  labs(
    title = "Blinking Probability Distributions with Predicted Probabilities",
    x = "Order of (scaled) Time (in 100ms bins)",
    y = "Blinking Density / Predicted Probability",
    fill = "Blink Logistic", # Legend for Blink_Logistic
    color = "Predicted Probabilities" # Title for Music legend
  ) +
  # Theme adjustments
  theme_minimal() +
  theme(
    strip.text = element_markdown(size = 12), # Use HTML-style facet labels
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  ) +
  # Custom colors for the Blink Logistic legend outline
  scale_fill_manual(values = c("0" = "coral", "1" = "turquoise")) +
  scale_color_manual(values = c("Happy" = "red", "Sad" = "darkgreen", "Silence" = "blue"))

#### C4b. Predicted probabilities only #########

facet_labels <- c(
  "Happy" = "<span style='color:#800020;'>Happy</span>",
  "Sad" = "<span style='color:#000080;'>Sad</span>",
  "Silence" = "<span style='color:#B8860B;'>Control</span>"
)

ggplot(final_data, aes(x = ORDER.s, y = predicted_prob, color = Music)) +
  # Predicted probabilities as dashed lines
  geom_line(linetype= "dotted", size = 0.85) +
  # Facet by Music condition with custom labels
  facet_wrap(~Music, labeller = labeller(Music = c(
    "Happy" = "Happy",
    "Sad" = "Sad",
    "Silence" = "Control"
  ))) +
  # Facet by Music condition with custom labels
  facet_wrap(~Music, labeller = labeller(Music = facet_labels)) +
  # Adjust titles, labels, and y-axis limits
  labs(
    title = "Smoothed Predicted Probabilities Across Music Conditions",
    x = "Order of (scaled) Time (in 100ms bins)",
    y = "Predicted Probability",
    color = "Predicted Probabilities" # Title for Music legend
  ) +
  scale_y_continuous(limits = c(0.0, 0.04)) + # Adjusted y-axis limits
  # Theme adjustments
  theme_minimal() +
  theme(
    strip.text = element_markdown(size = 12), # Use HTML-style facet labels
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  ) +
  # Custom colors for lines
  scale_color_manual(values = c("Happy" = "#800020", "Sad" = "#000080", "Silence" = "#B8860B"))

#### D. ANT1-2 ANALYSES --  ####
## Di. IVS & Ranodm Effects 
#1. Congruency = Congruent, Incongruent, Neutral
#2. Music_Cx = Happy, Sad, Silence
#3. CUE-TYPE = Double, No, Center, Valid, Invalid
#4. ANT_Cx = ANT_1, ANT_2
#5. Congruency_Cue = #1 x #3
#6. P = Participant

#7. BlinkLogistic = Blinking Probability (1,0)

##D.ii. load in & format data ## 
d <- read.csv("8-25_Database.csv") 

# Participants
d$P <- as.factor(d$P)

# DVs 
d$BlinkLogistic <- as.factor(d$BlinkLogistic)
d$DV <- as.factor(d$DV)

#ANT CUE/CONGRUENCY VARS
d$Congruency <- as.factor(d$Congruency)
d$Congruency_Cue <- as.factor(d$Congruency_Cue)
d$CUE_TYPE <- as.factor(d$CUE_TYPE)
d$EC_IncongruentNeutral <- as.factor(d$EC_IncongruentNeutral)
d$Orienting_InvalidCenter <- as.factor(d$Orienting_InvalidCenter)
d$Orienting_InvalidValid <- as.factor(d$Orienting_InvalidValid)
d$Orienting_ValidCenter <- as.factor(d$Orienting_ValidCenter)
d$EC_IncongruentCongruent <- as.factor(d$EC_IncongruentCongruent)

d$Music_Cx <- factor(d$Music_Cx, levels = c("SILENCE", "HAPPY", "SAD"))

d$Music_NoMusic <- as.factor(d$Music_NoMusic)

#ANT Phase VARS
d$TRIAL_CATwoRT <- as.factor(d$TRIAL_CATwoRT)
d$CUE_TARGET <- d$CUE_TARGET

#ANT Administration VARS
d$ANT_Cx <- as.factor(d$ANT_Cx)

##D.iii. Filter data for CUE WOIs ## 
filtered_data.2 <- d %>%
  filter(TimeStamp >= 2400 & TimeStamp <= 2700)

filtered_data.3 <- d %>%
  filter(TimeStamp >= 2750 & TimeStamp <= 3050)

# Filter the data EC WOIs 
filtered_data.4 <- d %>%
  filter(TimeStamp >= 2900 & TimeStamp <= 3200)

filtered_data.5 <- d %>%
  filter(TimeStamp >= 3350 & TimeStamp <= 4000)

## D.iv. filtered data by cue comparison

## Window 2 (2400-2700ms Cue Onset) 
filtered_data.2a <- filter(filtered_data.2, ALERTING_ORIENTING == "ALERTING")
filtered_data.2vcc <- filter(filtered_data.2, Orienting_ValidCenter == "Valid_Center")

## Window 3 (2750-3050 Pre-Trial Onset) 
filtered_data.3a <- filter(filtered_data.3, ALERTING_ORIENTING == "ALERTING")
filtered_data.3vcc <- filter(filtered_data.3, Orienting_ValidCenter == "Valid_Center")

## Window 4 (2900-3200 Flanker Trial Onset) 
filtered_data.4icn <- filter(filtered_data.4, EC_IncongruentNeutral == "Endogenous_EC")

## Window 5 (3350-4000 RT Window) 
filtered_data.5icn <- filter(filtered_data.5, EC_IncongruentNeutral == "Endogenous_EC")

############# D. ALERTING - No Cue Double Cue #########

########## D1a. ALERTING - Cue Onset Window (2400-2700ms ) No Cue Double Cue #########
alerting.w2.a.base1a <- glmer(BlinkLogistic ~ (ANT_Cx * CUE_TYPE * Music_Cx) + 
                          (1 | P), 
                        data = filtered_data.2a, 
                        family = "binomial", 
                        control = glmerControl(optimizer = "bobyqa"))
#fit

alerting.w2.a1 <- glmer(BlinkLogistic ~ (ANT_Cx * CUE_TYPE * Music_Cx) + 
                           (1 + ANT_Cx | P), 
                         data = filtered_data.2a, 
                         family = "binomial", 
                         control = glmerControl(optimizer = "bobyqa"))
#fit

alerting.w2.a2 <- glmer(BlinkLogistic ~ (ANT_Cx * CUE_TYPE * Music_Cx) + 
                          (1 + CUE_TYPE | P), 
                        data = filtered_data.2a, 
                        family = "binomial", 
                        control = glmerControl(optimizer = "bobyqa"))
#dnc

alerting.w2.a3 <- glmer(BlinkLogistic ~ (ANT_Cx * CUE_TYPE * Music_Cx) + 
                          (1 + Music_Cx | P), 
                        data = filtered_data.2a, 
                        family = "binomial", 
                        control = glmerControl(optimizer = "bobyqa"))
#dnc

alerting.w2.a4 <- glmer(BlinkLogistic ~ (ANT_Cx * CUE_TYPE * Music_Cx) + 
                          (1 + ANT_Cx + CUE_TYPE | P), 
                        data = filtered_data.2a, 
                        family = "binomial", 
                        control = glmerControl(optimizer = "bobyqa"))
#fit

anova(alerting.w2.a4, alerting.w2.a1, alerting.w2.a2, alerting.w2.a3, alerting.w2.a.base)
#4 is best 

summary(alerting.w2.a4)

##r2 values
alerting.w2.a.r2 <- performance::r2(alerting.w2.a4)
print(alerting.w2.a.r2) # % of conditional variance 

#95% CIs
ci.alerting.model <- confint(alerting.w2.a4, parm = "beta_", method = "Wald")

# Exponentiate the CIs
ci.or.alerting <- exp(ci.alerting.model)

#Coefficients (log-odds)
fixed_effects <- fixef(alerting.w2.a4)

# Calculate odds ratios (exponentiate the coefficients)
odds_ratios <- exp(fixed_effects)

# Print odds ratios and confidence intervals
data.frame(Odds_Ratios = odds_ratios, 
           CI_Lower = ci.or.alerting[,1], 
           CI_Upper = ci.or.alerting[,2])

#comparisons ~re-leveling #
filtered_data.2a$Music_Cxsadhsi <- factor(filtered_data.2a$Music_Cx, levels = c("SAD", "HAPPY", "SILENCE"))

alerting.w2.a4b <- glmer(BlinkLogistic ~ (ANT_Cx * CUE_TYPE * Music_Cxsadhsi) + 
                          (1 + ANT_Cx + CUE_TYPE | P), 
                        data = filtered_data.2a, 
                        family = "binomial", 
                        control = glmerControl(optimizer = "bobyqa"))
#fit
summary(alerting.w2.a4b)

#95% CIs
ci.alerting.model <- confint(alerting.w2.a4b, parm = "beta_", method = "Wald")

# Exponentiate the CIs
ci.or.alerting <- exp(ci.alerting.model)

#Coefficients (log-odds)
fixed_effects <- fixef(alerting.w2.a4b)

# Calculate odds ratios (exponentiate the coefficients)
odds_ratios <- exp(fixed_effects)

# Print odds ratios and confidence intervals
data.frame(Odds_Ratios = odds_ratios, 
           CI_Lower = ci.or.alerting[,1], 
           CI_Upper = ci.or.alerting[,2])

########## D1b. ALERTING - Post Cue Window (2750-3050) - No Cue Double Cue ######### 
alerting.w3.a.base <- glmer(BlinkLogistic ~ (ANT_Cx * CUE_TYPE * Music_Cx) + 
                          (1 + ANT_Cx | P), 
                        data = filtered_data.3a, 
                        family = "binomial", 
                        control = glmerControl(optimizer = "bobyqa"))
#fit

alerting.w3.a.base2 <- glmer(BlinkLogistic ~ (ANT_Cx * CUE_TYPE * Music_Cx) + 
                              (1 + ANT_Cx | P)
                            + (1 | CUE_TYPE), 
                            data = filtered_data.3a, 
                            family = "binomial", 
                            control = glmerControl(optimizer = "bobyqa"))
#singular

alerting.w3.a1 <- glmer(BlinkLogistic ~ (ANT_Cx * CUE_TYPE * Music_Cx) + 
                          (1 + ANT_Cx | P), 
                        data = filtered_data.3a, 
                        family = "binomial", 
                        control = glmerControl(optimizer = "bobyqa"))
#fit

alerting.w3.a2 <- glmer(BlinkLogistic ~ (ANT_Cx * CUE_TYPE * Music_Cx) + 
                          (1 + CUE_TYPE | P), 
                        data = filtered_data.3a, 
                        family = "binomial", 
                        control = glmerControl(optimizer = "bobyqa"))
#fit

alerting.w3.a3 <- glmer(BlinkLogistic ~ (ANT_Cx * CUE_TYPE * Music_Cx) + 
                          (1 + Music_Cx | P), 
                        data = filtered_data.3a, 
                        family = "binomial", 
                        control = glmerControl(optimizer = "bobyqa"))
#fit

alerting.w3.a4 <- glmer(BlinkLogistic ~ (ANT_Cx * CUE_TYPE * Music_Cx) + 
                          (1 + ANT_Cx + CUE_TYPE | P), 
                        data = filtered_data.3a, 
                        family = "binomial", 
                        control = glmerControl(optimizer = "bobyqa"))
#fit

alerting.w3.a5 <- glmer(BlinkLogistic ~ (ANT_Cx * CUE_TYPE * Music_Cx) + 
                          (1 + ANT_Cx + CUE_TYPE + Music_Cx | P), 
                        data = filtered_data.3a, 
                        family = "binomial", 
                        control = glmerControl(optimizer = "bobyqa"))
#singular

anova(alerting.w3.a.base, alerting.w3.a1, alerting.w3.a2, 
      alerting.w3.a3, alerting.w3.a4)
#2 is best 

summary(alerting.w3.a2)

##r2 values
alerting.w3.a.r2 <- performance::r2(alerting.w3.a2)
print(alerting.w3.a.r2) # % of conditional variance 

#95% CIs
ci.alerting.model <- confint(alerting.w3.a2, parm = "beta_", method = "Wald")

# Exponentiate the CIs
ci.or.alerting <- exp(ci.alerting.model)

#Coefficients (log-odds)
fixed_effects <- fixef(alerting.w3.a2)

# Calculate odds ratios (exponentiate the coefficients)
odds_ratios <- exp(fixed_effects)

# Print odds ratios and confidence intervals
data.frame(Odds_Ratios = odds_ratios, 
           CI_Lower = ci.or.alerting[,1], 
           CI_Upper = ci.or.alerting[,2])

#comparisons ~re-leveling #
filtered_data.3a$Music_Cxsadhsi <- factor(filtered_data.3a$Music_Cx, levels = c("SAD", "HAPPY", "SILENCE"))

alerting.w3.a2a <- glmer(BlinkLogistic ~ (ANT_Cx * CUE_TYPE * Music_Cxsadhsi) + 
                           (1 + ANT_Cx + CUE_TYPE | P), 
                         data = filtered_data.3a, 
                         family = "binomial", 
                         control = glmerControl(optimizer = "bobyqa"))
#fit

summary(alerting.w3.a2a)

#95% CIs
ci.alerting.model <- confint(alerting.w3.a2a, parm = "beta_", method = "Wald")

# Exponentiate the CIs
ci.or.alerting <- exp(ci.alerting.model)

#Coefficients (log-odds)
fixed_effects <- fixef(alerting.w3.a2a)

# Calculate odds ratios (exponentiate the coefficients)
odds_ratios <- exp(fixed_effects)

# Print odds ratios and confidence intervals
data.frame(Odds_Ratios = odds_ratios, 
           CI_Lower = ci.or.alerting[,1], 
           CI_Upper = ci.or.alerting[,2])

########## D2: ORIENTING - Center v. Valid Cue ##########
#comparisons ~re-leveling #
filtered_data.2vcc$CUE_TYPE <- factor(filtered_data.2vcc$CUE_TYPE, levels = c("Valid", "Center"))

########## D2a: ORIENTING - Cue Onset Window (2400-2700ms Cue Onset) Center v. Valid ##########

orienting.w2o.vcc.base <- glmer(BlinkLogistic ~ (ANT_Cx * CUE_TYPE * Music_Cx) + 
                                  (1 | P), 
                                data = filtered_data.2vcc, 
                                family = "binomial", 
                                control = glmerControl(optimizer = "bobyqa"))
#fit 

orienting.w2o.vcc.base2 <- glmer(BlinkLogistic ~ (ANT_Cx * CUE_TYPE * Music_Cx) + 
                                   (1 | P)
                                 + (1 | CUE_TYPE), 
                                 data = filtered_data.2vcc, 
                                 family = "binomial", 
                                 control = glmerControl(optimizer = "bobyqa"))
#singular

orienting.w2o.vcc1 <- glmer(BlinkLogistic ~ (ANT_Cx * CUE_TYPE * Music_Cx) + 
                              (1 + ANT_Cx | P), 
                            data = filtered_data.2vcc, 
                            family = "binomial", 
                            control = glmerControl(optimizer = "bobyqa"))
#fit 
orienting.w2o.vcc2 <- glmer(BlinkLogistic ~ (ANT_Cx * CUE_TYPE * Music_Cx) + 
                              (1 + CUE_TYPE | P), 
                            data = filtered_data.2vcc, 
                            family = "binomial", 
                            control = glmerControl(optimizer = "bobyqa"))
#singular
orienting.w2o.vcc3 <- glmer(BlinkLogistic ~ (ANT_Cx * CUE_TYPE * Music_Cx) + 
                              (1 + Music_Cx | P), 
                            data = filtered_data.2vcc, 
                            family = "binomial", 
                            control = glmerControl(optimizer = "bobyqa"))
#fit
anova(orienting.w2o.vcc.base, 
      orienting.w2o.vcc1, orienting.w2o.vcc3)
#1 is best 

summary(orienting.w2o.vcc1)
o.w2.vcc <- performance::r2(orienting.w2o.vcc1)
print(o.w2.vcc) 

#95% CIs
ci.orienting.w2.model <- confint(orienting.w2o.vcc1, parm = "beta_", method = "Wald")

# Exponentiate the CIs
ci.or.orienting <- exp(ci.orienting.w2.model)

#Coefficients (log-odds)
fixed_effects <- fixef(orienting.w2o.vcc1)

# Calculate odds ratios (exponentiate the coefficients)
odds_ratios <- exp(fixed_effects)

# Print odds ratios and confidence intervals
data.frame(Odds_Ratios = odds_ratios, 
           CI_Lower = ci.or.orienting[,1], 
           CI_Upper = ci.or.orienting[,2])

#comparisons ~re-leveling #
filtered_data.2vcc$Music_Cxsadhsi <- factor(filtered_data.2vcc$Music_Cx, levels = c("SAD", "HAPPY", "SILENCE"))

orienting.w2o.vcc1a <- glmer(BlinkLogistic ~ (ANT_Cx * CUE_TYPE * Music_Cxsadhsi) + 
                              (1 + ANT_Cx | P), 
                            data = filtered_data.2vcc, 
                            family = "binomial", 
                            control = glmerControl(optimizer = "bobyqa"))
#fit 
summary(orienting.w2o.vcc1a)

#95% CIs
ci.orienting.w2.model <- confint(orienting.w2o.vcc1a, parm = "beta_", method = "Wald")

# Exponentiate the CIs
ci.or.orienting <- exp(ci.orienting.w2.model)

#Coefficients (log-odds)
fixed_effects <- fixef(orienting.w2o.vcc1a)

# Calculate odds ratios (exponentiate the coefficients)
odds_ratios <- exp(fixed_effects)

# Print odds ratios and confidence intervals
data.frame(Odds_Ratios = odds_ratios, 
           CI_Lower = ci.or.orienting[,1], 
           CI_Upper = ci.or.orienting[,2])

########## D2b: ORIENTING - Post Cue Window (2750-3050 Pre-Trial Onset) Center v. Valid ##########
#comparisons ~re-leveling #
filtered_data.3vcc$CUE_TYPE <- factor(filtered_data.3vcc$CUE_TYPE, levels = c("Valid", "Center"))

orienting.w3o.vcc.base <- glmer(BlinkLogistic ~ (ANT_Cx * CUE_TYPE * Music_Cx) + 
                              (1 | P), 
                            data = filtered_data.3vcc, 
                            family = "binomial", 
                            control = glmerControl(optimizer = "bobyqa"))
#fit 

orienting.w3o.vcc.base2 <- glmer(BlinkLogistic ~ (ANT_Cx * CUE_TYPE * Music_Cx) + 
                                  (1 | P)
                                 + (1 | CUE_TYPE), 
                                data = filtered_data.3vcc, 
                                family = "binomial", 
                                control = glmerControl(optimizer = "bobyqa"))
#singular 

orienting.w3o.vcc1 <- glmer(BlinkLogistic ~ (ANT_Cx * CUE_TYPE * Music_Cx) + 
                              (1 + ANT_Cx | P), 
                            data = filtered_data.3vcc, 
                            family = "binomial", 
                            control = glmerControl(optimizer = "bobyqa"))
#fit

orienting.w3o.vcc2 <- glmer(BlinkLogistic ~ (ANT_Cx * CUE_TYPE * Music_Cx) + 
                              (1 + CUE_TYPE | P), 
                            data = filtered_data.3vcc, 
                            family = "binomial", 
                            control = glmerControl(optimizer = "bobyqa"))
#fit 

orienting.w3o.vcc3 <- glmer(BlinkLogistic ~ (ANT_Cx * CUE_TYPE * Music_Cx) + 
                              (1 + Music_Cx | P), 
                            data = filtered_data.3vcc, 
                            family = "binomial", 
                            control = glmerControl(optimizer = "bobyqa"))
#fit 

anova(orienting.w3o.vcc3, orienting.w3o.vcc2, orienting.w3o.vcc1, orienting.w3o.vcc.base)
#1 & 2 were best 

orienting.w3o.vcc4 <- glmer(BlinkLogistic ~ (ANT_Cx * CUE_TYPE * Music_Cx) + 
                              (1 + ANT_Cx + CUE_TYPE | P), 
                            data = filtered_data.3vcc, 
                            family = "binomial", 
                            control = glmerControl(optimizer = "bobyqa"))
#fit 

anova(orienting.w3o.vcc4, orienting.w3o.vcc2, orienting.w3o.vcc1)
#one is still better 
summary(orienting.w3o.vcc1)
o.w3.vcc1 <- performance::r2(orienting.w3o.vcc1)
print(o.w3.vcc1)

summary(filtered_data.3vcc$BlinkLogistic)

#95% CIs
ci.orienting.w3.model <- confint(orienting.w3o.vcc1, parm = "beta_", method = "Wald")

# Exponentiate the CIs
ci.or.orienting <- exp(ci.orienting.w3.model)

#Coefficients (log-odds)
fixed_effects <- fixef(orienting.w3o.vcc1)

# Calculate odds ratios (exponentiate the coefficients)
odds_ratios <- exp(fixed_effects)

# Print odds ratios and confidence intervals
data.frame(Odds_Ratios = odds_ratios, 
           CI_Lower = ci.or.orienting[,1], 
           CI_Upper = ci.or.orienting[,2])

#comparisons ~re-leveling #
filtered_data.3vcc$Music_Cxsadhsi <- factor(filtered_data.3vcc$Music_Cx, levels = c("SAD", "HAPPY", "SILENCE"))

orienting.w3o.vcc1a <- glmer(BlinkLogistic ~ (ANT_Cx * CUE_TYPE * Music_Cxsadhsi) + 
                              (1 + ANT_Cx | P), 
                            data = filtered_data.3vcc, 
                            family = "binomial", 
                            control = glmerControl(optimizer = "bobyqa"))
#fit
summary(orienting.w3o.vcc1a)

#95% CIs
ci.orienting.w3.model <- confint(orienting.w3o.vcc1a, parm = "beta_", method = "Wald")

# Exponentiate the CIs
ci.or.orienting <- exp(ci.orienting.w3.model)

#Coefficients (log-odds)
fixed_effects <- fixef(orienting.w3o.vcc1a)

# Calculate odds ratios (exponentiate the coefficients)
odds_ratios <- exp(fixed_effects)

# Print odds ratios and confidence intervals
data.frame(Odds_Ratios = odds_ratios, 
           CI_Lower = ci.or.orienting[,1], 
           CI_Upper = ci.or.orienting[,2])

############# D3. EXECUTIVE CONTROL - Flanker Onset Window (2900ms - 3200ms Flanker Trial Onset) Incongruent - Neutral #########
ec.w4icn.base <- glmer(BlinkLogistic ~ (ANT_Cx * Congruency * Music_Cx) + 
                       (1 | P), 
                     data = filtered_data.4icn, 
                     family = "binomial", 
                     control = glmerControl(optimizer = "bobyqa"))
#fit

ec.w4icn.base2 <- glmer(BlinkLogistic ~ (ANT_Cx * Congruency * Music_Cx) + 
                         (1 | P)
                       + (1 | Congruency), 
                       data = filtered_data.4icn, 
                       family = "binomial", 
                       control = glmerControl(optimizer = "bobyqa"))
#singular

ec.w4icn.m1 <- glmer(BlinkLogistic ~ (ANT_Cx * Congruency * Music_Cx) + 
                       (1 + ANT_Cx | P), 
                     data = filtered_data.4icn, 
                     family = "binomial", 
                     control = glmerControl(optimizer = "bobyqa"))
#fit

ec.w4icn.m2 <- glmer(BlinkLogistic ~ (ANT_Cx * Congruency * Music_Cx) + 
                       (1 + Congruency | P), 
                     data = filtered_data.4icn, 
                     family = "binomial", 
                     control = glmerControl(optimizer = "bobyqa"))
#singular 

ec.w4icn.m3 <- glmer(BlinkLogistic ~ (ANT_Cx * Congruency * Music_Cx) + 
                       (1 + Music_Cx | P), 
                     data = filtered_data.4icn, 
                     family = "binomial", 
                     control = glmerControl(optimizer = "bobyqa"))
#fit

ec.w4icn.m4 <- glmer(BlinkLogistic ~ (ANT_Cx * Congruency * Music_Cx) + 
                       (1 + ANT_Cx + Music_Cx | P), 
                     data = filtered_data.4icn, 
                     family = "binomial", 
                     control = glmerControl(optimizer = "bobyqa"))
#singular

anova(ec.w4icn.m3, ec.w4icn.m1, ec.w4icn.base) #1 is best 
summary(ec.w4icn.m1)
ec.w4.icn1 <- performance::r2(ec.w4icn.m1)
print(ec.w4.icn1) #xxx % of conditional variance 

summary(filtered_data.4icn$BlinkLogistic)

#95% CIs
ci.ec.w4.model <- confint(ec.w4icn.m1, parm = "beta_", method = "Wald")

# Exponentiate the CIs
ci.or.ec <- exp(ci.ec.w4.model)

#Coefficients (log-odds)
fixed_effects <- fixef(ec.w4icn.m1)

# Calculate odds ratios (exponentiate the coefficients)
odds_ratios <- exp(fixed_effects)

# Print odds ratios and confidence intervals
data.frame(Odds_Ratios = odds_ratios, 
           CI_Lower = ci.or.ec[,1], 
           CI_Upper = ci.or.ec[,2])

#comparisons ~re-leveling #
filtered_data.4icn$Music_Cxsadhsi <- factor(filtered_data.4icn$Music_Cx, levels = c("SAD", "HAPPY", "SILENCE"))

ec.w4icn.m1a <- glmer(BlinkLogistic ~ (ANT_Cx * Congruency * Music_Cxsadhsi) + 
                        (1 + ANT_Cx | P), 
                      data = filtered_data.4icn, 
                      family = "binomial", 
                      control = glmerControl(optimizer = "bobyqa"))
#fit

summary(ec.w4icn.m1a)

#95% CIs
ci.ec.w4a.model <- confint(ec.w4icn.m1a, parm = "beta_", method = "Wald")

# Exponentiate the CIs
ci.or.ec <- exp(ci.ec.w4a.model)

#Coefficients (log-odds)
fixed_effects <- fixef(ec.w4icn.m1a)

#odds ratios (exponentiate the coefficients)
odds_ratios <- exp(fixed_effects)

data.frame(Odds_Ratios = odds_ratios, 
           CI_Lower = ci.or.ec[,1], 
           CI_Upper = ci.or.ec[,2])

############# D4a. EXECUTIVE CONTROL - RT Window (3350-4000ms) Incongruent - Neutral #########

############# D4b. EXECUTIVE CONTROL (BLINKING) - RT Window (3350-4000ms) Incongruent - Neutral #########

ec.w5icn.base <- glmer(BlinkLogistic ~ (ANT_Cx * Congruency * Music_Cx) + 
                       (1 | P), 
                     data = filtered_data.5icn, 
                     family = "binomial", 
                     control = glmerControl(optimizer = "bobyqa"))
#fit

ec.w5icn.base2 <- glmer(BlinkLogistic ~ (ANT_Cx * Congruency * Music_Cx) + 
                         (1 | P)
                       + (1 | Congruency), 
                       data = filtered_data.5icn, 
                       family = "binomial", 
                       control = glmerControl(optimizer = "bobyqa"))
#singular

ec.w5icn.m1 <- glmer(BlinkLogistic ~ (ANT_Cx * Congruency * Music_Cx) + 
                       (1 + ANT_Cx | P), 
                     data = filtered_data.5icn, 
                     family = "binomial", 
                     control = glmerControl(optimizer = "bobyqa"))
#fit

ec.w5icn.m2 <- glmer(BlinkLogistic ~ (ANT_Cx * Congruency * Music_Cx) + 
                       (1 + Congruency | P), 
                     data = filtered_data.5icn, 
                     family = "binomial", 
                     control = glmerControl(optimizer = "bobyqa"))
#fit 

ec.w5icn.m3 <- glmer(BlinkLogistic ~ (ANT_Cx * Congruency * Music_Cx) + 
                       (1 + Music_Cx | P), 
                     data = filtered_data.5icn, 
                     family = "binomial", 
                     control = glmerControl(optimizer = "bobyqa"))
#singular

ec.w5icn.m4 <- glmer(BlinkLogistic ~ (ANT_Cx * Congruency * Music_Cx) + 
                       (1 + ANT_Cx + Music_Cx | P), 
                     data = filtered_data.5icn, 
                     family = "binomial", 
                     control = glmerControl(optimizer = "bobyqa"))
#singular

anova(ec.w5icn.base, ec.w5icn.m1, ec.w5icn.m2) #m1 is best 
summary(ec.w5icn.m1)
ec.w5.icn1 <- performance::r2(ec.w5icn.m1)
print(ec.w5.icn1) #xxx % of conditional variance 

summary(filtered_data.5icn$BlinkLogistic)

#95% CIs
ci.ec.w5.model <- confint(ec.w5icn.m1, parm = "beta_", method = "Wald")

# Exponentiate the CIs
ci.or.ec <- exp(ci.ec.w5.model)

#Coefficients (log-odds)
fixed_effects <- fixef(ec.w5icn.m1)

#odds ratios (exponentiate the coefficients)
odds_ratios <- exp(fixed_effects)

data.frame(Odds_Ratios = odds_ratios, 
           CI_Lower = ci.or.ec[,1], 
           CI_Upper = ci.or.ec[,2])

#comparisons ~re-leveling #
filtered_data.5icn$Music_Cxsadhsi <- factor(filtered_data.5icn$Music_Cx, levels = c("SAD", "HAPPY", "SILENCE"))

ec.w5icn.m1a <- glmer(BlinkLogistic ~ (ANT_Cx * Congruency * Music_Cxsadhsi) + 
                       (1 + ANT_Cx | P), 
                     data = filtered_data.5icn, 
                     family = "binomial", 
                     control = glmerControl(optimizer = "bobyqa"))
#fit

summary(ec.w5icn.m1a)

#95% CIs
ci.ec.w5a.model <- confint(ec.w5icn.m1a, parm = "beta_", method = "Wald")

# Exponentiate the CIs
ci.or.ec <- exp(ci.ec.w5a.model)

#Coefficients (log-odds)
fixed_effects <- fixef(ec.w5icn.m1a)

#odds ratios (exponentiate the coefficients)
odds_ratios <- exp(fixed_effects)

data.frame(Odds_Ratios = odds_ratios, 
           CI_Lower = ci.or.ec[,1], 
           CI_Upper = ci.or.ec[,2])

############# D4c. EXECUTIVE CONTROL (RT & BLINKING) - RT Window (3350-4000ms) Incongruent - Neutral #########

## Filter data for BlinkStart and RT by ICN Data ##
blink_rt_icn <- filtered_data.5icn %>% 
  filter(DV %in% c("BlinkStart", "RT"))

# Convert DV to a binary variable for whether it is a BlinkStart
blink_rt_icn$BlinkStart <- ifelse(blink_rt_icn$DV == "BlinkStart", 1, 0)

## TimeStamp (RT) from Blink occurrence by EC - Incongruent - Neutral by Music ##
blink_rt.w5_icon_base <- lmer(TimeStamp ~ (BlinkStart * ANT_Cx * Congruency * Music_Cx) 
           + (1 | P), data = blink_rt_icn)
#fit 
blink_rt.w5_icon_base2 <- lmer(TimeStamp ~ (BlinkStart * ANT_Cx * Congruency * Music_Cx) 
                              + (1 | P)
                              + (1 | Congruency), data = blink_rt_icn)
#dnc 

blink_rt.w5_icon_m1 <- lmer(TimeStamp ~ (BlinkStart * ANT_Cx * Congruency * Music_Cx) 
           + (1 + Congruency | P), data = blink_rt_icn)
#fit
blink_rt.w5_icon_m2 <- lmer(TimeStamp ~ (BlinkStart * ANT_Cx * Congruency * Music_Cx) 
            + (1 + Music_Cx | P), data = blink_rt_icn)
#fit

blink_rt.w5_icon_m3 <- lmer(TimeStamp ~ (BlinkStart * ANT_Cx * Congruency * Music_Cx) 
                            + (1 + ANT_Cx | P), data = blink_rt_icn)
#fit 

anova(blink_rt.w5_icon_base, blink_rt.w5_icon_m1, blink_rt.w5_icon_m2, blink_rt.w5_icon_m3) 
#m1 & 3 improved fit 

blink_rt.w5_icon_m4 <- lmer(TimeStamp ~ (BlinkStart * ANT_Cx * Congruency * Music_Cx) 
                            + (1 + ANT_Cx + Congruency | P), data = blink_rt_icn)
#fit 

anova(blink_rt.w5_icon_m1, blink_rt.w5_icon_m3, blink_rt.w5_icon_m4) #4 is best 
summary(blink_rt.w5_icon_m4)
ec.w5.icn.rt <- performance::r2(blink_rt.w5_icon_m4)
print(ec.w5.icn.rt) #xxx % of conditional variance 

mean(blink_rt_icn$TimeStamp)
sd(blink_rt_icn$TimeStamp)

summary(blink_rt_icn$BlinkStart)

#comparisons ~re-leveling #
blink_rt_icn$Music_Cxsadhsi <- factor(blink_rt_icn$Music_Cx, levels = c("SAD", "HAPPY", "SILENCE"))

blink_rt.w5_icon_m5a <- lmer(TimeStamp ~ (BlinkStart * ANT_Cx * Congruency * Music_Cxsadhsi) 
                            + (1 + ANT_Cx + Congruency | P), data = blink_rt_icn)

summary(blink_rt.w5_icon_m5a)

############# D5. PLOT ################

#ANT_2 and ANT_1 subsets
ant_t2_data <- blink_rt_icn %>% filter(ANT_Cx == "ANT_2")
ant_t1_data <- blink_rt_icn %>% filter(ANT_Cx == "ANT_1")

#ANT_1 on top
ggplot() +
  #ANT-T2
  geom_boxplot(data = ant_t2_data, aes(x = factor(BlinkStart), y = TimeStamp, fill = factor(BlinkStart))) +
  # Overlay ANT-T1
  geom_boxplot(data = ant_t1_data, aes(x = factor(BlinkStart), y = TimeStamp, color = factor(BlinkStart)),
               fill = NA, linetype = "dotted", size = 0.6, alpha = 0.9, outlier.shape = NA) +
  facet_wrap(~ Congruency, nrow = 1) +  # Horizontal facets
  scale_fill_manual(values = c("0" = "skyblue", "1" = "orange")) +  # Fill colors for ANT_2
  scale_color_manual(values = c("0" = "darkblue", "1" = "darkorange")) +  # Darker colors for ANT_1
  labs(
    x = "Blink Start (0 = No, 1 = Yes)",
    y = "Reaction Time (ms)",
    title = "ANT T2 Reaction Times by Blink Start (solid) with ANT T1 Outline Overlay (dotted), by Congruency"
  ) +
  theme_minimal()

######### E. EXECUTIVE CONTROL (RT & BLINKING - Time Series) - RT Window (3350-4000ms) Incongruent - Neutral ###########
filtered_data <- filtered_data.5icn %>%
  filter(DV %in% c("RT", "BlinkStart"))

# DV as factor
filtered_data$DV <- factor(filtered_data$DV, levels = c("RT", "BlinkStart"))

#Model Development
ec_tsbase <- lmer(TimeStamp ~ DV * ANT_Cx * Music_Cx * Congruency + (1 | P), 
                 data = filtered_data, 
                 REML = FALSE)
#fit

ec_tsbase2 <- lmer(TimeStamp ~ DV * ANT_Cx * Music_Cx * Congruency + (1 | P) + (1 | Congruency), 
                   data = filtered_data, 
                   REML = FALSE)
#singular

ec_ts1 <- lmer(TimeStamp ~ DV * ANT_Cx * Music_Cx * Congruency + (1 + ANT_Cx | P), 
               data = filtered_data, 
               REML = FALSE)
#fit

ec_ts2 <- lmer(TimeStamp ~ DV * ANT_Cx * Music_Cx * Congruency + (1 + Music_Cx | P), 
               data = filtered_data, 
               REML = FALSE)
#did not converge

ec_ts3 <- lmer(TimeStamp ~ DV * ANT_Cx * Music_Cx * Congruency + (1 + Congruency | P), 
               data = filtered_data, 
               REML = FALSE)
#fit

anova(ec_tsbase, ec_ts1, ec_ts3) #1 significantly improved fit

ec_ts5 <- lmer(TimeStamp ~ DV * ANT_Cx * Music_Cx * Congruency + (1 + ANT_Cx + Congruency | P), 
               data = filtered_data, 
               REML = FALSE)
#fit

anova(ec_ts1, ec_ts5) #5 is best  
summary(ec_ts5)
ec.w5.icn.rt.ts <- performance::r2(ec_ts5)
print(ec.w5.icn.rt.ts) #xxx % of conditional variance 

#comparisons ~re-leveling #
filtered_data$Music_Cxsadhsi <- factor(filtered_data$Music_Cx, levels = c("SAD", "HAPPY", "SILENCE"))

ec_ts5a <- lmer(TimeStamp ~ DV * ANT_Cx * Music_Cxsadhsi * Congruency + (1 + ANT_Cx + Congruency | P), 
               data = filtered_data, 
               REML = FALSE)
#fit
summary(ec_ts5a)

mean(filtered_data$TimeStamp)
sd(filtered_data$TimeStamp)

summary(filtered_data$DV)

