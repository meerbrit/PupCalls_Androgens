#### Analysis script for Flutamide study: REP (REPEAT) calls ###################
############## BWalkenhorst 2024 ###############################################

#### SETUP ####

# clear everything and use garbage collector at start
rm(list = ls(all.names = TRUE)) #includes hidden objects.
gc() 

# load all necessary libraries
library(readxl) 
library(ggplot2) 
library(car)
library(tidyverse)
library(tidybayes) 
library(brms)
library(bayestestR) #e.g. diagnostic_posterior
library(bayesplot)
library(ggokabeito) # colour palette
library(emmeans) # emtrends
library(extrafont)# use font_import() on first use

set.seed(23)
REP_data <- readRDS('data/REP_data.rds')

#Generic weakly informative prior for call length and interval duration
priors <- c(set_prior("normal(0,1)", class = "Intercept"), set_prior("normal(0,1)", class='b'))

# Custom ggplot theme 
theme_clean <- function() {
  theme_minimal(base_family='Calibri') +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = rel(2), hjust = 0.5),
      axis.title = element_text(face = "bold", size = rel(2)),
      axis.text = element_text(face = "bold", size= rel(1.25)),
      strip.text = element_text(face = "bold", size = rel(2), color='white'),
      strip.background = element_rect(fill = "grey80", color = NA),
      legend.title = element_text(face = "bold", size = rel(2)),
      legend.text = element_text(face = 'italic', size = rel(1.5)))
}

#### FUNCTIONS ####
# Process model and calculate natural scale results
process_model_results <- function(model, data, log_scale = TRUE) {
  # Map predictors to raw data columns for SD calculation
  sd_mapping <- list(
    AGE_z = "REC_AGE_D",
    WEIGHT_z = "WEIGHT_DIFF_PER",
    COMP_NORM_z = "COMP_NORM",
    GS_z = "GROUPSIZE",
    RAIN_z = "MonthlyRainfall"
  )
  
  # Extract fixed effects
  fixed_effects <- fixef(model)
  
  # Ensure fixed effects have the expected structure
  if (!"Estimate" %in% colnames(fixed_effects)) {
    stop("fixef(model) does not have an 'Estimate' column.")
  }
  
  intercept <- fixed_effects["Intercept", "Estimate"]
  effects <- fixed_effects[rownames(fixed_effects) != "Intercept", "Estimate"]
  
  # Handle case where no effects exist
  if (length(effects) == 0) {
    warning("No effects found in the model. Returning intercept only.")
    return(data.frame(
      Effect = "Intercept",
      Effect_Value = round(intercept, 3),
      Percentage_Change = NA,
      Absolute_Change = NA,
      Adjusted_Call = if (log_scale) round(exp(intercept), 4) else round(intercept, 4),
      SD_Unit = "NA"
    ))
  }
  
  # Process based on scale
  if (log_scale) {
    # Log-scale model: transform to natural scale
    intercept_natural <- exp(intercept)
    absolute_changes <- intercept_natural * (exp(effects) - 1)
    percentage_changes <- (exp(effects) - 1) * 100
    adjusted_change <- intercept_natural + absolute_changes
  } else {
    # Natural-scale model: use raw values
    intercept_natural <- intercept
    absolute_changes <- effects
    percentage_changes <- (effects / abs(intercept)) * 100
    adjusted_change <- intercept_natural + absolute_changes
  }
  
  # Calculate SDs only for continuous predictors
  sd_units <- sapply(names(effects), function(var) {
    raw_column <- sd_mapping[[var]]
    if (!is.null(raw_column) && raw_column %in% colnames(data)) {
      round(sd(data[[raw_column]], na.rm = TRUE), 3)
    } else {
      NA  # Return NA for factor levels or interactions
    }
  })
  
  # Combine results
  results <- data.frame(
    Effect = names(effects),
    Effect_Value = round(effects, 3),
    Percentage_Change = round(percentage_changes, 2),
    Absolute_Change = round(absolute_changes, 4),
    Adjusted_Call = round(adjusted_change, 4),
    SD_Unit = ifelse(is.na(sd_units), "NA", sd_units)  # Handle missing SDs
  )
  
  # Add intercept row
  results <- rbind(
    data.frame(
      Effect = "Intercept",
      Effect_Value = round(intercept, 3),
      Percentage_Change = NA,
      Absolute_Change = NA,
      Adjusted_Call = round(intercept_natural, 4),
      SD_Unit = "NA"
    ),
    results
  )
  
  return(results)
}

# get age values for EMMs
get_age_vars <- function(){
  age_vars <- c((30- mean(REP_data$REC_AGE_D))/sd(REP_data$REC_AGE_D), 
                (75- mean(REP_data$REC_AGE_D))/sd(REP_data$REC_AGE_D), 
                (120- mean(REP_data$REC_AGE_D))/sd(REP_data$REC_AGE_D))
  return(age_vars)
}

# Calculate the difference in EMMs in %
calculate_percent_difference <- function(contrasts, emmeans_output, log_scale = TRUE) {
  # Convert contrasts and emmeans_output to data frames
  contrasts_df <- as.data.frame(contrasts)
  emmeans_df <- as.data.frame(emmeans_output)
  
  # Extract group1 and group2 from the contrast column
  contrasts_df <- contrasts_df %>%
    tidyr::separate(contrast, into = c("group1", "group2"), sep = " - ", remove = FALSE)
  
  # Add reference group columns based on the contrast type (SEX or TREATMENT)
  if ("SEX" %in% colnames(contrasts_df)) {
    # For sex-based contrasts
    contrasts_df <- contrasts_df %>%
      left_join(
        emmeans_df %>% rename(ref_group1 = emmean), 
        by = c("SEX" = "SEX", "group1" = "TREATMENT", "AGE_z" = "AGE_z")
      ) %>%
      left_join(
        emmeans_df %>% rename(ref_group2 = emmean), 
        by = c("SEX" = "SEX", "group2" = "TREATMENT", "AGE_z" = "AGE_z")
      )
  } else if ("TREATMENT" %in% colnames(contrasts_df)) {
    # For treatment-based contrasts
    contrasts_df <- contrasts_df %>%
      left_join(
        emmeans_df %>% rename(ref_group1 = emmean), 
        by = c("TREATMENT" = "TREATMENT", "group1" = "SEX", "AGE_z" = "AGE_z")
      ) %>%
      left_join(
        emmeans_df %>% rename(ref_group2 = emmean), 
        by = c("TREATMENT" = "TREATMENT", "group2" = "SEX", "AGE_z" = "AGE_z")
      )
  }
  
  # Calculate percentage differences using reference values
  if (log_scale) {
    # Log-scale: transform back to original scale and calculate percentage differences
    contrasts_df <- contrasts_df %>%
      mutate(
        percent_difference = (exp(estimate) - 1) * 100,
        response_estimate = exp(estimate),
        response_lower = exp(lower.HPD),
        response_upper = exp(upper.HPD)
      )
  } else {
    # Raw-scale: calculate percentage differences relative to ref_group1 or ref_group2
    contrasts_df <- contrasts_df %>%
      mutate(
        percent_difference = (estimate / abs(ref_group1)) * 100,
        response_estimate = estimate,
        response_lower = lower.HPD,
        response_upper = upper.HPD
      )
  }
  
  # Return relevant columns
  return(contrasts_df %>%
           select(
             contrast, group1, group2, estimate, percent_difference,
             response_estimate, response_lower, response_upper, ref_group1, ref_group2
           ))
}

################################################################################
######################## Call length ###########################################
################################################################################

ggplot(REP_data, aes(x = BEG_avg_Len)) + 
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  labs(title = "Histogram of Call length", x = "Call length (/s)")


# 1) Check for non-linearity ####
B_REP_len_TA <- brms::brm(formula = BEG_avg_Len ~ TREATMENT * AGE_z * SEX + (1|LITTER_CODE/ID),
                          data = REP_data, family = lognormal(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                          prior = priors, threads = threading(4), 
                          file="B_REP_len_TA")
B_REP_len_TA <- add_criterion(B_REP_len_TA, c("loo", "loo_R2"), moment_match = TRUE)
summary(B_REP_len_TA)
plot(B_REP_len_TA)
pp_check(B_REP_len_TA, ndraws=100)

B_REP_len_TA2 <- brms::brm(formula = BEG_avg_Len ~  TREATMENT * (AGE_z + I(AGE_z^2)) * SEX + (1|LITTER_CODE/ID),
                          data = REP_data, family = lognormal(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                          prior = priors, threads = threading(4),
                          file="B_REP_len_TA2")
B_REP_len_TA2 <- add_criterion(B_REP_len_TA2, c("loo", "loo_R2"), moment_match = TRUE)

loo(B_REP_len_TA, B_REP_len_TA2)
# elpd_diff se_diff
# B_REP_len_TA2  0.0       0.0   
# B_REP_len_TA  -8.9       5.8   

# 2) Define model for REP LENGTH ####
# B_REP_len <- brms::brm(formula = BEG_avg_Len ~ TREATMENT * (AGE_z + I(AGE_z^2)) * SEX + 
#                          WEIGHT_z + COMP_NORM_z + GS_z + I(GS_z^2) + RAIN_z + (1|LITTER_CODE/ID),
#                        data = REP_data, family = lognormal(link='identity'),
#                        chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
#                        save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
#                        prior = priors, threads = threading(4),
#                        file="B_REP_len")
# B_REP_len <- add_criterion(B_REP_len, c("loo", "loo_R2"), moment_match = TRUE)

B_REP_len <- readRDS('models/B_REP_len.rds')


#### 3 RESULTS: REP LENGTH ####
summary(B_REP_len)

plot(B_REP_len)
pp_check(B_REP_len, ndraws = 100)

describe_posterior(
  B_REP_len,
  effects = "all", #fixed vs all (for random effects)
  component = "all",
  rope_range = rope_range(B_REP_len),  
  test = c("p_direction", "p_significance", "rope"),
  centrality = "all",
  dispersion = TRUE
)

process_model_results(B_REP_len, data = REP_data, log_scale = T)

loo_R2(B_REP_len, moment_match=T) 

bayes_R2(B_REP_len)


performance::variance_decomposition(B_REP_len)

### EMMs at different ages TAS ####
(emm_REP_LEN <- emmeans(B_REP_len,  ~ TREATMENT:SEX | AGE_z, at = list('AGE_z' = get_age_vars())))

pd(emm_REP_LEN)

equivalence_test(emm_REP_LEN, rope = rope_range(B_REP_len))

p_significance(emm_REP_LEN, threshold = rope_range(B_REP_len))

(pairs_within_sex <- # Contrasts between treatments within sex at each age
    contrast(emm_REP_LEN, method = "pairwise", by = c("SEX", 'AGE_z')))

calculate_percent_difference(pairs_within_sex, emm_REP_LEN, log_scale = T)

pd(pairs_within_sex)

equivalence_test(pairs_within_sex, rope = rope_range(B_REP_len))

p_significance(pairs_within_sex, threshold = rope_range(B_REP_len))

(pairs_within_treatment <- # Contrasts between sex within treatment at each age
contrast(emm_REP_LEN, method = "pairwise", by = c("TREATMENT", 'AGE_z')))

calculate_percent_difference(pairs_within_treatment, emm_REP_LEN, log_scale = T)

p_direction(pairs_within_treatment)

equivalence_test(pairs_within_treatment, range = rope_range(B_REP_len))

p_significance(pairs_within_treatment, threshold = rope_range(B_REP_len))

rm(emm_REP_LEN, pairs_within_sex, pairs_within_treatment)

### PLOTS: REP LENGTH ####
# coefficients:
posterior_desc <- describe_posterior(
  B_REP_len,
  effects = "fixed",
  component = "all",
  rope_range = rope_range(B_REP_len),
  centrality = "all",
  dispersion = TRUE
)
# drop sigma row
posterior_desc <- posterior_desc[-c(1,24),]

# clean up labels:
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'b_', '')
posterior_desc$TREATMENT <- ifelse(grepl("TREATMENTSUB", posterior_desc$Parameter), "SC",
                             ifelse(grepl("TREATMENTFLUT", posterior_desc$Parameter), "DT", "DC"))
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTFLUT', 'DT')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTSUB', 'SC')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'IAGE_zE2', 'Age^2')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'AGE_z', 'Age')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'SEXM', 'Male')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'WEIGHT_z', 'Weight offset')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'COMP_NORM_z', 'Competition')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'IGS_zE2', 'Group size^2')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'GS_z', 'Group size')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'RAIN_z', 'Monthly rainfall')

custom_order <- c('Monthly rainfall',
                  'Group size^2',
                  'Group size',
                  'Competition', 
                  'Weight offset', 
                  'DT:Age^2:Male','SC:Age^2:Male','Age^2:Male',
                  'DT:Age:Male','SC:Age:Male','Age:Male', 
                  'DT:Male','SC:Male','Male', 
                  'DT:Age^2','SC:Age^2',"Age^2", 
                  'DT:Age','SC:Age',"Age", 
                  'DT','SC')#,'Intercept') 

posterior_desc$Parameter <- factor(posterior_desc$Parameter, levels = custom_order)
posterior_desc$TREATMENT <- factor(posterior_desc$TREATMENT, levels = c("DC", "SC", "DT"))
# Coeff_REP_LEN 700*800
ggplot(posterior_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
 # ggtitle('REP call length')+
  theme_clean()+
  theme(legend.position="none")

# Ontogeny 30 - 130 plot ####
# get all needed values
(sd_age <- sd(REP_data$REC_AGE_D))#
(mean_age <- mean(REP_data$REC_AGE_D))#

range(REP_data$REC_AGE_D)# 31, 130
rec_age_c <- seq(30, 130, by=1)
age_z_vals <- (rec_age_c - mean(REP_data$REC_AGE_D))/sd(REP_data$REC_AGE_D)
age_z_vals <- seq(min(age_z_vals), max(age_z_vals), length.out = 20)

REP_LEN_pred <- B_REP_len %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(REP_data$TREATMENT),
                                    SEX = levels(REP_data$SEX),
                                    AGE_z = age_z_vals,
                                    WEIGHT_z = mean(REP_data$WEIGHT_z),
                                    COMP_NORM_z = mean(REP_data$COMP_NORM_z),
                                    GS_z = mean(REP_data$GS_z),
                                    RAIN_z = mean(REP_data$RAIN_z)),
               re_formula = NA,  robust=T) 

#unscale AGE_z values:
REP_LEN_pred$REC_AGE_D <- REP_LEN_pred$AGE_z * sd_age + mean_age
# ensure right format
REP_LEN_pred$REP_len <- REP_LEN_pred$.epred
REP_LEN_pred$TREATMENT <- factor(REP_LEN_pred$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))
REP_data$TREATMENT <- factor(REP_data$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))

#800*500: REP_LEN_TAS
ggplot(REP_LEN_pred, aes(x = REC_AGE_D, y = REP_len, color = TREATMENT, fill = TREATMENT)) +  
  geom_point(data=REP_data, aes(y=BEG_avg_Len))+
  stat_lineribbon(.width = .95) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Repeat call length (s)\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Repeat call length (s) \n", n.breaks = 10) +
  theme_clean()+
  facet_wrap(~SEX)

rm(REP_LEN_pred,  B_REP_len)

################################################################################
######################## Call interval length ##################################
################################################################################

ggplot(REP_data, aes(x = BEG_avg_Int)) + 
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  labs(title = "Histogram of Call interval duration", x = "Call interval (/s)")

#### BAYES MODELS ##############################################################
# 1) Check for non-linearity ####
B_REP_int_TA <- brms::brm(formula = BEG_avg_Int ~ TREATMENT * AGE_z * SEX + (1|LITTER_CODE/ID),
                          data = REP_data, family = lognormal(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                          prior = priors, threads = threading(4),
                          file="B_REP_int_TA")
B_REP_int_TA <- add_criterion(B_REP_int_TA, c("loo", "loo_R2"), moment_match = TRUE)
summary(B_REP_int_A)
plot(B_REP_int_A)
pp_check(B_REP_int_A, ndraws=100)

B_REP_int_TA2 <- brms::brm(formula = BEG_avg_Int ~  TREATMENT * (AGE_z + I(AGE_z^2))* SEX  + (1|LITTER_CODE/ID), 
                           data = REP_data, family = lognormal(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                           prior = priors, threads = threading(4),
                           file="B_REP_int_TA2")
B_REP_int_TA2 <- add_criterion(B_REP_int_TA2, c("loo", "loo_R2"), moment_match = TRUE)
summary(B_REP_int_A2)
plot(B_REP_int_A2)
pp_check(B_REP_int_A2, ndraws=100)

loo(B_REP_int_TA, B_REP_int_TA2)


# 2) Define models for REP interval length ####
# B_REP_int<- brms::brm(formula = BEG_avg_Int ~ TTREATMENT * (AGE_z + I(AGE_z^2)) * SEX + 
#                         WEIGHT_z + COMP_NORM_z + GS_z + I(GS_z^2) + RAIN_z + (1|LITTER_CODE/ID),
#                              data = REP_data, family = lognormal(link='identity'),
#                              chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
#                              save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
#                              prior = priors, threads = threading(4),
#                              file="B_REP_int")
# B_REP_int <- add_criterion(B_REP_int, c("loo", "loo_R2"), moment_match = TRUE)

B_REP_int <- readRDS("models/B_REP_int.rds")

#### RESULTS: REP interval length ####
summary(B_REP_int)

plot(B_REP_int)
pp_check(B_REP_int, ndraws = 100)

describe_posterior(
  B_REP_int,
  effects = "all", #fixed vs all (for random effects)
  component = "all",
  rope_range = rope_range(B_REP_int),  
  test = c("p_direction", "p_significance", "rope"),
  centrality = "all",
  dispersion = TRUE
)

process_model_results(B_REP_int, data = REP_data, log_scale = T)

loo_R2(B_REP_int) 

bayes_R2(B_REP_int)

performance::variance_decomposition(B_REP_int)

### EMMs at different ages TAS ####
(emm_REP_INT <- emmeans(B_REP_int,  ~ TREATMENT:SEX | AGE_z, at = list('AGE_z' = get_age_vars())))

pd(emm_REP_INT)

equivalence_test(emm_REP_INT, rope = rope_range(B_REP_int))

p_significance(emm_REP_INT, threshold = rope_range(B_REP_int))

# Contrasts between treatments within sex at each age
(pairs_within_sex <- contrast(emm_REP_INT, method = "pairwise", by = c("SEX", 'AGE_z')))

calculate_percent_difference(pairs_within_sex, emm_REP_INT, log_scale = T)

pd(pairs_within_sex)

equivalence_test(pairs_within_sex, rope = rope_range(B_REP_int))

p_significance(pairs_within_sex, threshold = rope_range(B_REP_int))

(pairs_within_treatment <- # Contrasts between sex within treatment at each age
        contrast(emm_REP_INT, method = "pairwise", by = c("TREATMENT", 'AGE_z')))

calculate_percent_difference(pairs_within_treatment, emm_REP_INT, log_scale = T)

pd(pairs_within_treatment)

equivalence_test(pairs_within_treatment, rope = rope_range(B_REP_int))

p_significance(pairs_within_treatment, threshold = rope_range(B_REP_int))

rm(emm_REP_INT, pairs_within_sex, pairs_within_treatment)

### PLOTS: REP INTERVAL LENGTH ####
# coefficients:
posterior_desc <- describe_posterior(
  B_REP_int,
  effects = "fixed",
  component = "all",
  rope_range = rope_range(B_REP_int),
  test = c("p_direction", "pnificance", "rope"),
  centrality = "all",
  dispersion = TRUE
)
# drop sigma row
posterior_desc <- posterior_desc[-c(24, 1),]

# clean up labels:
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'b_', '')
posterior_desc$TREATMENT <- ifelse(grepl("TREATMENTSUB", posterior_desc$Parameter), "SC",
                                   ifelse(grepl("TREATMENTFLUT", posterior_desc$Parameter), "DT", "DC"))
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTSUB', 'SC')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTFLUT', 'DT')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'IAGE_zE2', 'Age^2')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'AGE_z', 'Age')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'SEXM', 'Male')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'WEIGHT_z', 'Weight offset')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'COMP_NORM_z', 'Competition')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'IGS_zE2', 'Group size^2')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'GS_z', 'Group size')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'RAIN_z', 'Monthly rainfall')

custom_order <- c('DT:Age^2:Male','SC:Age^2:Male','Age^2:Male',
                  'DT:Age:Male','SC:Age:Male','Age:Male', 
                  'DT:Monthly rainfall','SC:Monthly rainfall','Monthly rainfall',
                  'DT:Group size^2','SC:Group size^2','Group size^2',
                  'DT:Group size','SC:Group size','Group size',
                  'DT:Competition','SC:Competition','Competition', 
                  'DT:Weight offset','SC:Weight offset','Weight offset', 
                  'DT:Male','SC:Male','Male', 
                  'DT:Age^2','SC:Age^2',"Age^2", 
                  'DT:Age','SC:Age',"Age", 
                  'DT','SC') 
posterior_desc$Parameter <- factor(posterior_desc$Parameter, levels = custom_order)
posterior_desc$TREATMENT <- factor(posterior_desc$TREATMENT, levels = c("DC", "SC", "DT"))

# Coeff_REP_INT 700*800
ggplot(posterior_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
 # ggtitle('REP call interval length')+
  theme_clean()+
  theme(legend.position="none")

# Ontogeny 30 - 130 ####
# get all needed values
sd_age <- sd(REP_data$REC_AGE_D)#
mean_age <- mean(REP_data$REC_AGE_D)#

range(REP_data$REC_AGE_D)# 31, 130
rec_age_c <- seq(30, 130, by=1)
age_z_vals <- (rec_age_c - mean(REP_data$REC_AGE_D))/sd(REP_data$REC_AGE_D)
age_z_vals <- seq(min(age_z_vals), max(age_z_vals), length.out = 20)

REP_INT_pred <- B_REP_int %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(REP_data$TREATMENT),
                                    SEX = levels(REP_data$SEX),
                                    AGE_z = age_z_vals,
                                    WEIGHT_z = mean(REP_data$WEIGHT_z),
                                    COMP_NORM_z = mean(REP_data$COMP_NORM_z),
                                    GS_z = mean(REP_data$GS_z),
                                    RAIN_z = mean(REP_data$RAIN_z)),
               re_formula = NA,  robust=T) 


#unscale AGE_z values:
REP_INT_pred$REC_AGE_D <- REP_INT_pred$AGE_z * sd_age + mean_age
# ensure right format
REP_INT_pred$REP_int <- REP_INT_pred$.epred
REP_INT_pred$TREATMENT <- factor(REP_INT_pred$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))
REP_data$TREATMENT <- factor(REP_data$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))

#800*500: REP_INT_TAS
ggplot(REP_INT_pred, aes(x = REC_AGE_D, y = REP_int, color = TREATMENT, fill = TREATMENT)) +  
  geom_point(data=REP_data, aes(y=BEG_avg_Int))+
  stat_lineribbon(.width = .95) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Repeat call interval duration (s)\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Repeat call interval duration (s) \n", n.breaks = 10) +
  theme_clean()+
  facet_wrap(~ SEX)

rm(REP_INT_pred, B_REP_int)

################################################################################
######################## Call rate #############################################
################################################################################
ggplot(REP_data, aes(x = BEG_rate)) + 
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  labs(title = "Histogram of Call rate", x = "Call rate (calls/s)")

priors <- c(set_prior("normal(2,2)", class = "Intercept"), set_prior("normal(0,1)", class='b'))

#### BAYES MODELS ##############################################################
# 1) Check for non-linearity ####
B_REP_rat_TA <- brms::brm(formula = BEG_rate ~ TREATMENT * AGE_z * SEX  + (1|LITTER_CODE/ID),
                          data = REP_data, family = gaussian(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                          prior = priors, threads = threading(4),
                          file="B_REP_rat_TA")
B_REP_rat_TA <- add_criterion(B_REP_rat_TA, c("loo", "loo_R2"), moment_match = TRUE)
summary(B_REP_rat_TA)
plot(B_REP_rat_TA)
pp_check(B_REP_rat_TA, ndraws=100) 

B_REP_rat_TA2 <- brms::brm(formula = BEG_rate ~ TREATMENT * (AGE_z + I(AGE_z^2)) * SEX  + (1|LITTER_CODE/ID), 
                          data = REP_data, family = gaussian(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                          prior = priors, threads = threading(4),
                          file="B_REP_rat_TA2")
B_REP_rat_TA2 <- add_criterion(B_REP_rat_TA2, c("loo", "loo_R2"), moment_match = TRUE)
summary(B_REP_rat_A2)
plot(B_REP_rat_A2)
pp_check(B_REP_rat_A2, ndraws=100)

loo(B_REP_rat_TA, B_REP_rat_TA2)

# 2) Define models for REP rate ####
# priors <- c(set_prior("normal(2,2)", class = "Intercept"), set_prior("normal(0,1)", class='b'))
# B_REP_rat <- brms::brm(formula = BEG_rate ~ TREATMENT * (AGE_z + I(AGE_z^2)) * SEX + 
#                          WEIGHT_z + COMP_NORM_z + GS_z + I(GS_z^2) + RAIN_z + (1|LITTER_CODE/ID),,
#                        data = REP_data, family = gaussian(link='identity'),
#                        chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
#                        save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
#                        prior = priors, threads = threading(4),
#                        file="B_REP_rat")
# B_REP_rat <- add_criterion(B_REP_rat, c("loo", "loo_R2"), moment_match = TRUE)

B_REP_rat <- readRDS("models/B_REP_rat.rds")

#### RESULTS: REP rate ####
summary(B_REP_rat)

plot(B_REP_rat)
pp_check(B_REP_rat, ndraws = 100)

describe_posterior(
  B_REP_rat,
  effects = "all", #fixed vs all (for random effects)
  component = "all",
  rope_range = rope_range(B_REP_rat),  
  test = c("p_direction", "p_significance", "rope"),
  centrality = "all",
  dispersion = TRUE
)

process_model_results(B_REP_rat, REP_data, log_scale = F)

loo_R2(B_REP_rat, moment_match=T) 

bayes_R2(B_REP_rat)

performance::variance_decomposition(B_REP_rat)

### EMMs at different ages TAS ####
(emm_REP_RAT <- emmeans(B_REP_rat,  ~ TREATMENT:SEX | AGE_z, at = list('AGE_z' = get_age_vars())))

pd(emm_REP_RAT)

equivalence_test(emm_REP_RAT, rope = rope_range(B_REP_rat))

p_significance(emm_REP_RAT, threshold = rope_range(B_REP_rat))

(pairs_within_sex <- # Contrasts between treatments within sex at each age
        contrasts_treatments_within_sex <- contrast(emm_REP_RAT, method = "pairwise", by = c("SEX", 'AGE_z')))

calculate_percent_difference(pairs_within_sex, emm_REP_RAT, log_scale = F)

pd(pairs_within_sex)

equivalence_test(pairs_within_sex, rope = rope_range(B_REP_rat))

p_significance(pairs_within_sex, threshold = rope_range(B_REP_rat))

(pairs_within_treatment <- # Contrasts between sex within treatment at each age
         contrast(emm_REP_RAT, method = "pairwise", by = c("TREATMENT", 'AGE_z')))

calculate_percent_difference(pairs_within_treatment, emm_REP_RAT, log_scale = F)


pd(pairs_within_treatment)

equivalence_test(pairs_within_treatment, rope = rope_range(B_REP_rat))

p_significance(pairs_within_treatment, threshold = rope_range(B_REP_rat))

rm(emm_REP_RAT, pairs_within_sex, pairs_within_treatment)

### PLOTS: REP RATE ####
# coefficients:
posterior_desc <- describe_posterior(
  B_REP_rat,
  effects = "fixed",
  component = "all",
  rope_range = rope_range(B_REP_rat),
  centrality = "all",
  dispersion = TRUE
)
# drop sigma row
posterior_desc <- posterior_desc[-c(24, 1),]

# clean up labels:
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'b_', '')
posterior_desc$TREATMENT <- ifelse(grepl("TREATMENTSUB", posterior_desc$Parameter), "SC",
                                   ifelse(grepl("TREATMENTFLUT", posterior_desc$Parameter), "DT", "DC"))
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTSUB', 'SC')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTFLUT', 'DT')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'IAGE_zE2', 'Age^2')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'AGE_z', 'Age')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'SEXM', 'Male')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'WEIGHT_z', 'Weight offset')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'COMP_NORM_z', 'Competition')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'IGS_zE2', 'Group size^2')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'GS_z', 'Group size')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'RAIN_z', 'Monthly rainfall')

custom_order <- c('Monthly rainfall',
                  'Group size^2',
                  'Group size',
                  'Competition', 
                  'Weight offset', 
                  'DT:Age^2:Male','SC:Age^2:Male','Age^2:Male',
                  'DT:Age:Male','SC:Age:Male','Age:Male', 
                  'DT:Male','SC:Male','Male', 
                  'DT:Age^2','SC:Age^2',"Age^2", 
                  'DT:Age','SC:Age',"Age", 
                  'DT','SC') 

posterior_desc$Parameter <- factor(posterior_desc$Parameter, levels = custom_order)
posterior_desc$TREATMENT <- factor(posterior_desc$TREATMENT, levels = c("DC", "SC", "DT"))

# Coeff_REP_RAT 700*800
ggplot(posterior_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
 # ggtitle('REP call rate')+
  theme_clean()+
  theme(legend.position="none")

rm(posterior_desc)

# Ontogeny 30 - 130 
# get all needed values
(sd_age <- sd(REP_data$REC_AGE_D))#
(mean_age <- mean(REP_data$REC_AGE_D))#

range(REP_data$REC_AGE_D)# 31, 130
rec_age_c <- seq(30, 130, by=1)
age_z_vals <- (rec_age_c - mean(REP_data$REC_AGE_D))/sd(REP_data$REC_AGE_D)
age_z_vals <- seq(min(age_z_vals), max(age_z_vals), length.out = 20)

REP_RAT_pred <- B_REP_rat %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(REP_data$TREATMENT),
                                    SEX = levels(REP_data$SEX),
                                    AGE_z = age_z_vals,
                                    WEIGHT_z = mean(REP_data$WEIGHT_z),
                                    COMP_NORM_z = mean(REP_data$COMP_NORM_z),
                                    GS_z = mean(REP_data$GS_z),
                                    RAIN_z = mean(REP_data$RAIN_z)),
               re_formula = NA,  robust=T) 

#unscale AGE_z values:
REP_RAT_pred$REC_AGE_D <- REP_RAT_pred$AGE_z * sd_age + mean_age
# ensure right format
REP_RAT_pred$REP_rate <- REP_RAT_pred$.epred
REP_RAT_pred$TREATMENT <- factor(REP_RAT_pred$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))
REP_data$TREATMENT <- factor(REP_data$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))

#800*500: REP_RAT_TAS
ggplot(REP_RAT_pred, aes(x = REC_AGE_D, y = REP_rate, color = TREATMENT, fill = TREATMENT)) +  
  geom_point(data=REP_data, aes(y=BEG_rate))+
  stat_lineribbon(.width = .95) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Repeat call rate (calls/s)\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Repeat call rate (calls/s) \n", n.breaks = 10) +
  theme_clean()+
  facet_wrap(~SEX)

rm(REP_RAT_pred, B_REP_rat)

# Cleanup ####
rm(priors, REP_data)
