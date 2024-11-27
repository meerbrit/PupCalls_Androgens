#### Analysis script for Flutamide study: Digging (DIG) calls ###################
############## BWalkenhorst 2024 ###############################################

#### SETUP ####
# clear everything and use garbage collector at start
rm(list = ls(all.names = TRUE)) #includes hidden objects.
gc() 

# load all necessary libraries
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
library(cowplot)

set.seed(23)

DIG_data <- readRDS('data/DIG_data.rds')

#Generic weakly informative prior: normal(0, 1) for call length and interval duration
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

# Calculate the % difference for emtrends/slopes
calculate_slope_percent_difference <- function(data) {
  data <- as.data.frame(data)
  
  # Separate slope data and contrast data
  slope_data <- data %>% filter(contrast == ".") # Rows without contrasts
  contrast_data <- data %>% filter(contrast != ".") # Rows with contrasts
  
  # Extract group1 and group2 from the contrast column
  contrast_data <- contrast_data %>%
    tidyr::separate(contrast, into = c("group1", "group2"), sep = " - ", remove = FALSE)
  
  # Extract TREATMENT and SEX from group1 and group2
  contrast_data <- contrast_data %>%
    mutate(
      group1_treatment = gsub(" .*", "", group1), # Extract TREATMENT (e.g., CTRL)
      group1_sex = gsub(".* ", "", group1),       # Extract SEX (e.g., F)
      group2_treatment = gsub(" .*", "", group2),
      group2_sex = gsub(".* ", "", group2)
    )
  
  # Match slopes for group1 and group2
  contrast_data <- contrast_data %>%
    left_join(
      slope_data %>% rename(group1_trend = AGE_z.trend, group1_lower = lower.HPD, group1_upper = upper.HPD),
      by = c("group1_treatment" = "TREATMENT", "group1_sex" = "SEX")
    ) %>%
    left_join(
      slope_data %>% rename(group2_trend = AGE_z.trend, group2_lower = lower.HPD, group2_upper = upper.HPD),
      by = c("group2_treatment" = "TREATMENT", "group2_sex" = "SEX")
    )
  
  # Calculate percentage differences
  contrast_data <- contrast_data %>%
    mutate(
      percent_difference = (AGE_z.trend / abs(group1_trend)) * 100
    )
  
  # Return only relevant columns
  return(contrast_data %>%
           select(
             contrast, group1, group2, TREATMENT, AGE_z.trend,
             lower.HPD, upper.HPD, group1_trend, group2_trend, percent_difference
           ))
}

# get key ages for EMMs
get_age_vars <- function(){
  age_vars <- c((30- mean(DIG_data$REC_AGE_D))/sd(DIG_data$REC_AGE_D), 
                (75- mean(DIG_data$REC_AGE_D))/sd(DIG_data$REC_AGE_D), 
                (120- mean(DIG_data$REC_AGE_D))/sd(DIG_data$REC_AGE_D))
  return(age_vars)
}

# calculate EMM contrast differences in %
calculate_emm_percent_difference <- function(contrasts, emmeans_output, log_scale = TRUE) {
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
ggplot(DIG_data, aes(x = DIG_avg_Len)) + 
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  labs(title = "Histogram of Call length", x = "Call length (/s)")

#### BAYES MODELS ##############################################################
# 1) Check for non-linearity ####
B_DIG_len_TA <- brms::brm(formula = DIG_avg_Len ~ TREATMENT * AGE_z * SEX + (1|LITTER_CODE/ID),
                           data = DIG_data, family = lognormal(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                           prior = priors, threads = threading(4),
                           file="B_DIG_len_TA")
B_DIG_len_TA <- add_criterion(B_DIG_len_TA, c("loo", "loo_R2"), moment_match = TRUE)
summary(B_DIG_len_A)
plot(B_DIG_len_A)
pp_check(B_DIG_len_A, ndraws=100)

B_DIG_len_TA2 <- brms::brm(formula = DIG_avg_Len ~  TREATMENT * (AGE_z + I(AGE_z^2)) * SEX + (1|LITTER_CODE/ID),
                          data = DIG_data, family = lognormal(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                          prior = priors, threads = threading(4),
                          file="B_DIG_len_TA2")
B_DIG_len_TA2 <- add_criterion(B_DIG_len_TA2, c("loo", "loo_R2"), moment_match = TRUE)
summary(B_DIG_len_A2)
plot(B_DIG_len_A2)
pp_check(B_DIG_len_A2, ndraws=100)

loo(B_DIG_len_TA, B_DIG_len_TA2)


# 2) Define models for DIG LENGTH ####
# B_DIG_len <- brms::brm(formula = DIG_avg_Len ~ TREATMENT * AGE_z * SEX + 
#                          WEIGHT_z + COMP_NORM_z + GS_z + I(GS_z^2) + RAIN_z +  (1|LITTER_CODE/ID),
#                              data = DIG_data, family = lognormal(link='identity'),
#                              chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
#                              save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
#                              prior = priors, threads = threading(4),
#                              file="B_DIG_len")
# B_DIG_len <- add_criterion(B_DIG_len, c("loo", "loo_R2"), moment_match = TRUE)

B_DIG_len <- readRDS("models/B_DIG_len.rds")

#### RESULTS: DIG LENGTH ####
summary(B_DIG_len)

plot(B_DIG_len)
pp_check(B_DIG_len, ndraws = 100)

describe_posterior(
  B_DIG_len,
  effects = "all", #fixed vs all (for random effects)
  component = "all",
  rope_range = rope_range(B_DIG_len),  
  test = c("p_direction", "p_significance", "rope"),
  centrality = "all",
  dispersion = TRUE
)

process_model_results(B_DIG_len, DIG_data, log_scale = T)

loo_R2(B_DIG_len) 

bayes_R2(B_DIG_len)

performance::variance_decomposition(B_DIG_len)


#### EMMs: TAS ####
(emm_DIG_LEN <- emtrends(B_DIG_len, pairwise ~ TREATMENT:SEX, var='AGE_z'))

pd(emm_DIG_LEN) 

equivalence_test(emm_DIG_LEN, range = rope_range(B_DIG_len)) 

p_significance(emm_DIG_LEN, threshold = rope_range(B_DIG_len))

# get the percentage differences:
calculate_slope_percent_difference(emm_DIG_LEN)

rm(emm_DIG_LEN)

### EMMs at different ages TAS ####
(emm_DIG_LEN <- emmeans(B_DIG_len,  ~ TREATMENT:SEX | AGE_z, at = list('AGE_z' = get_age_vars())))

pd(emm_DIG_LEN)

equivalence_test(emm_DIG_LEN, rope = rope_range(B_DIG_len))

p_significance(emm_DIG_LEN, threshold = rope_range(B_DIG_len))

(pairs_within_sex <- # Contrasts between treatments within sex at each age
    contrast(emm_DIG_LEN, method = "pairwise", by = c("SEX", 'AGE_z')))
calculate_emm_percent_difference(pairs_within_sex, emm_DIG_LEN, log_scale = T)

pd(pairs_within_sex)

equivalence_test(pairs_within_sex, rope = rope_range(B_DIG_len))

p_significance(pairs_within_sex, threshold = rope_range(B_DIG_len))

(pairs_within_treatment <- # Contrasts between sex within treatment at each age
    contrast(emm_DIG_LEN, method = "pairwise", by = c("TREATMENT", 'AGE_z')))

calculate_emm_percent_difference(pairs_within_treatment, emm_DIG_LEN, log_scale = T)

p_direction(pairs_within_treatment)

equivalence_test(pairs_within_treatment, range = rope_range(B_DIG_len))

p_significance(pairs_within_treatment, threshold = rope_range(B_DIG_len))

rm(emm_DIG_LEN, pairs_within_sex, pairs_within_treatment)

### PLOTS: DIG LENGTH ####
# coefficients:
posterior_desc <- describe_posterior(
  B_DIG_len,
  effects = "fixed",
  component = "all",
  rope_range = rope_range(B_DIG_len),
  test = c("p_direction", "p_significance", "rope"),
  centrality = "all",
  dispersion = TRUE
)
# drop sigma row
posterior_desc <- posterior_desc[-c(18, 1),]

# clean up labels:
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'b_', '')
posterior_desc$TREATMENT <- ifelse(grepl("TREATMENTSUB", posterior_desc$Parameter), "SC",
                             ifelse(grepl("TREATMENTFLUT", posterior_desc$Parameter), "DT", "DC"))
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTFLUT', 'DT')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTSUB', 'SC')
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
                  'DT:Age:Male','SC:Age:Male','Age:Male', 
                  'DT:Male','SC:Male','Male', 
                  'DT:Age','SC:Age',"Age", 
                  'DT','SC') 

posterior_desc$Parameter <- factor(posterior_desc$Parameter, levels = custom_order)
posterior_desc$TREATMENT <- factor(posterior_desc$TREATMENT, levels = c("DC", "SC", "DT"))

# Coeff_DIG_LEN 700*800
ggplot(posterior_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
 # ggtitle('DIG call length')+
  theme_clean()+
  theme(legend.position="none")

# Ontogeny 30 - 130 
# get all needed values
(sd_age <- sd(DIG_data$REC_AGE_D))#25.35114
(mean_age <- mean(DIG_data$REC_AGE_D))#86.1989

range(DIG_data$REC_AGE_D)# 31, 130
rec_age_c <- seq(30, 130, by=1)
age_z_vals <- (rec_age_c - mean(DIG_data$REC_AGE_D))/sd(DIG_data$REC_AGE_D)
age_z_vals <- seq(min(age_z_vals), max(age_z_vals), length.out = 20)

DIG_LEN_pred <- B_DIG_len %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(DIG_data$TREATMENT),
                                    SEX = levels(DIG_data$SEX),
                                    AGE_z = age_z_vals,
                                    WEIGHT_z = mean(DIG_data$WEIGHT_z),
                                    COMP_NORM_z = mean(DIG_data$COMP_NORM_z),
                                    GS_z = mean(DIG_data$GS_z),
                                    RAIN_z = mean(DIG_data$RAIN_z)),
            re_formula = NA,  robust=T)

#unscale AGE_z values:
DIG_LEN_pred$REC_AGE_D <- DIG_LEN_pred$AGE_z * sd_age + mean_age
# ensure right format
DIG_LEN_pred$DIG_len <- DIG_LEN_pred$.epred
DIG_LEN_pred$TREATMENT <- factor(DIG_LEN_pred$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))
DIG_data$TREATMENT <- factor(DIG_data$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))

#800*500: DIG_LEN_TAS
ggplot(DIG_LEN_pred, aes(x = REC_AGE_D, y = DIG_len, color = TREATMENT, fill = TREATMENT)) +  
  geom_point(data=DIG_data, aes(y=DIG_avg_Len))+
  stat_lineribbon(.width = .95) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Digging call length (s)\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Digging call length (s) \n", n.breaks = 10) +
  theme_clean()+ 
  facet_wrap(~SEX)

rm(DIG_LEN_pred, B_DIG_len)

################################################################################
######################## Call interval length ##################################
################################################################################
ggplot(DIG_data, aes(x = DIG_avg_Int)) + 
     geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
     labs(title = "Histogram of Call interval duration", x = "Call interval (/s)")

#### BAYES MODELS ##############################################################
# 1) Check for non-linearity ####
B_DIG_int_TA <- brms::brm(formula = DIG_avg_Int ~ TREATMENT * AGE_z * SEX  + (1|LITTER_CODE/ID),
                          data = DIG_data, family = lognormal(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                          prior = priors, threads = threading(4),
                          file="B_DIG_int_TA")
B_DIG_int_TA <- add_criterion(B_DIG_int_TA, c("loo", "loo_R2"), moment_match = TRUE)
summary(B_DIG_int_TA)
plot(B_DIG_int_TA)
pp_check(B_DIG_int_TA, ndraws=100)

B_DIG_int_TA2 <- brms::brm(formula = DIG_avg_Int ~ TREATMENT * (AGE_z + I(AGE_z^2)) * SEX  + (1|LITTER_CODE/ID), 
                          data = DIG_data, family = lognormal(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                          prior = priors, threads = threading(4),
                          file="B_DIG_int_TA2")
B_DIG_int_TA2 <- add_criterion(B_DIG_int_TA2, c("loo", "loo_R2"), moment_match = TRUE)
summary(B_DIG_int_TA2)
plot(B_DIG_int_TA2)
pp_check(B_DIG_int_TA2, ndraws=100)

loo(B_DIG_int_TA, B_DIG_int_TA2)


# 2) Define models for DIG interval length ####
# B_DIG_int <- brms::brm(formula = DIG_avg_Int ~ TREATMENT * AGE_z * SEX + 
#                          WEIGHT_z + COMP_NORM_z + GS_z + I(GS_z^2) + RAIN_z + (1|LITTER_CODE/ID),
#                              data = DIG_data, family = lognormal(link='identity'),
#                              chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
#                              save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
#                              prior = priors, threads = threading(4),
#                              file="B_DIG_int")
# B_DIG_int <- add_criterion(B_DIG_int, c("loo", "loo_R2"), moment_match = TRUE)

B_DIG_int <- readRDS("models/B_DIG_int.rds")

#### RESULTS: DIG interval length ####
summary(B_DIG_int)

plot(B_DIG_int)
pp_check(B_DIG_int, ndraws = 100)

describe_posterior(
  B_DIG_int,
  effects = "all", #fixed vs all (for random effects)
  component = "all",
  rope_range = rope_range(B_DIG_int),  
  test = c("p_direction", "p_significance", "rope"),
  centrality = "all",
  dispersion = TRUE
)

process_model_results(B_DIG_int, DIG_data, log_scale = T)

loo_R2(B_DIG_int) 

bayes_R2(B_DIG_int)

performance::variance_decomposition(B_DIG_int)

### EMMS: DIG interval length: TAS ####
(emm_DIG_INT <- emtrends(B_DIG_int, pairwise ~ TREATMENT:SEX, var = 'AGE_z'))

pd(emm_DIG_INT) 

equivalence_test(emm_DIG_INT, range = rope_range(B_DIG_int)) 


p_significance(emm_DIG_INT, threshold = rope_range(B_DIG_int)) 

# get the percentage differences:
calculate_slope_percent_difference(emm_DIG_INT)

rm(emm_DIG_INT)

### EMMs at different ages TAS ####
(emm_DIG_INT <- emmeans(B_DIG_int,  ~ TREATMENT:SEX | AGE_z, at = list('AGE_z' = get_age_vars())))

pd(emm_DIG_INT)

equivalence_test(emm_DIG_INT, rope = rope_range(B_DIG_int))

p_significance(emm_DIG_INT, threshold = rope_range(B_DIG_int))

(pairs_within_sex <- # Contrasts between treatments within sex at each age
   contrast(emm_DIG_INT, method = "pairwise", by = c("SEX", 'AGE_z')))

calculate_emm_percent_difference(pairs_within_sex, emm_DIG_INT, log_scale = T)

pd(pairs_within_sex)

equivalence_test(pairs_within_sex, rope = rope_range(B_DIG_int))

p_significance(pairs_within_sex, threshold = rope_range(B_DIG_int))

(pairs_within_treatment <- # Contrasts between sex within treatment at each age
  contrast(emm_DIG_INT, method = "pairwise", by = c("TREATMENT", 'AGE_z')))

calculate_emm_percent_difference(pairs_within_treatment, emm_DIG_INT, log_scale = T)

p_direction(pairs_within_treatment)

equivalence_test(pairs_within_treatment, range = rope_range(B_DIG_int))

p_significance(pairs_within_treatment, threshold = rope_range(B_DIG_int))

rm(emm_DIG_INT, pairs_within_sex, pairs_within_treatment)

### PLOTS: DIG INTERVAL DURATION ####
# coefficients:
posterior_desc <- describe_posterior(
  B_DIG_int,
  effects = "fixed",
  component = "all",
  rope_range = rope_range(B_DIG_int),
  test = c("p_direction", "p_significance", "rope"),
  centrality = "all",
  dispersion = TRUE
)
# drop sigma row
posterior_desc <- posterior_desc[-c(18,1),]

# clean up labels:
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'b_', '')
posterior_desc$TREATMENT <- ifelse(grepl("TREATMENTSUB", posterior_desc$Parameter), "SC",
                                   ifelse(grepl("TREATMENTFLUT", posterior_desc$Parameter), "DT", "DC"))
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTFLUT', 'DT')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTSUB', 'SC')
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
                  'DT:Age:Male','SC:Age:Male','Age:Male', 
                  'DT:Male','SC:Male','Male', 
                  'DT:Age','SC:Age',"Age", 
                  'DT','SC') 
posterior_desc$Parameter <- factor(posterior_desc$Parameter, levels = custom_order)
posterior_desc$TREATMENT <- factor(posterior_desc$TREATMENT, levels = c("DC", "SC", "DT"))

# Coeff_DIG_INT 700*800
ggplot(posterior_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
 # ggtitle('DIG call interval duration')+
  theme_clean()+
  theme(legend.position="none")

# Ontogeny 30 - 130 
# get all needed values
(sd_age <- sd(DIG_data$REC_AGE_D))#
(mean_age <- mean(DIG_data$REC_AGE_D))#

range(DIG_data$REC_AGE_D)# 31, 130
rec_age_c <- seq(30, 130, by=1)
age_z_vals <- (rec_age_c - mean(DIG_data$REC_AGE_D))/sd(DIG_data$REC_AGE_D)
age_z_vals <- seq(min(age_z_vals), max(age_z_vals), length.out = 20)

DIG_INT_pred <- B_DIG_int %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(DIG_data$TREATMENT),
                                    SEX = levels(DIG_data$SEX),
                                    AGE_z = age_z_vals,
                                    WEIGHT_z = mean(DIG_data$WEIGHT_z),
                                    COMP_NORM_z = mean(DIG_data$COMP_NORM_z),
                                    GS_z = mean(DIG_data$GS_z),
                                    RAIN_z = mean(DIG_data$RAIN_z)),
              re_formula = NA,  robust=T)

#unscale AGE_z values:
DIG_INT_pred$REC_AGE_D <- DIG_INT_pred$AGE_z * sd_age + mean_age
# ensure right format
DIG_INT_pred$DIG_int <- DIG_INT_pred$.epred
DIG_INT_pred$TREATMENT <- factor(DIG_INT_pred$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))
DIG_data$TREATMENT <- factor(DIG_data$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))

#800*500: DIG_INT_TAS
ggplot(DIG_INT_pred, aes(x = REC_AGE_D, y = DIG_int, color = TREATMENT, fill = TREATMENT)) +  
  geom_point(data=DIG_data, aes(y=DIG_avg_Int))+
  stat_lineribbon(.width = .95) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Digging call interval duration (s)\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Digging call interval duration (s) \n", n.breaks = 10) +
  theme_clean()+
  facet_wrap(~SEX)

rm(DIG_INT_pred, B_DIG_int)

################################################################################
######################## Call rate #############################################
################################################################################
ggplot(DIG_data, aes(x = DIG_rate)) + 
       geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
       labs(title = "Histogram of Call Rates", x = "Call Rate (calls/s)")


# define prior for DIG rate
priors <- c(set_prior("normal(2,2)", class = "Intercept"), set_prior("normal(0,1)", class='b'))

#### BAYES MODELS ##############################################################
# 1) Check for non-linearity
B_DIG_rat_TA <- brms::brm(formula = DIG_rate ~ TREATMENT * AGE_z * SEX  + (1|LITTER_CODE/ID),
                          data = DIG_data, family = gaussian(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                          prior = priors, threads = threading(4),
                          file="B_DIG_rat_TA")
B_DIG_rat_TA <- add_criterion(B_DIG_rat_TA, c("loo", "loo_R2"), moment_match = TRUE)
summary(B_DIG_rat_A)
plot(B_DIG_rat_A)
pp_check(B_DIG_rat_A, ndraws=100)

B_DIG_rat_TA2 <- brms::brm(formula = DIG_rate ~ TREATMENT * (AGE_z + I(AGE_z^2)) * SEX + (1|LITTER_CODE/ID), 
                          data = DIG_data, family = gaussian(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                          prior = priors, threads = threading(4),
                          file="B_DIG_rat_TA2")
B_DIG_rat_TA2 <- add_criterion(B_DIG_rat_TA2, c("loo", "loo_R2"), moment_match = TRUE)
summary(B_DIG_rat_A2)
plot(B_DIG_rat_A2)
pp_check(B_DIG_rat_A2, ndraws=100)

loo(B_DIG_rat_TA, B_DIG_rat_TA2)

# 2) Define models for DIG rate ####
# priors <- c(set_prior("normal(2,2)", class = "Intercept"), set_prior("normal(0,1)", class='b'))
# B_DIG_rat <- brms::brm(formula = DIG_rate ~ TREATMENT * AGE_z * SEX + 
#                          WEIGHT_z + COMP_NORM_z + GS_z + I(GS_z^2) + RAIN_z +(1|LITTER_CODE/ID),
#                        data = DIG_data, family = gaussian(link='identity'),
#                        chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
#                        save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
#                        prior = priors, threads = threading(4),
#                        file="B_DIG_rat")
# B_DIG_rat <- add_criterion(B_DIG_rat, c("loo", "loo_R2"), moment_match = TRUE)

B_DIG_rat <- readRDS("models/B_DIG_rat.rds")

#### RESULTS: DIG rate ####
summary(B_DIG_rat)

plot(B_DIG_rat)
pp_check(B_DIG_rat, ndraws = 100) 

describe_posterior(
  B_DIG_rat,
  effects = "all", #fixed vs all (for random effects)
  component = "all",
  rope_range = rope_range(B_DIG_rat),  
  test = c("p_direction", "p_significance", "rope"),
  centrality = "all",
  dispersion = TRUE
)

process_model_results(B_DIG_rat, DIG_data, log_scale = F)

loo_R2(B_DIG_rat) 
bayes_R2(B_DIG_rat)

performance::variance_decomposition(B_DIG_rat)

#### EMMs DIG RAT TAS ####
(emm_DIG_RAT <- emtrends(B_DIG_rat, pairwise ~ TREATMENT:SEX, var = 'AGE_z'))

pd(emm_DIG_RAT) 

equivalence_test(emm_DIG_RAT, range = rope_range(B_DIG_rat)) 

p_significance(emm_DIG_RAT, threshold = rope_range(B_DIG_rat))

## get the percentage differences:
calculate_slope_percent_difference(emm_DIG_RAT)

rm(emm_DIG_RAT)

### EMMs at different ages TAS ####
(emm_DIG_RAT <- emmeans(B_DIG_rat,  ~ TREATMENT:SEX | AGE_z, at = list('AGE_z' = get_age_vars())))

pd(emm_DIG_RAT)
equivalence_test(emm_DIG_RAT, rope = rope_range(B_DIG_rat))

p_significance(emm_DIG_RAT, threshold = rope_range(B_DIG_rat))

(pairs_within_sex <- # Contrasts between treatments within sex at each age
  contrast(emm_DIG_RAT, method = "pairwise", by = c("SEX", 'AGE_z')))

calculate_emm_percent_difference(pairs_within_sex, emm_DIG_RAT, log_scale = T)

pd(pairs_within_sex)

equivalence_test(pairs_within_sex, rope = rope_range(B_DIG_rat))

p_significance(pairs_within_sex, threshold = rope_range(B_DIG_rat))

(pairs_within_treatment <- # Contrasts between sex within treatment at each age
 contrast(emm_DIG_RAT, method = "pairwise", by = c("TREATMENT", 'AGE_z')))

calculate_emm_percent_difference(pairs_within_treatment, emm_DIG_RAT, log_scale = T)

p_direction(pairs_within_treatment)

equivalence_test(pairs_within_treatment, range = rope_range(B_DIG_rat))

p_significance(pairs_within_treatment, threshold = rope_range(B_DIG_rat))

rm(emm_DIG_RAT, pairs_within_sex, pairs_within_treatment)

### PLOTS: DIG RATE ####
# coefficients:
posterior_desc <- describe_posterior(
  B_DIG_rat,
  effects = "fixed",
  component = "all",
  rope_range = rope_range(B_DIG_rat),
  centrality = "all",
  dispersion = TRUE
)
# drop sigma row
posterior_desc <- posterior_desc[-c(18, 1),]

# clean up labels:
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'b_', '')
posterior_desc$TREATMENT <- ifelse(grepl("TREATMENTSUB", posterior_desc$Parameter), "SC",
                                   ifelse(grepl("TREATMENTFLUT", posterior_desc$Parameter), "DT", "DC"))
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTFLUT', 'DT')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTSUB', 'SC')
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
                  'DT:Age:Male','SC:Age:Male','Age:Male', 
                  'DT:Male','SC:Male','Male', 
                  'DT:Age','SC:Age',"Age", 
                  'DT','SC') 
posterior_desc$Parameter <- factor(posterior_desc$Parameter, levels = custom_order)
posterior_desc$TREATMENT <- factor(posterior_desc$TREATMENT, levels = c("DC", "SC", "DT"))

# Coeff_DIG_RAT 700*800
ggplot(posterior_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
 # ggtitle('DIG call rate')+
  theme_clean()+
  theme(legend.position="none")

rm(posterior_desc)

# Ontogeny 30 - 130 ####
# get all needed values
sd_age <- sd(DIG_data$REC_AGE_D)#
mean_age <- mean(DIG_data$REC_AGE_D)#

range(DIG_data$REC_AGE_D)# 31, 130
rec_age_c <- seq(30, 130, by=1)
age_z_vals <- (rec_age_c - mean(DIG_data$REC_AGE_D))/sd(DIG_data$REC_AGE_D)
age_z_vals <- seq(min(age_z_vals), max(age_z_vals), length.out = 20)

DIG_RAT_pred <- B_DIG_rat %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(DIG_data$TREATMENT),
                                    SEX = levels(DIG_data$SEX),
                                    AGE_z = age_z_vals,
                                    WEIGHT_z = mean(DIG_data$WEIGHT_z),
                                    COMP_NORM_z = mean(DIG_data$COMP_NORM_z),
                                    GS_z = mean(DIG_data$GS_z),
                                    RAIN_z = mean(DIG_data$RAIN_z)),
             re_formula = NA,  robust=T)

#unscale AGE_z values:
DIG_RAT_pred$REC_AGE_D <- DIG_RAT_pred$AGE_z * sd_age + mean_age
# ensure right format
DIG_RAT_pred$DIG_rate <- DIG_RAT_pred$.epred
DIG_RAT_pred$TREATMENT <- factor(DIG_RAT_pred$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))
DIG_data$TREATMENT <- factor(DIG_data$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))

# TAS
ggplot(DIG_RAT_pred, aes(x = REC_AGE_D, y = DIG_rate, color = TREATMENT, fill = TREATMENT)) +  
  geom_point(data=DIG_data, aes(y=DIG_rate))+
  stat_lineribbon(.width = .95) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Digging call rate (calls/s)\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Digging call rate (calls/s) \n", n.breaks = 10) +
  theme_clean()+
  facet_wrap(~SEX)

rm(DIG_RAT_pred, B_DIG_rat)

# Cleanup ####
rm(DIG_data, priors)
