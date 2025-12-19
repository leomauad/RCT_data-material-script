# ==============================================================================
# RCT Analysis Script - Mixed Models & ART ANOVA
# Outcomes: Pain (NPS), Handgrip Strength, Isokinetic Strength
# ==============================================================================

# ---- 1. Setup ----
library(tidyverse)
library(lme4)
library(lmerTest)
library(ARTool)
library(effectsize)
library(ggpubr)

setwd("~/Documents/rct/manuscript#1/dataset")
data <- read.csv(file.choose())

data <- data %>%
  mutate(
    id = row_number(),
    intervencao = factor(intervencao, levels = c(1, 2), labels = c("Intervention", "Control")),
    sexo = factor(sexo),
    IMC_Cat = factor(IMC_Cat),
    IPAQ_CAT = factor(IPAQ_CAT)
  )

# ---- 2. Analysis Functions ----

prepare_long_data <- function(data, outcome_pre, outcome_post) {
  data %>%
    select(id, intervencao, pre = all_of(outcome_pre), post = all_of(outcome_post)) %>%
    pivot_longer(cols = c(pre, post), names_to = "time", values_to = "outcome") %>%
    mutate(time = factor(time, levels = c("pre", "post"), labels = c("Pre", "Post"))) %>%
    filter(complete.cases(.))
}

run_analysis <- function(long_data) {
  # Mixed model
  mixed_model <- lmer(outcome ~ intervencao * time + (1|id), data = long_data)
  
  # ART model (non-parametric alternative)
  art_model <- art(outcome ~ intervencao * time, data = long_data)
  
  # Post-hoc comparisons
  posthoc <- art.con(art_model, "intervencao:time", adjust = "holm") %>%
    as.data.frame() %>%
    separate(contrast, into = c("group1", "group2"), sep = " - ") %>%
    mutate(p.adj = p.adjust(p.value, method = "holm"),
           sig = ifelse(p.adj < 0.05, "*", ""))
  
  simple_effects <- art.con(art_model, "intervencao:time", adjust = "bonferroni")
  effect <- eta_squared(art_model)
  
  list(mixed = mixed_model, art = art_model, posthoc = posthoc, 
       simple_effects = simple_effects, eta_sq = effect)
}

create_plot <- function(long_data, outcome_name, y_label = outcome_name) {
  plot_data <- long_data %>%
    group_by(intervencao, time) %>%
    summarise(mean = mean(outcome), se = sd(outcome)/sqrt(n()),
              ci_lower = mean - 1.96*se, ci_upper = mean + 1.96*se, .groups = "drop")
  
  ggplot(plot_data, aes(x = time, y = mean, color = intervencao)) +
    geom_line(aes(group = intervencao), linewidth = 1.5) +
    geom_point(size = 4, shape = 15) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.1, linewidth = 0.8) +
    scale_color_manual(values = c("orange", "tan4")) +
    labs(title = paste("Intervention Effects on", outcome_name),
         subtitle = "Aligned Rank Transformation ANOVA",
         x = "Time Point", y = y_label, color = "Group") +
    theme_pubr() +
    theme(text = element_text(family = "Times"),
          plot.title = element_text(size = 16, face = "bold"),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 10))
}

check_baseline <- function(data, var_pre, var_post) {
  list(
    pre = t.test(as.formula(paste(var_pre, "~ intervencao")), data = data),
    post = t.test(as.formula(paste(var_post, "~ intervencao")), data = data)
  )
}

# ---- 3. Pain (NPS) Analysis ----

long_pain <- prepare_long_data(data, "END1", "END2")
results_pain <- run_analysis(long_pain)
baseline_pain <- check_baseline(data, "END1", "END2")

summary(results_pain$mixed)
anova(results_pain$mixed)
summary(results_pain$art)
results_pain$posthoc
results_pain$simple_effects
results_pain$eta_sq

plot_pain <- create_plot(long_pain, "Pain", "Pain Scale")
print(plot_pain)

# ---- 4. Handgrip Strength Analysis (Unaffected Side) ----

long_handgrip <- prepare_long_data(data, "PPnA1", "PPnA2")
results_handgrip <- run_analysis(long_handgrip)
baseline_handgrip <- check_baseline(data, "PPnA1", "PPnA2")

summary(results_handgrip$mixed)
anova(results_handgrip$mixed)
summary(results_handgrip$art)
results_handgrip$posthoc
results_handgrip$simple_effects
results_handgrip$eta_sq

plot_handgrip <- create_plot(long_handgrip, "Handgrip Strength (Unaffected)", "Handgrip Strength (kg.f)")
print(plot_handgrip)

# ---- 5. Isokinetic Strength Analysis (Abduction) ----

long_iso <- prepare_long_data(data, "ISO_in_abd1", "ISO_in_abd2")
results_iso <- run_analysis(long_iso)
baseline_iso <- check_baseline(data, "ISO_in_abd1", "ISO_in_abd2")

summary(results_iso$mixed)
anova(results_iso$mixed)
summary(results_iso$art)
results_iso$posthoc
results_iso$simple_effects
results_iso$eta_sq

plot_iso <- create_plot(long_iso, "Isokinetic Abduction", "Peak Torque (N.m)")
print(plot_iso)

# ---- 6. Assumption Checking (Optional) ----

check_assumptions <- function(mixed_model) {
  resids <- residuals(mixed_model)
  par(mfrow = c(1, 2))
  qqnorm(resids); qqline(resids)
  plot(fitted(mixed_model), resids, main = "Residuals vs Fitted"); abline(h = 0, col = "red")
  shapiro.test(resids)
}

# Uncomment to run:
# check_assumptions(results_pain$mixed)
# check_assumptions(results_handgrip$mixed)
# check_assumptions(results_iso$mixed)

# ---- 7. Export Results (Optional) ----

# Save all plots to PDF
# pdf("rct_results_plots.pdf", width = 10, height = 6)
# print(plot_pain)
# print(plot_handgrip)
# print(plot_iso)
# dev.off()