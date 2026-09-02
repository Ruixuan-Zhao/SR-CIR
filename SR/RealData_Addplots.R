################################################################################
################ 1.0 mg vs Placebo: RR curves with 95% CI ######################
################################################################################

# Construct the fitted survival ratio curve and its pointwise Wald confidence interval.
make_rr_ci_data <- function(point, cov_alpha, estimator,
                            time_grid = every_years,
                            conf_level = 0.95) {
  # The first two entries are the coefficients of t and t^2.
  alpha <- as.numeric(point[1:2])
  
  # Each row is x(t) = (t, t^2)
  x_time <- cbind(time_grid, time_grid^2)
  
  # log SR(t) = alpha_1 * t + alpha_2 * t^2
  log_rr <- drop(x_time %*% alpha)
  
  # Var{log SR(t)} = x(t)^T Var(alpha) x(t)
  var_log_rr <- rowSums((x_time %*% cov_alpha) * x_time)
  se_log_rr <- sqrt(pmax(var_log_rr, 0))

  z_value <- stats::qnorm(1 - (1 - conf_level) / 2)

  # Exponentiate the log-scale estimate and limits to return to the SR scale.
  data.frame(
    time = time_grid,
    RR = exp(log_rr),
    lower = exp(log_rr - z_value * se_log_rr),
    upper = exp(log_rr + z_value * se_log_rr),
    estimator = estimator
  )
}

# Reconstruct the full MLE sandwich covariance matrix. 
mle2_cov_full <- res2$cov %*%
  crossprod(res2$S_mat) %*%
  t(res2$cov)

# Obtain the covariance block for (alpha_1, alpha_2) only.
mle2_cov_alpha <- mle2_cov_full[1:2, 1:2, drop = FALSE]

# The DR sandwich covariance for alpha.
dr2_cov_alpha <- dr2_out$sand$cov.mat[1:2, 1:2, drop = FALSE]

# Generate the MLE SR curve and pointwise confidence limits.
data_RR_1mg_MLE_CI <- make_rr_ci_data(
  point = res2$point.est,
  cov_alpha = mle2_cov_alpha,
  estimator = "MLE"
)


data_RR_1mg_DR_CI <- make_rr_ci_data(
  point = dr2_out$dr$dr.point,
  cov_alpha = dr2_cov_alpha,
  estimator = "DR"
)


rr_ci_range <- range(
  c(
    data_RR_1mg_MLE_CI$lower,
    data_RR_1mg_MLE_CI$upper,
    data_RR_1mg_DR_CI$lower,
    data_RR_1mg_DR_CI$upper
  ),
  finite = TRUE
)

# Add 5% vertical space so the confidence limits do not touch the panel edge.
rr_ci_padding <- 0.05 * diff(rr_ci_range)

# Keep the lower y-axis limit nonnegative because a survival ratio is positive.
rr_ci_limits <- c(
  max(0, rr_ci_range[1] - rr_ci_padding),
  rr_ci_range[2] + rr_ci_padding
)

# Stack the MLE and DR results into one long-format plotting data set.
data_RR_1mg_CI <- rbind(
  data_RR_1mg_MLE_CI,
  data_RR_1mg_DR_CI
)

# Fix the estimator order used in the legend and line-type scale.
data_RR_1mg_CI$estimator <- factor(
  data_RR_1mg_CI$estimator,
  levels = c("DR", "MLE")
)

# Plot both estimators in one panel. 
RR_1mg_MLE_DR_CI_plot <- ggplot(
  data_RR_1mg_CI,
  aes(
    x = time,
    y = RR,
    linetype = estimator,
    group = estimator
  )
) +
  # Draw the two pointwise confidence bands.
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    fill = "#D55E00",
    alpha = 0.10,
    color = NA,
    show.legend = FALSE
  ) +
  # SR = 1 represents equal survival probabilities in the two treatment arms.
  geom_hline(
    yintercept = 1,
    linetype = "dotted",
    color = "grey45"
  ) +
  # Draw lower CI boundaries using the estimator-specific line types.
  geom_line(
    aes(y = lower),
    color = "#D55E00",
    linewidth = 0.4,
    alpha = 0.55,
    show.legend = FALSE
  ) +
  # Draw upper CI boundaries using the estimator-specific line types.
  geom_line(
    aes(y = upper),
    color = "#D55E00",
    linewidth = 0.4,
    alpha = 0.55,
    show.legend = FALSE
  ) +
  geom_line(
    color = "#D55E00",
    linewidth = 1.1
  ) +
  scale_linetype_manual(
    name = "Estimator",
    values = c("DR" = "dashed", "MLE" = "solid"),
    breaks = c("DR", "MLE")
  ) +
  scale_x_continuous(
    breaks = 0:6,
    limits = c(0, 6)
  ) +
  coord_cartesian(ylim = rr_ci_limits) +
  labs(
    x = "Years",
    y = "Survival Ratio"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 16),
    legend.key.width = grid::unit(1.6, "cm"),
    legend.position = "right",
    panel.grid.minor.y = element_line(
      color = "grey90",
      linewidth = 0.3
    )
)


RR_1mg_MLE_DR_CI_plot

ggsave(
  "Figures/RR_1mg_MLE_DR_CI_plot.png",
  plot = RR_1mg_MLE_DR_CI_plot,
  width = 9,
  height = 6,
  dpi = 300
)


# Extract and print the survival-ratio estimates and pointwise
# 95% confidence intervals at 6 years.
six_year_rr_ci <- subset(
  data_RR_1mg_CI,
  abs(time - 6) < 1e-8,
  select = c(estimator, RR, lower, upper)
)

names(six_year_rr_ci) <- c(
  "Estimator",
  "Estimate",
  "CI_lower",
  "CI_upper"
)

rownames(six_year_rr_ci) <- NULL

# Print the numerical results.
print(six_year_rr_ci, digits = 10)

# Print the results rounded to two decimal places.
for (i in seq_len(nrow(six_year_rr_ci))) {
  cat(sprintf(
    "%s: %.2f (pointwise 95%% CI: %.2f--%.2f)\n",
    as.character(six_year_rr_ci$Estimator[i]),
    six_year_rr_ci$Estimate[i],
    six_year_rr_ci$CI_lower[i],
    six_year_rr_ci$CI_upper[i]
  ))
}


################################################################################
################ Covariate-specific fitted survival curves #####################
################################################################################

# Use the first and third age quartiles as representative prediction values;
profile_age <- unname(
  quantile(prost2$age, probs = c(0.25, 0.75), na.rm = TRUE)
)

# Hold hemoglobin at its sample median in both profiles for comparability.
profile_hg <- median(prost2$hg, na.rm = TRUE)

# Use exactly the same centering constants as in model fitting.
age_center <- unname(attr(age2, "scaled:center"))
hg_center <- unname(attr(hg2, "scaled:center"))

# Define two hypothetical profiles.
survival_profiles <- data.frame(
  profile = c(
    "Lower-risk profile",
    "Higher-risk profile"
  ),
  age = c(profile_age[1], profile_age[2]),
  normal_activity = c(1, 0),
  cardiovascular_history = c(0, 1),
  hg = c(profile_hg, profile_hg)
)

# Construct X in the same order as X2:
# intercept, centered age, normal activity, hx, centered hemoglobin
X_survival_profiles <- cbind(
  intercept = 1,
  age = survival_profiles$age - age_center,
  normal_activity = survival_profiles$normal_activity,
  hx = survival_profiles$cardiovascular_history,
  hg = survival_profiles$hg - hg_center
)

# Print the raw and model-scale profiles
print(survival_profiles)
print(X_survival_profiles)

# Extract alpha and beta
alpha2_mle <- unname(
  res2$point.est[
    grepl("^alpha", names(res2$point.est))
  ]
)
beta2_mle <- unname(
  res2$point.est[
    grepl("^beta", names(res2$point.est))
  ]
)

# Use the same annual grid as in model fitting
profile_timepoints <- timepoints_fit2

# Build f(t) = (t, t^2) for both representative profiles at every time point.
fun_tx_profiles <- make_fun_tx_t_t2(
  time_grid = profile_timepoints,
  n = nrow(X_survival_profiles)
)

# Compute the survival probabilities under placebo and 1.0 mg.
profile_interval_prob <- getProb.sequent(
  alpha = alpha2_mle,
  beta = beta2_mle,
  x_b = X_survival_profiles,
  timepoints = profile_timepoints,
  fun_tx = fun_tx_profiles
)

# Extract probabilities under placebo group (column 1).
interval_prob_placebo <- do.call(
  cbind,
  lapply(profile_interval_prob, function(prob) prob[, 1])
)

# Extract probabilities under 1.0 mg group (column 2).
interval_prob_1mg <- do.call(
  cbind,
  lapply(profile_interval_prob, function(prob) prob[, 2])
)


cumprod_by_row <- function(prob_matrix) {
  t(
    apply(
      prob_matrix,
      MARGIN = 1,
      FUN = cumprod
    )
  )
}

survival_placebo <- cbind(
  1,
  cumprod_by_row(interval_prob_placebo)
)

survival_1mg <- cbind(
  1,
  cumprod_by_row(interval_prob_1mg)
)

# Reshape the profile-by-time matrices into long format for ggplot2.
profile_survival_data <- do.call(
  rbind,
  lapply(seq_len(nrow(survival_profiles)), function(i) {
    rbind(
      data.frame(
        time = profile_timepoints,
        survival = survival_placebo[i, ],
        profile = survival_profiles$profile[i],
        treatment = "Placebo"
      ),
      data.frame(
        time = profile_timepoints,
        survival = survival_1mg[i, ],
        profile = survival_profiles$profile[i],
        treatment = "1.0 mg Estrogen"
      )
    )
  })
)

# Set the left-to-right order of the two profile panels.
profile_survival_data$profile <- factor(
  profile_survival_data$profile,
  levels = c(
    "Lower-risk profile",
    "Higher-risk profile"
  )
)

# Set the treatment order used in the color legend.
profile_survival_data$treatment <- factor(
  profile_survival_data$treatment,
  levels = c(
    "Placebo",
    "1.0 mg Estrogen"
  )
)

# Plot treatment-specific survival curves within each profile.
profile_survival_plot <- ggplot(
  profile_survival_data,
  aes(
    x = time,
    y = survival,
    color = treatment,
    group = treatment
  )
) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 1.8) +
  # Place the lower- and higher-risk profiles in separate panels.
  facet_wrap(
    ~ profile,
    nrow = 1
  ) +
  scale_color_manual(
    name = "Group",
    values = c(
      "Placebo" = "#999999",
      "1.0 mg Estrogen" = "#E41A1C"
    ),
    labels = c(
      "Placebo",
      "1.0 mg Estrogen"
    )
  ) +
  scale_x_continuous(
    breaks = 0:6,
    limits = c(0, 6)
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2)
  ) +
  coord_cartesian(
    ylim = c(0, 1)
  ) +
  labs(
    x = "Years",
    y = "Fitted Survival Probability"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 16),
    legend.key.width = grid::unit(1.6, "cm"),
    legend.position = "right",
    panel.grid.minor.y = element_line(
      color = "grey90",
      linewidth = 0.3
    )
)


profile_survival_plot

ggsave(
  "Figures/Fitted_survival_1mg_profiles.png",
  plot = profile_survival_plot,
  width = 12,
  height = 6,
  dpi = 300
)

