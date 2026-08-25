library(survival)
library(survminer)
library(ggplot2)

####################### Step 1: Data Pre-processing ############################
load("prostate.rda")
table(prostate$rx) # treatment
table(prostate$status) # Status on the end of follow-up

prostate$eventType = prostate$status
levels(prostate$eventType) = list(alive="alive",pdeath="dead - prostatic ca",odeath = c(levels(prostate$eventType)[c(3:10)]))
table(prostate$eventType)

summary(prostate$dtime) # Observed months of follow-up
# covariates age, pf, hx, hg
prost = prostate[!is.na(prostate$age), c(2,3:6,8,9,13,19)] # patno == 42, age == NA
prost$dtime = prost$dtime/12 # convert months to years
dim(prost)


################### Step 2: Kaplan-Meier Survival Curve ########################

data_tot = data.frame(
  time = prost$dtime,
  status = 1 - as.numeric(prost$status=="alive"),
  group = prost$rx
)

km_fit_tot = survfit(Surv(time, status) ~ group, data = data_tot)
summary(km_fit_tot)

# Kaplan-Meier Survival Curve
KM_sur_plot = ggsurvplot(km_fit_tot, data = data_tot, 
                         conf.int = TRUE,   # Show confidence intervals
                         conf.int.alpha = 0.2,
                         risk.table = TRUE, # Add risk table below the plot
                         pval = FALSE,       # Show p-value
                         legend.labs = c("Placebo", "0.2 mg Estrogen", "1.0 mg Estrogen", "5.0 mg Estrogen"),  # Labels for groups
                         legend.title = "Group", 
                         xlab = "Years", 
                         ylab = "Survival Probability", # title = "Kaplan-Meier Survival Curve for Prostate Cancer",
                         palette = c("#999999", "#377EB8", "#E41A1C", "#4DAF4A"),# Colors for different curves
                         # palette = "Set1",
                         linewidth = 1.2,
                         alpha = 1,
                         xlim = c(0, 6),
                         break.time.by = 1,
                         ggtheme = theme_minimal())

KM_sur_plot$plot = KM_sur_plot$plot + theme(plot.title = element_text(hjust = 0.5), 
                                            axis.text = element_text(size=16),
                                            axis.title = element_text(size=18),
                                            legend.text = element_text(size=16),
                                            legend.title = element_text(size=18)) 

KM_sur_plot$plot

ggsave("Figures/KM_plot.png", plot = KM_sur_plot$plot, width = 9, height = 6, dpi = 300)

