################################################################################
######################### f(t;X) = (t, t^2)^\top ###############################
################################################################################
library(ggplot2)
library(ggpubr)
library(numDeriv)

source("getProb_sequent.R")
source("MLE_point.R")
source("CallMLE.R")
source("WrapResults.R")
source("DR_point_var.R")
source("CallDR.R")
source("DR_var_sandwich.R")

run_dr_sand <- function(a, fun_tx, x_b_nuis, x_c_nuis, y, delta_c, res_mle,
                        max.step = 1000, thres = 1e-6) {
  p_alpha <- ncol(fun_tx[[1]])
  p_beta <- ncol(x_b_nuis)
  
  alpha_nuis <- res_mle$point.est[seq_len(p_alpha)]
  beta_nuis <- res_mle$point.est[p_alpha + seq_len(p_beta)]
  
  res_ps <- fit_ps(a = a, x_c_nuis = x_c_nuis)
  gamma_nuis <- res_ps$point.est
  
  res_dr <- DREst(a = a, fun_tx = fun_tx,
                  x_b_nuis = x_b_nuis, x_c_nuis = x_c_nuis,
                  y = y, delta_c = delta_c,
                  alpha.nuis = alpha_nuis,
                  beta.nuis = beta_nuis,
                  gamma.nuis = gamma_nuis,
                  max.step = max.step, thres = thres,
                  alpha.start = alpha_nuis)
  
  res_hess <- hessian.DR(a = a, fun_tx = fun_tx,
                         x_b_nuis = x_b_nuis, x_c_nuis = x_c_nuis,
                         y = y, delta_c = delta_c,
                         alpha.dr = res_dr$dr.point,
                         alpha.nuis = alpha_nuis,
                         beta.nuis = beta_nuis,
                         gamma.nuis = gamma_nuis)
  
  res_sand <- var.DR.sandwich(dr.point = res_dr$dr.point,
                              S_mat_psc = res_hess$S_mat_psc,
                              S_mat_lop = res_mle$S_mat,
                              S_mat_dr = res_dr$S_mat,
                              cov.psc = res_hess$cov.psc,
                              cov.lop = res_mle$cov,
                              dU.by.dalpha = res_hess$dU.by.dalpha,
                              dU.by.dab = res_hess$dU.by.dab,
                              dU.by.dgamma = res_hess$dU.by.dgamma)
  
  dr_table <- rbind(res_dr$dr.point,
                    res_sand$sd.est,
                    res_sand$p.value,
                    res_sand$conf.lower,
                    res_sand$conf.upper)
  rownames(dr_table) <- c("point.est", "se.est", "p.value", "conf.lower", "conf.upper")
  colnames(dr_table) <- names(res_dr$dr.point)
  
  list(dr = res_dr, ps = res_ps, hessian = res_hess,
       sand = res_sand, table = dr_table)
}

make_fun_tx_t_t2 <- function(time_grid, n) {
  lapply(time_grid, function(t_i) {
    cbind(rep(t_i, n), rep(t_i^2, n))
  })
}

################## placebo vs. 0.2 mg estrogen #################################
# Variables
prost1 = prost[prost$rx %in% levels(prost$rx)[c(1,2)],]
A1 = as.numeric(prost1$rx == "0.2 mg estrogen")
# table(A1) # 0 = "placebo" ; 1 = "0.2 mg estrogen"
Y1 = prost1$dtime
# summary(Y1)
# tab1 = table(A1, prost1$eventType)
# print(tab1)
# print(prop.table(tab1, 1))

Delta1 = 1 - as.numeric(prost1$eventType == "alive")
# table(Delta1)

# Covariates
age1 = scale(prost1$age, scale = FALSE)
pf1 = as.numeric(prost1$pf == "normal activity")
hx1 = prost1$hx
hg1 = scale(prost1$hg, scale = FALSE)
n1 = length(A1)
X1 = cbind(rep(1,n1),age1, pf1, hx1, hg1)
# dim(X1)

## Fit model
timepoints = (0:6)

m = length(timepoints) - 1
Y_re1 = Y1
Y_re1[(Y1 <= timepoints[2]) & (Y1 >= timepoints[1])] = timepoints[2]
for (i in 1:m){
  Y_re1[(Y1 <= timepoints[i+1]) & (Y1 > timepoints[i])] = timepoints[i+1]
}

Y_re1[Y1 > timepoints[m+1]] = timepoints[m+1]
Delta1[Y1 > timepoints[m+1]] = 0
# table(Y_re1)

### MLE
timepoints_fit1 = sort(unique(c(0, Y_re1)))
fun_tx1 = make_fun_tx_t_t2(timepoints_fit1, length(A1))
res1 = MLEst(a = A1, fun_tx = fun_tx1, x_b = X1, 
             y = Y_re1, delta_c = Delta1, 
             alpha.start = c(0,0), beta.start = rep(0,dim(X1)[2]),
             max.step = 1000, thres = 1e-6, method = 1)


res_table1 = rbind(res1$point.est, res1$se.est, res1$p.value,res1$conf.lower,res1$conf.upper)
rownames(res_table1) = c("point.est", "se.est", "p.value", "conf.lower", "conf.upper")
print(res_table1)


### DR
dr1_out <- run_dr_sand(a = A1, fun_tx = fun_tx1,
                       x_b_nuis = X1, x_c_nuis = X1,
                       y = Y_re1, delta_c = Delta1,
                       res_mle = res1)
dr_table1 <- dr1_out$table
print(dr_table1)


################## placebo vs. 1.0 mg estrogen #################################
# Variables
prost2 = prost[prost$rx %in% levels(prost$rx)[c(1,3)],]
A2 = as.numeric(prost2$rx == "1.0 mg estrogen")
# table(A2) # 0 = "placebo" ; 1 = "1.0 mg estrogen"

# Observed time and censoring indicator
Y2 = prost2$dtime
# summary(Y2)
# tab2 = table(A2, prost2$eventType)
# print(tab2)
# print(prop.table(tab2, 1))

Delta2 = 1 - as.numeric(prost2$eventType == "alive")
# table(Delta2)


# Covariates
age2 = scale(prost2$age, scale = FALSE)
pf2 = as.numeric(prost2$pf == "normal activity")
hx2 = prost2$hx
hg2 = scale(prost2$hg, scale = FALSE)
n2 = length(A2)
X2 = cbind(rep(1,n2),age2, pf2, hx2, hg2)
# dim(X2)


## Fit model

timepoints = (0:6)
m = length(timepoints) - 1
Y_re2 = Y2
Y_re2[(Y2 <= timepoints[2]) & (Y2 >= timepoints[1])] = timepoints[2]
for (i in 1:m){
  Y_re2[(Y2 <= timepoints[i+1]) & (Y2 > timepoints[i])] = timepoints[i+1]
}

Y_re2[Y2 > timepoints[m+1]] = timepoints[m+1]
Delta2[Y2 > timepoints[m+1]] = 0
# table(Y_re2)


### MLE
timepoints_fit2 = sort(unique(c(0, Y_re2)))
fun_tx2 = make_fun_tx_t_t2(timepoints_fit2, length(A2))


res2 = MLEst(a = A2, fun_tx = fun_tx2, x_b = X2,
             y = Y_re2, delta_c = Delta2,
             alpha.start = c(0,0), beta.start = rep(0,dim(X2)[2]),
             max.step = 1000, thres = 1e-6, method = 1)


res_table2 = rbind(res2$point.est,res2$se.est,res2$p.value,res2$conf.lower,res2$conf.upper)
rownames(res_table2) = c("point.est", "se.est", "p.value", "conf.lower", "conf.upper")
print(res_table2)

### DR
dr2_out <- run_dr_sand(a = A2, fun_tx = fun_tx2,
                       x_b_nuis = X2, x_c_nuis = X2,
                       y = Y_re2, delta_c = Delta2,
                       res_mle = res2)
dr_table2 <- dr2_out$table
print(dr_table2)




################## placebo vs. 5.0 mg estrogen #################################
# Variable
prost3 = prost[prost$rx %in% levels(prost$rx)[c(1,4)],]
A3 = as.numeric(prost3$rx == "5.0 mg estrogen")
# table(A3) # 0 = "placebo" ; 1 = "5.0 mg estrogen"

# Observed time and censoring indicator
Y3 = prost3$dtime
# summary(Y3)

# tab3 = table(A3, prost3$eventType)
# print(tab3)
# print(prop.table(tab3, 1))

Delta3 = 1 - as.numeric(prost3$eventType == "alive")
# table(Delta3)

# Covariates
age3 = scale(prost3$age, scale = FALSE)
pf3 = as.numeric(prost3$pf == "normal activity")
hx3 = prost3$hx
hg3 = scale(prost3$hg, scale = FALSE)
n3 = length(A3)
X3 = cbind(rep(1,n3),age3, pf3, hx3, hg3)
# dim(X3)

## Fit model

timepoints = (0:6)

m = length(timepoints) - 1
Y_re3 = Y3
Y_re3[(Y3 <= timepoints[2]) & (Y3 >= timepoints[1])] = timepoints[2]
for (i in 1:m){
  Y_re3[(Y3 <= timepoints[i+1]) & (Y3 > timepoints[i])] = timepoints[i+1]
}

Y_re3[Y3 > timepoints[m+1]] = timepoints[m+1]
Delta3[Y3 > timepoints[m+1]] = 0
# table(Y_re3)

### MLE
timepoints_fit3 = sort(unique(c(0, Y_re3)))
fun_tx3 = make_fun_tx_t_t2(timepoints_fit3, length(A3))


res3 = MLEst(a = A3, fun_tx = fun_tx3, x_b = X3, 
             y = Y_re3, delta_c = Delta3, 
             alpha.start = c(0,0), beta.start = rep(0,dim(X3)[2]),
             max.step = 1000, thres = 1e-6, method = 1)


res_table3 = rbind(res3$point.est,res3$se.est,res3$p.value,res3$conf.lower,res3$conf.upper)
rownames(res_table3) = c("point.est", "se.est", "p.value", "conf.lower", "conf.upper")
print(res_table3)

### DR
dr3_out <- run_dr_sand(a = A3, fun_tx = fun_tx3,
                       x_b_nuis = X3, x_c_nuis = X3,
                       y = Y_re3, delta_c = Delta3,
                       res_mle = res3)
dr_table3 <- dr3_out$table
print(dr_table3)

table_to_pdf_data <- function(res_table, contrast, estimator) {
  parameter <- colnames(res_table)
  if (is.null(parameter) || any(is.na(parameter)) || any(parameter == "")) {
    parameter <- paste0("param", seq_len(ncol(res_table)))
  }
  
  data.frame(
    Contrast = contrast,
    Estimator = estimator,
    Parameter = parameter,
    Point = round(as.numeric(res_table["point.est", ]), 4),
    SE = round(as.numeric(res_table["se.est", ]), 4),
    P_value = signif(as.numeric(res_table["p.value", ]), 3),
    CI_lower = round(as.numeric(res_table["conf.lower", ]), 4),
    CI_upper = round(as.numeric(res_table["conf.upper", ]), 4),
    row.names = NULL
  )
}

mle_pdf_table <- rbind(
  table_to_pdf_data(res_table1, "0.2 mg Estrogen vs Placebo", "MLE"),
  table_to_pdf_data(res_table2, "1.0 mg Estrogen vs Placebo", "MLE"),
  table_to_pdf_data(res_table3, "5.0 mg Estrogen vs Placebo", "MLE")
)

dr_pdf_table <- rbind(
  table_to_pdf_data(dr_table1, "0.2 mg Estrogen vs Placebo", "DR"),
  table_to_pdf_data(dr_table2, "1.0 mg Estrogen vs Placebo", "DR"),
  table_to_pdf_data(dr_table3, "5.0 mg Estrogen vs Placebo", "DR")
)

table_pdf_file <- "MLE_DR_tables_t_t2.pdf"
grDevices::pdf(table_pdf_file, width = 12, height = 8)
print(ggtexttable(mle_pdf_table, rows = NULL, theme = ttheme("light")))
print(ggtexttable(dr_pdf_table, rows = NULL, theme = ttheme("light")))
grDevices::dev.off()



################################## Plots #######################################


every_years = (0:72)/12
make_rr_data <- function(point, group, estimator) {
  alpha <- as.numeric(point[1:2])
  
  data.frame(
    time = every_years,
    RR = exp(alpha[1] * every_years + alpha[2] * every_years^2),
    group = group,
    estimator = estimator
  )
}

data_RR <- rbind(
  make_rr_data(res1$point.est,
               "0.2 mg Estrogen vs Placebo", "MLE"),
  make_rr_data(res2$point.est,
               "1.0 mg Estrogen vs Placebo", "MLE"),
  make_rr_data(res3$point.est,
               "5.0 mg Estrogen vs Placebo", "MLE"),
  make_rr_data(dr1_out$dr$dr.point,
               "0.2 mg Estrogen vs Placebo", "DR"),
  make_rr_data(dr2_out$dr$dr.point,
               "1.0 mg Estrogen vs Placebo", "DR"),
  make_rr_data(dr3_out$dr$dr.point,
               "5.0 mg Estrogen vs Placebo", "DR")
)


RR_compare_all_plot <- ggplot(
  data_RR,
  aes(
    x = time,
    y = RR,
    color = group,
    linetype = estimator,
    group = interaction(group, estimator)
  )
) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "grey45") +
  geom_line(linewidth = 1.1) +
  scale_color_manual(
    name = "Contrast",
    values = c(
      "0.2 mg Estrogen vs Placebo" = "#0072B2",
      "1.0 mg Estrogen vs Placebo" = "#D55E00",
      "5.0 mg Estrogen vs Placebo" = "#009E73"
    ),
    labels = c(
      "0.2 mg vs Placebo",
      "1.0 mg vs Placebo",
      "5.0 mg vs Placebo"
    )
  ) +
  scale_linetype_manual(
    name = "Estimator",
    values = c("MLE" = "solid", "DR" = "dashed")
  ) +
  scale_x_continuous(
    breaks = 0:6,
    limits = c(0, 6)
  ) +
  scale_y_continuous(
    breaks = seq(0.9, 1.9, by = 0.2),
    minor_breaks = seq(0.9, 1.9, by = 0.1)
  ) +
  coord_cartesian(ylim = c(0.85, 1.95)) +
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

RR_compare_all_plot

ggsave(
  "Figures/RR_compare_all_plot.png",
  plot = RR_compare_all_plot,
  width = 9,
  height = 6,
  dpi = 300
)
