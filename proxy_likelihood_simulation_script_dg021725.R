##########
# This script performs the simulation studies from the proxy-lhood chapter in
#     Wixson's Dissertation
##########

library(tidyverse)
source("./proxy_likelihood_fitting_scripts_dg021725.R")

#####
# I took the data from the second simulation for the plots as this was how I 
#     checked that the loop was working and I wanted to make the plots at the 
#     same time. 
# 
# For the full simulation study I changed num_sims to 100 and changed other 
#     relevant variables in lines 26-34
#
# Each simulation needs different parameter values 
# 
# The user will need to manually adjust which order MA models are plotted 
# 
# For the MA(15) from the model fitted to present climate FWI you need to
#   1. load the innovations parameter values (commented in lines 50-53)
#   2. change the generating model below (commented line 76)

num_sims     <- 2
length_ts    <- 10000
gen_model    <- "ARMA11"
max_MA_order <- 20
plot_name <- 
  paste0("./tpd_plot_arma11sim_t", theta, "_p", phi, "_version.pdf")
results_name <- 
  paste0("./results_diss_arma11sim_t", theta, "_p", phi, "_version.pdf")
plot_max_lag <- 15

###
# ARMA(1,1)
###

set.seed(18723872)
theta <- 0.9
phi <- 0.6
phi_text <- paste("\u03D5=", phi)
theta_text <- paste("\u03B8=", theta)
gen_model_name <- paste0("Model TL-ARMA(", phi_text, ", ", theta_text, ")")

###
# MA(15)
###

# set.seed(18723872)
# innovs <- readRDS("./ERA5_CO_innov_thetas.rds")
# innovs <- innovs$theta_late[15, 1:15]
# gen_model_name <- "Model TL-MA(15)"

###
# simulations
results_arma11 <- tibble(chosen_model_a = rep("?", num_sims), 
                         correct_a = rep(FALSE, num_sims),
                         chosen_model_b = rep("?", num_sims), 
                         correct_b = rep(FALSE, num_sims), 
                         approx_a = rep(FALSE, num_sims), 
                         approx_b = rep(FALSE, num_sims), 
                         theta_lower = rep(0, num_sims), 
                         theta_upper = rep(0, num_sims),
                         theta_in = rep(FALSE, num_sims),
                         phi_lower = rep(0, num_sims), 
                         phi_upper = rep(0, num_sims), 
                         phi_in = rep(FALSE, num_sims))

scores_AIC <- scores_BIC <- matrix(0, ncol = max_MA_order+2, nrow = num_sims)

for(i in 1:num_sims){
  tryCatch(
    {
      temp_data <- gen_arma11(n = length_ts, phi = phi, theta = theta)
      # temp_data <- gen_maq(n = length_ts, thetas = innovs)
      temp_pairs <- get_all_large_pairs(temp_data, thresh = 0.95)
      arma11_fit <- optim_ARMA11_from_list(theta_phi = c(0.1, 0.2),
                                           points_list = temp_pairs$points)
      fitted_arma11[i, ] <- arma11_fit$par
      ar1_fit <- optim_AR1_from_list(phi = 0.1, points_list = temp_pairs$points)
      fitted_ar1[i] <- ar1_fit$par
      ma_fits <- optim_MAq_sequentially(max_q = max_MA_order, 
                                        points_list = temp_pairs$points)
      
      fitted_models <- ma_fits
      fitted_models$AR1 <- ar1_fit
      fitted_models$ARMA11 <- arma11_fit
      
      scores_temp <- get_scores(library_fits = fitted_models, ts_list = temp_pairs, 
                                split_length = 500, ts_length = length_ts)
      scores_AIC[i,] <- scores_temp$CLAIC
      scores_BIC[i,] <- scores_temp$CLBIC
      
      results_arma11$chosen_model_a[i] <- scores_temp$Min_CLAIC
      if(scores_temp$Min_CLAIC == gen_model){
        results_arma11$correct_a[i] <- TRUE
        results_arma11$approx_a[i]  <- TRUE 
      } else if(
        abs(scores_temp$CLAIC[which(scores_temp$CLAIC == 
                                    min(na.omit(scores_temp$CLAIC)))] - 
            scores_temp$CLAIC[7]) < 2){
        results_arma11$approx_a[i]  <- TRUE 
      }
      results_arma11$chosen_model_b[i] <- scores_temp$Min_CLBIC
      if(scores_temp$Min_CLBIC == gen_model){
        results_arma11$correct_b[i] <- TRUE
        results_arma11$approx_b[i]  <- TRUE 
      } else if(
        abs(scores_temp$CLBIC[which(scores_temp$CLBIC == 
                                    min(na.omit(scores_temp$CLBIC)))] - 
            scores_temp$CLBIC[7]) < 2){
        results_arma11$approx_b[i]  <- TRUE 
      }
      
      results_arma11$theta_lower[i] = cis[1,1]
      results_arma11$theta_upper[i] = cis[1,2]
      results_arma11$theta_in[i] = (theta < cis[1,2]) && (theta > cis[1,1])
      results_arma11$phi_lower[i] =  cis[2,1]
      results_arma11$phi_upper[i] = cis[2,2]
      results_arma11$phi_in[i] = (phi < cis[2,2]) && (phi > cis[2,1])
      
    }, error=function(e){cat("ERROR: optimization error, skipping iteration \n")}
  )
  
  if(i %% 10 == 0){ 
    print(paste0("Simulation ", i, " complete"))}
}

###
# get plots
###.

tpdf_tibble <- tibble(lag = 0:plot_max_lag,
                      tpd = c(1, arma_11_tpdf(params = c(theta, phi),
                                              max_lag. = plot_max_lag)),
                      tpd_est_arma11 = 
                        c(1, arma_11_tpdf(params = arma11_fit$par,
                                          max_lag. = plot_max_lag)),
                      tpd_est_ar1 = 
                        c(1, ar_1_tpdf(ar1_fit$par, max_lag. = plot_max_lag)),
                      tpd_est_ma5 = 
                        c(1, ma_q_tpdf(ma_fits$MA5$par,
                                       max_lag. = plot_max_lag)),
                      tpd_est_ma10 = 
                        c(1, ma_q_tpdf(ma_fits$MA10$par,
                                       max_lag. = plot_max_lag)),
                      tpd_est_ma15 = 
                        c(1, ma_q_tpdf(ma_fits$MA15$par,
                                       max_lag. = plot_max_lag)),
                      tpd_est_ma20 = 
                        c(1, ma_q_tpdf(ma_fits$MA20$par,
                                       max_lag. = plot_max_lag)),
                      emp_tpd = TPDF(temp_data, maxlag = plot_max_lag))

p2 <- ggplot(data = tpdf_tibble, aes(x = lag, y = tpd)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0), lwd = 0.75, col = "black") +
  geom_segment(mapping = aes(x = lag  + 0.1, y = emp_tpd,
                             xend = lag + 0.1, yend = 0),
               col = "orange", lwd= 1) +
  geom_segment(mapping = aes(x = lag  + 0.2, y = tpd_est_arma11,
                             xend = lag + 0.2, yend = 0),
               col = "#009E73", lwd= 1) +
  geom_segment(mapping = aes(x = lag  + 0.3, y = tpd_est_ar1,
                             xend = lag + 0.3, yend = 0),
               col = "#CC79A7", lwd= 1) +
  geom_segment(mapping = aes(x = lag  + 0.4, y = tpd_est_ma10,
                             xend = lag + 0.4, yend = 0),
               col = "#0072B2", lwd= 1) +
  # geom_segment(mapping = aes(x = lag  + 0.5, y = tpd_est_ma20,
  #                            xend = lag + 0.5, yend = 0),
  #              col = "grey62", lwd = 0.75) +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 2,
           label= gen_model_name) +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 4,
           label = "Empirical", col = "orange") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 6,
           label = "fitted ARMA(1,1)", color = "#009E73") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 8,
           label = "fitted AR(1)", color = "#CC79A7") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 10,
           label = "fitted MA(10)", color = "#0072B2") +
  # annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 12,
  #          label = "fitted MA(20)", color = "grey62") +
  labs(x = "Lag", y = "TPD") +
  scale_x_continuous(breaks = 0:15) +
  theme_minimal()
p2
ggsave(p2, filename = plot_name, device = cairo_pdf,
       width = 7, height = 2.25, units = "in")

#####
# Results
##### 
results_list <- list(results = results_arma11, 
                     scores_AIC = scores_AIC, 
                     scores_BIC = scores_BIC)
resname <- paste0("./results_diss_arma11sim_t", theta, "_p", phi, "_version.pdf")
saveRDS(results_list, results_name)

