############ 
# This script performs the simulations in the model selection section 
#### 

set.seed(287) 
length_ts    <- 10000
source("./proxy_likelihood_fitting.R")
library(tidyverse)

# To demonstrate that two different models can have similar TPDs

# 1: Generate example ts
par_sim <- c(0.5, 0.1) # theta=0.5, phi=0.1
temp_data  <- gen_arma11(n = length_ts, theta = par_sim[1], phi = par_sim[2]) 
temp_pairs <- get_all_large_pairs(temp_data, thresh = 0.95) # get pairs for lhood

# 2: fit with proxy lhood
out_fitted <- optim_ARMA11_from_list(theta_phi = c(0.1, 0.1), 
                                     points_list = temp_pairs$points)
out_fitted$par
min(eigen(out_fitted$hessian)$value)

# 3: compute tpdf for fitted lambda values
out_fitted$tpdf <- arma_11_tpdf(out_fitted$par)
est_tpdf        <- TPDF(temp_data, maxlag = 20)[-1]

# 4: get innovations fit for comparison
plot_lag <- 5
out_fitted$innov <- innovations(c(1, est_tpdf), max_q = plot_lag)
round(out_fitted$innov[[1]][plot_lag, ], 3)
out_fitted$innov$tpd <- ma_q_tpdf(out_fitted$innov[[1]][2, 1:plot_lag])

# 6: Create a nicer looking plot for the paper
tpdf_tibble_1  <- tibble(lag = 0:plot_lag,
                         tpd = c(1, arma_11_tpdf(par_sim)[1:plot_lag]),
                         tpd_fitted = c(1, out_fitted$tpdf[1:plot_lag]),
                         tpd_est_natural = c(1, est_tpdf[1:plot_lag]), 
                         tpd_innov = c(1, out_fitted$innov$tpd[1:plot_lag]))

p_1 <- ggplot(data = tpdf_tibble_1, aes(x = lag, y = tpd)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0), lwd = 0.75, col = "black") +
  geom_segment(mapping = aes(x = lag  + 0.15, y = tpd_est_natural,
                             xend = lag + 0.15, yend = 0),
               col = "orange", lwd= 1) +
  geom_segment(mapping = aes(x = lag  + 0.3, y = tpd_fitted,
                             xend = lag + 0.3, yend = 0),
               col = "#009E73", lwd= 1) +
  geom_segment(mapping = aes(x = lag  + 0.45, y = tpd_innov,
                             xend = lag + 0.45, yend = 0),
               col = "#CC79A7", lwd= 1) +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 2,
           label= "Model TL-ARMA(\u03D5=0.1, \u03B8=0.5)") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 4,
           label = "Natural est.", col = "orange") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 6,
           label = "Fitted TL-ARMA(1,1)", color = "#009E73") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 8,
           label = "Fitted TL-MA(2)", color = "#CC79A7") +
  labs(x = "Lag", y = "TPD") +
  scale_x_continuous(breaks = 0:20) +
  theme_minimal()
p_1
ggsave("./proxy_paper/fitted_tpd_comparison_t01p05_081826.png", plot = p_1)

# 7: model selection
fitted_models <- optim_MAq_multiple(theta_init = 0.1, max_q = 2, 
                                    points_list = temp_pairs$points)
fitted_models$ARMA11 <- out_fitted
fitted_models <- fitted_models[-1]
scores_temp <- get_scores(library_fits = fitted_models, ts_list = temp_pairs, 
                          split_length = 500, ts_length = length_ts)
scores_temp$Min_CLAIC
# [1] "ARMA11"
scores_temp$CLAIC[1] - scores_temp$CLAIC[2]
# [1] 5.207453


# To test this first example 100 times
set.seed(287)
best_CLAIC <- best_CLBIC <- score_diff <- rep(NA, 100)
for(i in 1:100){
  # generate time series
  temp_data  <- gen_arma11(n = length_ts, theta = par_sim[1], phi = par_sim[2]) 
  # get pairs for lhood
  temp_pairs <- get_all_large_pairs(temp_data, thresh = 0.95) 
  # fit with proxy lhood
  out_fitted <- optim_ARMA11_from_list(theta_phi = c(0.1, 0.1), 
                                       points_list = temp_pairs$points)
  fitted_models <- optim_MAq_multiple(theta_init = 0.1, max_q = 2, 
                                      points_list = temp_pairs$points)
  fitted_models$ARMA11 <- out_fitted
  # remove MA1
  fitted_models <- fitted_models[-1]
  # compute scores and select best model
  scores_temp <- get_scores(library_fits = fitted_models, ts_list = temp_pairs, 
                            split_length = 500, ts_length = length_ts)
  best_CLAIC[i] <- scores_temp$Min_CLAIC
  best_CLBIC[i] <- scores_temp$Min_CLBIC
  score_diff[i] <- scores_temp$CLAIC[1] - scores_temp$CLAIC[2]
  if(i%%10 == 0){
    print(paste0("Compeleted sim ", i, " at ", Sys.time()))
  }
}


arma0501_results <- list(best_CLAIC = best_CLAIC, 
                         best_CLBIC = best_CLBIC, 
                         score_diff = score_diff)
saveRDS(arma0501_results, "./simulation_results_arma0906_081816.rds")




# To test a model with stronger, longer dependence
set.seed(876)
round(arma_11_tpdf(c(0.5, 0.1))[c(3, 11)], 3)
#[1] 0.005 0.000
round(arma_11_tpdf(c(0.9, 0.6))[c(3, 11)], 3)
#[1] 0.288 0.005

# 1: Generate example ts
par_sim <- c(0.9, 0.6) # theta=0.9, phi=0.6
temp_data  <- gen_arma11(n = length_ts, theta = par_sim[1], phi = par_sim[2]) 
temp_pairs <- get_all_large_pairs(temp_data, thresh = 0.95) # get pairs for lhood

# 2: fit with proxy lhood
out_fitted <- optim_ARMA11_from_list(theta_phi = c(0.1, 0.1), 
                                     points_list = temp_pairs$points)
out_fitted$par
min(eigen(out_fitted$hessian)$value)

# 3: compute tpdf for fitted lambda values
out_fitted$tpdf <- arma_11_tpdf(out_fitted$par)
est_tpdf        <- TPDF(temp_data, maxlag = 20)[-1]

# 4: get competing fit for comparison
fitted_models <- optim_MAq_multiple(theta_init = 0.1, 
                                        points_list = temp_pairs$points)
tpds <- lapply(fitted_models, function(x){ma_q_tpdf(x$par)})
for(i in 1:20){fitted_models[[i]]$tpd = tpds[[i]]}
fitted_models$AR1 <- optim_AR1_from_list(phi = 0.1, 
                                         points_list = temp_pairs$points)
fitted_models$AR1$tpd <- ar_1_tpdf(fitted_models$AR1$par)
fitted_models$ARMA11 <- out_fitted

# 6: Create a nicer looking plot for the paper
plot_lag <- 15
tpdf_tibble_1  <- tibble(lag = 0:plot_lag,
                         tpd = c(1, arma_11_tpdf(par_sim)[1:plot_lag]),
                         tpd_fitted = c(1, fitted_models$ARMA11$tpdf[1:plot_lag]),
                         tpd_est_natural = c(1, est_tpdf[1:plot_lag]),
                         tpd_ar1 = c(1, fitted_models$AR1$tpd[1:plot_lag]),
                         tpd_ma15 = c(1, fitted_models$MA15$tpd[1:plot_lag]))

p_2 <- ggplot(data = tpdf_tibble_1, aes(x = lag, y = tpd)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0), lwd = 0.75, col = "black") +
  geom_segment(mapping = aes(x = lag  + 0.15, y = tpd_est_natural,
                             xend = lag + 0.15, yend = 0),
               col = "orange", lwd= 1) +
  geom_segment(mapping = aes(x = lag  + 0.3, y = tpd_fitted,
                             xend = lag + 0.3, yend = 0),
               col = "#009E73", lwd= 1) +
  geom_segment(mapping = aes(x = lag  + 0.45, y = tpd_ar1,
                             xend = lag + 0.45, yend = 0),
               col = "#CC79A7", lwd= 1) +
  geom_segment(mapping = aes(x = lag  + 0.6, y = tpd_ma15,
                             xend = lag + 0.6, yend = 0),
               col = "#0072B2", lwd= 1) +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 2,
           label= "Model TL-ARMA(\u03D5=0.6, \u03B8=0.9)") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 4,
           label = "Natural est.", col = "orange") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 6,
           label = "Fitted TL-ARMA(1,1)", color = "#009E73") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 8,
           label = "Fitted TL-AR(1)", color = "#CC79A7") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 10,
           label = "Fitted TL-MA(15)", color = "#0072B2") +
  labs(x = "Lag", y = "TPD") +
  scale_x_continuous(breaks = 0:20) +
  theme_minimal()
p_2
ggsave("./proxy_paper/fitted_tpd_comparison_t019p06_081826.png", plot = p_2)

# 7: model selection
scores_temp <- get_scores(library_fits = fitted_models, ts_list = temp_pairs, 
                          split_length = 500, ts_length = length_ts)
scores_temp$AIC <- 2*scores_temp$lhood + 2*c(1:20, 1, 2)

scores_temp$Min_CLAIC
# [1] "ARMA11"



# To test this second example 100 times
set.seed(287)
best_CLAIC <- best_CLBIC <- best_CLAIC2 <- best_AIC <- best_AIC2 <- 
  score_diffA <- score_diffB <- rep(NA, 100)
for(i in 1:100){
  # generate time series
  temp_data  <- gen_arma11(n = length_ts, theta = par_sim[1], phi = par_sim[2]) 
  # get pairs for lhood
  temp_pairs <- get_all_large_pairs(temp_data, thresh = 0.95) 
  # fit with proxy lhood
  out_fitted <- optim_ARMA11_from_list(theta_phi = c(0.1, 0.1), 
                                       points_list = temp_pairs$points)
  fitted_models <- optim_MAq_multiple(theta_init = 0.1, max_q = 15, 
                                      points_list = temp_pairs$points)
  fitted_models$AR1 <- optim_AR1_from_list(phi = 0.1, 
                                           points_list = temp_pairs$points)
  fitted_models$ARMA11 <- out_fitted
  
  # compute scores and select best model
  scores_temp <- get_scores(library_fits = fitted_models, ts_list = temp_pairs, 
                            split_length = 500, ts_length = length_ts)
  scores_temp$AIC <- 2*scores_temp$lhood + 2*c(1:15, 1, 2)
  best_CLAIC[i]  <- scores_temp$Min_CLAIC
  best_CLBIC[i]  <- scores_temp$Min_CLBIC
  name_ind <- which(names(scores_temp$lhood) == "ARMA11")
  score_diffA[i] <- scores_temp$CLAIC[name_ind] - min(na.omit(scores_temp$CLAIC))
  score_diffB[i] <- scores_temp$CLBIC[name_ind] - min(na.omit(scores_temp$CLBIC))
  best_CLAIC2[i] <- ifelse(score_diffA[i] < 2, "ARMA11", 
                           names(scores_temp$lhood[which.min(scores_temp$CLAIC)]))
  best_AIC[i]    <- names(scores_temp$lhood)[which.min(scores_temp$AIC)]
  best_AIC2[i]   <- ifelse(unname(scores_temp$AIC[name_ind]) - 
                             min(scores_temp$AIC) < 2, 
                           "ARMA11",   best_AIC[i])
  if(i%%10 == 0){
    print(paste0("Compeleted sim ", i, " at ", Sys.time()))
  }
}

arma0906_results <- list(best_CLAIC = best_CLAIC, 
                           best_CLAIC2 = best_CLAIC2, 
                           best_CLBIC = best_CLBIC, 
                           best_AIC = best_AIC, 
                           best_AIC2 = best_AIC2, 
                           score_diffA = score_diffA, 
                           score_diffB = score_diffB)
saveRDS(arma0906_results, "./simulation_results_arma0906_081816.rds")






# To test a model with more complex dependence
set.seed(49291)
innovs <- readRDS("./ERA5_CO_innov_thetas.rds")
innovs <- innovs$theta_late[15, 1:15]
round(arma_11_tpdf(c(0.5, 0.1))[c(3, 11)], 3)
#[1] 0.005 0.000
round(arma_11_tpdf(c(0.9, 0.6))[c(3, 11)], 3)
#[1] 0.288 0.005
round(ma_q_tpdf(innovs)[c(3, 11)], 3)
#[1] 0.320 0.161

# 1: Generate example ts
par_sim <- innovs #thetas from wixson and cooley 2023
temp_data  <- gen_maq(n = length_ts, thetas = par_sim) 
temp_pairs <- get_all_large_pairs(temp_data, thresh = 0.95) # get pairs for lhood

# 2: fit with proxy lhood
out_ma   <- optim_MAq_multiple(theta_init = 0.1, 
                                   points_list = temp_pairs$points)
out_arma <- optim_ARMA11_from_list(theta_phi = c(0.1, 0.1), 
                                   points_list = temp_pairs$points)
out_ar   <- optim_AR1_from_list(phi = 0.1, points_list = temp_pairs$points)
out_ma$MA15$par
min(eigen(out_ma$MA15$hessian)$value)

# 3: compute tpdf for fitted lambda values
out_ma$MA15$tpdf <- ma_q_tpdf(out_ma$MA15$par)
est_tpdf         <- TPDF(temp_data, maxlag = 20)[-1]

# 4: get competing fit for comparison
fitted_models        <- out_ma

tpds <- lapply(fitted_models, function(x){ma_q_tpdf(x$par)})
for(i in 1:20){fitted_models[[i]]$tpd = tpds[[i]]}
fitted_models$ARMA11     <- out_arma
fitted_models$AR1        <- out_ar
fitted_models$AR1$tpd    <- ar_1_tpdf(fitted_models$AR1$par)
fitted_models$ARMA11$tpd <- arma_11_tpdf(fitted_models$ARMA11$par)


# 6: Create a nicer looking plot for the paper
plot_lag <- 20
tpdf_tibble_3  <- tibble(lag = 0:plot_lag,
                         tpd = c(1, ma_q_tpdf(par_sim)[1:plot_lag]),
                         tpd_fitted = c(1, fitted_models$MA15$tpd[1:plot_lag]),
                         tpd_arma11 = c(1, fitted_models$ARMA11$tpd[1:plot_lag]),
                         tpd_est_natural = c(1, est_tpdf[1:plot_lag]),
                         tpd_ar1 = c(1, fitted_models$AR1$tpd[1:plot_lag]))

p_3 <- ggplot(data = tpdf_tibble_3, aes(x = lag, y = tpd)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0), lwd = 0.75, col = "black") +
  geom_segment(mapping = aes(x = lag  + 0.15, y = tpd_est_natural,
                             xend = lag + 0.15, yend = 0),
               col = "orange", lwd= 1) +
  geom_segment(mapping = aes(x = lag  + 0.3, y = tpd_arma11,
                             xend = lag + 0.3, yend = 0),
               col = "#009E73", lwd= 1) +
  geom_segment(mapping = aes(x = lag  + 0.45, y = tpd_ar1,
                             xend = lag + 0.45, yend = 0),
               col = "#CC79A7", lwd= 1) +
  geom_segment(mapping = aes(x = lag  + 0.6, y = tpd_fitted,
                             xend = lag + 0.6, yend = 0),
               col = "#0072B2", lwd= 1) +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 2,
           label= "Model TL-MA(15)") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 4,
           label = "Natural est.", col = "orange") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 6,
           label = "Fitted TL-ARMA(1,1)", color = "#009E73") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 8,
           label = "Fitted TL-AR(1)", color = "#CC79A7") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 10,
           label = "Fitted TL-MA(15)", color = "#0072B2") +
  labs(x = "Lag", y = "TPD") +
  scale_x_continuous(breaks = 0:20) +
  theme_minimal()
p_3
ggsave("./proxy_paper/fitted_tpd_comparison_wildfireInnovThetas_081826.png", 
       plot = p_1)

# 7: model selection
scores_temp <- get_scores(library_fits = fitted_models, ts_list = temp_pairs, 
                          split_length = 500, ts_length = length_ts)
scores_temp$AIC <- 2*scores_temp$lhood + 2*c(1:20, 1, 2)
scores_temp$Min_CLAIC
# [1] "MA15"
names(scores_temp$lhood)[which.min(scores_temp$AIC)]
# [1] "MA15"


# To test this second example 100 times
set.seed(88324)
best_CLAIC <- best_CLBIC <- best_CLAIC2 <- best_AIC <- best_AIC2 <- 
  score_diffA <- score_diffB <- rep(NA, 100)
for(i in 1:100){
  # generate time series
  temp_data  <- gen_maq(n = length_ts, thetas = par_sim) 
  # get pairs for lhood
  temp_pairs <- get_all_large_pairs(temp_data, thresh = 0.95) 
  # fit with proxy lhood
  out_fitted <- optim_ARMA11_from_list(theta_phi = c(0.1, 0.1), 
                                       points_list = temp_pairs$points)
  fitted_models <- optim_MAq_multiple(theta_init = 0.1, 
                                      points_list = temp_pairs$points)
  fitted_models$AR1 <- optim_AR1_from_list(phi = 0.1, 
                                           points_list = temp_pairs$points)
  fitted_models$ARMA11 <- out_fitted
  
  # compute scores and select best model
  scores_temp <- get_scores(library_fits = fitted_models, ts_list = temp_pairs, 
                            split_length = 500, ts_length = length_ts)
  scores_temp$AIC <- 2*scores_temp$lhood + 2*c(1:20, 1, 2)
  best_CLAIC[i]  <- scores_temp$Min_CLAIC
  best_CLBIC[i]  <- scores_temp$Min_CLBIC
  name_ind <- which(names(scores_temp$lhood) == "MA15")
  score_diffA[i] <- scores_temp$CLAIC[name_ind] - min(na.omit(scores_temp$CLAIC))
  score_diffB[i] <- scores_temp$CLBIC[name_ind] - min(na.omit(scores_temp$CLBIC))
  best_CLAIC2[i] <- ifelse(score_diffA[i] < 2, "MA15", 
                           names(scores_temp$lhood[which.min(scores_temp$CLAIC)]))
  best_AIC[i]    <- names(scores_temp$lhood)[which.min(scores_temp$AIC)]
  best_AIC2[i]   <- ifelse(unname(scores_temp$AIC[name_ind]) - 
                             min(scores_temp$AIC) < 2, 
                           "MA15",   best_AIC[i])
  if(i%%10 == 0){
    print(paste0("Compeleted sim ", i, " at ", Sys.time()))
  }
}

ma15_innov_results <- list(best_CLAIC = best_CLAIC, 
                           best_CLAIC2 = best_CLAIC2, 
                           best_CLBIC = best_CLBIC, 
                           best_AIC = best_AIC, 
                           best_AIC2 = best_AIC2, 
                           score_diffA = score_diffA, 
                           score_diffB = score_diffB)
saveRDS(ma15_innov_results, "./simulation_results_ma15innov_081816.rds")

       