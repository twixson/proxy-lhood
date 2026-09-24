############ 
# This script performs the simulations in the model fitting section
#### 

# To use the proxy-likelihood to fit models:
source("./proxy_likelihood_fitting.R")
library(tidyverse)
library(cowplot)
set.seed(873)
length_ts    <- 10000

# 0: get some data to estimate the TPD of: 
q_sim <- 5  # order of MA model to simulate from
thetas_sim <- runif(q_sim)
while (min(abs(polyroot(c(1,thetas_sim)))) < 1){ # check if valid model
  thetas_sim <- runif(q_sim)
}
temp_data  <- gen_maq(n = length_ts, thetas = thetas_sim) # generate ts
temp_pairs <- get_all_large_pairs(temp_data, thresh = 0.95) # get pairs for lhood

# 1: get initial values for optimization:
init_tpdf    <- TPDF(temp_data, maxlag = 20)[-1]
init_lambdas <- get_lambda(init_tpdf)

# 2: fit with proxy lhood
out_fitted <- optim_MAq_single(thetas = rep(0.1, 20), 
                               points_list = temp_pairs$points)
out_fitted$par
min(eigen(out_fitted$hessian)$value)

# 3: compute tpdf for fitted lambda values
out_fitted$tpdf <- ma_q_tpdf(out_fitted$par)

plot(0:20, c(1, ma_q_tpdf(thetas_sim)), type = "h", xlab = "lag", ylab = "TPD")
lines(0:20 + 0.1, c(1, out_fitted$tpdf), col = 2, type = "h")
lines(0:20 + 0.2, c(1, init_tpdf), col = 4, type = "h")
legend("topright", bty = "n", lty = 1, col = c(1, 2, 4),
       legend = c("model", "likelihood", "standard"), 
       text.col = c(1, 2, 4))

# 4: get innovations fit for comparison
out_fitted$innov <- innovations(c(1, init_tpdf), max_q = 20)
round(out_fitted$innov[[1]][20, ], 3)
lines(0:20 + 0.3, c(1, ma_q_tpdf(out_fitted$innov[[1]][20, 1:20])), col = 6, type = "h")

# 6: Create a nicer looking plot for the paper
plot_max_lag <- length(out_fitted$tpdf)
tpdf_tibble11  <- tibble(lag = 0:plot_max_lag,
                       tpd = c(1, ma_q_tpdf(thetas_sim)),
                       tpd_fitted = c(1, out_fitted$tpdf),
                       tpd_est_natural = c(1, init_tpdf), 
                       tpd_innov = c(1, ma_q_tpdf(out_fitted$innov[[1]][20, 1:20])))

p1 <- ggplot(data = tpdf_tibble11, aes(x = lag, y = tpd)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0), lwd = 0.75, col = "black") +
  geom_segment(mapping = aes(x = lag  + 0.15, y = tpd_fitted,
                             xend = lag + 0.15, yend = 0),
               col = "orange", lwd= 1) +
  geom_segment(mapping = aes(x = lag  + 0.3, y = tpd_est_natural,
                             xend = lag + 0.3, yend = 0),
               col = "#009E73", lwd= 1) +
  geom_segment(mapping = aes(x = lag  + 0.45, y = tpd_innov,
                             xend = lag + 0.45, yend = 0),
               col = "#CC79A7", lwd= 1) +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 2,
           label= "Model") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 4,
           label = "HR fit", col = "orange") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 6,
           label = "Empirical est.", color = "#009E73") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 8,
           label = "Innov. fit", color = "#CC79A7") +
  labs(x = "Lag", y = "TPD") +
  scale_x_continuous(breaks = 0:20) +
  theme_minimal()
p1





########
# Try again with different order
########


# 0: get some data to estimate the TPD of: 
q_sim <- 10  # order of MA model to simulate from
thetas_sim <- runif(q_sim)
while (min(abs(polyroot(c(1,thetas_sim)))) < 1){ # check if valid model
  thetas_sim <- runif(q_sim)
}
temp_data  <- gen_maq(n = length_ts, thetas = thetas_sim) # generate ts
temp_pairs <- get_all_large_pairs(temp_data, thresh = 0.95) # get pairs for lhood

# 1: get initial values for optimization:
init_tpdf    <- TPDF(temp_data, maxlag = 20)[-1]
init_lambdas <- get_lambda(init_tpdf)

# 2: fit with proxy lhood
out_fitted <- optim_MAq_single(thetas = rep(0.1, 20), 
                               points_list = temp_pairs$points)
out_fitted$par
min(eigen(out_fitted$hessian)$value)

# 3: compute tpdf for fitted lambda values
out_fitted$tpdf <- ma_q_tpdf(out_fitted$par)

plot(0:20, c(1, ma_q_tpdf(thetas_sim)), type = "h", xlab = "lag", ylab = "TPD")
lines(0:20 + 0.1, c(1, out_fitted$tpdf), col = 2, type = "h")
lines(0:20 + 0.2, c(1, init_tpdf), col = 4, type = "h")
legend("topright", bty = "n", lty = 1, col = c(1, 2, 4),
       legend = c("model", "likelihood", "standard"), 
       text.col = c(1, 2, 4))

# 4: get innovations fit for comparison
out_fitted$innov <- innovations(c(1, init_tpdf), max_q = 20)
round(out_fitted$innov[[1]][20, ], 3)
lines(0:20 + 0.3, c(1, ma_q_tpdf(out_fitted$innov[[1]][20, 1:20])), col = 6, type = "h")

# 6: Create a nicer looking plot for the paper
plot_max_lag <- length(out_fitted$tpdf)
tpdf_tibble22  <- tibble(lag = 0:plot_max_lag,
                        tpd = c(1, ma_q_tpdf(thetas_sim)),
                        tpd_fitted = c(1, out_fitted$tpdf),
                        tpd_est_natural = c(1, init_tpdf), 
                        tpd_innov = c(1, ma_q_tpdf(out_fitted$innov[[1]][20, 1:20])))

p2 <- ggplot(data = tpdf_tibble22, aes(x = lag, y = tpd)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0), lwd = 0.75, col = "black") +
  geom_segment(mapping = aes(x = lag  + 0.15, y = tpd_fitted,
                             xend = lag + 0.15, yend = 0),
               col = "orange", lwd= 1) +
  geom_segment(mapping = aes(x = lag  + 0.3, y = tpd_est_natural,
                             xend = lag + 0.3, yend = 0),
               col = "#009E73", lwd= 1) +
  geom_segment(mapping = aes(x = lag  + 0.45, y = tpd_innov,
                             xend = lag + 0.45, yend = 0),
               col = "#CC79A7", lwd= 1) +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 2,
           label= "Model") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 4,
           label = "HR fit", col = "orange") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 6,
           label = "Empirical est.", color = "#009E73") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 8,
           label = "Innov. fit", color = "#CC79A7") +
  labs(x = "Lag", y = "TPD") +
  scale_x_continuous(breaks = 0:20) +
  theme_minimal()
p2




########
# Try again with different order (again)
########


# 0: get some data to estimate the TPD of: 
q_sim <- 15  # order of MA model to simulate from
thetas_sim <- runif(q_sim)
while (min(abs(polyroot(c(1,thetas_sim)))) < 1){ # check if valid model
  thetas_sim <- runif(q_sim)
}
temp_data  <- gen_maq(n = length_ts, thetas = thetas_sim) # generate ts
temp_pairs <- get_all_large_pairs(temp_data, thresh = 0.95) # get pairs for lhood

# 1: get initial values for optimization:
init_tpdf    <- TPDF(temp_data, maxlag = 20)[-1]
init_lambdas <- get_lambda(init_tpdf)

# 2: fit with proxy lhood
out_fitted <- optim_MAq_single(thetas = rep(0.1, 20), 
                               points_list = temp_pairs$points)
out_fitted$par
min(eigen(out_fitted$hessian)$value)

# 3: compute tpdf for fitted lambda values
out_fitted$tpdf <- ma_q_tpdf(out_fitted$par)

plot(0:20, c(1, ma_q_tpdf(thetas_sim)), type = "h", xlab = "lag", ylab = "TPD")
lines(0:20 + 0.1, c(1, out_fitted$tpdf), col = 2, type = "h")
lines(0:20 + 0.2, c(1, init_tpdf), col = 4, type = "h")
legend("topright", bty = "n", lty = 1, col = c(1, 2, 4),
       legend = c("model", "likelihood", "standard"), 
       text.col = c(1, 2, 4))

# 4: get innovations fit for comparison
out_fitted$innov <- innovations(c(1, init_tpdf), max_q = 20)
round(out_fitted$innov[[1]][20, ], 3)
lines(0:20 + 0.3, c(1, ma_q_tpdf(out_fitted$innov[[1]][20, 1:20])), col = 6, type = "h")

# 6: Create a nicer looking plot for the paper
plot_max_lag <- length(out_fitted$tpdf)
tpdf_tibble33  <- tibble(lag = 0:plot_max_lag,
                         tpd = c(1, ma_q_tpdf(thetas_sim)),
                         tpd_fitted = c(1, out_fitted$tpdf),
                         tpd_est_natural = c(1, init_tpdf), 
                         tpd_innov = c(1, ma_q_tpdf(out_fitted$innov[[1]][20, 1:20])))

p3 <- ggplot(data = tpdf_tibble33, aes(x = lag, y = tpd)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0), lwd = 0.75, col = "black") +
  geom_segment(mapping = aes(x = lag  + 0.15, y = tpd_fitted,
                             xend = lag + 0.15, yend = 0),
               col = "orange", lwd= 1) +
  geom_segment(mapping = aes(x = lag  + 0.3, y = tpd_est_natural,
                             xend = lag + 0.3, yend = 0),
               col = "#009E73", lwd= 1) +
  geom_segment(mapping = aes(x = lag  + 0.45, y = tpd_innov,
                             xend = lag + 0.45, yend = 0),
               col = "#CC79A7", lwd= 1) +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 2,
           label= "Model") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 4,
           label = "HR fit", col = "orange") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 6,
           label = "Empirical est.", color = "#009E73") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 8,
           label = "Innov. fit", color = "#CC79A7") +
  labs(x = "Lag", y = "TPD") +
  scale_x_continuous(breaks = 0:20) +
  theme_minimal()
p3




########
# Try again with ARMA model
########


# 0: get some data to estimate the TPD of: 
par_sim <- runif(2) # generate ts
temp_data  <- gen_arma11(n = length_ts, theta = par_sim[1], phi = par_sim[2]) 
temp_pairs <- get_all_large_pairs(temp_data, thresh = 0.95) # get pairs for lhood

# 1: get initial values for optimization:
init_tpdf    <- TPDF(temp_data, maxlag = 20)[-1]
init_lambdas <- get_lambda(init_tpdf)

# 2: fit with proxy lhood
out_fitted <- optim_ARMA11_from_list(theta_phi = c(0.1, 0.1), 
                                     points_list = temp_pairs$points)
out_fitted$par
min(eigen(out_fitted$hessian)$value)

# 3: compute tpdf for fitted lambda values
out_fitted$tpdf <- arma_11_tpdf(out_fitted$par)

plot(0:20, c(1, arma_11_tpdf(par_sim)), type = "h", xlab = "lag", ylab = "TPD")
lines(0:20 + 0.1, c(1, out_fitted$tpdf), col = 2, type = "h")
lines(0:20 + 0.2, c(1, init_tpdf), col = 4, type = "h")
legend("topright", bty = "n", lty = 1, col = c(1, 2, 4),
       legend = c("model", "likelihood", "standard"), 
       text.col = c(1, 2, 4))

# 4: get innovations fit for comparison
out_fitted$innov <- innovations(c(1, init_tpdf), max_q = 20)
round(out_fitted$innov[[1]][20, ], 3)
lines(0:20 + 0.3, c(1, ma_q_tpdf(out_fitted$innov[[1]][20, 1:20])), col = 6, type = "h")

# 6: Create a nicer looking plot for the paper
plot_max_lag <- length(out_fitted$tpdf)
tpdf_tibble44  <- tibble(lag = 0:plot_max_lag,
                         tpd = c(1, arma_11_tpdf(par_sim)),
                         tpd_fitted = c(1, out_fitted$tpdf),
                         tpd_est_natural = c(1, init_tpdf), 
                         tpd_innov = c(1, ma_q_tpdf(out_fitted$innov[[1]][20, 1:20])))

p4 <- ggplot(data = tpdf_tibble44, aes(x = lag, y = tpd)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0), lwd = 0.75, col = "black") +
  geom_segment(mapping = aes(x = lag  + 0.15, y = tpd_fitted,
                             xend = lag + 0.15, yend = 0),
               col = "orange", lwd= 1) +
  geom_segment(mapping = aes(x = lag  + 0.3, y = tpd_est_natural,
                             xend = lag + 0.3, yend = 0),
               col = "#009E73", lwd= 1) +
  geom_segment(mapping = aes(x = lag  + 0.45, y = tpd_innov,
                             xend = lag + 0.45, yend = 0),
               col = "#CC79A7", lwd= 1) +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 2,
           label= "Model") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 4,
           label = "HR fit", col = "orange") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 6,
           label = "Empirical est.", color = "#009E73") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 8,
           label = "Innov. fit", color = "#CC79A7") +
  labs(x = "Lag", y = "TPD") +
  scale_x_continuous(breaks = 0:20) +
  theme_minimal()
p4



#######
# Combine plots for paper
#######
a <- cowplot::plot_grid(p1, p2, p3, p4)
save_plot("./plots/model_fitting_4model_plot.png", a, base_height = 3.4, base_asp = 3)








#########################################################################

#########################################################################

#########################################################################




######
# larger simulation to get idea of estimation accuracy and variablility
######
num_sims <- 100



errors_ma5 <- data.frame(hr_full_sq = rep(0, num_sims),
                         natural_full_sq = rep(0, num_sims),
                         hr_short_sq = rep(0, num_sims), 
                         natural_short_sq = rep(0, num_sims), 
                         hr_full_ab = rep(0, num_sims),
                         natural_full_ab = rep(0, num_sims),
                         hr_short_ab = rep(0, num_sims), 
                         natural_short_ab = rep(0, num_sims))
errors_ma10 <- data.frame(hr_full_sq = rep(0, num_sims),
                          natural_full_sq = rep(0, num_sims),
                          hr_short_sq = rep(0, num_sims), 
                          natural_short_sq = rep(0, num_sims), 
                          hr_full_ab = rep(0, num_sims),
                          natural_full_ab = rep(0, num_sims),
                          hr_short_ab = rep(0, num_sims), 
                          natural_short_ab = rep(0, num_sims))
errors_ma15 <- data.frame(hr_full_sq = rep(0, num_sims),
                          natural_full_sq = rep(0, num_sims),
                          hr_short_sq = rep(0, num_sims), 
                          natural_short_sq = rep(0, num_sims), 
                          hr_full_ab = rep(0, num_sims),
                          natural_full_ab = rep(0, num_sims),
                          hr_short_ab = rep(0, num_sims), 
                          natural_short_ab = rep(0, num_sims))
errors_arma11 <- data.frame(hr_full_sq = rep(0, num_sims),
                            ls_full_sq = rep(0, num_sims),
                            la_full_sq = rep(0, num_sims), 
                            hr_full_ab = rep(0, num_sims),
                            ls_full_ab = rep(0, num_sims),
                            la_full_ab = rep(0, num_sims))

set.seed(982)
for(i in 1:num_sims){
  # ma(5)
  q_sim <- 5
  thetas_sim <- runif(q_sim)
  while (min(abs(polyroot(c(1,thetas_sim)))) < 1){ # check if valid model
    thetas_sim <- runif(q_sim)
  }
  temp_data  <- gen_maq(n = length_ts, thetas = thetas_sim) # generate ts
  temp_pairs <- get_all_large_pairs(temp_data, thresh = 0.95) # get pairs for lhood
  init_tpdf  <- TPDF(temp_data, maxlag = 20)[-1]
  out_fitted <- optim_MAq_multiple(theta_init = 0.1, 
                                       points_list = temp_pairs$points)
  out_fitted$tpdf_short <- ma_q_tpdf(out_fitted[[q_sim]]$par)
  out_fitted$tpdf_long  <- ma_q_tpdf(out_fitted[[20]]$par)
  out_fitted$innov <- innovations(c(1, init_tpdf), max_q = 20)
  out_fitted$innov_short <- ma_q_tpdf(out_fitted$innov[[1]][q_sim, ])
  out_fitted$innov_long  <- ma_q_tpdf(out_fitted$innov[[1]][20, ])
  model_tpd        <- ma_q_tpdf(thetas_sim)
  errors_ma5[i, 1] <- sum((model_tpd - out_fitted$tpdf_long)^2)
  errors_ma5[i, 2] <- sum((model_tpd - out_fitted$innov_long)^2)
  errors_ma5[i, 3] <- sum((model_tpd - out_fitted$tpdf_short)[1:q_sim]^2)
  errors_ma5[i, 4] <- sum((model_tpd - out_fitted$innov_short)[1:q_sim]^2)
  errors_ma5[i, 5] <- sum(abs(model_tpd - out_fitted$tpdf_long))
  errors_ma5[i, 6] <- sum(abs(model_tpd - out_fitted$innov_long))
  errors_ma5[i, 7] <- sum(abs(model_tpd - out_fitted$tpdf_short)[1:q_sim])
  errors_ma5[i, 8] <- sum(abs(model_tpd - out_fitted$innov_short)[1:q_sim])
  
  if(i%%10 == 0){print(paste0("Iteration ", i, " complete"))}
}

get_ab_results <- function(x){
  list(mean_hr_short = mean(x$hr_short_ab), 
       sd_HR_short = sd(x$hr_short_ab),
       mean_innov_short = mean(x$natural_short_ab), 
       sd_innov_short = sd(x$natural_short_ab),
       num_HR_better_short = sum(x$hr_short_ab - x$natural_short_ab < 0),
       mean_hr_full = mean(x$hr_full_ab), 
       sd_HR_full = sd(x$hr_full_ab),
       mean_innov_full = mean(x$natural_full_ab), 
       sd_innov_full = sd(x$natural_full_ab), 
       num_HR_better_full = sum(x$hr_full_ab - x$natural_full_ab < 0))
}
results1 <- get_ab_results(errors_ma5)


set.seed(982)
for(i in 1:num_sims){
  # ma(10)
  q_sim <- 10
  thetas_sim <- runif(q_sim)
  while (min(abs(polyroot(c(1,thetas_sim)))) < 1){ # check if valid model
    thetas_sim <- runif(q_sim)
  }
  temp_data  <- gen_maq(n = length_ts, thetas = thetas_sim) # generate ts
  temp_pairs <- get_all_large_pairs(temp_data, thresh = 0.95) # get pairs for lhood
  init_tpdf  <- TPDF(temp_data, maxlag = 20)[-1]
  out_fitted <- optim_MAq_multiple(theta_init = 0.1, 
                                       points_list = temp_pairs$points)
  out_fitted$tpdf_short <- ma_q_tpdf(out_fitted[[q_sim]]$par)
  out_fitted$tpdf_long  <- ma_q_tpdf(out_fitted[[20]]$par)
  out_fitted$innov <- innovations(c(1, init_tpdf), max_q = 20)
  out_fitted$innov_short <- ma_q_tpdf(out_fitted$innov[[1]][q_sim, ])
  out_fitted$innov_long  <- ma_q_tpdf(out_fitted$innov[[1]][20, ])
  model_tpd        <- ma_q_tpdf(thetas_sim)
  errors_ma10[i, 1] <- sum((model_tpd - out_fitted$tpdf_long)^2)
  errors_ma10[i, 2] <- sum((model_tpd - out_fitted$innov_long)^2)
  errors_ma10[i, 3] <- sum((model_tpd - out_fitted$tpdf_short)[1:q_sim]^2)
  errors_ma10[i, 4] <- sum((model_tpd - out_fitted$innov_short)[1:q_sim]^2)
  errors_ma10[i, 5] <- sum(abs(model_tpd - out_fitted$tpdf_long))
  errors_ma10[i, 6] <- sum(abs(model_tpd - out_fitted$innov_long))
  errors_ma10[i, 7] <- sum(abs(model_tpd - out_fitted$tpdf_short)[1:q_sim])
  errors_ma10[i, 8] <- sum(abs(model_tpd - out_fitted$innov_short)[1:q_sim])
  if(i%%10 == 0){print(paste0("Iteration ", i, " complete"))}
}

results2 <- get_ab_results(errors_ma10)


set.seed(982)
for(i in 1:num_sims){
  # ma(15)
  q_sim <- 15
  thetas_sim <- runif(q_sim)
  while (min(abs(polyroot(c(1,thetas_sim)))) < 1){ # check if valid model
    thetas_sim <- runif(q_sim)
  }
  temp_data  <- gen_maq(n = length_ts, thetas = thetas_sim) # generate ts
  temp_pairs <- get_all_large_pairs(temp_data, thresh = 0.95) # get pairs for lhood
  init_tpdf  <- TPDF(temp_data, maxlag = 20)[-1]
  out_fitted <- optim_MAq_multiple(theta_init = 0.1,  
                                       points_list = temp_pairs$points)
  out_fitted$tpdf_short <- ma_q_tpdf(out_fitted[[q_sim]]$par)
  out_fitted$tpdf_long  <- ma_q_tpdf(out_fitted[[20]]$par)
  out_fitted$innov <- innovations(c(1, init_tpdf), max_q = 20)
  out_fitted$innov_short <- ma_q_tpdf(out_fitted$innov[[1]][q_sim, ])
  out_fitted$innov_long  <- ma_q_tpdf(out_fitted$innov[[1]][20, ])
  model_tpd        <- ma_q_tpdf(thetas_sim)
  errors_ma15[i, 1] <- sum((model_tpd - out_fitted$tpdf_long)^2)
  errors_ma15[i, 2] <- sum((model_tpd - out_fitted$innov_long)^2)
  errors_ma15[i, 3] <- sum((model_tpd - out_fitted$tpdf_short)[1:q_sim]^2)
  errors_ma15[i, 4] <- sum((model_tpd - out_fitted$innov_short)[1:q_sim]^2)
  errors_ma15[i, 5] <- sum(abs(model_tpd - out_fitted$tpdf_long))
  errors_ma15[i, 6] <- sum(abs(model_tpd - out_fitted$innov_long))
  errors_ma15[i, 7] <- sum(abs(model_tpd - out_fitted$tpdf_short)[1:q_sim])
  errors_ma15[i, 8] <- sum(abs(model_tpd - out_fitted$innov_short)[1:q_sim])
  if(i%%10 == 0){print(paste0("Iteration ", i, " complete"))}
}

results3 <- get_ab_results(errors_ma15)


### 
# arma(1,1) fit
###

arma_least_squares <- function(theta_phi, tpd){
  model_tpd <- arma_11_tpdf(params = theta_phi)
  sum((model_tpd - tpd)^2)
}
arma_least_absolute <- function(theta_phi, tpd){
  model_tpd <- arma_11_tpdf(params = theta_phi)
  sum(abs(model_tpd - tpd))
}

set.seed(982)
for(i in 1:num_sims){
  # arma(1,1)
  par_sim <- runif(2) # generate ts
  temp_data  <- gen_arma11(n = length_ts, theta = par_sim[1], phi = par_sim[2]) 
  temp_pairs <- get_all_large_pairs(temp_data, thresh = 0.95) # get pairs for lhood
  init_tpdf  <- TPDF(temp_data, maxlag = 20)[-1]
  out_fitted <- optim_ARMA11_from_list(theta_phi = c(0.1, 0.1),
                                       points_list = temp_pairs$points)
  out_fitted$tpdf   <- arma_11_tpdf(out_fitted$par)
  model_tpd         <- arma_11_tpdf(par_sim)
  out_fitted$ls_par <- optim(out_fitted$par, fn = arma_least_squares, 
                             tpd = init_tpdf)$par
  out_fitted$la_par <- optim(out_fitted$par, fn = arma_least_absolute, 
                             tpd = init_tpdf)$par
  q_sim <- max(which(model_tpd >= 0.01))
  errors_arma11[i, 1] <- sum((model_tpd - out_fitted$tpdf)^2)
  errors_arma11[i, 2] <- arma_least_squares(out_fitted$ls_par, model_tpd)
  errors_arma11[i, 3] <- arma_least_squares(out_fitted$la_par, model_tpd)
  errors_arma11[i, 4] <- sum(abs(model_tpd - out_fitted$tpdf))
  errors_arma11[i, 5] <- arma_least_absolute(out_fitted$ls_par, model_tpd)
  errors_arma11[i, 6] <- arma_least_absolute(out_fitted$la_par, model_tpd)
  if(i%%10 == 0){print(paste0("Iteration ", i, " complete"))}
}


get_ab_results <- function(x){
  list(mean_hr_full = mean(x$hr_full_ab), 
       sd_HR_full = sd(x$hr_full_ab),
       mean_ls_full = mean(x$ls_full_ab), 
       sd_ls_full = sd(x$ls_full_ab), 
       mean_la_full = mean(x$la_full_ab), 
       sd_la_full = sd(x$la_full_ab), 
       num_HR_better_ls = sum(x$hr_full_ab - x$ls_full_ab < 0),
       num_HR_better_la = sum(x$hr_full_ab - x$la_full_ab < 0))
}

results4 <- get_ab_results(errors_arma11)

saveRDS(errors_arma11, "./data/proxy_sim_model_fitting_arma_results.rds")


full_results <- list(errors_ma5 = errors_ma5, 
                     errors_ma10 = errors_ma10, 
                     errors_ma15 = errors_ma15, 
                     errors_arma11 = errors_arma11, 
                     results_ma5 = results1, 
                     results_ma10 = results2,
                     results_ma15 = results3, 
                     results_arma11 = results4)
saveRDS(full_results, "./data/proxy_sim_model_fitting_results_all_081926.rds")
