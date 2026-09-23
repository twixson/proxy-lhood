############ 
# This script performs the simulations in the TPD estimation section
#### 

set.seed(873)
length_ts    <- 10000
source("./code/proxy_likelihood_fitting.R")
library(cowplot)
library(tidyverse)

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
# check if matrix of TPD values is positive definite, if not find a nearby PD
init_tpdm    <- toeplitz(c(1, init_tpdf))
if(min(eigen(init_tpdm)$val) < 0){
  init_tpdm  <- nearPD(init_tpdm)
  init_tpdf  <- init_tpdm[1, 2:21] 
  print("~~~ FOUND A NON-PD EMPIRICAL TPDF ~~~")
}
init_lambdas <- get_lambda(init_tpdf)

# 2: optimize the negative log-likelihood
out_llhood <- optim(par = init_lambdas, # initial values
                    ll_from_list,       # function (negative llhood)
                    points_list = temp_pairs$points, # pairwise data
                    control = list(maxit = 10000),   # need lots of iterations
                    hessian = TRUE)     # return the hessian to check PD
out_llhood$par
min(eigen(out_llhood$hessian)$value)

# 3: compute tpdf for fitted lambda values
out_llhood$tpdf <- sapply(out_llhood$par, hr_tpdm)
  # check if matrix of TPD values is positive definite, if not find a nearby PD
out_llhood$tpdm <- toeplitz(c(1, out_llhood$tpdf))
if(min(eigen(init_tpdm)$val) < 0){
  out_llhood$tpdm  <- nearPD(out_llhood$tpdm)
  out_llhood$tpdf  <- out_llhood$tpdf[1, 2:21] 
  print("~~~ FOUND A NON-PD PROXY TPDF ~~~")
}

plot(0:20, c(1, ma_q_tpdf(thetas_sim)), type = "h", xlab = "lag", ylab = "TPD")
lines(0:20 + 0.1, c(1, out_llhood$tpdf), col = 2, type = "h")
lines(0:20 + 0.2, c(1, init_tpdf), col = 4, type = "h")
legend("topright", bty = "n", lty = 1, col = c(1, 2, 4),
       legend = c("model", "likelihood", "standard"), 
       text.col = c(1, 2, 4))

# 4: Try innovations on likelihood fitted tpd values
out_llhood$innov <- innovations(c(1, out_llhood$tpdf), max_q = 20)
out_llhood$innov_tpdf <- ma_q_tpdf(out_llhood$innov[[1]][20, 1:20])
# check if matrix of TPD values is positive definite, if not find a nearby PD
out_llhood$innov_tpdm <- toeplitz(c(1, out_llhood$innov_tpdf))
if(min(eigen(init_tpdm)$val) < 0){
  out_llhood$innov_tpdm  <- nearPD(out_llhood$innov_tpdm)
  out_llhood$innov_tpdf  <- out_llhood$innov_tpdf[1, 2:21] 
  print("~~~ FOUND A NON-PD LHOOD INNOVATIONS TPDF ~~~")
}
round(out_llhood$innov[[1]][20, ], 3)
lines(0:20 + 0.3, c(1, out_llhood$innov_tpdf), col = 6, type = "h")

# 5: Compare to basic innovations fit
innov_emp  <- innovations(c(1, init_tpdf), max_q = 20)
innov_tpdf <- ma_q_tpdf(innov_emp[[1]][20, ])
# check if matrix of TPD values is positive definite, if not find a nearby PD
innov_tdpm <- toeplitz(c(1, innov_tpdf))
if(min(eigen(init_tpdm)$val) < 0){
  innov_tdpm  <- nearPD(innov_tpdm)
  innov_tpdf  <- innov_tpdf[1, 2:21] 
  print("~~~ FOUND A NON-PD INNOVATIONS TPDF ~~~")
}

plot(0:20, c(1, ma_q_tpdf(thetas_sim)), type = "h", xlab = "lag", ylab = "TPD")
lines(0:20 + 0.1, c(1, out_llhood$innov_tpdf), type = "h", col = 2)
lines(0:20 + 0.2, c(1, innov_tpdf), type = "h", col = 4)
legend("topright", bty = "n", lty = 1, col = c(1, 2, 4),
       legend = c("model", "likelihood", "standard"), 
       text.col = c(1, 2, 4))
      # it looks the same as the above plot, unsurprising I suppose... 

# 6: Create a nicer looking plot for the paper
plot_max_lag <- length(out_llhood$tpdf)
tpdf_tibble1  <- tibble(lag = 0:plot_max_lag,
                       tpd = c(1, ma_q_tpdf(thetas_sim)),
                       tpd_est_hr = c(1, out_llhood$tpdf),
                       tpd_est_natural = c(1, init_tpdf))

p1 <- ggplot(data = tpdf_tibble1, aes(x = lag, y = tpd)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0), lwd = 0.75, col = "black") +
  geom_segment(mapping = aes(x = lag  + 0.15, y = tpd_est_hr,
                             xend = lag + 0.15, yend = 0),
               col = "orange", lwd= 1) +
  geom_segment(mapping = aes(x = lag  + 0.3, y = tpd_est_natural,
                             xend = lag + 0.3, yend = 0),
               col = "#009E73", lwd= 1) +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 2,
           label= "Model") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 4,
           label = "HR est.", col = "orange") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 6,
           label = "Empirical est.", color = "#009E73") +
  annotate(geom = "text", x=0,y=Inf, hjust = -0.15, vjust = 1,
           label = paste0("TL-MA(5)")) + 
  labs(x = "Lag", y = "TPD") +
  scale_x_continuous(breaks = 0:20) +
  theme_minimal()
p1

# 7: Compute the square error between the estimates and the model TPD 
errors1 <- 
  data.frame(
    hr_full = sum((tpdf_tibble1$tpd - tpdf_tibble1$tpd_est_hr)^2), 
    natural_full = sum((tpdf_tibble1$tpd - tpdf_tibble1$tpd_est_natural)^2),
    hr = sum((tpdf_tibble1$tpd[1:q_sim] - 
                      tpdf_tibble1$tpd_est_hr[1:q_sim])^2), 
    natural_short = sum((tpdf_tibble1$tpd[1:q_sim] - 
                           tpdf_tibble1$tpd_est_natural[1:q_sim])^2))





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
# check if matrix of TPD values is positive definite, if not find a nearby PD
init_tpdm    <- toeplitz(c(1, init_tpdf))
if(min(eigen(init_tpdm)$val) < 0){
  init_tpdm  <- nearPD(init_tpdm)
  init_tpdf  <- init_tpdm[1, 2:21] 
  print("~~~ FOUND A NON-PD EMPIRICAL TPDF ~~~")
}
init_lambdas <- get_lambda(init_tpdf)

# 2: optimize the negative log-likelihood
out_llhood <- optim(par = init_lambdas, # initial values
                    ll_from_list,       # function (negative llhood)
                    points_list = temp_pairs$points, # pairwise data
                    control = list(maxit = 10000),   # need lots of iterations
                    hessian = TRUE)     # return the hessian to check PD
out_llhood$par
min(eigen(out_llhood$hessian)$value)

# 3: compute tpdf for fitted lambda values
out_llhood$tpdf <- sapply(out_llhood$par, hr_tpdm)
# check if matrix of TPD values is positive definite, if not find a nearby PD
out_llhood$tpdm <- toeplitz(c(1, out_llhood$tpdf))
if(min(eigen(init_tpdm)$val) < 0){
  out_llhood$tpdm  <- nearPD(out_llhood$tpdm)
  out_llhood$tpdf  <- out_llhood$tpdf[1, 2:21] 
  print("~~~ FOUND A NON-PD PROXY TPDF ~~~")
}

plot(0:20, c(1, ma_q_tpdf(thetas_sim)), type = "h", xlab = "lag", ylab = "TPD")
lines(0:20 + 0.1, c(1, out_llhood$tpdf), col = 2, type = "h")
lines(0:20 + 0.2, c(1, init_tpdf), col = 4, type = "h")
legend("topright", bty = "n", lty = 1, col = c(1, 2, 4),
       legend = c("model", "likelihood", "standard"), 
       text.col = c(1, 2, 4))

# 4: Try innovations on likelihood fitted tpd values
out_llhood$innov <- innovations(c(1, out_llhood$tpdf), max_q = 20)
out_llhood$innov_tpdf <- ma_q_tpdf(out_llhood$innov[[1]][20, 1:20])
# check if matrix of TPD values is positive definite, if not find a nearby PD
out_llhood$innov_tpdm <- toeplitz(c(1, out_llhood$innov_tpdf))
if(min(eigen(init_tpdm)$val) < 0){
  out_llhood$innov_tpdm  <- nearPD(out_llhood$innov_tpdm)
  out_llhood$innov_tpdf  <- out_llhood$innov_tpdf[1, 2:21] 
  print("~~~ FOUND A NON-PD LHOOD INNOVATIONS TPDF ~~~")
}
round(out_llhood$innov[[1]][20, ], 3)
lines(0:20 + 0.3, c(1, out_llhood$innov_tpdf), col = 6, type = "h")

# 5: Compare to basic innovations fit
innov_emp  <- innovations(c(1, init_tpdf), max_q = 20)
innov_tpdf <- ma_q_tpdf(innov_emp[[1]][20, ])
# check if matrix of TPD values is positive definite, if not find a nearby PD
innov_tdpm <- toeplitz(c(1, innov_tpdf))
if(min(eigen(init_tpdm)$val) < 0){
  innov_tdpm  <- nearPD(innov_tpdm)
  innov_tpdf  <- innov_tpdf[1, 2:21] 
  print("~~~ FOUND A NON-PD INNOVATIONS TPDF ~~~")
}

plot(0:20, c(1, ma_q_tpdf(thetas_sim)), type = "h", xlab = "lag", ylab = "TPD")
lines(0:20 + 0.1, c(1, out_llhood$innov_tpdf), type = "h", col = 2)
lines(0:20 + 0.2, c(1, innov_tpdf), type = "h", col = 4)
legend("topright", bty = "n", lty = 1, col = c(1, 2, 4),
       legend = c("model", "likelihood", "standard"), 
       text.col = c(1, 2, 4))
# it looks the same as the above plot, unsurprising I suppose... 

# 6: Create a nicer looking plot for the paper
plot_max_lag <- length(out_llhood$tpdf)
tpdf_tibble2  <- tibble(lag = 0:plot_max_lag,
                       tpd = c(1, ma_q_tpdf(thetas_sim)),
                       tpd_est_hr = c(1, out_llhood$tpdf),
                       tpd_est_natural = c(1, init_tpdf))

p2 <- ggplot(data = tpdf_tibble2, aes(x = lag, y = tpd)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0), lwd = 0.75, col = "black") +
  geom_segment(mapping = aes(x = lag  + 0.15, y = tpd_est_hr,
                             xend = lag + 0.15, yend = 0),
               col = "orange", lwd= 1) +
  geom_segment(mapping = aes(x = lag  + 0.3, y = tpd_est_natural,
                             xend = lag + 0.3, yend = 0),
               col = "#009E73", lwd= 1) +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 2,
           label= "Model") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 4,
           label = "HR est.", col = "orange") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 6,
           label = "Empirical est.", color = "#009E73") +
  annotate(geom = "text", x=0,y=Inf, hjust = -0.15, vjust = 1,
           label = paste0("TL-MA(10)")) + 
  labs(x = "Lag", y = "TPD") +
  scale_x_continuous(breaks = 0:20) +
  theme_minimal()
p2

# 7: Compute the square error between the estimates and the model TPD 
errors2 <- 
  data.frame(
    hr_full = sum((tpdf_tibble2$tpd - tpdf_tibble2$tpd_est_hr)^2), 
    natural_full = sum((tpdf_tibble2$tpd - tpdf_tibble2$tpd_est_natural)^2),
    hr_short = sum((tpdf_tibble2$tpd[1:q_sim] - 
                      tpdf_tibble2$tpd_est_hr[1:q_sim])^2), 
    natural_short = sum((tpdf_tibble2$tpd[1:q_sim] - 
                           tpdf_tibble2$tpd_est_natural[1:q_sim])^2))

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
# check if matrix of TPD values is positive definite, if not find a nearby PD
init_tpdm    <- toeplitz(c(1, init_tpdf))
if(min(eigen(init_tpdm)$val) < 0){
  init_tpdm  <- nearPD(init_tpdm)
  init_tpdf  <- init_tpdm[1, 2:21] 
  print("~~~ FOUND A NON-PD EMPIRICAL TPDF ~~~")
}
init_lambdas <- get_lambda(init_tpdf)

# 2: optimize the negative log-likelihood
out_llhood <- optim(par = init_lambdas, # initial values
                    ll_from_list,       # function (negative llhood)
                    points_list = temp_pairs$points, # pairwise data
                    control = list(maxit = 10000),   # need lots of iterations
                    hessian = TRUE)     # return the hessian to check PD
out_llhood$par
min(eigen(out_llhood$hessian)$value)

# 3: compute tpdf for fitted lambda values
out_llhood$tpdf <- sapply(out_llhood$par, hr_tpdm)
# check if matrix of TPD values is positive definite, if not find a nearby PD
out_llhood$tpdm <- toeplitz(c(1, out_llhood$tpdf))
if(min(eigen(init_tpdm)$val) < 0){
  out_llhood$tpdm  <- nearPD(out_llhood$tpdm)
  out_llhood$tpdf  <- out_llhood$tpdf[1, 2:21] 
  print("~~~ FOUND A NON-PD PROXY TPDF ~~~")
}

plot(0:20, c(1, ma_q_tpdf(thetas_sim)), type = "h", xlab = "lag", ylab = "TPD")
lines(0:20 + 0.1, c(1, out_llhood$tpdf), col = 2, type = "h")
lines(0:20 + 0.2, c(1, init_tpdf), col = 4, type = "h")
legend("topright", bty = "n", lty = 1, col = c(1, 2, 4),
       legend = c("model", "likelihood", "standard"), 
       text.col = c(1, 2, 4))

# 4: Try innovations on likelihood fitted tpd values
out_llhood$innov <- innovations(c(1, out_llhood$tpdf), max_q = 20)
out_llhood$innov_tpdf <- ma_q_tpdf(out_llhood$innov[[1]][20, 1:20])
# check if matrix of TPD values is positive definite, if not find a nearby PD
out_llhood$innov_tpdm <- toeplitz(c(1, out_llhood$innov_tpdf))
if(min(eigen(init_tpdm)$val) < 0){
  out_llhood$innov_tpdm  <- nearPD(out_llhood$innov_tpdm)
  out_llhood$innov_tpdf  <- out_llhood$innov_tpdf[1, 2:21] 
  print("~~~ FOUND A NON-PD LHOOD INNOVATIONS TPDF ~~~")
}
round(out_llhood$innov[[1]][20, ], 3)
lines(0:20 + 0.3, c(1, out_llhood$innov_tpdf), col = 6, type = "h")

# 5: Compare to basic innovations fit
innov_emp  <- innovations(c(1, init_tpdf), max_q = 20)
innov_tpdf <- ma_q_tpdf(innov_emp[[1]][20, ])
# check if matrix of TPD values is positive definite, if not find a nearby PD
innov_tdpm <- toeplitz(c(1, innov_tpdf))
if(min(eigen(init_tpdm)$val) < 0){
  innov_tdpm  <- nearPD(innov_tpdm)
  innov_tpdf  <- innov_tpdf[1, 2:21] 
  print("~~~ FOUND A NON-PD INNOVATIONS TPDF ~~~")
}

plot(0:20, c(1, ma_q_tpdf(thetas_sim)), type = "h", xlab = "lag", ylab = "TPD")
lines(0:20 + 0.1, c(1, out_llhood$innov_tpdf), type = "h", col = 2)
lines(0:20 + 0.2, c(1, innov_tpdf), type = "h", col = 4)
legend("topright", bty = "n", lty = 1, col = c(1, 2, 4),
       legend = c("model", "likelihood", "standard"), 
       text.col = c(1, 2, 4))
# it looks the same as the above plot, unsurprising I suppose... 

# 6: Create a nicer looking plot for the paper
plot_max_lag <- length(out_llhood$tpdf)
tpdf_tibble3  <- tibble(lag = 0:plot_max_lag,
                       tpd = c(1, ma_q_tpdf(thetas_sim)),
                       tpd_est_hr = c(1, out_llhood$tpdf),
                       tpd_est_natural = c(1, init_tpdf))

p3 <- ggplot(data = tpdf_tibble3, aes(x = lag, y = tpd)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0), lwd = 0.75, col = "black") +
  geom_segment(mapping = aes(x = lag  + 0.15, y = tpd_est_hr,
                             xend = lag + 0.15, yend = 0),
               col = "orange", lwd= 1) +
  geom_segment(mapping = aes(x = lag  + 0.3, y = tpd_est_natural,
                             xend = lag + 0.3, yend = 0),
               col = "#009E73", lwd= 1) +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 2,
           label= "Model ") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 4,
           label = "HR est.", col = "orange") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 6,
           label = "Empirical est.", color = "#009E73") +
  annotate(geom = "text", x=0,y=Inf, hjust = -0.15, vjust = 1,
           label = paste0("TL-MA(15)")) + 
  labs(x = "Lag", y = "TPD") +
  scale_x_continuous(breaks = 0:20) +
  theme_minimal()
p3

# 7: Compute the square error between the estimates and the model TPD 
errors3 <- 
  data.frame(
    hr_full = sum((tpdf_tibble3$tpd - tpdf_tibble3$tpd_est_hr)^2), 
    natural_full = sum((tpdf_tibble3$tpd - tpdf_tibble3$tpd_est_natural)^2),
    hr_short = sum((tpdf_tibble3$tpd[1:q_sim] - 
                      tpdf_tibble3$tpd_est_hr[1:q_sim])^2), 
    natural_short = sum((tpdf_tibble3$tpd[1:q_sim] - 
                           tpdf_tibble3$tpd_est_natural[1:q_sim])^2))





########
# Try again with ARMA model
########


# 0: get some data to estimate the TPD of: 
par_sim <- runif(2) # generate ts
temp_data  <- gen_arma11(n = length_ts, theta = par_sim[1], phi = par_sim[2]) 
temp_pairs <- get_all_large_pairs(temp_data, thresh = 0.95) # get pairs for lhood

# 1: get initial values for optimization:
init_tpdf    <- TPDF(temp_data, maxlag = 20)[-1]
# check if matrix of TPD values is positive definite, if not find a nearby PD
init_tpdm    <- toeplitz(c(1, init_tpdf))
if(min(eigen(init_tpdm)$val) < 0){
  init_tpdm  <- nearPD(init_tpdm)
  init_tpdf  <- init_tpdm[1, 2:21] 
  print("~~~ FOUND A NON-PD EMPIRICAL TPDF ~~~")
}
init_lambdas <- get_lambda(init_tpdf)

# 2: optimize the negative log-likelihood
out_llhood <- optim(par = init_lambdas, # initial values
                    ll_from_list,       # function (negative llhood)
                    points_list = temp_pairs$points, # pairwise data
                    control = list(maxit = 10000),   # need lots of iterations
                    hessian = TRUE)     # return the hessian to check PD
out_llhood$par
min(eigen(out_llhood$hessian)$value)

# 3: compute tpdf for fitted lambda values
out_llhood$tpdf <- sapply(out_llhood$par, hr_tpdm)
# check if matrix of TPD values is positive definite, if not find a nearby PD
out_llhood$tpdm <- toeplitz(c(1, out_llhood$tpdf))
if(min(eigen(init_tpdm)$val) < 0){
  out_llhood$tpdm  <- nearPD(out_llhood$tpdm)
  out_llhood$tpdf  <- out_llhood$tpdf[1, 2:21] 
  print("~~~ FOUND A NON-PD PROXY TPDF ~~~")
}

plot(0:20, c(1, arma_11_tpdf(par_sim)), type = "h", xlab = "lag", ylab = "TPD")
lines(0:20 + 0.1, c(1, out_llhood$tpdf), col = 2, type = "h")
lines(0:20 + 0.2, c(1, init_tpdf), col = 4, type = "h")
legend("topright", bty = "n", lty = 1, col = c(1, 2, 4),
       legend = c("model", "likelihood", "standard"), 
       text.col = c(1, 2, 4))

# 4: Try innovations on likelihood fitted tpd values
out_llhood$innov <- innovations(c(1, out_llhood$tpdf), max_q = 20)
out_llhood$innov_tpdf <- ma_q_tpdf(out_llhood$innov[[1]][20, 1:20])
# check if matrix of TPD values is positive definite, if not find a nearby PD
out_llhood$innov_tpdm <- toeplitz(c(1, out_llhood$innov_tpdf))
if(min(eigen(init_tpdm)$val) < 0){
  out_llhood$innov_tpdm  <- nearPD(out_llhood$innov_tpdm)
  out_llhood$innov_tpdf  <- out_llhood$innov_tpdf[1, 2:21] 
  print("~~~ FOUND A NON-PD LHOOD INNOVATIONS TPDF ~~~")
}
round(out_llhood$innov[[1]][20, ], 3)
lines(0:20 + 0.3, c(1, out_llhood$innov_tpdf), col = 6, type = "h")

# 6: Create a nicer looking plot for the paper
plot_max_lag <- length(out_llhood$tpdf)
tpdf_tibble4  <- tibble(lag = 0:plot_max_lag,
                       tpd = c(1, arma_11_tpdf(par_sim)),
                       tpd_est_hr = c(1, out_llhood$tpdf),
                       tpd_est_natural = c(1, init_tpdf))

p4 <- ggplot(data = tpdf_tibble4, aes(x = lag, y = tpd)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0), lwd = 0.75, col = "black") +
  geom_segment(mapping = aes(x = lag  + 0.15, y = tpd_est_hr,
                             xend = lag + 0.15, yend = 0),
               col = "orange", lwd= 1) +
  geom_segment(mapping = aes(x = lag  + 0.3, y = tpd_est_natural,
                             xend = lag + 0.3, yend = 0),
               col = "#009E73", lwd= 1) +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 2,
           label= "Model") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 4,
           label = "HR est.", col = "orange") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 6,
           label = "Empirical est.", color = "#009E73") +
  annotate(geom = "text", x=0,y=Inf, hjust = -0.15, vjust = 1,
           label = paste0("TL-ARMA(1,1)")) + 
  labs(x = "Lag", y = "TPD") + 
       # title = paste0("TL-ARMA(", round(par_sim, 2), ", ", 
       #              round(par_sim[2], 2), ")")) +
  scale_x_continuous(breaks = 0:20) +
  theme_minimal()
p4

# 7: Compute the square error between the estimates and the model TPD 
errors4 <- 
  data.frame(
    hr_full = sum((tpdf_tibble4$tpd - tpdf_tibble4$tpd_est_hr)^2), 
    natural_full = sum((tpdf_tibble4$tpd - tpdf_tibble4$tpd_est_natural)^2),
    hr_short = sum((tpdf_tibble4$tpd[1:10] - 
                      tpdf_tibble4$tpd_est_hr[1:10])^2), 
    natural_short = sum((tpdf_tibble4$tpd[1:10] - 
                           tpdf_tibble4$tpd_est_natural[1:10])^2))



#######
# Combine plots for paper
#######
a <- cowplot::plot_grid(p1, p2, p3, p4)
save_plot("./plots/tpd_est_4model_plot.png", a, base_height = 3.4, base_asp = 3)

b <- cowplot::plot_grid(p1, p2)
save_plot("./plots/tpd_est_2model_plot.png", b, base_height = 3.4, base_asp = 3.5)




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
                         natural_full_sq = rep(0, num_sims),
                         hr_short_sq = rep(0, num_sims), 
                         natural_short_sq = rep(0, num_sims), 
                         hr_full_ab = rep(0, num_sims),
                         natural_full_ab = rep(0, num_sims),
                         hr_short_ab = rep(0, num_sims), 
                         natural_short_ab = rep(0, num_sims))

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
  init_tpdf    <- TPDF(temp_data, maxlag = 20)[-1]
  # check if matrix of TPD values is positive definite, if not find a nearby PD
  init_tpdm    <- toeplitz(c(1, init_tpdf))
  if(min(eigen(init_tpdm)$val) < 0){
    init_tpdm  <- nearPD(init_tpdm)
    init_tpdf  <- init_tpdm[1, 2:21] 
    print("~~~ FOUND A NON-PD EMPIRICAL TPDF ~~~")
  }
  init_lambdas <- get_lambda(init_tpdf)
  out_llhood <- optim(par = init_lambdas, # initial values
                      ll_from_list,       # function (negative llhood)
                      points_list = temp_pairs$points, # pairwise data
                      control = list(maxit = 10000))   # need lots of iterations
  out_llhood$tpdf <- sapply(out_llhood$par, hr_tpdm)
  # check if matrix of TPD values is positive definite, if not find a nearby PD
  out_llhood$tpdm <- toeplitz(c(1, out_llhood$tpdf))
  if(min(eigen(init_tpdm)$val) < 0){
    out_llhood$tpdm  <- nearPD(out_llhood$tpdm)
    out_llhood$tpdf  <- out_llhood$tpdf[1, 2:21] 
    print("~~~ FOUND A NON-PD PROXY TPDF ~~~")
  }
  model_tpd <- ma_q_tpdf(thetas_sim)
  errors_ma5[i, 1] <- sum((model_tpd - out_llhood$tpdf)^2)
  errors_ma5[i, 2] <- sum((model_tpd - init_tpdf)^2)
  errors_ma5[i, 3] <- sum((model_tpd - out_llhood$tpdf)[1:q_sim]^2)
  errors_ma5[i, 4] <- sum((model_tpd - init_tpdf)[1:q_sim]^2)
  errors_ma5[i, 5] <- sum(abs(model_tpd - out_llhood$tpdf))
  errors_ma5[i, 6] <- sum(abs(model_tpd - init_tpdf))
  errors_ma5[i, 7] <- sum(abs(model_tpd - out_llhood$tpdf)[1:q_sim])
  errors_ma5[i, 8] <- sum(abs(model_tpd - init_tpdf)[1:q_sim])
  
  if(i%%10 == 0){print(paste0("Iteration ", i, " complete"))}
}


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
  init_tpdf    <- TPDF(temp_data, maxlag = 20)[-1]
  # check if matrix of TPD values is positive definite, if not find a nearby PD
  init_tpdm    <- toeplitz(c(1, init_tpdf))
  if(min(eigen(init_tpdm)$val) < 0){
    init_tpdm  <- nearPD(init_tpdm)
    init_tpdf  <- init_tpdm[1, 2:21] 
    print("~~~ FOUND A NON-PD EMPIRICAL TPDF ~~~")
  }
  init_lambdas <- get_lambda(init_tpdf)
  out_llhood <- optim(par = init_lambdas, # initial values
                      ll_from_list,       # function (negative llhood)
                      points_list = temp_pairs$points, # pairwise data
                      control = list(maxit = 10000))   # need lots of iterations
  out_llhood$tpdf <- sapply(out_llhood$par, hr_tpdm)
  # check if matrix of TPD values is positive definite, if not find a nearby PD
  out_llhood$tpdm <- toeplitz(c(1, out_llhood$tpdf))
  if(min(eigen(init_tpdm)$val) < 0){
    out_llhood$tpdm  <- nearPD(out_llhood$tpdm)
    out_llhood$tpdf  <- out_llhood$tpdf[1, 2:21] 
    print("~~~ FOUND A NON-PD PROXY TPDF ~~~")
  }
  model_tpd <- ma_q_tpdf(thetas_sim)
  errors_ma10[i, 1] <- sum((model_tpd - out_llhood$tpdf)^2)
  errors_ma10[i, 2] <- sum((model_tpd - init_tpdf)^2)
  errors_ma10[i, 3] <- sum((model_tpd - out_llhood$tpdf)[1:q_sim]^2)
  errors_ma10[i, 4] <- sum((model_tpd - init_tpdf)[1:q_sim]^2)
  errors_ma10[i, 5] <- sum(abs(model_tpd - out_llhood$tpdf))
  errors_ma10[i, 6] <- sum(abs(model_tpd - init_tpdf))
  errors_ma10[i, 7] <- sum(abs(model_tpd - out_llhood$tpdf)[1:q_sim])
  errors_ma10[i, 8] <- sum(abs(model_tpd - init_tpdf)[1:q_sim])
  if(i%%10 == 0){print(paste0("Iteration ", i, " complete"))}
}



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
  init_tpdf    <- TPDF(temp_data, maxlag = 20)[-1]
  # check if matrix of TPD values is positive definite, if not find a nearby PD
  init_tpdm    <- toeplitz(c(1, init_tpdf))
  if(min(eigen(init_tpdm)$val) < 0){
    init_tpdm  <- nearPD(init_tpdm)
    init_tpdf  <- init_tpdm[1, 2:21] 
    print("~~~ FOUND A NON-PD EMPIRICAL TPDF ~~~")
  }
  init_lambdas <- get_lambda(init_tpdf)
  out_llhood <- optim(par = init_lambdas, # initial values
                      ll_from_list,       # function (negative llhood)
                      points_list = temp_pairs$points, # pairwise data
                      control = list(maxit = 10000))   # need lots of iterations
  out_llhood$tpdf <- sapply(out_llhood$par, hr_tpdm)
  # check if matrix of TPD values is positive definite, if not find a nearby PD
  out_llhood$tpdm <- toeplitz(c(1, out_llhood$tpdf))
  if(min(eigen(init_tpdm)$val) < 0){
    out_llhood$tpdm  <- nearPD(out_llhood$tpdm)
    out_llhood$tpdf  <- out_llhood$tpdf[1, 2:21] 
    print("~~~ FOUND A NON-PD PROXY TPDF ~~~")
  }
  model_tpd <- ma_q_tpdf(thetas_sim)
  errors_ma15[i, 1] <- sum((model_tpd - out_llhood$tpdf)^2)
  errors_ma15[i, 2] <- sum((model_tpd - init_tpdf)^2)
  errors_ma15[i, 3] <- sum((model_tpd - out_llhood$tpdf)[1:q_sim]^2)
  errors_ma15[i, 4] <- sum((model_tpd - init_tpdf)[1:q_sim]^2)
  errors_ma15[i, 5] <- sum(abs(model_tpd - out_llhood$tpdf))
  errors_ma15[i, 6] <- sum(abs(model_tpd - init_tpdf))
  errors_ma15[i, 7] <- sum(abs(model_tpd - out_llhood$tpdf)[1:q_sim])
  errors_ma15[i, 8] <- sum(abs(model_tpd - init_tpdf)[1:q_sim])
  if(i%%10 == 0){print(paste0("Iteration ", i, " complete"))}
}



set.seed(982)
for(i in 1:num_sims){
  # arma(1,1)
  par_sim <- runif(2) # generate ts
  temp_data  <- gen_arma11(n = length_ts, theta = par_sim[1], phi = par_sim[2]) 
  temp_pairs <- get_all_large_pairs(temp_data, thresh = 0.95) # get pairs for lhood
  init_tpdf    <- TPDF(temp_data, maxlag = 20)[-1]
  # check if matrix of TPD values is positive definite, if not find a nearby PD
  init_tpdm    <- toeplitz(c(1, init_tpdf))
  if(min(eigen(init_tpdm)$val) < 0){
    init_tpdm  <- nearPD(init_tpdm)
    init_tpdf  <- init_tpdm[1, 2:21] 
    print("~~~ FOUND A NON-PD EMPIRICAL TPDF ~~~")
  }
  init_lambdas <- get_lambda(init_tpdf)
  out_llhood <- optim(par = init_lambdas, # initial values
                      ll_from_list,       # function (negative llhood)
                      points_list = temp_pairs$points, # pairwise data
                      control = list(maxit = 10000))   # need lots of iterations
  out_llhood$tpdf <- sapply(out_llhood$par, hr_tpdm)
  # check if matrix of TPD values is positive definite, if not find a nearby PD
  out_llhood$tpdm <- toeplitz(c(1, out_llhood$tpdf))
  if(min(eigen(init_tpdm)$val) < 0){
    out_llhood$tpdm  <- nearPD(out_llhood$tpdm)
    out_llhood$tpdf  <- out_llhood$tpdf[1, 2:21] 
    print("~~~ FOUND A NON-PD PROXY TPDF ~~~")
  }
  model_tpd <- arma_11_tpdf(params = par_sim)
  q_sim <- max(which(model_tpd >= 0.01))
  errors_arma11[i, 1] <- sum((model_tpd - out_llhood$tpdf)^2)
  errors_arma11[i, 2] <- sum((model_tpd - init_tpdf)^2)
  errors_arma11[i, 3] <- sum((model_tpd - out_llhood$tpdf)[1:q_sim]^2)
  errors_arma11[i, 4] <- sum((model_tpd - init_tpdf)[1:q_sim]^2)
  errors_arma11[i, 5] <- sum(abs(model_tpd - out_llhood$tpdf))
  errors_arma11[i, 6] <- sum(abs(model_tpd - init_tpdf))
  errors_arma11[i, 7] <- sum(abs(model_tpd - out_llhood$tpdf)[1:q_sim])
  errors_arma11[i, 8] <- sum(abs(model_tpd - init_tpdf)[1:q_sim])
  if(i%%10 == 0){print(paste0("Iteration ", i, " complete"))}
}


get_ab_results <- function(x){
  list(mean_hr_short = mean(x$hr_short_ab), 
       sd_HR_short = sd(x$hr_short_ab),
       mean_nat_short = mean(x$natural_short_ab), 
       sd_nat_short = sd(x$natural_short_ab), 
       num_HR_short_better = sum(x$hr_short_ab - x$natural_short_ab < 0),
       mean_diff_short = mean(x$hr_short_ab - x$natural_short_ab), 
       sd_diff_short = sd(x$hr_short_ab - x$natural_short_ab),
       mean_hr_full = mean(x$hr_full_ab), 
       sd_HR_full = sd(x$hr_full_ab),
       mean_nat_full = mean(x$natural_full_ab), 
       sd_nat_full = sd(x$natural_full_ab), 
       num_HR_full_better = sum(x$hr_full_ab - x$natural_full_ab < 0),
       mean_diff_full = mean(x$hr_full_ab - x$natural_full_ab), 
       sd_diff_full = sd(x$hr_full_ab - x$natural_full_ab))
}

error_results <- list(errors_ma5 = errors_ma5, 
                      errors_ma10 = errors_ma10, 
                      errors_ma15 = errors_ma15, 
                      errors_arma11 = errors_arma11)
saveRDS(error_results, "./data/proxy_tpd_estimator_simulation_results.rds")
