library(tidyverse)
source("./proxy_likelihood_fitting_scripts_dg021725.R")
fire_data <- readRDS("./ERA5_CO_timeseries_data_dg021825.rds")

#### EARLY (past climate)
early_ts  <- fire_data$early
early_tpdf <- TPDF(early_ts, maxlag = 20, thresh = 0.95)

out_early <- get_all_large_pairs(early_ts, thresh = 0.95)
early_ARMA11_fit <- optim_ARMA11_from_list(theta_phi = c(0.1, 0.1),
                                           points_list = out_early$points)
early_ARMA11_tpdf <- arma_11_tpdf(params = early_ARMA11_fit$par)
early_AR1_fit <- optim_AR1_from_list(phi = 0.1,
                                     points_list = out_early$points)
early_AR1_tpdf <- ar_1_tpdf(phi = early_AR1_fit$par)
early_MA_fits <- optim_MAq_sequentially(max_q = 20, 
                                        points_list = out_early$points)

fitted_models_e <- early_MA_fits
fitted_models_e$AR1 <- early_AR1_fit
fitted_models_e$ARMA11 <- early_ARMA11_fit
scores_early <- get_scores(library_fits = fitted_models_e, ts_list = out_early, 
                           split_length = 153, ts_length = 3060)

tpdf_tibble_early <- tibble(lag = 0:20,
                            emp_tpd = early_tpdf,
                            tpd_est_ar1 = c(1, ar_1_tpdf(early_AR1_fit$par, max_lag. = 20)),
                            tpd_est_arma11 = c(1, arma_11_tpdf(params = early_ARMA11_fit$par,
                                                               max_lag. = 20)),
                            tpd_est_ma5 = c(1, ma_q_tpdf(early_MA_fits$MA5$par,
                                                         max_lag. = 20)), 
                            tpd_est_ma10 = c(1, ma_q_tpdf(early_MA_fits$MA10$par,
                                                          max_lag. = 20)),
                            tpd_est_ma15 = c(1, ma_q_tpdf(early_MA_fits$MA15$par,
                                                          max_lag. = 20)), 
                            tpd_est_ma20 = c(1, ma_q_tpdf(early_MA_fits$MA20$par,
                                                          max_lag. = 20)))

p_early <- ggplot(data = tpdf_tibble_early, aes(x = lag, y = emp_tpd)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0), lwd = 0.75, col = "black") +
  geom_segment(mapping = aes(x = lag  + 0.1, y = tpd_est_arma11,
                             xend = lag + 0.1, yend = 0),
               col = "orange", lwd = 0.75) +
  geom_segment(mapping = aes(x = lag  + 0.2, y = tpd_est_ar1,
                             xend = lag + 0.2, yend = 0),
               col = "#009E73", lwd = 0.75) +
  geom_segment(mapping = aes(x = lag  + 0.3, y = tpd_est_ma10,
                             xend = lag + 0.3, yend = 0),
               col = "#CC79A7", lwd = 0.75) +
  geom_segment(mapping = aes(x = lag  + 0.4, y = tpd_est_ma15,
                             xend = lag + 0.4, yend = 0),
               col = "#0072B2", lwd = 0.75) +
  geom_segment(mapping = aes(x = lag  + 0.5, y = tpd_est_ma20,
                             xend = lag + 0.5, yend = 0),
               col = "grey62", lwd = 0.75) +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 2,
           label= "Empirical") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 4,
           label = "fitted ARMA(1,1)", col = "orange") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 6,
           label = "fitted AR(1)", color = "#009E73") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 8,
           label = "fitted MA(10)", color = "#CC79A7") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 10,
           label = "fitted MA(15)", color = "#0072B2") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 12,
           label = "fitted MA(20)", color = "grey62") + # #D55E00
  labs(x = "Lag", y = "TPD") +
  scale_x_continuous(breaks = 0:20) +
  theme_minimal()
p_early
ggsave(p_early, filename = "./early_ERA5CO_tpd_plot.pdf", 
       device = cairo_pdf,
       width = 7, height = 3, units = "in")


######## LATE
late_ts   <- fire_data$late
late_tpdf  <- TPDF(late_ts, maxlag = 20, thresh = 0.95)

out_late <- get_all_large_pairs(late_ts, thresh = 0.95)
late_ARMA11_fit <- optim_ARMA11_from_list(theta_phi = c(0.1, 0.1),
                                          points_list = out_late$points)
late_ARMA11_tpdf <- arma_11_tpdf(params = late_ARMA11_fit$par)
late_AR1_fit <- optim_AR1_from_list(phi = 0.1, weights = -2, 
                                    points_list = out_late$points)
late_AR1_tpdf <- ar_1_tpdf(phi = late_AR1_fit$par)
late_MA_fits <- optim_MAq_sequentially(max_q = 20, 
                                       points_list = out_late$points)

fitted_models <- late_MA_fits
fitted_models$AR1 <- late_AR1_fit
fitted_models$ARMA11 <- late_ARMA11_fit

scores_late <- get_scores(library_fits = fitted_models, ts_list = out_late, 
                          split_length = 153, ts_length = 3060)

tpdf_tibble_late <- tibble(lag = 0:20,
                           emp_tpd = late_tpdf,
                           tpd_est_ar1 = c(1, ar_1_tpdf(late_AR1_fit$par, max_lag. = 20)),
                           tpd_est_arma11 = c(1, arma_11_tpdf(params = late_ARMA11_fit$par,
                                                              max_lag. = 20)),
                           tpd_est_ma5 = c(1, ma_q_tpdf(late_MA_fits$MA5$par,
                                                        max_lag. = 20)), 
                           tpd_est_ma10 = c(1, ma_q_tpdf(late_MA_fits$MA10$par,
                                                         max_lag. = 20)),
                           tpd_est_ma15 = c(1, ma_q_tpdf(late_MA_fits$MA15$par,
                                                         max_lag. = 20)), 
                           tpd_est_ma20 = c(1, ma_q_tpdf(late_MA_fits$MA20$par,
                                                         max_lag. = 20)))

p_late <- ggplot(data = tpdf_tibble_late, aes(x = lag, y = emp_tpd)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0), lwd = 0.75, col = "black") +
  geom_segment(mapping = aes(x = lag  + 0.1, y = tpd_est_arma11,
                             xend = lag + 0.1, yend = 0),
               col = "orange", lwd = 0.75) +
  geom_segment(mapping = aes(x = lag  + 0.2, y = tpd_est_ar1,
                             xend = lag + 0.2, yend = 0),
               col = "#009E73", lwd = 0.75) +
  geom_segment(mapping = aes(x = lag  + 0.3, y = tpd_est_ma10,
                             xend = lag + 0.3, yend = 0),
               col = "#CC79A7", lwd = 0.75) +
  geom_segment(mapping = aes(x = lag  + 0.4, y = tpd_est_ma15,
                             xend = lag + 0.4, yend = 0),
               col = "#0072B2", lwd = 0.75) +
  geom_segment(mapping = aes(x = lag  + 0.5, y = tpd_est_ma20,
                             xend = lag + 0.5, yend = 0),
               col = "grey62", lwd = 0.75) +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 2,
           label= "Empirical") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 4,
           label = "fitted ARMA(1,1)", col = "orange") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 6,
           label = "fitted AR(1)", color = "#009E73") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 8,
           label = "fitted MA(10)", color = "#CC79A7") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 10,
           label = "fitted MA(15)", color = "#0072B2") +
  annotate(geom = "text", x=Inf,y=Inf, hjust = 1, vjust = 12,
           label = "fitted MA(20)", color = "grey62") +
  labs(x = "Lag", y = "TPD") +
  scale_x_continuous(breaks = 0:20) +
  theme_minimal()
p_late
ggsave(p_late, filename = "./late_ERA5CO_tpd_plot.pdf", 
       device = cairo_pdf,
       width = 7, height = 3, units = "in")

