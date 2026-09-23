############ 
# This script compiles the necessary functions to do a TLETS proxy-likelihood 
#       simulation study. 
#### 

source("./generate_TLETS.R")
lam_vec <- readRDS("./HR_lam_from_sigma_vector_by_10e-4.rds")
# grad_vec <- readRDS("./HR_grad_from_sigma_vector_5e-4_15e-4_25e-4.rds")

#####
# Functions to get list of large lag-h points
get_large_lag_h_pairs <- function(ts, h, thresh){
  n         <- length(ts)
  temp_mat  <- matrix(c(ts[1:(n-h)], ts[(1+h):n]), ncol = 2, byrow = F)
  temp_rad  <- apply(temp_mat, 1, function(x){sqrt(sum(x^2))})
  large_ind <- which(temp_rad > unname(quantile(temp_rad, prob = thresh)))
  return(list(indices = large_ind, points = temp_mat[large_ind, ]))
}

get_all_large_pairs <- function(ts, max_h=20, thresh = 0.975, bias_corr = F){
  if(bias_corr){
    temp_mean <- mean(ts)
    n_temp <- length(ts)
    ts <- pmax(rep(0, n_temp), x - mean(x))
  }
  temp_out <- lapply(1:max_h, get_large_lag_h_pairs, 
                     ts = ts, thresh = thresh)
  indices  <- lapply(temp_out, function(x){x$indices})
  points   <- lapply(temp_out, function(x){x$points})
  return(list(indices = indices, points = points))
}

#####
# ma(q)-to-sigma functions
ma_q_tpdf <- function(thetas, max_lag. = 20){
  q <- length(thetas)
  sigmas <- rep(0, max_lag.)
  thetas <- c(1, thetas, rep(0, max_lag.))
  for(i in 1:q){
    sigmas[i] <- sum(pmax(thetas[1:(q+1)], rep(0, q+1)) * 
                       pmax(thetas[(i+1):(i+q+1)], rep(0, q+1))) / 
      sum(thetas^2) # this is sigma_0 and thus forces the scale to be 1
  }
  return(sigmas)
}

#####
# AR(1)-to-sigma functions
ar_1_tpdf <- function(phi, max_lag. = 20){
  pmax(phi^(1:max_lag.), rep(0, max_lag.)) # scale 1 implies just phi^h
}

#####
# ARMA(1,1)-to-sigma functions 
get_arma_11_h <- function(params, h){
  theta <- params[1]
  phi <- params[2]
  out <- 0
  if(phi > 0){
    if(phi + theta > 0){
      num <- (phi + theta) * phi^(h-1) * (1 + phi*theta)
      den <- 1 + 2 * phi * theta + theta^2
      out <- num / den
    }
  } else if(phi + theta > 0){
    if(h %% 2 == 0){
      num <- (phi + theta)^2 * phi ^ h
      den <- 1 - phi^4 + (phi + theta)^2
      out <- num / den
    } else {
      num <- (phi + theta) * phi^(h - 1) * (1 - phi^4)
      den <- 1 - phi^4 + (phi + theta)^2
      out <- num / den
    }
  } else if(phi + theta < 0){
    if(h %% 2 == 0){
      num <- (phi + theta) * phi^(h - 1) * (1 + theta * phi^3)
      den <- 1 + phi^2 * theta^2 + 2 * phi^3 * theta
      out <- num / den
    }
  }
  return(out)
}

arma_11_tpdf <- function(params, max_lag. = 20){
  temp_out <- rep(0, max_lag.)
  for(h in 1:max_lag.){
    temp_out[h] <- get_arma_11_h(params, h)
  }
  return(temp_out)
}

#####
# functions to convert tpdf to HR lambdas and get gradient

get_lambda <- function(sigma, lam_vec. = lam_vec){
  p_max <- length(sigma)
  out <- numeric(p_max)
  if(p_max == 1){
    temp_ind <- floor(sigma * 1000)
    if(temp_ind < 1){
      out <- 100 - ((1000 * sigma) %% 1) * (100 - lam_vec.[1])
    } else {
      out <- lam_vec.[temp_ind] - 
        ((1000 * sigma) %% 1) * (lam_vec.[temp_ind] - lam_vec.[temp_ind + 1])
    }
  } else{
    for(p in 1:p_max){
      temp_ind <- floor(sigma[p] * 1000)
      if(temp_ind < 1){
        out[p] <- 100 - ((1000 * sigma[p]) %% 1) * (100 - lam_vec.[1])
      }else { 
        out[p] <- lam_vec.[temp_ind] - 
          ((1000 * sigma[p]) %% 1) * 
          (lam_vec.[temp_ind] - lam_vec.[temp_ind + 1])
      }
    }
  }
  return(out)
}

#####
# This function is emperical gradient 

# get_grad <- function(sigma, grad_vec. = grad_vec){
#   p_max <- length(sigma)
#   out <- numeric(p_max)
#   if(p_max == 1){
#     temp_ind <- floor((sigma + 0.0005) * 1000)
#     step_size <- round((1000 * sigma + 0.5), 10) %% 1
#     if(temp_ind <= 1){
#       out <- step_size * (-5000 - grad_vec.[1])
#     } else {
#       out <- grad_vec.[temp_ind] - 
#         step_size * (grad_vec.[temp_ind] - grad_vec.[temp_ind + 1])
#     }
#   } else{
#     for(p in 1:p_max){
#       temp_ind <- floor((sigma[p] + 0.0005) * 1000)
#       step_size <- round((1000 * sigma[p] + 0.5), 10) %% 1
#       if(temp_ind <= 1){
#         out[p] <- step_size * (-10000 - grad_vec.[1])
#       } else { 
#         out[p] <- grad_vec.[temp_ind] - 
#           step_size * (grad_vec.[temp_ind] - grad_vec.[temp_ind + 1])
#       }
#     }
#   }
#   return(out)
# }

# smoothing splines version of lambda and gradient of lambda functions
sigma <- seq(0.001, 1, length.out = 1000)
temp_curve <- smooth.spline(sigma, lam_vec)
get_lambda <- function(sigma, splines = temp_curve){
  predict(splines, sigma)$y
}
get_grad <- function(sigma, grad_splines = temp_curve){
  predict(grad_splines, sigma)$y
}

#####
# composite HR-function
HR_biv_density <- function(point, l){
  if(is.null(dim(point))){
    x1 <- max(point[1], 0.0001)
    x2 <- max(point[2], 0.0001)
  } else { 
    x1 <- pmax(point[,1], rep(0.0001, length(point[,1])))
    x2 <- pmax(point[,2], rep(0.0001, length(point[,1])))
  }
  
  a1  <- l - (1/l)*log(x1/x2)
  a2  <- l - (1/l)*log(x2/x1)
  cf1 <- pnorm(a1)
  cf2 <- pnorm(a2)
  df1 <- dnorm(a1)
  df2 <- dnorm(a2)
  G      <- exp(- cf1 / x1^2 - cf2 / x2^2)
  dvx1   <- - 2 * cf1/(x1^3) - df1 / (l * x1^3) + df2 /  (l * x1 * x2^2)
  dvx2   <- - 2 * cf2/(x2^3) - df2 / (l * x2^3) + df1 /  (l * x2 * x1^2)
  dvx1x2 <- df1 * (-l^2 - log(x1/x2)) / (l^3 * x1^3 * x2) + 
    df2 * (-l^2 - log(x2/x1)) / (l^3 * x1 * x2^3)
  
  return(G * (dvx1 * dvx2 - dvx1x2))
}

#####
# log-likelihood contributions
ll_one_h <- function(X, l){
  temp_contribs <<- HR_biv_density(X, l = l)
  if(is.null(is.na(temp_contribs))){
    if(min(temp_contribs) < 1e-10){
      print("At least one point has zero HR density, replaced with 1e-10")
    }
  }
  ll_contrib <- log(pmax(temp_contribs, rep(1e-10, length(temp_contribs))))
  sum(-ll_contrib)
}

ll_from_list <- function(points_list, lambdas, weights. = 0, print_contribs = F){
  temp <- numeric(length(lambdas))
  for(i in 1:(length(lambdas))){
    #print(paste0("lag ", i))
    temp[i] <- ll_one_h(X = points_list[[i]], l = lambdas[i])
  }
  if(print_contribs == TRUE){
    print("Negative log likelihood contributions by lag (before weights):")
    print(temp)
  }
  weighted_out <- (1/(1:length(lambdas))^(weights.)) * temp
  if(print_contribs == TRUE){
    print("Negative log likelihood contributions by lag (after weights):")
    print(weighted_out)
  }
  sum(weighted_out)
}

##### 
# log-likelihood function wrapper for implementation in optim
get_AR1_ll_from_list <- function(phi, points_list, max_lag=20, weights = 0){
  temp_tpdf <- ar_1_tpdf(phi, max_lag. = max_lag) 
  lambdas   <- get_lambda(temp_tpdf)
  ll_from_list(points_list = points_list, lambdas = lambdas, weights. = weights)
}

optim_AR1_from_list <- function(phi, points_list, max_lag=20, weights = 0){
  temp_out <- optim(get_AR1_ll_from_list, 
                    par = phi, 
                    points_list = points_list, 
                    max_lag = max_lag, 
                    weights = weights, 
                    method = "Brent", 
                    lower = 0, upper = 1, 
                    hessian = T)
  
  if(temp_out$hessian <= 0){
    k = 0
    while(temp_out$hessian < 0 && k < 5){
      print(paste0("ERROR: negative hessian in AR1 model. Try again, k = ", k))
      k = k+1
      temp_out <- optim(get_AR1_ll_from_list, 
                        par = runif(1), 
                        points_list = points_list, 
                        max_lag = max_lag, 
                        weights = weights, 
                        method = "Brent", 
                        lower = 0, upper = 1, 
                        hessian = T)
    }
    if(temp_out$hessian <= 0){
      temp_out$hessian <- 0.00001
      print(paste0("ERROR: negative hessian in AR1 model. Tried k = ", k, 
                   " times, setting hessian to be flat"))
    }
  }
  return(temp_out)
  
}

get_ARMA11_ll_from_list <- function(theta_phi, points_list, 
                                    max_lag=20, weights = 0){
  temp_tpdf <- arma_11_tpdf(params = theta_phi,  max_lag. = max_lag) 
  lambdas   <- get_lambda(temp_tpdf)
  ll_from_list(points_list = points_list, lambdas = lambdas, weights. = weights)
}

optim_ARMA11_from_list <- function(theta_phi, points_list, 
                                   max_lag = 20, weights = 0, 
                                   omethod = "CG"){
  temp_out <- optim(par = theta_phi, 
                    get_ARMA11_ll_from_list,
                    points_list = points_list, 
                    max_lag = max_lag, 
                    weights = weights, 
                    hessian = T, 
                    method = omethod)
  
  if(det(temp_out$hessian) <= 0){
    k = 0
    while(det(temp_out$hessian) <= 0 && k < 5){
      print(paste0("ERROR: negative hessian in ARMA11 model. Try again, k = ", k))
      k = k+1
      arma11_fit <- optim(par = runif(2), 
                          get_ARMA11_ll_from_list,
                          points_list = points_list, 
                          max_lag = max_lag, 
                          weights = weights, 
                          hessian = T, 
                          method = "L-BFGS-B", 
                          lower = c(-1, 0), upper = c(1,1))
    }
    if(det(temp_out$hessian)<= 0){
      temp_out$hessian <- matrix(c(0.001, 0, 0, 0.001), nrow = 2)
      print(paste0("ERROR: negative hessian in ARMA11 model. Tried k = ", k, 
                   " times, setting hessian to be flat"))
    }
  }
  return(temp_out)
}

get_MAq_from_list <- function(thetas, points_list, max_lag=20, weights = 0, 
                              print_optim_steps = F){
  temp_tpdf <- ma_q_tpdf(thetas, max_lag. = max_lag)
  lambdas   <- get_lambda(temp_tpdf)
  if(print_optim_steps == TRUE){
    print("thetas: ")
    print(round(thetas[1:length(thetas)], 3))
  }
  ll_from_list(points_list = points_list, lambdas = lambdas, weights. = weights)
}

optim_MAq_single <- function(thetas, points_list, max_lag=20, weights = 0){
  if(length(thetas)==1){
    optim(get_MAq_from_list, 
          par = thetas, 
          points_list = points_list, 
          max_lag = max_lag, 
          weights = weights, 
          method = "Brent", 
          lower = 0, upper = 1, 
          hessian = T)
  } else {
    optim(get_MAq_from_list, 
          par = thetas, 
          points_list = points_list, 
          max_lag = max_lag, 
          weights = weights, 
          hessian = T, 
          method = "L-BFGS-B", 
          lower = 0, upper = 1)
  }
}


optim_MAq_sequentially <- function(max_q, points_list, max_lag = 20, weights = 0, 
                                   omethod = "CG"){
  ma_seq_fit <- list()
  ma_seq_fit$MA1 <- optim(get_MAq_from_list, 
                          par = 0.1, 
                          points_list = points_list, 
                          max_lag = max_lag, 
                          weights = weights, 
                          method = "Brent", 
                          lower = 0, upper = 1, 
                          hessian = T)
  for(q in 2:max_q){
    temp_out <- optim(get_MAq_from_list,
                      par = c(ma_seq_fit[[q-1]]$par, 0.5),
                      points_list = points_list,
                      max_lag = max_lag,
                      weights = weights, 
                      hessian = T, 
                      method = "L-BFGS-B", 
                      lower = 0, upper = 1)
    if(temp_out$convergence != 0){
      print("ERROR: No convergence, trying again")
      temp_out <- optim(get_MAq_from_list,
                        par = c(ma_seq_fit[[q-1]]$par, runif(1)),
                        points_list = points_list,
                        max_lag = max_lag,
                        weights = weights, 
                        hessian = T, 
                        method = "L-BFGS-B", 
                        lower = 0, upper = 1)
    }
    if(det(temp_out$hessian) <= 0){
      k = 0
      while(det(temp_out$hessian) <= 0 && k < 5){
        print(paste0("ERROR: negative hessian in MA(",q, 
                     ") model. Try again, k = ", k))
        k = k+1
        temp_out <- optim(get_MAq_from_list,
                          par = c(ma_seq_fit[[q-1]]$par, runif(1)),
                          points_list = points_list,
                          max_lag = max_lag,
                          weights = weights, 
                          hessian = T, 
                          method = "L-BFGS-B", 
                          lower = -1, upper = 1)
      }
      if(det(temp_out$hessian) <= 0){
        temp_out$hessian <- diag(0.001, nrow = q)
        print(paste0("ERROR: negative hessian in MA model. Tried k = ", k, 
                     " times, setting hessian to be flat"))
      }
    }
    ma_seq_fit[[paste0("MA", q)]] <- temp_out
    if(q %% 5 == 0){
      print(paste0("~~~ Order ", q, " fitting complete ~~~"))
    }
  }
  return(ma_seq_fit)
}

get_innov_MAq_thetas <- function(ts, max_lag=20){
  tpdf <- TPDF(ts, maxlag = max_lag)
  innovations(tpdf, max_q = max_lag)[[1]]
}


#####
# Functions to fit the whole library of models and to extract info from list
extract_info <- function(list1){
  output <- list()
  output$parameters  <- lapply(list1, function(x){x$par})
  output$lhood_score <- sapply(list1, function(x){x$value})
  output$hessian     <- lapply(list1, function(x){x$hessian})
  return(output)
}

fit_library <- function(data_list, max_lag = 20, weights = 0){
  out <- optim_MAq_sequentially(max_q = max_lag, 
                                points_list = data_list, 
                                max_lag = max_lag, 
                                weights = weights)
  out$AR1 <- optim_AR1_from_list(phi = 0.2, 
                                 points_list = data_list, 
                                 max_lag = max_lag, 
                                 weights = weights)
  out$ARMA11 <- optim_ARMA11_from_list(theta_phi = c(0.2, 0.3), 
                                       points_list = data_list, 
                                       max_lag = max_lag, 
                                       weights = weights)
  return(extract_info(out))
}

#####
# score contribution of one point
HR_score_contribution <- function(x, lambda.){ # eq (6) in manuscript
  xi   <- x[1]
  xj   <- x[2]
  lij  <- log(xi/xj)
  lji  <- log(xj/xi)
  aij  <- lambda. - lij/lambda.
  aji  <- lambda. - lji/lambda.
  
  dP1  <- (1 + lij/(lambda.^2)) * dnorm(aij)
  dP2  <- (1 + lji/(lambda.^2)) * dnorm(aji)
  dVi  <- 
    -2/(xi^3) * pnorm(aij) - 
    1/(lambda. * (xi^3)) * dnorm(aij) + 
    1/(lambda. * xi * (xj^2)) * dnorm(aji)
  dVj  <- 
    -2/(xj^3) * pnorm(aji) - 
    1/(lambda. * (xj^3)) * dnorm(aji) + 
    1/(lambda. * (xi^2) * xj) * dnorm(aij)
  dVidVj <- 
    (-lambda.^2 - lij)/(lambda.^3 * (xi^3) * xj) * dnorm(aij) + 
    (-lambda.^2 - lji)/(lambda.^3 * xi * (xj^3)) * dnorm(aji)
  dldVi <-
    -2/(xi^3) * (1 + lij/(lambda.^2)) * dnorm(aij) + 
    (lambda.^2 + 1 - lij^2/(lambda.^4))/(xi^3) * dnorm(aij) - 
    (lambda.^2 + 1 - lji^2/(lambda.^4))/(xi * (xj^2)) * dnorm(aji)
  dldVj <- 
    -2/(xj^3) * (1 + lji/(lambda.^2)) * dnorm(aji) + 
    (lambda.^2 + 1 - lji^2/(lambda.^4))/(xj^3) * dnorm(aji) - 
    (lambda.^2 + 1 - lij^2/(lambda.^4))/((xi^2) * xj) * dnorm(aij)
  dldVidVj <- 
    (1 + 1/(lambda.^2) + (1/(lambda.^2) + 3/(lambda.^4))*lij - 
       lij^2/(lambda.^4) - lij^3/(lambda.^6)) * dnorm(aij)/(xi^3*xj) + 
    (1 + 1/(lambda.^2) + (1/(lambda.^2) + 3/(lambda.^4))*lji - 
       lji^2/(lambda.^4) - lji^3/(lambda.^6)) * dnorm(aji)/(xi*xj^3)
  
  t3_num <- dVj * dldVi + dVi * dldVj - dldVidVj
  t3_den <- dVi * dVj - dVidVj
  return(-dP1/(xi^2) -dP2/(xj^2) + t3_num/t3_den)
}

HR_score_pairwise <- function(data, lambda){
  data <- ifelse(data > 0.000001, data, 0.000001)
  scores <- apply(data, 1, HR_score_contribution, lambda. = lambda)
  sum(scores)
}

##### 
# d/dtheta sigma(theta)
#####
# AR(1)
ddtheta_ar1 <- function(param, max_lag.=20){
  return((1:max_lag.)*param^(0:(max_lag.-1)))
}
# MA(q)
ddtheta_maq <- function(params, dindex, max_lag.=20){
  params_ext <- c(1, params)
  params_0s  <- c(rep(0, max_lag.), params_ext, rep(0, max_lag.))
  params_max <- pmax(rep(0, length(params_0s)), params_0s)
  num1       <- params_max[(dindex + max_lag. + 1) + (1:max_lag.)] + 
    params_max[(dindex + max_lag. + 1) - (1:max_lag.)]
  den1       <- sum(params_ext^2)
  num2       <- apply(matrix(1:max_lag., nrow = 1), 2, function(h){
    2 * params_ext[dindex + 1] * 
      sum(params_max[(max_lag. + 1):(max_lag. + length(params_ext))] * 
            params_max[(max_lag. + 1 + h):(max_lag. + length(params_ext) + h)])
  })
  den2       <- den1^2
  return(num1/den1 + num2/den2)
}
# ARMA(1,1)
ddtheta_arma11_theta <- function(phi, theta, max_lag.=20){
  dval <- num1 <- num2 <- numeric(max_lag.)
  if(phi > 0){
    if(phi + theta > 0){ # case 1
      num1 <- theta * phi^((1:max_lag.)-1) * (1 + phi * theta) * (1 + phi)
      den1 <- 1 + 2*theta*phi + theta^2
      num2 <- 2*(phi + theta)^2 * phi^((1:max_lag.)-1) * (1 + phi * theta)
      den2 <- den1^2
      dval <- num1/den1 - num2/den2
    } else { # case 2
      dval <- rep(0, max_lag)
    }
  } else if(phi + theta > 0){ # case 3/4
    evens <- seq(2, max_lag., by = 2)
    odds  <- seq(1, max_lag., by = 2)
    num1[evens] <- 2 * (phi + theta) * phi^evens
    den1        <- 1 - phi^4 + (phi + theta)^2
    num2[evens] <- 2 * (phi + theta)^3 * phi^evens
    den2        <- den1^2
    num1[odds]  <- phi^(odds - 1) * (1 - phi^4)
    num2[odds]  <- 2 * (phi + theta)^2 * phi^(odds - 1) * (1 - phi^4)
    dval        <- num1/den1 - num2/den2
  } else { # case 4/5
    evens <- seq(2, max_lag., by = 2)
    odds  <- seq(1, max_lag., by = 2)
    num1[evens] <- phi^(evens - 1) * (1 + 2 * theta * phi^3 + phi^4) * 
      (1 + theta^2 * phi^2 + 2 * theta * phi^3)
    den1        <- 1 + theta^2 * phi^2 + 2 * theta * phi^3
    num2[evens] <- 2 * (theta * phi^2 + phi^3) * (phi + theta) * 
      phi^(evens - 1) * (1 + theta * phi^3)
    den2        <- den1^2
    num1[odds]  <- 0
    num2[odds]  <- 0
    dval <- num1/den1 + num2/den2
  }
  return(dval)
}
ddtheta_arma11_phi <- function(phi, theta, max_lag.=20){
  dval <- num1 <- num2 <- numeric(max_lag.)
  if(phi > 0){
    if(phi + theta > 0){ # case 1
      num1 <- phi^((1:max_lag.) - 2) * 
        ((1:max_lag.) * phi + (1:max_lag.) * theta - theta) * (1 + phi * theta) + 
        phi^((1:max_lag.) - 1) * theta * (phi + theta)
      den1 <- 1 + 2*theta*phi + theta^2
      num2 <- 2 * theta * (phi + theta) * phi^((1:max_lag.) - 1) * (1 + phi * theta)
      den2 <- den1^2
      dval <- num1/den1 - num2/den2
    } else { # case 2
      dval <- rep(0, max_lag)
    }
  } else if(phi + theta > 0){ # case 3/4
    evens <- seq(2, max_lag., by = 2)
    odds  <- seq(1, max_lag., by = 2)
    num1[evens] <- 2 * (phi + theta) * phi^evens + 
      (phi + theta)^2 * evens * phi^(evens - 1)
    den1        <- 1 - phi^4 + (phi + theta)^2
    num2[evens] <- (phi + theta)^2 * phi^evens * 
      (4 * phi^3 + 2 * phi + 2 * theta)
    den2        <- den1^2
    num1[odds]  <- phi^(odds - 2) * 
      ((odds * phi + odds * theta - theta) * (1 - phi^4) + 
         4 * phi^5 + 4 * theta * phi ^4)
    num2[odds]  <- phi^(odds - 1) * (phi + theta) * 
      (1 - phi^4) * (2 * phi ^3 + phi + theta)
    dval        <- num1/den1 - num2/den2
  } else { # case 4/5
    evens <- seq(2, max_lag., by = 2)
    odds  <- seq(1, max_lag., by = 2)
    num1[evens] <- phi^(evens - 2) * 
      ((evens * phi + evens * theta - theta) * 
         (1 + theta * phi^3) + 3 * theta * phi ^3 * (phi + theta))
    den1        <- 1 + theta^2 * phi^2 + 2 * theta * phi^3
    num2[evens] <- 2 * theta * phi^evens * (phi + theta) * 
      (1 + theta * phi^3) * (theta + 3 * phi)
    den2        <- den1^2
    num1[odds]  <- 0
    num2[odds]  <- 0
    dval <- num1/den1 + num2/den2
  }
  return(dval)
}

#####
# d/dSigma Lambda(Sigma)
ddsigma <- function(tpdf_vec){
  return(get_grad(tpdf_vec))
}

##### 
# J-matrix computation: start by getting list of M scores (pseudo-replicates)
#           
scores_list <- function(points_list, 
                        indices, 
                        split_length, 
                        ts_length, 
                        model, 
                        Theta_hat, 
                        max_lag=20){
  if(!(model %in% c("ar1", "arma11", "maq"))){
    stop("Invalid model, try 'ar1', 'arma11', or 'maq'")
  } 
  # ddtheta is the derivative of the model parameter to tpdf function
  # ddsigmas is the derivative of the tpd value to lambdas function
  # ddt_dds is product of ddtheta and ddsigmas
  # ddll is the derivative of the log likelihood
  if(model == "ar1"){
    sigmas    <- ar_1_tpdf(Theta_hat, max_lag. = max_lag)
    lambdas   <- get_lambda(sigmas)
    ddtheta   <- ddtheta_ar1(param = Theta_hat, max_lag. = max_lag)
  } else if(model == "arma11"){
    sigmas    <- arma_11_tpdf(Theta_hat, max_lag. = max_lag)
    lambdas   <- get_lambda(sigmas)
    ddtheta_t <- ddtheta_arma11_theta(theta = Theta_hat[1], 
                                      phi = Theta_hat[2], 
                                      max_lag. = max_lag)
    ddtheta_p <- ddtheta_arma11_phi(theta = Theta_hat[1], 
                                    phi = Theta_hat[2], 
                                    max_lag. = max_lag)
    ddtheta   <- matrix(c(ddtheta_t, ddtheta_p), ncol = 2, byrow = F)
  } else if(model == "maq"){
    q_max     <- length(Theta_hat)
    sigmas    <- ma_q_tpdf(Theta_hat, max_lag. = max_lag)
    lambdas   <- get_lambda(sigmas)
    ddtheta   <- matrix(NA, ncol = length(Theta_hat), nrow = max_lag)
    for(i in 1:length(Theta_hat)){
      ddtheta[, i] <- ddtheta_maq(params = Theta_hat, 
                                  dindex = i, 
                                  max_lag. = max_lag)
    }
  }
  ddsigmas <- ddsigma(tpdf_vec = sigmas)
  
  ddts     <- ddtheta * ddsigmas
  
  M        <- ts_length/split_length
  M_splits <- c(seq(1, ts_length, by = split_length), ts_length+1)
  ddll     <- matrix(NA, nrow = M, ncol = max_lag)
  for(h in 1:max_lag){
    for(s in 1:M){
      if(length(which(
        indices[[h]] %in% M_splits[s]:(M_splits[s+1]-1))) == 1){
        ddll[s, h] <- 
          HR_score_contribution(x = 
                                  points_list[[h]][which(
                                    indices[[h]] %in% 
                                      M_splits[s]:(M_splits[s+1]-1)), ], 
                                lambda = lambdas[h])
      }else if(length(which(
        indices[[h]] %in% M_splits[s]:(M_splits[s+1]-1))) > 1){
        ddll[s, h] <- 
          HR_score_pairwise(data = 
                              points_list[[h]][which(
                                indices[[h]] %in% 
                                  M_splits[s]:(M_splits[s+1]-1)), ], 
                            lambda = lambdas[h])
      } else {ddll[s, h] <- 0}
    }
  }
  ddll %*% ddts
}

get_penalty <- function(scores_matrix, hess_mat){
  basic_penalty <- sum(diag(cov(scores_matrix) %*% solve(hess_mat)))
  return(basic_penalty)
}

get_scores <- function(library_fits, ts_list, split_length, ts_length){
  model_names       <- names(library_fits)
  no_penalty <- CLAIC <- CLBIC <- numeric(length(library_fits))
  names(no_penalty) <- model_names
  models            <- c(rep("maq", length(model_names)-2), 
                         "ar1", "arma11")
  penalty <- rep(0, length(models))
  for(i in 1:length(library_fits)){
    no_penalty[i]    <- library_fits[[i]]$value
    temp_scores      <- scores_list(points_list = ts_list$points, 
                                    indices = ts_list$indices, 
                                    split_length = split_length, 
                                    ts_length = ts_length, 
                                    model = models[i], 
                                    Theta_hat = library_fits[[i]]$par)
    if(is.null(library_fits[[i]]$hessian)){
      print(paste0("ERROR: Model ", model_names[i], 
                   " has no hessian, used identity"))
      library_fits[[i]]$hessian <- diag(rep(1, length(library_fits[[i]]$par)))
    }
    if(min(eigen(library_fits[[i]]$hessian)$values) < 0.00001){
      print(paste0("ERROR: Optimization failed in model ", model_names[i], 
                   " (hessian is not PD) unknown uncertainty."))
      # penalty[[i]] <- list("basic" = NA, "new" = NA)
      penalty[i] <- NA
      CLAIC[i] <- NA
      CLBIC[i] <- NA
    } else {
      penalty[i] <- get_penalty(scores_matrix = temp_scores, 
                                hess_mat = library_fits[[i]]$hessian)
      CLAIC[i] <- get_CLAIC(no_penalty[i], penalty[i])
      CLBIC[i] <- get_CLBIC(no_penalty[i], penalty[i], n = ts_length)
    }
  }
  no_pen_min <- model_names[which.min(no_penalty)]
  not_NA     <- !is.na(CLAIC)
  CLAIC_min  <- model_names[not_NA][which.min(na.omit(CLAIC))]
  CLBIC_min  <- model_names[not_NA][which.min(na.omit(CLBIC))]
  return(list("Min_lhood" = no_pen_min, 
              "Min_CLAIC" = CLAIC_min, 
              "Min_CLBIC" = CLBIC_min, 
              "CLAIC"     = CLAIC, 
              "CLBIC"     = CLBIC, 
              "lhood"     = no_penalty, 
              "basic_penalty" = round(penalty, 2)))
}

#####
# Scripts to get CLAIC and CLBIC 
get_CLAIC <- function(neg_llhood, penalty){
  2*neg_llhood + 2*penalty
}

get_CLBIC <- function(neg_llhood, penalty, n){
  2*neg_llhood + log(n)*penalty
}

#####
# Scripts to get AN CI's 
get_cov_matrix <- function(scores_matrix, hess_mat){
  hess_inv <- solve(hess_mat)
  hess_inv %*% cov(scores_matrix) %*% t(hess_inv)
}

get_CIs <- function(params, cov_mat, t_df){
  crit_val <- qt(0.975, df = t_df)
  lower <- params - crit_val*sqrt(diag(cov_mat))
  upper <- params + crit_val*sqrt(diag(cov_mat))
  out   <- matrix(c(lower, upper), ncol = 2, byrow = F)
  colnames(out) <- c("lower", "upper")
  return(out)
}







