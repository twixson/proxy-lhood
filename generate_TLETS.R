########
# Generate TLETS data
########
library(evd)

###
# Transformed linear operations
f <- function(y){
  if(y < -36){
    y = -36
    print("*!*!* mass at machine zero *!*!*")
  } 
  return(ifelse(y > 100, y, log(1+exp(y))))
}

finv <- function(x){
  return(ifelse(x > 100, x, log(exp(x)-1)))
}

tadd <- function(a,b){
  return(f(finv(a)+finv(b)))
}

tmult <- function(a,b){
  return(f(a*finv(b)))
}

innovations <- function(tpdf, max_q =50){
  #initialize variables nu and theta_hat
  nu <- rep(NA, max_q + 1) 
  theta_hat <- matrix(0, nrow = max_q, ncol = max_q)
  #set nu_0 equal to tpdf(0)
  nu[1] <- tpdf[1]
  #compute n theta_hat's and nu_n up to n=max_q iterations
  for(n in 1:max_q){
    for(k in 0:(n-1)){
      if(k==0){
        temp <- tpdf[n+1]
      } else {temp <- 0
        for(j in 0:(k-1)){
          temp <- temp + theta_hat[k, k-j]*theta_hat[n, n-j]*nu[j+1]}
        temp <- tpdf[n-k+1] - temp}
      theta_hat[n, n-k] <- nu[k+1]^(-1)*temp}
    nu[n+1] <- tpdf[1]
    for(l in 0:(n-1)){nu[n+1] <- nu[n+1] - theta_hat[n, n-l]^2*nu[l+1]}}
  results <- list()
  results[[1]] <- theta_hat
  results[[2]] <- nu
  return(results)
}


### 
# tpdf function - from ts vector
TPDF <- function(ts, thresh = 0.975, maxlag = 50){
  tpdf <- rep(0, maxlag+1)
  if(class(ts) == "list"){
    tpdf[1] <- 1
    for(i in 2:(length(ts)+1)){
      temp_mat <- ts[[i-1]]
      temp_r   <- rowSums(temp_mat^2)
      tpdf[i]  <- 2*mean(temp_mat[, 1]*temp_mat[, 2]/temp_r)
    }
  } else {
    rt0 <- sqrt(ts^2 + ts^2)
    r00 <- sort(rt0)[floor(thresh*(length(rt0)+1))]
    tpdf <- rep(0, maxlag+1)
    tpdf[1] <- 2/length(which(rt0>r00)) * sum((ts^2)[which(rt0>r00)] /
                                                (rt0[which(rt0>r00)]^2))
    n <- length(ts)
    rt <- matrix(NA, nrow = n-1, ncol = maxlag)
    r0 <- rep(NA, maxlag)
    for(i in 1:maxlag){
      temp_means <- c(mean(ts[1:(n-i)]), mean(ts[(1+i):n]))
      temp_mat <- matrix(
        c(pmax(rep(0, (n-i)), ts[1:(n-i)] - temp_means[1]), 
          pmax(rep(0, (n-i)), ts[(1+i):n] - temp_means[2])), 
        ncol = 2, byrow = F)
      rt[,i] <- c(sqrt(temp_mat[,1]^2 + temp_mat[,2]^2), rep(NA, i-1))
      temp <- na.omit(rt[,i])
      if(length(temp)){ # ignore zero-length columns
        k <- floor(thresh*(length(temp)+1))
        r0[i] <- sort(temp)[k]
      } else {r0[i] <- 0}
      indices <- which(rt[,i]>r0[i])
      tpdf[i+1] <- 2/length(indices) * sum((temp_mat[,1]*temp_mat[,2])[indices] /
                                             (rt[indices, i]^2))
    }
  }
  return(tpdf)
}

transform_marginal <- function(x, q = 0.975){
  quant_q <- unname(quantile(prob = q, x))
  large_x <- which(x >= quant_q)
  params  <- fpot(x, threshold = quant_q, shape = 1/2, std.err = F)$estimate
  x_new   <- ecdf(x)(x)
  for(i in 1:length(x)){
    if(i %in% large_x){
      x_new[i] <- q + (1-q)*pgpd(x[i], shape = 1/2, loc = quant_q, scale = params)
    }
  }
  
  return(qfrechet(x_new, shape = 2))
}

gen_ar1 <- function(n, phi){
  RVnoise   <- rfrechet(1000+n, shape = 2)
  ar1_ts    <- numeric(1000 + n)
  ar1_ts[1] <- RVnoise[1]
  for(i in 2:(1000+n)){
    ar1_ts[i] <- tadd(tmult(phi, ar1_ts[i-1]), RVnoise[i])
  }
  ar1_ts <- ar1_ts[1001:(1000+n)]
  transform_marginal(ar1_ts)
}

# temp_phi <- 0.6
# temp <- gen_ar1(n = 10000, phi = temp_phi)
# hist(temp, 1000, xlim = c(0, 20), freq = F)
# curve(dfrechet(x, shape = 2), from = 0, to = 20, add = T, col = 2, 1000)
# plot(temp[1:300], type = "l")
# temp_tpdf <- TPDF(temp, maxlag = 20)
# plot(temp_tpdf, type = "h", lwd = 2, ylim = c(0, 1))
# lines((1:21) + 0.1, temp_phi^(0:20), type = "h", col = 2, lwd = 2)

gen_arma11 <- function(n, phi, theta){
  RVnoise   <- rfrechet(1000+n, shape = 2)
  arma11_ts    <- numeric(1000 + n)
  arma11_ts[1] <- RVnoise[1]
  for(i in 2:(1000+n)){
    arma11_ts[i] <- f(phi * finv(arma11_ts[i-1]) + 
                        finv(RVnoise[i]) + 
                        theta * finv(RVnoise[i-1]))
  }
  arma11_ts <- arma11_ts[1001:(1000+n)]
  a <- transform_marginal(arma11_ts)
}

# temp_theta <- 0.9
# temp_phi   <- 0.5
# temp <- gen_arma11(10000, phi = temp_phi, theta = temp_theta)
# hist(temp, 1000, xlim = c(0, 20), freq = F)
# curve(dfrechet(x, shape = 2), from = 0, to = 20, add = T, col = 2, 1000)
# plot(temp[1:300], type = "l")
# temp_tpdf <- TPDF(temp, maxlag = 20)
# plot(temp_tpdf, type = "h", lwd = 2, ylim = c(0, 1))
# # plot AR(1) tpdf for comparison for now
# lines((1:21) + 0.1, temp_phi^(0:20), type = "h", col = 2, lwd = 2)

gen_maq <- function(n, thetas){
  q         <- length(thetas)
  RVnoise   <- rfrechet(q+n+1, shape = 2)
  maq_ts    <- numeric(n)
  temp_vals <- numeric(q+1)
  for(i in (q+1):(n+q+1)){
    temp_vals[1] <- RVnoise[i]
    for(j in 1:q){
      temp_vals[j+1] <- tmult(thetas[j], RVnoise[i-j])
    }
    temp_val <- temp_vals[1]
    for(j in 1:q){
      temp_val <- tadd(temp_val, temp_vals[j+1])
    }
    maq_ts[i-q] <- temp_val
  }
  transform_marginal(maq_ts)
}

# temp_thetas <- c(0.7, 0.5, 0.44, 0.23, 0.1)
# temp <- gen_maq(10000, thetas = temp_thetas)
# hist(temp, 1000, xlim = c(0, 20), freq = F)
# curve(dfrechet(x, shape = 2), from = 0, to = 20, add = T, col = 2, 1000)
# plot(temp[1:300], type = "l")
# temp_tpdf <- TPDF(temp, maxlag = 20)
# plot(temp_tpdf, type = "h", lwd = 2, ylim = c(0, 1))
# thetas_tpdf <- c(1, temp_thetas, rep(0, 10))
# theor_tpdf  <- rep(0, 21)
# theor_tpdf[1] <- 1
# for(i in 2:(length(temp_thetas)+1)){
#   theor_tpdf[i] <- sum(thetas_tpdf[1:length(temp_thetas)] * 
#                          thetas_tpdf[i:(length(temp_thetas) + i - 1)]) / 
#     sum(thetas_tpdf^2)
# }
# lines((1:21) + 0.1, theor_tpdf, type = "h", col = 2, lwd = 2)

