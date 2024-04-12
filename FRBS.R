library(changepoints)


# Functions
## Bases functions (eigenfunctions of covariance function of Xi)
phi_fct = function(k, t){
  if(k == 1){
    phi = 1
  }else{
    phi = sqrt(2)*cos((k-1)*pi*t)
  }
  return(phi)
}

## generate Xi
x_fct = function(i,t){
  return(sum(sapply(1:K, function(x){v_mat[x,i]*u_vec[x]*phi_fct(x,t)})))
}

## generate slope functions (before and after change)
beta1_fct = function(t){
  return(sum(sapply(1:K, function(x){3*(-1)^(x+1)*x^(-2)*phi_fct(x,t)})))
}
beta2_fct = function(t){
  return(sum(sapply(1:K, function(x){4*(-1)^(x+1)*exp(-x)*phi_fct(x,t)})))
}


## Functional linear regression in RKHS
rkhs_reg = function(yc, Xc_mat, lambda){
  if(is.vector(Xc_mat) && is.atomic(Xc_mat)){
    Xc_mat = matrix(Xc_mat, length(Xc_mat), 1)
  }
  n = length(yc)
  p = dim(Xc_mat)[2]
  grid_for_x = seq(0,1,length=p) # assume evenly spacing
  K_grid = matrix(0, p, p)
  for(j in 1:p){
    for(i in 1:j){
      K_grid[i,j] = cosh(grid_for_x[i])*cosh(1-grid_for_x[j])/sinh(1)
    }
  }
  K_grid = K_grid + t(K_grid) - diag(diag(K_grid))
  B = (1/p^2) * Xc_mat%*%K_grid%*%t(Xc_mat) # for fitting nxn
  C = (1/p) * K_grid%*%t(Xc_mat) # for evaluating estimate pxn
  design_mat = B + lambda*diag(rep(1,n))
  ahat <- lm(yc~design_mat-1)$coefficients# nx1
  beta_hat_grid = C %*% ahat # px1
  yhat = B %*% ahat
  return(list(RSS = sum((yc - yhat)^2), yhat = yhat, beta_hat_grid = beta_hat_grid))
}

## CUSUM statistics
CUSUM_rkhs = function(y, X_mat, lambda, s, e, t){
  if(is.vector(X_mat) && is.atomic(X_mat)){
    X_mat = matrix(X_mat, length(X_mat), 1)
  }
  yc12 = scale(y[s:e], center = TRUE, scale = FALSE)
  Xc12_mat = scale(X_mat[s:e,], center = TRUE, scale = FALSE)
  yc1 = scale(y[s:t], center = TRUE, scale = FALSE)
  Xc1_mat = scale(X_mat[s:t,], center = TRUE, scale = FALSE)
  yc2 = scale(y[(t+1):e], center = TRUE, scale = FALSE)
  Xc2_mat = scale(X_mat[(t+1):e,], center = TRUE, scale = FALSE)
  result = rkhs_reg(yc12, Xc12_mat, lambda)$RSS - rkhs_reg(yc1, Xc1_mat, lambda)$RSS - rkhs_reg(yc2, Xc2_mat, lambda)$RSS
  return(result)
}

## seeded intervals
seeded.intervals = function(n, delta){
  M = ceiling(log2(n/delta))+1
  n_vec = c(2^c(1:M) - 1)
  s_vec = n/2^c(1:M)
  l_vec = 2*s_vec
  output = NULL
  for(k in 1:M){
    for(i in 1:n_vec[k]){
      output = rbind(output, c(ceiling((i-1)*s_vec[k])+1, floor((i-1)*s_vec[k]+l_vec[k])))
    }
  }
  return(output)
}

## FRBS
SBS.rkhs = function(y, X_mat, lambda, s, e, Alpha, Beta, delta, level = 0){
  print(paste0("SBS at level: ", level))
  Alpha_new = pmax(Alpha, s)
  Beta_new = pmin(Beta, e)
  idx = which(Beta_new - Alpha_new > 2*delta)
  Alpha_new = Alpha_new[idx]
  Beta_new = Beta_new[idx]
  M = length(Alpha_new)
  S = NULL
  Dval = NULL
  Level = NULL
  Parent = NULL
  if(M == 0){
    return(list(S = S, Dval = Dval, Level = Level, Parent = Parent))
  }else{
    level = level + 1
    parent = matrix(c(s, e), nrow = 2)
    a = rep(0, M)
    b = rep(0, M)
    for(m in 1:M){
      s_star = Alpha_new[m] + delta
      e_star = Beta_new[m] - delta
      temp = rep(0, e_star - s_star + 1)
      for(t in s_star:e_star){
        temp[t-s_star+1] = CUSUM_rkhs(y, X_mat, lambda, Alpha_new[m], Beta_new[m], t)
      }
      best_value = max(temp)
      best_t = which.max(temp) + s_star - 1
      a[m] = best_value
      b[m] = best_t
    }
    m_star = which.max(a)
  }
  temp1 = SBS.rkhs(y, X_mat, lambda, s, b[m_star]-1, Alpha, Beta, delta, level)
  temp2 = SBS.rkhs(y, X_mat, lambda, b[m_star], e, Alpha, Beta, delta, level)
  S = c(temp1$S, b[m_star], temp2$S)
  Dval = c(temp1$Dval, a[m_star], temp2$Dval)
  Level = c(temp1$Level, level, temp2$Level)
  Parent = cbind(temp1$Parent, parent, temp2$Parent)
  result = list(S = S, Dval = Dval, Level = Level, Parent = Parent)
  class(result) = "BS"
  return(result)
}

## Compute training and testing errors
CV.rkhs.regression = function(y, X, lambda, zeta, delta = 20){
  N = dim(X)[1]
  p = dim(X)[2]
  even_indexes = seq(2, N, 2)
  odd_indexes = seq(1, N, 2)
  train.X = X[odd_indexes,]
  train.y = y[odd_indexes]
  validation.X = X[even_indexes,]
  validation.y = y[even_indexes]
  s_intervals = seeded.intervals(length(train.y), 20)
  init_train = SBS.rkhs(train.y, train.X, lambda, 1, length(train.y), s_intervals[,1], s_intervals[,2], delta)
  init_train_cpt = thresholdBS(init_train, zeta)$cpt_hat[,1]
  if(length(init_train_cpt) >= 1){
    init_train_cpt_long = c(0, init_train_cpt, nrow(train.X))
    train_error = 0
    test_error = 0
    init_train_beta = NULL
    for(k in 1:(length(init_train_cpt)+1)){
      train.yc_temp = scale(train.y[(init_train_cpt_long[k]+1):(init_train_cpt_long[k+1])], center = TRUE, scale = FALSE)
      train.Xc_temp = scale(train.X[(init_train_cpt_long[k]+1):(init_train_cpt_long[k+1]),], center = TRUE, scale = FALSE)
      test.yc_temp = scale(validation.y[(init_train_cpt_long[k]+1):(init_train_cpt_long[k+1])], center = TRUE, scale = FALSE)
      test.Xc_temp = scale(validation.X[(init_train_cpt_long[k]+1):(init_train_cpt_long[k+1]),], center = TRUE, scale = FALSE)
      train_temp = rkhs_reg(train.yc_temp, train.Xc_temp, lambda)
      train_error = train_error + train_temp$RSS
      test_error = test_error + sum((test.yc_temp - test.Xc_temp%*%train_temp$beta_hat_grid/p)^2)
      init_train_beta = cbind(init_train_beta, train_temp$beta_hat_grid)
    }
    init_cpt = odd_indexes[init_train_cpt]
    K_hat = length(init_train_cpt)
  }else{
    init_cpt = init_train_cpt
    K_hat = 0
    train.yc_temp = scale(train.y, center = TRUE, scale = FALSE)
    train.Xc_temp = scale(train.X, center = TRUE, scale = FALSE)
    test.yc_temp = scale(validation.y, center = TRUE, scale = FALSE)
    test.Xc_temp = scale(validation.X, center = TRUE, scale = FALSE)
    train_temp = rkhs_reg(train.yc_temp, train.Xc_temp, lambda)
    train_error = train_temp$RSS
    test_error = sum((test.yc_temp - test.Xc_temp%*%train_temp$beta_hat_grid/p)^2)
    init_train_beta = train_temp$beta_hat_grid
  }
  result = list(cpt_hat = init_cpt, K_hat, test_error = test_error, train_error = train_error, beta_hat = init_train_beta)
  return(result)
}

## Cross validation to select the tuning parameters (lambda and zeta)
CV.search.rkhs.regression = function(y, X, lambda_set, zeta_set, delta = 20){
  output = sapply(1:length(lambda_set), function(i) sapply(1:length(zeta_set), 
                                                           function(j) CV.rkhs.regression(y, X, lambda_set[i], zeta_set[j])))
  cpt_hat = output[seq(1,5*length(zeta_set),5),]## estimated change points
  K_hat = output[seq(2,5*length(zeta_set),5),]## number of estimated change points
  test_error = output[seq(3,5*length(zeta_set),5),]## validation loss
  train_error = output[seq(4,5*length(zeta_set),5),]## training loss
  beta_hat = output[seq(5,5*length(zeta_set),5),]
  result = list(cpt_hat = cpt_hat, K_hat = K_hat, test_error = test_error, train_error = train_error, beta_hat = beta_hat)
  return(result)
}

## Local refinement based on FRBS
local.refine.RKHS.regression = function(cpt_init, beta_hat, y, X, lambda, w = 0.9){
  n = nrow(X)
  cpt_init_long = c(0, cpt_init, n)
  cpt_init_numb = length(cpt_init)
  cpt_refined = rep(0, cpt_init_numb+1)
  for (k in 1:cpt_init_numb){
    s = w*cpt_init_long[k] + (1-w)*cpt_init_long[k+1]
    e = (1-w)*cpt_init_long[k+1] + w*cpt_init_long[k+2]
    lower = ceiling(s) + 2
    upper = floor(e) - 2
    
    
    b = sapply(lower:upper, function(eta)(lassoDPDU_error(y[ceiling(s):eta], cbind(rep(1, n), X)[ceiling(s):eta,], beta_hat[,k]) + lassoDPDU_error(y[(eta+1):floor(e)], cbind(rep(1, n), X)[(eta+1):floor(e),], beta_hat[,k+1])))
    cpt_refined[k+1] = ceiling(s) + which.min(b)
  }
  return(cpt_refined[-1])
}


## estimate the jump size
kappa2_rkhs_est = function(y_vec, x_mat, cpt_rkhs_hat){
  n = length(y_vec)
  p = dim(x_mat)[2]
  Khat_rkhs = length(cpt_rkhs_hat)
  cpt_rkhs_hat_long = c(0, cpt_rkhs_hat, n)
  K_mat = matrix(0, p, p)
  for(j in 1:p){
    for(i in 1:j){
      K_mat[i,j] = cosh(grid[i])*cosh(1-grid[j])/sinh(1)
    }
  }
  K_mat = K_mat + t(K_mat) - diag(diag(K_mat))
  beta_hat_mat = matrix(NA, p, Khat_rkhs+1)
  kappa2_hat_vec = rep(NA, Khat_rkhs)
  for(i in 1:(Khat_rkhs+1)){
    y_vec_c = scale(y_vec[(cpt_rkhs_hat_long[i]+1):cpt_rkhs_hat_long[i+1]], center = TRUE, scale = FALSE)
    x_mat_c = scale(x_mat[(cpt_rkhs_hat_long[i]+1):cpt_rkhs_hat_long[i+1],], center = TRUE, scale = FALSE)
    beta_hat_mat[,i] = rkhs_reg(y_vec_c, x_mat_c, lambda_rkhs_CV)$beta_hat_grid
  }
  for(i in 1:Khat_rkhs){
    kappa2_hat_vec[i] = (1/p^2)*t(beta_hat_mat[,(i+1)] - beta_hat_mat[,i])%*%K_mat%*%(beta_hat_mat[,(i+1)] - beta_hat_mat[,i])
  }
  return(list(kappa2=kappa2_hat_vec, beta=beta_hat_mat))
}

## Simulate two sided Brownian Motion with drift
simu.2BM_Drift = function(n, drift, LRV){
  z_vec = rnorm(2*n)
  w_vec = c(rev(cumsum(z_vec[n:1])/sqrt(1:n)), 0, cumsum(z_vec[(n+1):(2*n)])/sqrt(1:n))
  v_vec = drift*abs(seq(-n, n)) + sqrt(LRV)*w_vec
  return(v_vec)
}

## LRV estimation
LRV.rkhs = function(cpt_init, beta_hat_mat, y_vec, x_mat, w = 0.9, block_size){
  n = nrow(x_mat)
  p = ncol(x_mat)
  cpt_init_long = c(0, cpt_init, n)
  cpt_init_numb = length(cpt_init)
  lrv_hat = rep(NA, cpt_init_numb)
  xc_mat = scale(x_mat, center = TRUE, scale = FALSE)
  cov_mat_hat = t(xc_mat)%*%xc_mat/n
  for (k in 1:cpt_init_numb){
    kappa2_hat = as.numeric((1/p^2)*t(beta_hat_mat[,(k+1)] - beta_hat_mat[,k])%*%cov_mat_hat%*%(beta_hat_mat[,(k+1)] - beta_hat_mat[,k]))
    s = w*cpt_init_long[k] + (1-w)*cpt_init_long[k+1]
    e = (1-w)*cpt_init_long[k+1] + w*cpt_init_long[k+2]
    z_vec = rep(NA, floor(e)-ceiling(s)+1)
    for (t in ceiling(s):floor(e)){
      z_vec[t-ceiling(s)+1] = as.numeric(2*y_vec[t] - (1/p)*crossprod(x_mat[t,], beta_hat_mat[,k]+beta_hat_mat[,k+1])) * (1/p)* crossprod(x_mat[t,], beta_hat_mat[,k+1]-beta_hat_mat[,k])
    }
    pair_numb = floor((floor(e)-ceiling(s)+1)/(2*block_size))
    z_mat1 = matrix(z_vec[1:(block_size*pair_numb*2)], nrow = block_size)
    z_mat1_colsum = apply(z_mat1, 2, sum)
    z_mat2 = matrix(rev(z_vec)[1:(block_size*pair_numb*2)], nrow = block_size)
    z_mat2_colsum = apply(z_mat2, 2, sum)
    lrv_hat[k] = (mean((z_mat1_colsum[2*(1:pair_numb)-1] - z_mat1_colsum[2*(1:pair_numb)])^2/(2*block_size)) + mean((z_mat2_colsum[2*(1:pair_numb)-1] - z_mat2_colsum[2*(1:pair_numb)])^2/(2*block_size)))/(2*kappa2_hat)
  }
  return(list(kappa2_hat = kappa2_hat, lrv_hat = lrv_hat))
}


# Simulations
## one changepoint
# Observations (n = 200, p = 100)
# generate responses (with errors N(0,1))
n = 200 # number of time points
p = 100 # number of observations on each curve
K = 50  # number of bases
m = 10000 # number of grids (approximate the integrals)

u_vec = sapply(1:K, function(x){(-1)^(x+1)*x^(-2/2)}) # vector of sqrt eigenvalues of covariance function of Xi (up to a sign, K largest)
# impose temporal dependence on Xi (each row of v_mat is a AR1 process with mean zero and marginal variance 1)
v_mat = matrix(NA, K, n)
set.seed(123)
for(k in 1:K){
  v_mat[k,] = arima.sim(list(ar = 0.7), sd = sqrt(1-0.7^2), n = n)
}
# approximate the integrals
grid = seq(0,1,length=m) # assume evenly spacing
phi_mat = matrix(NA, K, m)
beta1_vec = rep(NA, m)
beta2_vec = rep(NA, m)
for (k in 1:K){
  phi_mat[k,] = sapply(grid, function(x){phi_fct(k,x)})
}
beta1_vec = sapply(grid, function(x){beta1_fct(x)})
beta2_vec = sapply(grid, function(x){beta2_fct(x)})

par(mfrow = c(2,1))
plot(grid, beta1_vec, type = "l", col = "blue") # 
plot(grid, beta2_vec, type = "l", col = "red")
par(mfrow = c(1,1))

integ1_vec = (1/m) * phi_mat %*% beta1_vec
integ2_vec = (1/m) * phi_mat %*% beta2_vec

set.seed(1234)
y1_vec = sapply(1:(n/2), function(i){sum(sapply(1:K, function(x){v_mat[x,i]*u_vec[x]*integ1_vec[x]}))}) + rnorm(n/2)
y2_vec = sapply((n/2+1):n, function(i){sum(sapply(1:K, function(x){v_mat[x,i]*u_vec[x]*integ2_vec[x]}))}) + rnorm(n/2)
y_vec = c(y1_vec, y2_vec)
# generate discretized Xi
x_mat = matrix(NA, n, p)
sample_mat = matrix(rep(seq(0,1,len=p),each=n),nrow=n) # fixed evenly spaced grids

for (i in 1:n){
  x_mat[i,] = sapply(sample_mat[i,], function(x){x_fct(i,x)})
}

lambda_rkhs_set = c(0.05, 0.1, 0.3, 0.5, 1, 2)
zeta_rkhs_set = c(5,10,15,20)
temp_rkhs_CV = CV.search.rkhs.regression(y_vec, x_mat, lambda_rkhs_set, zeta_rkhs_set)
s_intervals = seeded.intervals(n, 20)
min_idx = as.vector(arrayInd(which.min(temp_rkhs_CV$test_error), dim(temp_rkhs_CV$test_error))) 
lambda_rkhs_CV = lambda_rkhs_set[min_idx[2]]
zeta_rkhs_CV = zeta_rkhs_set[min_idx[1]]
temp_rkhs = SBS.rkhs(y_vec, x_mat, lambda_rkhs_CV, 1, n, s_intervals[,1], s_intervals[,2], 20)
cpt_flr_hat = thresholdBS(temp_rkhs, zeta_rkhs_CV)

temp_reg = CV.search.DPDU.regression(y_vec, x_mat, c(0.05, 0.1, 0.3, 0.5, 1, 2), c(1,3,5))
min_idx = as.vector(arrayInd(which.min(temp_reg$test_error), dim(temp_reg$test_error)))
cpt_reg_hat = unlist(temp_reg$cpt_hat[min_idx[1], min_idx[2]])



## two changepoint
# Observations (n = 300, p = 100)
# generate responses (with errors N(0,1))
n = 300 # number of time points
p = 100 # number of observations on each curve
K = 50  # number of bases
m = 10000 # number of grids (approximate the integrals)

u_vec = sapply(1:K, function(x){(-1)^(x+1)*x^(-2/2)}) # vector of sqrt eigenvalues of covariance function of Xi (up to a sign, K largest)
# impose temporal dependence on Xi (each row of v_mat is a AR1 process with mean zero and marginal variance 1)
v_mat = matrix(NA, K, n)
set.seed(123)
for(k in 1:K){
  v_mat[k,] = arima.sim(list(ar = 0.7), sd = sqrt(1-0.7^2), n = n)
}
# approximate the integrals
grid = seq(0,1,length=m) # assume evenly spacing
phi_mat = matrix(NA, K, m)
beta1_vec = rep(NA, m)
beta2_vec = rep(NA, m)
for (k in 1:K){
  phi_mat[k,] = sapply(grid, function(x){phi_fct(k,x)})
}
beta1_vec = sapply(grid, function(x){beta1_fct(x)})
beta2_vec = sapply(grid, function(x){beta2_fct(x)})

par(mfrow = c(2,1))
plot(grid, beta1_vec, type = "l", col = "blue") # 
plot(grid, beta2_vec, type = "l", col = "red")
par(mfrow = c(1,1))

integ1_vec = (1/m) * phi_mat %*% beta1_vec
integ2_vec = (1/m) * phi_mat %*% beta2_vec
set.seed(123)
y1_vec = sapply(1:(n/3), function(i){sum(sapply(1:K, function(x){v_mat[x,i]*u_vec[x]*integ1_vec[x]}))}) + rnorm(n/3)
y2_vec = sapply((n/3+1):(2*n/3), function(i){sum(sapply(1:K, function(x){v_mat[x,i]*u_vec[x]*integ2_vec[x]}))}) + rnorm(n/3)
y3_vec = sapply((2*n/3+1):n, function(i){sum(sapply(1:K, function(x){v_mat[x,i]*u_vec[x]*integ1_vec[x]}))}) + rnorm(n/3)
y_vec = c(y1_vec, y2_vec, y3_vec)
# generate discretized Xi
x_mat = matrix(NA, n, p)
sample_mat = matrix(rep(seq(0,1,len=p),each=n),nrow=n) # fixed evenly spaced grids

for (i in 1:n){
  x_mat[i,] = sapply(sample_mat[i,], function(x){x_fct(i,x)})
}

lambda_rkhs_set = c(0.05, 0.1, 0.3, 0.5, 1, 2)
zeta_rkhs_set = c(10,15,20,25)
temp_rkhs_CV = CV.search.rkhs.regression(y_vec, x_mat, lambda_rkhs_set, zeta_rkhs_set)
s_intervals = seeded.intervals(n, 20)
min_idx = as.vector(arrayInd(which.min(temp_rkhs_CV$test_error), dim(temp_rkhs_CV$test_error))) 
lambda_rkhs_CV = lambda_rkhs_set[min_idx[2]]
zeta_rkhs_CV = zeta_rkhs_set[min_idx[1]]
cpt_rkhs_hat = unlist(temp_rkhs_CV$cpt_hat[min_idx[1], min_idx[2]])
temp_rkhs = SBS.rkhs(y_vec, x_mat, lambda_rkhs_CV, 1, n, s_intervals[,1], s_intervals[,2], 20)
cpt_flr_hat = thresholdBS(temp_rkhs, zeta_rkhs_CV)
Khat_rkhs = length(cpt_rkhs_hat)

temp_reg = CV.search.DPDU.regression(y_vec, x_mat, c(0.05, 0.1, 0.3, 0.5, 1, 2), c(1,3,5))
min_idx = as.vector(arrayInd(which.min(temp_reg$test_error), dim(temp_reg$test_error)))
cpt_reg_hat = unlist(temp_reg$cpt_hat[min_idx[1], min_idx[2]])

