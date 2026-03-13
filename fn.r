arcoeff_fun <-
function(coef = 1) {
  function(s1, s2) {
    s1n <- .norm01(s1); s2n <- .norm01(s2)
    
    raw <- switch(as.integer(coef),
                  # 1) smooth linear plane (keep)
                  `1` = 0.8*(s1n - 0.5) + 0.6*(s2n - 0.5),
                  
                  # 2) smooth radial hill (Gaussian bump)
                  `2` = 1.5 * exp(-((s1n-0.5)^2 + (s2n-0.5)^2)/0.20) - 0.5,
                  
                  # 3) smooth sinusoidal surface (low frequency)
                  `3` = 0.9*sin(2*pi*s1n)*cos(2*pi*s2n),
                  
                  # 4) smooth polynomial surface
                  `4` = 0.7*(s1n^2 + s2n^2) - 0.8*(s1n*s2n) - 0.15,
                  
                  # 5) smooth mixture of sin/cos + interaction
                  `5` = 0.6*sin(pi*s1n) + 0.5*cos(2*pi*s2n) - 0.35*(s1n-0.5)*(s2n-0.5)
    )
    
    0.9 * tanh(raw)  # keep in (-0.9,0.9)
  }
}
arcoeff_fun_2 <-
function(coef = 1) {
  function(s1, s2) {
    s1 <- as.numeric(s1); s2 <- as.numeric(s2)
    s1n <- .norm01(s1); s2n <- .norm01(s2)
    
    # ---- radius helpers ----
    r_center <- sqrt((s1n - 0.5)^2 + (s2n - 0.5)^2)
    r_center_sc <- r_center / sqrt(0.5^2 + 0.5^2)   # [0,1]
    
    # case 7: quarter concentric circles with origin at lower-left corner (0,0)
    r_ll <- sqrt((s1n - 0.0)^2 + (s2n - 0.0)^2)      # [0, sqrt(2)]
    r_ll_sc <- r_ll / sqrt(2)                        # [0,1]
    
    raw1 <- switch(as.integer(coef),
             `1` =  2*(s1n - 0.5),
             `2` =  2.0*sin(pi*s1n) + 0.4*cos(2*pi*s2n),
             `3` =  0.9*sin(2*pi*s1n)*cos(pi*s2n),
            `4` =  2*(s1n^2 - 0.5*s1n*s2n) - 0.25,
                   `5` =  0.8*cos(2*pi*s1n) + 0.25*sin(2*pi*s2n),
                   `6` =  cos(2*pi*r_center_sc),              
                   `7` =  1.1*cos(2*pi*r_ll_sc^1.5),
                   `8` =  2*(0.9*(s1n - 0.5) + 0.9*(s2n - 0.5) +
                               0.35*sin(2*pi*(s1n + 0.35*s2n)) +
                               0.25*cos(2*pi*(0.6*s1n - s2n))),
                   `9` =  1.2*(s1n - 0.5) - 1.0*(s2n - 0.5) +
                     0.90*sin(6*pi*(s1n + 0.23*s2n)) +
                     0.70*cos(8*pi*(0.55*s1n - 0.85*s2n)) +
                     0.55*sin(10*pi*(s1n*s2n + 0.15)),
                   
                   `10` =  1.4*(s1n - 0.5) - 1.2*(s2n - 0.5) +
                     1.10*sin(10*pi*(s1n + 0.27*s2n + 0.12*s1n*s2n)) +
                     0.95*cos(12*pi*(0.62*s1n - 0.91*s2n + 0.08*s1n^2)) +
                     0.80*sin(14*pi*(s1n*s2n + 0.10)) +
                     0.55*exp(-((s1n-0.25)^2 + (s2n-0.75)^2)/0.010) -
                     0.50*exp(-((s1n-0.78)^2 + (s2n-0.30)^2)/0.008)
    )
    
    raw2 <- switch(as.integer(coef),
                   `1` =  1*(s2n - 0.5),
                   `2` =  0.9*cos(pi*s2n) + 0.2*sin(2*pi*s1n),
                   `3` =  0.8*cos(2*pi*s2n)*sin(pi*s1n),
                   `4` =  1.0*(s2n^2 - 0.5*s1n*s2n) - 0.20,
                   `5` =  0.7*sin(2*pi*s2n) - 0.25*cos(2*pi*s1n),
                   `6` =  1*(s1n - 0.5) + 0.6*(s2n - 0.5),      
                   `7` =  1.1*(s1n - 0.5) + 1.1*(s2n - 0.5),
                   `8` =  1.0*(s2n - 0.5) - 0.8*(s1n - 0.5) +
                     0.30*cos(2*pi*(s2n + 0.25*s1n)) +
                     0.20*sin(2*pi*(0.8*s1n + 0.7*s2n)),
                   `9` =  0.8*(s2n - 0.5) + 0.6*(s1n - 0.5) +
                     0.85*cos(6*pi*(s2n + 0.31*s1n)) +
                     0.65*sin(8*pi*(0.75*s1n + 0.40*s2n)) +
                     0.50*cos(10*pi*(s1n*s2n + 0.10)),
                   `10` =  1.0*(s2n - 0.5) + 0.9*(s1n - 0.5) +
                     1.05*cos(10*pi*(s2n + 0.33*s1n + 0.10*s1n*s2n)) +
                     0.90*sin(12*pi*(0.80*s1n + 0.45*s2n + 0.06*s2n^2)) +
                     0.75*cos(14*pi*(s1n*s2n + 0.12)) -
                     0.55*exp(-((s1n-0.30)^2 + (s2n-0.65)^2)/0.012) +
                     0.50*exp(-((s1n-0.82)^2 + (s2n-0.22)^2)/0.009)
    )
    
    k1 <- 0.85 * tanh(raw1)
    k2 <- 0.75 * tanh(raw2)
    
    a2 <- k2
    a1 <- k1 * (1 - k2)
    cbind(a1, a2)
  }
}
B_spline_1d <-
function(omega, K = 30, degree = 3) {
  splines::bs(omega, df = K, degree = degree, intercept = TRUE)
}
B_spline_2d <-
function(omega_mat, K1 = 12, K2 = 12, degree = 3) {
  w1 <- omega_mat[,1]; w2 <- omega_mat[,2]
  B1 <- splines::bs(w1, df = K1, degree = degree, intercept = TRUE)
  B2 <- splines::bs(w2, df = K2, degree = degree, intercept = TRUE)
  n <- nrow(omega_mat)
  K <- ncol(B1) * ncol(B2)
  B <- matrix(0, n, K)
  for (i in 1:n) B[i,] <- as.vector(B1[i,] %o% B2[i,])
  B
}
build_W <-
function(location, dia, gamma) {
  m <- nrow(location)
  D <- as.matrix(dist(location))
  W <- matrix(0, m, m)
  idx <- which(D <= dia & D > 0, arr.ind = TRUE)
  if (nrow(idx) > 0) {
    W[idx] <- exp(-(D[idx]^2) / (2*gamma^2))
  }
  diag(W) <- 0
  W
}
CalculateGVC <-
function(f, f.hat, bandwidth){
  M = length(f)
  sum = 0
  q = c(0.5,rep(1,M-2),0.5)
  for(i in 1:M){
    num = -log(f[i]/f.hat[i])+(f[i]-f.hat[i])/f.hat[i]
    dem = (1 - (1/(2*bandwidth + 1)))^2
    sum = sum + q[i]*(num/dem)
  }
  return(sum/M)
}
choose_lambda_cv <-
function(I_mat, B, R, location,
                             r, dia, gamma,
                             lambda1_range, lambda2_range,
                             n_folds = 5,
                             cd_maxit = 5, cd_tol = 1e-2,
                             n_iter_fit = 25, step_A = 1, step_Theta = 1e-3,
                             Ucap = 30, verbose = TRUE) {
  

  m <- ncol(I_mat)
  folds <- sample(rep(1:n_folds, length.out = m))
  
  cv_mean <- function(lambda1, lambda2) {
    lambda1 <- max(lambda1_range[1], min(lambda1, lambda1_range[2]))
    lambda2 <- max(lambda2_range[1], min(lambda2, lambda2_range[2]))
    mean(sapply(1:n_folds, function(k)
      cv_score(I_mat, B, R, location, folds, k,
               r=r, dia=dia, gamma=gamma,
               lambda1=lambda1, lambda2=lambda2,
               n_iter_fit=n_iter_fit, step_A=step_A, step_Theta=step_Theta, Ucap=Ucap)
    ))
  }
  
  # coordinate descent with 1D optimize in transformed scale
  eps0 <- 1e-12
  lambda1 <- sqrt(lambda1_range[1]*lambda1_range[2])
  lambda2 <- sqrt(max(lambda2_range[1],0)+eps0) * sqrt(lambda2_range[2]+eps0) - eps0
  trace <- data.frame(iter=integer(), lambda1=double(), lambda2=double(), cv=double())
  
  opt_lambda1 <- function(l2) {
    f <- function(t) cv_mean(10^t, l2)
    optimize(f, interval = c(log10(lambda1_range[1]), log10(lambda1_range[2])))
  }
  opt_lambda2 <- function(l1) {
    f <- function(t) cv_mean(l1, exp(t)-eps0)
    optimize(f, interval = c(log(lambda2_range[1]+eps0), log(lambda2_range[2]+eps0)))
  }
  
  for (it in 1:cd_maxit) {
    r1 <- opt_lambda1(lambda2); l1_new <- 10^(r1$minimum)
    r2 <- opt_lambda2(l1_new);  l2_new <- exp(r2$minimum) - eps0
    l2_new <- max(lambda2_range[1], min(l2_new, lambda2_range[2]))
    cvv <- cv_mean(l1_new, l2_new)
    
    trace <- rbind(trace, data.frame(iter=it, lambda1=l1_new, lambda2=l2_new, cv=cvv))
    if (verbose) cat(sprintf("[CV] it=%d lambda1=%.3e lambda2=%.3e CV=%.6f\n", it, l1_new, l2_new, cvv))
    
    rel <- max(abs(log(l1_new+eps0)-log(lambda1+eps0)), abs(log(l2_new+eps0)-log(lambda2+eps0)))
    lambda1 <- l1_new; lambda2 <- l2_new
    if (rel < cd_tol) break
  }
  list(lambda1=lambda1, lambda2=lambda2, cv_trace=trace)
}
clamp <-
function(x, lo, hi) pmin(pmax(x, lo), hi)
cv_score <-
function(I_mat, B, R, location,
                     folds, fold_id,
                     r, dia, gamma,
                     lambda1, lambda2,
                     n_iter_fit, step_A, step_Theta, Ucap) {
  test_idx <- which(folds == fold_id)
  train_idx <- setdiff(seq_len(ncol(I_mat)), test_idx)
  
  I_train <- I_mat[,train_idx, drop=FALSE]
  loc_train <- location[train_idx,, drop=FALSE]
  loc_test  <- location[test_idx,, drop=FALSE]
  
  W_train <- build_W(loc_train, dia = dia, gamma = gamma)
  
  fit <- fit_core(I_train, B, R, W_train,
                  r = r,
                  lambda1 = lambda1, lambda2 = lambda2,
                  n_iter = n_iter_fit,
                  step_A = step_A,
                  step_Theta = step_Theta,
                  Ucap = Ucap,
                  verbose = FALSE)
  
  Theta_hat <- fit$Theta
  A_train_hat <- fit$A
  G <- B %*% Theta_hat
  
  score <- 0
  for (k in seq_along(test_idx)) {
    idx <- test_idx[k]
    Ivec <- I_mat[, idx]
    a_hat <- fit_a_test(Ivec, B, Theta_hat,
                        A_train_hat, loc_train, loc_test[k,],
                        dia = dia, gamma = gamma, lambda2 = lambda2,
                        n_iter = 30, step = 1, Ucap = Ucap)
    u <- clamp(as.vector(G %*% a_hat), -Ucap, Ucap)
    score <- score + sum(u + Ivec * exp(-u))
  }
  score
}
deg_vec <-
function(W) rowSums(W)
detect_dim <-
function(Xs) {
  x1 <- Xs[[1]]
  if (is.matrix(x1) && length(dim(x1)) == 2) return(2L)
  if (is.array(x1)  && length(dim(x1)) == 2) return(2L)
  1L
}
fftshift_mat <-
function(M) {
  # shift zero-frequency to center for 2D matrix
  n1 <- nrow(M); n2 <- ncol(M)
  i1 <- c((floor(n1/2)+1):n1, 1:floor(n1/2))
  i2 <- c((floor(n2/2)+1):n2, 1:floor(n2/2))
  M[i1, i2, drop = FALSE]
}
fit_a_test <-
function(I_vec, B, Theta_hat,
                       A_train, loc_train, loc_test,
                       dia, gamma, lambda2,
                       n_iter = 30, step = 1, Ucap = 30) {
  # neighbor weights from test point to train points
  d <- sqrt((loc_train[,1]-loc_test[1])^2 + (loc_train[,2]-loc_test[2])^2)
  w <- ifelse(d <= dia, exp(-(d^2)/(2*gamma^2)), 0)
  dsum <- sum(w)
  
  G <- B %*% Theta_hat
  r <- ncol(Theta_hat)
  
  a <- if (dsum > 0) colSums(A_train * w)/dsum else rep(0, r)
  
  for (it in 1:n_iter) {
    u <- clamp(as.vector(G %*% a), -Ucap, Ucap)
    s <- as.vector(I_vec * exp(-u))
    grad_w <- as.vector(crossprod(G, (1 - s)))
    
    if (lambda2 > 0 && dsum > 0) {
      sum_w_A <- colSums(A_train * w)
      grad_sp <- 4*lambda2*(dsum*a - sum_w_A)
      H_sp <- (4*lambda2*dsum) * diag(r)
    } else {
      grad_sp <- rep(0, r)
      H_sp <- matrix(0, r, r)
    }
    
    Gs <- G * sqrt(pmax(s, 0))
    H <- crossprod(Gs) + H_sp + ridge_eps()*diag(r)
    delta <- solve(H, grad_w + grad_sp)
    a_new <- a - step*as.vector(delta)
    if (any(!is.finite(a_new))) break
    a <- a_new
  }
  a
}
fit_core <-
function(I_mat, B, R, W,
                     r = 2,
                     lambda1, lambda2,
                     n_iter = 30,
                     step_A = 1,
                     step_Theta = 1e-3,
                     reid_every = 5,
                     Ucap = 30,
                     verbose = FALSE) {
  logI <- log(I_mat)
  init <- init_svd(B, logI, r = r)
  Theta <- init$Theta
  A <- init$A
  
  obj_trace <- numeric(n_iter)
  for (it in 1:n_iter) {
    A <- update_A(B, Theta, A, I_mat, W, lambda2, step_A = step_A, n_inner = 1, Ucap = Ucap)
    Theta <- update_Theta(B, Theta, A, I_mat, R, lambda1, step_Theta = step_Theta, Ucap = Ucap)
    
    if (reid_every > 0 && it %% reid_every == 0) {
      tmp <- reidentify(Theta, A)
      Theta <- tmp$Theta; A <- tmp$A
    }
    
    obj_trace[it] <- obj_srcse(B, Theta, A, I_mat, R, W, lambda1, lambda2, Ucap = Ucap)
    if (verbose) cat(sprintf("iter=%d obj=%.6f\n", it, obj_trace[it]))
    if (!is.finite(obj_trace[it])) break
    if (it >= 2 && abs(obj_trace[it]-obj_trace[it-1])/(abs(obj_trace[it-1])+1e-8) < 1e-6) {
      obj_trace <- obj_trace[1:it]; break
    }
  }
  
  Uhat <- clamp(B %*% Theta %*% t(A), -Ucap, Ucap)
  sdfhat <- exp(Uhat)
  list(sdfhat = sdfhat, A = A, Theta = Theta, obj_trace = obj_trace)
}
fit_sdf <-
function(Xs, location,
                    r = 2,
                    # basis params
                    K = 30, K1 = 12, K2 = 12,
                    degree = 3, diff_order = 2,
                    # spatial weights
                    dia = 3, gamma = 1,
                    # optimization controls
                    n_iter = 30, step_A = 1, step_Theta = 1e-3, Ucap = 30,
                    # lambda selection
                    cv = TRUE,
                    lambda1_range = c(1e-6, 1e0),
                    lambda2_range = c(0,   1e1),
                    lambda1 = 1e-3, lambda2 = 1,
                    n_folds = 5,
                    # plotting A
                    plotA = TRUE,
                    interp_lambda = NULL,
                    verbose = TRUE) {
  
  stopifnot(is.list(Xs))
  m <- length(Xs)
  stopifnot(nrow(location) == m, ncol(location) == 2)
  
  d <- detect_dim(Xs)
  
  # Build I_mat and omega + basis
  if (d == 1L) {
    p1 <- periodogram_1d(Xs[[1]])
    omega <- p1$omega
    n_omega <- length(omega)
    I_mat <- matrix(NA_real_, n_omega, m)
    I_mat[,1] <- p1$I
    for (i in 2:m) I_mat[,i] <- periodogram_1d(Xs[[i]])$I
    
    B <- B_spline_1d(omega, K = K, degree = degree)
    R <- R_diff_1d(ncol(B), diff_order = diff_order)
    
  } else {
    p1 <- periodogram_2d(Xs[[1]], drop_zero = TRUE)
    omega <- p1$omega
    n_omega <- nrow(omega)
    I_mat <- matrix(NA_real_, n_omega, m)
    I_mat[,1] <- p1$I
    for (i in 2:m) I_mat[,i] <- periodogram_2d(Xs[[i]], drop_zero = TRUE)$I
    
    B <- B_spline_2d(omega, K1 = K1, K2 = K2, degree = degree)
    # effective sizes for R
    K1_eff <- ncol(splines::bs(omega[,1], df=K1, degree=degree, intercept=TRUE))
    K2_eff <- ncol(splines::bs(omega[,2], df=K2, degree=degree, intercept=TRUE))
    R <- R_diff_2d(K1_eff, K2_eff, diff_order = diff_order)
  }
  
  W <- build_W(location, dia = dia, gamma = gamma)
  
  # choose lambdas
  cv_trace <- NULL
  if (cv) {
    if (verbose) cat("[fit_sdf] CV=TRUE: choosing lambda1/lambda2 by replicate-holdout CV...\n")
    sel <- choose_lambda_cv(I_mat, B, R, location,
                            r=r, dia=dia, gamma=gamma,
                            lambda1_range=lambda1_range,
                            lambda2_range=lambda2_range,
                            n_folds=n_folds,
                            n_iter_fit=max(10, round(n_iter*0.8)),
                            step_A=step_A, step_Theta=step_Theta,
                            Ucap=Ucap, verbose=verbose)
    lambda1 <- sel$lambda1
    lambda2 <- sel$lambda2
    cv_trace <- sel$cv_trace
  } else {
    if (verbose) cat(sprintf("[fit_sdf] CV=FALSE: using lambda1=%.3e lambda2=%.3e\n", lambda1, lambda2))
  }
  
  # final fit on full
  fit <- fit_core(I_mat, B, R, W,
                  r=r, lambda1=lambda1, lambda2=lambda2,
                  n_iter=n_iter, step_A=step_A, step_Theta=step_Theta,
                  Ucap=Ucap, verbose=FALSE)
  
  if (plotA) {
    plotA_imageplot(location, fit$A, interp_lambda = interp_lambda, main_prefix = "A_hat")
  }
  
  list(
    d = d,
    omega = omega,
    location = location,
    sdfhat = fit$sdfhat,
    A = fit$A,
    Theta = fit$Theta,
    lambda1 = lambda1,
    lambda2 = lambda2,
    cv_trace = cv_trace,
    B = B,
    I = I_mat
  )
}
fit_sdf2 <-
function(Xs, location,
                    r = 2,
                    # basis params
                    K = 30, K1 = 12, K2 = 12,
                    degree = 3, diff_order = 2,
                    # spatial weights
                    dia = 3, gamma = 1,
                    # optimization controls
                    n_iter = 30, step_A = 1, step_Theta = 1e-3, Ucap = 30,
                    # lambda selection
                    cv = TRUE,
                    lambda1_range = c(1e-6, 1e0),
                    lambda2_range = c(0,   1e1),
                    lambda1 = 1e-3, lambda2 = 1,
                    n_folds = 5,
                    # plotting A
                    plotA = TRUE,
                    interp_lambda = NULL,
                    verbose = TRUE) {
  lambda2 = 0
  stopifnot(is.list(Xs))
  m <- length(Xs)
  stopifnot(nrow(location) == m, ncol(location) == 2)
  
  d <- detect_dim(Xs)
  
  # Build I_mat and omega + basis
  if (d == 1L) {
    p1 <- periodogram_1d(Xs[[1]])
    omega <- p1$omega
    n_omega <- length(omega)
    I_mat <- matrix(NA_real_, n_omega, m)
    I_mat[,1] <- p1$I
    for (i in 2:m) I_mat[,i] <- periodogram_1d(Xs[[i]])$I
    
    B <- B_spline_1d(omega, K = K, degree = degree)
    R <- R_diff_1d(ncol(B), diff_order = diff_order)
    
  } else {
    p1 <- periodogram_2d(Xs[[1]], drop_zero = TRUE)
    omega <- p1$omega
    n_omega <- nrow(omega)
    I_mat <- matrix(NA_real_, n_omega, m)
    I_mat[,1] <- p1$I
    for (i in 2:m) I_mat[,i] <- periodogram_2d(Xs[[i]], drop_zero = TRUE)$I
    
    B <- B_spline_2d(omega, K1 = K1, K2 = K2, degree = degree)
    # effective sizes for R
    K1_eff <- ncol(splines::bs(omega[,1], df=K1, degree=degree, intercept=TRUE))
    K2_eff <- ncol(splines::bs(omega[,2], df=K2, degree=degree, intercept=TRUE))
    R <- R_diff_2d(K1_eff, K2_eff, diff_order = diff_order)
  }
  
  W <- build_W(location, dia = dia, gamma = gamma)
  
  # choose lambdas
  cv_trace <- NULL
  if (cv) {
    if (verbose) cat("[fit_sdf] CV=TRUE: choosing lambda1/lambda2 by replicate-holdout CV...\n")
    sel <- choose_lambda_cv(I_mat, B, R, location,
                            r=r, dia=dia, gamma=gamma,
                            lambda1_range=lambda1_range,
                            lambda2_range=lambda2_range,
                            n_folds=n_folds, 
                            n_iter_fit=max(10, round(n_iter*0.8)),
                            step_A=step_A, step_Theta=step_Theta,
                            Ucap=Ucap, verbose=verbose)
    lambda1 <- sel$lambda1
    lambda2 <- sel$lambda2
    cv_trace <- sel$cv_trace
  } else {
    if (verbose) cat(sprintf("[fit_sdf] CV=FALSE: using lambda1=%.3e lambda2=%.3e\n", lambda1, lambda2))
  }
  
  # final fit on full
  fit <- fit_core(I_mat, B, R, W,
                  r=r, lambda1=lambda1, lambda2=lambda2,
                  n_iter=n_iter, step_A=step_A, step_Theta=step_Theta,
                  Ucap=Ucap, verbose=FALSE)
  
  if (plotA) {
    plotA_imageplot(location, fit$A, interp_lambda = interp_lambda, main_prefix = "A_hat")
  }
  
  list(
    d = d,
    omega = omega,
    location = location,
    sdfhat = fit$sdfhat,
    A = fit$A,
    Theta = fit$Theta,
    lambda1 = lambda1,
    lambda2 = lambda2,
    cv_trace = cv_trace,
    B = B,
    I = I_mat
  )
}
init_svd <-
function(B, logI_mat, r) {
  BtB <- crossprod(B)
  Psi <- solve(BtB + ridge_eps()*diag(ncol(B))) %*% crossprod(B, logI_mat)
  sv <- svd(Psi)
  r_use <- min(r, length(sv$d))
  Theta0 <- sv$u[,1:r_use, drop=FALSE]
  A0 <- sv$v[,1:r_use, drop=FALSE] %*% diag(sv$d[1:r_use], r_use, r_use)
  list(Theta = Theta0, A = A0)
}
localEstimates <-
function(xn,yn,bandwith){
  #boxcar smoother W = 1/(2p+1)
  w=rep(1/(2*bandwith + 1),length(xn))
  f.hat = as.numeric((t(w)%*%yn)/sum(w) )
  #f.hat = sum(yn)/(2*bandwith + 1)
  return(f.hat)
}
make_lattice_location <-
function(m, n, xlim = c(0, 1), ylim = c(0, 1)) {
  stopifnot(m >= 2, n >= 2, length(xlim) == 2, length(ylim) == 2)
  x <- seq(xlim[1], xlim[2], length.out = m)
  y <- seq(ylim[1], ylim[2], length.out = n)
  as.matrix(expand.grid(s1 = x, s2 = y))  # (m*n) x 2
}
mt_smooth_periodogram_2d <-
function(I_mat, nw1 = 3, nw2 = 3, k1 = NULL, k2 = NULL,
                                     circular = TRUE,
                                     eps = .Machine$double.eps) {
  I_mat <- as.matrix(I_mat)
  n1 <- nrow(I_mat); n2 <- ncol(I_mat)
  N <- n1 * n2

  if (is.null(k1)) k1 <- max(1, floor(2 * nw1))
  if (is.null(k2)) k2 <- max(1, floor(2 * nw2))

  # DPSS tapers (multitaper preferred; fallback to sapa)
  get_dpss <- function(n, nw, k) {
    if (requireNamespace("multitaper", quietly = TRUE)) {
      multitaper::dpss(n, nw = nw, k = k)$v
    } else if (requireNamespace("sapa", quietly = TRUE)) {
      sapa::dpss(n, nw = nw, k = k)$v
    } else {
      stop("Install 'multitaper' or 'sapa' to get DPSS tapers.")
    }
  }

  V1 <- get_dpss(n1, nw1, k1)  # n1 x k1
  V2 <- get_dpss(n2, nw2, k2)  # n2 x k2

  # FFT helpers
  fft2 <- function(M) fft(M)
  ifft2 <- function(M) fft(M, inverse = TRUE) / length(M)

  # For circular convolution on the same grid:
  # conv(I, W) = IFFT( FFT(I) * FFT(W) )
  FI <- fft2(I_mat)

  out <- matrix(0, n1, n2)
  K <- k1 * k2

  for (a in 1:k1) {
    for (b in 1:k2) {
      taper2d <- V1[, a] %o% V2[, b]   # n1 x n2
      # spectral window for this taper: |FFT(taper)|^2 / N
      H <- fft2(taper2d)
      W <- (Mod(H)^2) / N

      if (!circular) {
        stop("non-circular convolution not implemented in this short version; set circular=TRUE.")
      } else {
        tmp <- ifft2(FI * fft2(W))
        out <- out + Re(tmp)
      }
    }
  }

  out <- out / K
  out[out < eps] <- eps
  out
}
neighborhood <-
function(x,x.star,y,bandwidth){
  index=c(1:length(x))
  indexn = index[(x<(x.star+bandwidth+1)) & (x>(x.star-bandwidth-1))]
  xn = x[indexn]
  yn = y[indexn]
  out=list(xn,yn)
  names(out) = c("xn", "yn")
  return(out)
}
obj_srcse <-
function(B, Theta, A, I_mat, R, W, lambda1, lambda2, Ucap = 30) {
  U <- B %*% Theta %*% t(A)
  whittle_loss(U, I_mat, Ucap = Ucap) +
    pen1_val(Theta, R, lambda1) +
    pen2_val(A, W, lambda2)
}
omega_signed <-
function(n) {
  j <- 0:(n-1)
  ifelse(j <= floor(n/2), 2*pi*j/n, 2*pi*(j-n)/n)
}
pacf_to_ar2 <-
function(k1, k2) {
  a2 <- k2
  a1 <- k1 * (1 - k2)
  cbind(a1, a2)
}
pen1_val <-
function(Theta, R, lambda1) lambda1 * sum(diag(t(Theta) %*% R %*% Theta))
pen2_val <-
function(A, W, lambda2) {
  if (lambda2 <= 0) return(0)
  m <- nrow(A); val <- 0
  for (i in 1:m) {
    nbr <- which(W[i,] > 0)
    if (length(nbr) > 0) {
      dif <- sweep(A[nbr,,drop=FALSE], 2, A[i,], "-")
      val <- val + sum(W[i,nbr] * rowSums(dif^2))
    }
  }
  lambda2 * val
}
periodogram_1d <-
function(x) {
  x <- as.numeric(x)
  n <- length(x)
  x <- x - mean(x)
  D <- fft(x) / sqrt(n)
  jmax <- floor(n/2)
  # drop frequency 0
  Dpos <- D[2:(jmax+1)]
  omega <- 2*pi*(1:jmax)/n
  I <- (Mod(Dpos)^2)/(2*pi)
  I[I < .Machine$double.eps] <- .Machine$double.eps
  list(omega = omega, I = I)
}
periodogram_2d <-
function(X, drop_zero = TRUE) {
  X <- as.matrix(X)
  n1 <- nrow(X); n2 <- ncol(X)
  N <- n1*n2
  X <- X - mean(X)
  D <- fft(X)/sqrt(N)
  I <- (Mod(D)^2)/((2*pi)^2)
  I[I < .Machine$double.eps] <- .Machine$double.eps
  
  w1 <- omega_signed(n1)
  w2 <- omega_signed(n2)
  omega_mat <- do.call(rbind, lapply(1:n2, function(j2) cbind(w1, rep(w2[j2], n1))))
  I_vec <- as.vector(I)
  
  if (drop_zero) {
    keep <- !(abs(omega_mat[,1]) < 1e-12 & abs(omega_mat[,2]) < 1e-12)
    omega_mat <- omega_mat[keep, , drop = FALSE]
    I_vec <- I_vec[keep]
  }
  list(omega = omega_mat, I = I_vec)
}
plotA_imageplot <-
function(location, A, nx = 80, ny = 80, interp_lambda = NULL,
                            main_prefix = "A") {
  if (!requireNamespace("fields", quietly = TRUE)) {
    stop("Package 'fields' is required. install.packages('fields').")
  }
  x <- location[,1]; y <- location[,2]
  gx <- seq(min(x), max(x), length.out = nx)
  gy <- seq(min(y), max(y), length.out = ny)
  grid <- expand.grid(gx, gy)
  names(grid) <- c("x","y")
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mfrow = c(1, ncol(A)), mar = c(3,3,3,5))
  
  for (k in 1:ncol(A)) {
    fit <- fields::Tps(location, A[,k], lambda = interp_lambda)
    Z <- matrix(predict(fit, grid), nrow = length(gx), ncol = length(gy))
    fields::image.plot(gx, gy, Z,
                       xlab = expression(s[1]), ylab = expression(s[2]),
                       main = paste0(main_prefix, "_", k))
    points(x, y, pch = 16, cex = 0.6)
  }
  invisible(NULL)
}
ProcessRawPdg <-
function(pdg,bandwidth){
  temp = rep(0,bandwidth)
  end = length(pdg)
  temp = rev(pdg[2:(bandwidth+1)])
  pdg.final = c(temp,pdg)
  temp = rev(pdg[(end-bandwidth-1):(end-1)])
  pdg.final = c(pdg.final,temp)
  return(pdg.final)
}
R_diff_1d <-
function(K, diff_order = 2) {
  if (diff_order <= 0) return(diag(K))
  Dm <- diff(diag(K), differences = diff_order)
  t(Dm) %*% Dm
}
R_diff_2d <-
function(K1_eff, K2_eff, diff_order = 2) {
  R1 <- R_diff_1d(K1_eff, diff_order)
  R2 <- R_diff_1d(K2_eff, diff_order)
  kronecker(R1, diag(K2_eff)) + kronecker(diag(K1_eff), R2)
}
reidentify <-
function(Theta, A) {
  M <- Theta %*% t(A)
  sv <- svd(M)
  r <- ncol(Theta)
  Theta_new <- sv$u[,1:r, drop=FALSE]
  A_new <- sv$v[,1:r, drop=FALSE] %*% diag(sv$d[1:r], r, r)
  sgn <- sign(Theta_new[1,]); sgn[sgn==0] <- 1
  Theta_new <- Theta_new %*% diag(sgn, r, r)
  A_new <- A_new %*% diag(sgn, r, r)
  list(Theta = Theta_new, A = A_new)
}
ridge_eps <-
function() 1e-8
sdfhat_vec1599_to_mat <-
function(sdfhat_vec, n1 = 40, n2 = 40, fill0 = NULL) {
  sdfhat_vec <- as.numeric(sdfhat_vec)
  stopifnot(length(sdfhat_vec) == n1 * n2 - 1)
  
  # reproduce omega_mat construction exactly as in periodogram_2d()
  w1 <- omega_signed(n1)
  w2 <- omega_signed(n2)
  omega_mat <- do.call(rbind, lapply(1:n2, function(j2) cbind(w1, rep(w2[j2], n1))))
  
  # locate the (0,0) frequency index in the SAME ordering
  idx0 <- which(abs(omega_mat[,1]) < 1e-12 & abs(omega_mat[,2]) < 1e-12)
  if (length(idx0) != 1) stop("Cannot uniquely identify the zero-frequency index.")
  
  # build full length n1*n2 vector in the ORIGINAL I_vec ordering
  full <- rep(NA_real_, n1 * n2)
  if (!is.null(fill0)) full[idx0] <- as.numeric(fill0)  # fill DC if desired
  
  full[-idx0] <- sdfhat_vec
  
  # reshape back to matrix consistent with as.vector(I) in periodogram_2d()
  matrix(full, nrow = n1, ncol = n2)
}
sim_ar_sdf <-
function(m = 40, n = 256,location = NULL,ar = 1,coef = 1,grid_n = 80,
                       ...) {
 
  if (is.null(location)) location <- cbind(runif(m, 0, 1), runif(m, 0, 1))
  stopifnot(nrow(location) == m, ncol(location) == 2)
  stopifnot(ar %in% c(1,2))

  s1 <- location[,1]; s2 <- location[,2]
  
  # ---- Generate coefficients on locations
  if (ar == 1) {
    fcoef1 <- arcoeff_fun(coef = coef)
    a1_loc <- as.numeric(fcoef1(s1, s2))
    a2_loc <- NULL
    
    Xs <- lapply(1:m, function(i) sim_ar1(n = n, a1 = a1_loc[i], sigma = 1))
    
    # True spectral density for AR(1) (closed form) on Fourier grid
    omega <- 2*pi*(1:floor(n/2))/n
    f_true <- sapply(1:m, function(i) {
      den <- (1 + a1_loc[i]^2 - 2*a1_loc[i]*cos(omega))
      (1/(2*pi)) / den
    })
    f_true <- as.matrix(f_true)
    
  } else {
    fcoef2 <- arcoeff_fun_2(coef = coef)
    A12 <- fcoef2(s1, s2)      # m x 2
    a1_loc <- as.numeric(A12[,1])
    a2_loc <- as.numeric(A12[,2])
    
    Xs <- lapply(1:m, function(i) sim_ar2(n = n, a1 = a1_loc[i], a2 = a2_loc[i], sigma = 1))
    
    # True AR(2) spectrum has closed form too, but it's a bit longer.
    # To keep this sim lightweight/robust, we will evaluate performance via
    # mean MSE between fitted log-spectrum and a "plug-in" truth computed on the same omega grid:
    omega <- 2*pi*(1:floor(n/2))/n
    # f(omega) = (sigma^2/(2pi)) / |1 - a1 e^{-iw} - a2 e^{-i2w}|^2
    f_true <- sapply(1:m, function(i) {
      z1 <- exp(-1i*omega)
      denom <- Mod(1 - a1_loc[i]*z1 - a2_loc[i]*(z1^2))^2
      (1/(2*pi)) / denom
    })
    f_true <- as.matrix(Re(f_true))
  }
  
  # ---- Fit
  fit <- fit_sdf(Xs = Xs, location = location, ...)
  fit0 = fit_sdf2(Xs = Xs, location = location, ...)
  # ---- Error metric
  fhat <- fit$sdfhat
  # both are n_omega x m
  if (all(dim(fhat) == dim(f_true))) {
    mean_ise <- mean((log(fhat) - log(f_true))^2)
  } else {
    mean_ise <- NA_real_
  }
  
  # ---- Coefficients on uniform grid covering location (return matrix or list of matrices)
  x <- location[,1]; y <- location[,2]
  gx <- seq(min(x), max(x), length.out = grid_n)
  gy <- seq(min(y), max(y), length.out = grid_n)
  grid <- expand.grid(gx, gy)
  names(grid) <- c("x","y")
  
  if (ar == 1) {
    a_grid <- as.numeric(arcoeff_fun(coef = coef)(grid[,1], grid[,2]))
    arcoeff_grid <- matrix(a_grid, nrow = grid_n, ncol = grid_n)
  } else {
    A12g <- arcoeff_fun_2(coef = coef)(grid[,1], grid[,2]) # (grid_n^2) x 2
    arcoeff_grid <- list(
      a1 = matrix(as.numeric(A12g[,1]), nrow = grid_n, ncol = grid_n),
      a2 = matrix(as.numeric(A12g[,2]), nrow = grid_n, ncol = grid_n)
    )
  }
  
  list(
    mean_ise = mean_ise,
    lambda1 = fit$lambda1,
    lambda2 = fit$lambda2,
    arcoeff_grid = arcoeff_grid,
    location = location,
    fit = fit,
    ftrue = f_true,
    fhat = fhat,
    Xs = Xs,
    ar = ar,
    a1_loc = a1_loc,
    a2_loc = a2_loc,
    f_true = f_true,
    B = fit$B,
    I = fit$I,
    fhat0 =  fit0$sdfhat,
    mean_ise0 = mean((log(fit0$sdfhat) - log(f_true))^2)
  )
}
sim_ar1 <-
function(n, a1, sigma = 1) {
  e <- rnorm(n, 0, sigma)
  x <- numeric(n)
  x[1] <- e[1] / sqrt(1 - a1^2 + 1e-8)
  for (t in 2:n) x[t] <- a1*x[t-1] + e[t]
  x
}
sim_ar2 <-
function(n, a1, a2, sigma = 1) {
  e <- rnorm(n, 0, sigma)
  x <- numeric(n)
  # initialize with noise (safe and simple)
  x[1] <- e[1]
  x[2] <- a1*x[1] + e[2]
  for (t in 3:n) x[t] <- a1*x[t-1] + a2*x[t-2] + e[t]
  x
}
sim2_2DAR <-
function(
    grid_nr = 400, grid_nc = 1600,      # overall grid size (rows, cols)
    sample_nr = 40, sample_nc = 40,     # each replicate size
    n_block_col = 4,                    # 4 blocks left->right
    ar_list = NULL,                     # length n_block_col; each list(ax, ay, axy, sd)
    contam_prob = 0.10,                 # contamination prob
    innov = c("gaussian","t"),
    t_df = 5,
    # NEW: physical domain extents
    y_range = c(0, 1),                  # rows -> [0,1]
    x_range = c(0, 4),                  # cols -> [0,4]
    do_fit = TRUE,
    fit_args = list(),
    verbose = TRUE
) {

  innov <- match.arg(innov)
  
  stopifnot(grid_nr %% sample_nr == 0, grid_nc %% sample_nc == 0)
  n_regions_r <- grid_nr / sample_nr
  n_regions_c <- grid_nc / sample_nc
  m <- n_regions_r * n_regions_c
  
  stopifnot(grid_nc %% n_block_col == 0)
  block_nc <- grid_nc / n_block_col
  stopifnot(block_nc %% sample_nc == 0)
  samples_per_block_col <- block_nc / sample_nc
  
  # default AR params
  if (is.null(ar_list)) {
    if (n_block_col != 4) stop("Please provide ar_list when n_block_col != 4.")
    ar_list <- list(
      list(ax = 0.6, ay = 0.10, axy = 0.05, sd = 1.00),
      list(ax = 0.55, ay = 0.20, axy = 0.06, sd = 1.00),
      list(ax = 0.45, ay = 0.30, axy = 0.07, sd = 1.00),
      list(ax = 0.3, ay = 0.40, axy = 0.08, sd = 1.00)
    )
  }
  stopifnot(length(ar_list) == n_block_col)
  
  # map sample (rr,cc) -> true block id (left->right)
  get_block_id <- function(rr, cc) {
    (cc - 1) %/% samples_per_block_col + 1
  }
  
  # adjacent-block contamination helper (fixed)
  pick_adj_block <- function(b0, B) {
    b0 <- as.integer(b0); B <- as.integer(B)
    if (B <= 1L) return(b0)
    if (b0 <= 1L) return(2L)
    if (b0 >= B)  return(B - 1L)
    sample(c(b0 - 1L, b0 + 1L), size = 1L)
  }
  
  # 2D-AR simulator for one sample
  sim_one_2dar <- function(nr, nc, ax, ay, axy, sd) {
    eps <- if (innov == "gaussian") {
      matrix(rnorm(nr * nc, sd = sd), nr, nc)
    } else {
      sc <- sd / sqrt(t_df / (t_df - 2))
      matrix(rt(nr * nc, df = t_df) * sc, nr, nc)
    }
    X <- matrix(0, nr, nc)
    for (i in 1:nr) {
      for (j in 1:nc) {
        x1 <- if (i > 1) ax  * X[i-1, j] else 0
        x2 <- if (j > 1) ay  * X[i, j-1] else 0
        x3 <- if (i > 1 && j > 1) axy * X[i-1, j-1] else 0
        X[i, j] <- x1 + x2 + x3 + eps[i, j]
      }
    }
    X
  }
  
  # ---- NEW: map grid indices to physical coordinates ----
  # rows -> y in [y_range[1], y_range[2]]
  # cols -> x in [x_range[1], x_range[2]]
  dx <- (x_range[2] - x_range[1]) / (grid_nc - 1)
  dy <- (y_range[2] - y_range[1]) / (grid_nr - 1)
  
  Xs <- vector("list", m)
  location <- matrix(NA_real_, m, 2)  # (x, y)
  true_block <- integer(m)
  used_block <- integer(m)
  contaminated <- logical(m)
  
  idx <- 0L
  for (rr in 1:n_regions_r) {
    for (cc in 1:n_regions_c) {
      idx <- idx + 1L
      b0 <- get_block_id(rr, cc)
      b_use <- b0
      contaminated[idx] <- FALSE
      
      if (runif(1) < contam_prob && n_block_col > 1) {
        b_use <- pick_adj_block(b0, n_block_col)
        contaminated[idx] <- TRUE
      }
      
      par <- ar_list[[b_use]]
      Xs[[idx]] <- sim_one_2dar(sample_nr, sample_nc, par$ax, par$ay, par$axy, par$sd)
      
      true_block[idx] <- b0
      used_block[idx] <- b_use
      
      # center location of this sample in physical coords:
      # columns correspond to x, rows correspond to y
      x0 <- x_range[1] + ((cc - 1) * sample_nc) * dx
      x1 <- x_range[1] + ((cc * sample_nc - 1)) * dx
      y0 <- y_range[1] + ((rr - 1) * sample_nr) * dy
      y1 <- y_range[1] + ((rr * sample_nr - 1)) * dy
      location[idx, ] <- c((x0 + x1)/2, (y0 + y1)/2)
    }
  }
  
  if (verbose) {
    cat(sprintf("[sim2_2DAR] domain: y in [%.2f,%.2f], x in [%.2f,%.2f]\n",
                y_range[1], y_range[2], x_range[1], x_range[2]))
    cat(sprintf("[sim2_2DAR] grid=%dx%d, sample=%dx%d, m=%d\n",
                grid_nr, grid_nc, sample_nr, sample_nc, m))
    cat(sprintf("[sim2_2DAR] blocks=%d (left->right), samples/block=%d\n",
                n_block_col, n_regions_r * samples_per_block_col))
    cat(sprintf("[sim2_2DAR] contaminated: %d/%d (%.1f%%)\n",
                sum(contaminated), m, 100*mean(contaminated)))
    cat(sprintf("[sim2_2DAR] innov=%s\n", innov))
  }
  
  out <- list(
    Xs = Xs,
    location = location,
    true_block = true_block,
    used_block = used_block,
    contaminated = contaminated,
    ar_list = ar_list,
    grid_nr = grid_nr, grid_nc = grid_nc,
    sample_nr = sample_nr, sample_nc = sample_nc,
    x_range = x_range, y_range = y_range
  )
  
  if (do_fit) {
    if (verbose) cat("[sim2_2DAR] fitting SRCSE via fit_sdf()...\n")
    out$fit <- do.call(fit_sdf, c(list(Xs = Xs, location = location), fit_args))

    out$fit0 = fit_sdf2(Xs = Xs, location = location, r = fit_args$r, K1 = fit_args$K1, K2 = fit_args$K2,
                                                                     cv = fit_args$cv, lambda1 = fit_args$lambda1, lambda2=0,
                                                                     n_iter=30, dia = fit_args$dia,gamma=fit_args$gamma,plotA=F)
    }
  
  out
}
sim2_rect4 <-
function(
    grid_nr = 400, grid_nc = 1600,      # overall grid
    sample_nr = 40, sample_nc = 40,     # each replicate size (40x40)
    n_block_col = 4,                    # 4 blocks left->right, each 400x400
    matern_list = NULL,                 # length 4, left->right
    contam_prob = 0.10,
    seed = 1,
    do_fit = TRUE,
    fit_args = list(),
    verbose = TRUE
) {
  if (!requireNamespace("fields", quietly = TRUE)) {
    stop("Package 'fields' is required. install.packages('fields')")
  }
  set.seed(seed)
  
  stopifnot(grid_nr %% sample_nr == 0, grid_nc %% sample_nc == 0)
  n_regions_r <- grid_nr / sample_nr   # 10
  n_regions_c <- grid_nc / sample_nc   # 40
  m <- n_regions_r * n_regions_c       # 400
  
  # block width = 400, so block columns = 1600/400 = 4
  stopifnot(grid_nc %% 400 == 0, grid_nc / 400 == n_block_col)
  samples_per_block_col <- 400 / sample_nc  # 10
  stopifnot(samples_per_block_col == 10)
  
  # default matern params: ordered change from left->right
  if (is.null(matern_list)) {
    matern_list <- list(
      list(range = 0.05, nu = 0.5, sigma2 = 1.0, nugget = 0.05),  # block 1
      list(range = 0.07, nu = 1.0, sigma2 = 1.0, nugget = 0.05),  # block 2
      list(range = 0.10, nu = 1.5, sigma2 = 1.0, nugget = 0.05),  # block 3
      list(range = 0.14, nu = 2.0, sigma2 = 1.0, nugget = 0.05)   # block 4
    )
  }
  stopifnot(length(matern_list) == n_block_col)
  
  # region (rr,cc) -> block id (1..4), left->right
  get_block_id <- function(rr, cc) {
    bc <- (cc - 1) %/% samples_per_block_col + 1  # 1..4
    bc
  }
  
  # local coords within one sample (40x40) in physical domain [0,1]^2
  dx <- 1 / (grid_nc - 1)
  dy <- 1 / (grid_nr - 1)
  x_local <- (0:(sample_nc - 1)) * dx
  y_local <- (0:(sample_nr - 1)) * dy
  coords_local <- as.matrix(expand.grid(x = x_local, y = y_local))  # p x 2
  D <- fields::rdist(coords_local)
  p <- nrow(coords_local)
  
  # precompute Cholesky per block
  chol_list <- vector("list", n_block_col)
  for (b in 1:n_block_col) {
    par <- matern_list[[b]]
    C <- par$sigma2 * fields::Matern(D, range = par$range, smoothness = par$nu)
    diag(C) <- diag(C) + par$nugget
    C <- (C + t(C)) / 2
    chol_list[[b]] <- chol(C)
  }
  
  # simulate all samples
  Xs <- vector("list", m)
  location <- matrix(NA_real_, m, 2)
  true_block <- integer(m)
  used_block <- integer(m)
  contaminated <- logical(m)
  
  idx <- 0L
  for (rr in 1:n_regions_r) {
    for (cc in 1:n_regions_c) {
      idx <- idx + 1L
      b0 <- get_block_id(rr, cc)
      b_use <- b0
      
      # contamination: with prob contam_prob, use another block's params
      if (runif(1) < contam_prob && n_block_col > 1) {
        b_use <- sample(setdiff(1:n_block_col, b0), 1)
        contaminated[idx] <- TRUE
      } else contaminated[idx] <- FALSE
      
      z <- drop(t(chol_list[[b_use]]) %*% rnorm(p))
      Xs[[idx]] <- matrix(z, nrow = sample_nr, ncol = sample_nc)
      
      true_block[idx] <- b0
      used_block[idx] <- b_use
      
      # center location of this sample in [0,1]^2
      x0 <- ((cc - 1) * sample_nc) * dx
      x1 <- ((cc * sample_nc - 1)) * dx
      y0 <- ((rr - 1) * sample_nr) * dy
      y1 <- ((rr * sample_nr - 1)) * dy
      location[idx, ] <- c((x0 + x1)/2, (y0 + y1)/2)
    }
  }
  
  if (verbose) {
    cat(sprintf("[sim2_rect4] grid=%dx%d, sample=%dx%d, m=%d\n",
                grid_nr, grid_nc, sample_nr, sample_nc, m))
    cat(sprintf("[sim2_rect4] blocks=%d (left->right), samples/block=%d\n",
                n_block_col, n_regions_r * samples_per_block_col))
    cat(sprintf("[sim2_rect4] contaminated: %d/%d (%.1f%%)\n",
                sum(contaminated), m, 100*mean(contaminated)))
  }
  
  out <- list(
    Xs = Xs,
    location = location,
    true_block = true_block,
    used_block = used_block,
    contaminated = contaminated,
    matern_list = matern_list,
    grid_nr = grid_nr, grid_nc = grid_nc,
    sample_nr = sample_nr, sample_nc = sample_nc
  )
  
  if (do_fit) {
    if (verbose) cat("[sim2_rect4] fitting SRCSE via fit_sdf()...\n")
    out$fit <- do.call(fit_sdf, c(list(Xs = Xs, location = location), fit_args))
  }
  
  out
}
SmoothGVC <-
function(matrix.Rawpdg,span){
  GVCp = rep(0,length(span))
  span.min = 0
  n.trial = dim(matrix.Rawpdg)[1]
  matrix.smooth = matrix(0,dim(matrix.Rawpdg)[1],dim(matrix.Rawpdg)[2])
  y.estimate = matrix(0,length(span),dim(matrix.Rawpdg)[2])
  for(k in 1:n.trial){
    for(i in 1:(length(span))){
      bandwidth = span[i]
      first.trial.process = ProcessRawPdg(matrix.Rawpdg[k,],bandwidth)
      x.feq = seq(1,length(first.trial.process),1)
      for(j in 1:length(matrix.Rawpdg[k,])){
        window = neighborhood(x=x.feq,x.star=x.feq[j]+bandwidth,y=first.trial.process,bandwidth)
        y.estimate[i,j]=localEstimates(window$xn,window$yn,bandwidth)
      }
      GVCp[i] = CalculateGVC(matrix.Rawpdg[k,], y.estimate[i,], bandwidth)
    }
    matrix.smooth[k,] = y.estimate[which(GVCp == min(GVCp)),]
  }  
  return(matrix.smooth)
}
true_spectrum_2dar <-
function(ax, ay, axy, sd = 1, n1 = 40, n2 = 40,
                               shifted = TRUE) {
  # omega grids (signed)
  w1 <- omega_signed(n1)  # length n1
  w2 <- omega_signed(n2)  # length n2
  
  # build f(w1,w2) on grid
  F <- matrix(NA_real_, nrow = n1, ncol = n2)
  
  for (i in 1:n1) {
    e1 <- exp(-1i * w1[i])
    for (j in 1:n2) {
      e2 <- exp(-1i * w2[j])
      denom <- 1 - ax*e1 - ay*e2 - axy*(e1*e2)
      F[i, j] <- (sd^2) / ((2*pi)^2) / (Mod(denom)^2)
    }
  }
  
  # If you want FFT-native ordering (unshifted), undo the shift.
  # Here F is already in signed frequency ordering (centered). To match raw fft() order,
  # you can "ifftshift" it. Most of your pipeline uses fftshift, so keep shifted=TRUE.
  if (!shifted) {
    # inverse shift: apply fftshift again (since shift is its own inverse for even sizes)
    F <- fftshift_mat(F)
  }
  
  list(f = F, w1 = w1, w2 = w2)
}
update_A <-
function(B, Theta, A, I_mat, W, lambda2,
                     step_A = 1, n_inner = 1, Ucap = 30) {
  m <- nrow(A); r <- ncol(A)
  G <- B %*% Theta
  dvec <- deg_vec(W)
  
  for (inner in 1:n_inner) {
    for (i in 1:m) {
      ai <- A[i,]
      u <- as.vector(G %*% ai)
      u <- clamp(u, -Ucap, Ucap)
      s <- as.vector(I_mat[,i] * exp(-u))
      
      grad_w <- as.vector(crossprod(G, (1 - s)))
      
      nbr <- which(W[i,] > 0)
      if (lambda2 > 0 && length(nbr) > 0) {
        sum_w_ar <- colSums(A[nbr,,drop=FALSE] * W[i,nbr])
        grad_sp <- 4*lambda2*(dvec[i]*ai - sum_w_ar)
        H_sp <- (4*lambda2*dvec[i]) * diag(r)
      } else {
        grad_sp <- rep(0, r)
        H_sp <- matrix(0, r, r)
      }
      
      Gs <- G * sqrt(pmax(s, 0))
      H <- crossprod(Gs) + H_sp + ridge_eps()*diag(r)
      delta <- solve(H, grad_w + grad_sp)
      A_new <- ai - step_A * as.vector(delta)
      if (all(is.finite(A_new))) A[i,] <- A_new
    }
  }
  A
}
update_Theta <-
function(B, Theta, A, I_mat, R, lambda1,
                         step_Theta = 1e-3, backtrack = 0.5, max_bt = 20,
                         Ucap = 30) {
  U <- B %*% Theta %*% t(A)
  Uc <- clamp(U, -Ucap, Ucap)
  S <- I_mat * exp(-Uc)
  Rmat <- 1 - S
  grad <- crossprod(B, Rmat %*% A) + 2*lambda1*(R %*% Theta)
  
  obj0 <- whittle_loss(U, I_mat, Ucap = Ucap) + pen1_val(Theta, R, lambda1)
  tstep <- step_Theta
  
  for (k in 1:max_bt) {
    Theta_new <- Theta - tstep * grad
    U_new <- B %*% Theta_new %*% t(A)
    obj_new <- whittle_loss(U_new, I_mat, Ucap = Ucap) + pen1_val(Theta_new, R, lambda1)
    if (is.finite(obj_new) && obj_new <= obj0) return(Theta_new)
    tstep <- tstep * backtrack
  }
  Theta
}
VecToMatrix <-
function(inputData, ch, T){
  matrixp=matrix(0,T/2-1,ch)
  J=T/2
  for(i in 1:ch){
    matrixp[,i]= inputData[((i-1)*J+1):(i*J)]
  }
  return(matrixp)
}
whittle_loss <-
function(U, I_mat, Ucap = 30) {
  Uc <- clamp(U, -Ucap, Ucap)
  sum(Uc + I_mat * exp(-Uc))
}



block_mean_drop0 <- function(x, block_len = 6, drop_tail = TRUE, all_zero_value = NA_real_) {
  x <- as.numeric(x)
  stopifnot(block_len >= 1, is.finite(block_len))
  
  n <- length(x)
  if (drop_tail) {
    k <- floor(n / block_len)
    if (k == 0) return(numeric(0))
    x <- x[1:(k * block_len)]
  } else {
    k <- ceiling(n / block_len)
    x <- c(x, rep(NA_real_, k * block_len - n))
  }
  
  X <- matrix(x, nrow = block_len, ncol = k)
  
  # compute mean of nonzero and non-NA entries within each block
  out <- apply(X, 2, function(v) {
    v <- v[is.finite(v) & v != 0]
    if (length(v) == 0) all_zero_value else mean(v)
  })
  
  as.numeric(out)
}


block_mean_drop0_list <- function(Xs, block_len = 6, drop_tail = TRUE, all_zero_value = 0) {
  stopifnot(is.list(Xs))
  lapply(Xs, block_mean_drop0,
         block_len = block_len,
         drop_tail = drop_tail,
         all_zero_value = all_zero_value)
}

scale_location <- function(location) {
  rangex = max(location[,1]) - min(location[,1])
  rangey = max(location[,2]) - min(location[,2])
  if(rangex >= rangey){
    location[,2] = (location[,2] - min(location[,2]))/rangey
    location[,1] = (location[,1] - min(location[,1]))/rangey
  }
  if(rangex < rangey){
    location[,2] = (location[,2] - min(location[,2]))/rangex
    location[,1] = (location[,1] - min(location[,1]))/rangex
  }
  return(location)
}