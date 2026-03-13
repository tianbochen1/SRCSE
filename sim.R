library(fields)
source('fn.R')
library(clusteval)
library(TSclust)

par(mfcol = c(2,4), mar=c(3.5,3.5,2,1.5), mgp=c(3,0.5,0)) 
cases = c(1,4,6,7)

for(i in 1:4){
  out1 <- sim_ar_sdf(
    m=220, n=128, ar=2, coef=cases[i],
    cv=FALSE, lambda1=0.1, lambda2=0,
    dia=0.3, gamma=0.1,
    plotA=F
  )
  s = out1$arcoeff_grid
  
  image(s$a1, xaxt='n', yaxt='n', xlab='', ylab='', 
        col = tim.colors(),
        main = bquote(alpha[1](s[1],s[2]): bold(C)*bold(.(i))  ))
  axis(side = 1, tck = -0.02)
  axis(side = 2, tck = -0.02)
  title(xlab=expression(s[1]), ylab=expression(s[2]), 
        line=2, cex.lab=1.2)
  
  image(s$a2, xaxt='n', yaxt='n', xlab='', ylab='', 
        col = tim.colors(),
        main = bquote(alpha[2](s[1],s[2]): bold(C)*bold(.(i))))
  axis(side = 1, tck = -0.02)
  axis(side = 2, tck = -0.02)
  title(xlab=expression(s[1]), ylab=expression(s[2]), 
        line=2, cex.lab=1.2)
}


########################################sim1#############################
source('fn.R')
set.seed(1)
res = matrix(0,4,5)
nrep =  200
mn = matrix(c(200,200,400,400,128,256,128,256),4,2)
KK = c(10,16,10,16)
RR = c(3,5,3,5)
for(k in 1:4){
  c1 = rep(0,5)
  for(i in 1:nrep){
    out1 <- sim_ar_sdf(
      m=mn[k,1], n=mn[k,2], ar=2, coef=1,
      cv=FALSE, lambda1=1, lambda2=0.1,
      dia=0.3, gamma=0.1,K=KK[k],r=RR[k],
      plotA=F
    )
    I = out1$I
    B = out1$B
    ftrue = out1$ftrue
    
    c1[1] = c1[1] + out1$mean_ise/nrep
    
    c1[2] = c1[2] + out1$mean_ise0/nrep
    
    smo1 =  B %*% solve(t(B)%*%B) %*% t(B) %*% I
    for(j in 1:length(smo1)){
      if(smo1[j]<=0){smo1[j]=0.000001}
    }
    c1[3] =  c1[3]+  mean((log(smo1) - log(ftrue))^2)/nrep
    
    smo2 = I
    for(j in 1:dim(I)[2]){
      smo2[,j] = smooth.spline(I[,j],cv=T)$y
    }
    for(j in 1:length(smo2)){
      if(smo2[j]<=0){smo2[j]=0.000001}
    }
    c1[4] =c1[4] + mean((log(smo2) - log(ftrue))^2)/nrep
    
    span = seq(from = 3, to = 15, by = 1)
    sgcv = t(SmoothGVC(t(I),span))
    c1[5] =c1[5] + mean((log(sgcv) - log(ftrue))^2)/nrep
    
    print(i)
  }
  res[k,] = c1
}
res

############clu
set.panel(2,5)
par(mar=c(3.5,3.5,1.5,1.5),mgp=c(2,0.5,0))
d = dist(t(srcse))
h = hclust(d, method = 'ward.D')
res = cutree(h,4)
plot(out1$location,col=res,pch=20,cex=1.5
     ,xlab='',ylab='',xaxt='n',yaxt='n', main = 'SRCSE')
title(xlab=expression(s[1]), ylab=expression(s[2]), 
      line=2, cex.lab=1.2)
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)

d = dist(t(se))
h = hclust(d, method = 'ward.D')
res = cutree(h,4)
plot(out1$location,col=res,pch=20,cex=1.5
     ,xlab='',ylab='',xaxt='n',yaxt='n', main='SE')
title(xlab=expression(s[1]), ylab=expression(s[2]), 
      line=2, cex.lab=1.2)
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)

d = dist(t(smo1))
h = hclust(d, method = 'ward.D')
res = cutree(h,4)
plot(out1$location,col=res,pch=20,cex=1.5
     ,xlab='',ylab='',xaxt='n',yaxt='n', main='SMO1')
title(xlab=expression(s[1]), ylab=expression(s[2]), 
      line=2, cex.lab=1.2)
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)

d = dist(t(smo2))
h = hclust(d, method = 'ward.D')
res = cutree(h,4)
plot(out1$location,col=res,pch=20,cex=1.5
     ,xlab='',ylab='',xaxt='n',yaxt='n', main='SMO2')
title(xlab=expression(s[1]), ylab=expression(s[2]), 
      line=2, cex.lab=1.2)
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)

d = dist(t(sgcv))
h = hclust(d, method = 'ward.D')
res = cutree(h,4)
plot(out1$location,col=res,pch=20,cex=1.5
     ,xlab='',ylab='',xaxt='n',yaxt='n', main='GGCV')
title(xlab=expression(s[1]), ylab=expression(s[2]), 
      line=2, cex.lab=1.2)
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)


out1 <- sim_ar_sdf(
  m=400, n= 256, ar=2, coef=4,
  cv=FALSE, lambda1=1, lambda2=0.2,
  dia=0.3, gamma=0.1,K=10,r=5,
  plotA=F
)

I = out1$I
B = out1$B
ftrue = out1$ftrue
srcse = out1$fhat
se = out1$fhat0
smo1 =  B %*% solve(t(B)%*%B) %*% t(B) %*% I
smo2 = I
for(j in 1:dim(I)[2]){
  smo2[,j] = smooth.spline(I[,j],cv=T)$y
}
span = seq(from = 5, to = 10, by = 1)
sgcv = t(SmoothGVC(t(I),span))



d = dist(t(srcse))
h = hclust(d, method = 'ward.D')
res = cutree(h,4)
res[which(res==2)] = 5
res[which(res==3)] = 2
res[which(res==5)] = 3
plot(out1$location,col=res,pch=20,cex=1.5
     ,xlab='',ylab='',xaxt='n',yaxt='n', main='SRCSE')
title(xlab=expression(s[1]), ylab=expression(s[2]), 
      line=2, cex.lab=1.2)
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)

d = dist(t(se))
h = hclust(d, method = 'ward.D')
res = cutree(h,4)
res[which(res==2)] = 5
res[which(res==3)] = 2
res[which(res==5)] = 3

plot(out1$location,col=res,pch=20,cex=1.5
     ,xlab='',ylab='',xaxt='n',yaxt='n', main='SE')
title(xlab=expression(s[1]), ylab=expression(s[2]), 
      line=2, cex.lab=1.2)
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)

d = dist(t(smo1))
h = hclust(d, method = 'ward.D')
res = cutree(h,4)
plot(out1$location,col=res,pch=20,cex=1.5
     ,xlab='',ylab='',xaxt='n',yaxt='n', main='SMO1')
title(xlab=expression(s[1]), ylab=expression(s[2]), 
      line=2, cex.lab=1.2)
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)

d = dist(t(smo2))
h = hclust(d, method = 'ward.D')
res = cutree(h,4)
plot(out1$location,col=res,pch=20,cex=1.5
     ,xlab='',ylab='',xaxt='n',yaxt='n', main='SMO2')
title(xlab=expression(s[1]), ylab=expression(s[2]), 
      line=2, cex.lab=1.2)
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)

d = dist(t(sgcv))
h = hclust(d, method = 'ward.D')
res = cutree(h,4)
plot(out1$location,col=res,pch=20,cex=1.5
     ,xlab='',ylab='',xaxt='n',yaxt='n', main='GGCV')
title(xlab=expression(s[1]), ylab=expression(s[2]), 
      line=2, cex.lab=1.2)
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)

####################sim2
source('fn.r')
set.seed(1)
nrep = 200
sim = rep(0,5)
ari = rep(0,5)
for(i in 1:nrep){
  res2 <- sim2_2DAR(
    grid_nr= 800, grid_nc=3200,
    sample_nr= 80, sample_nc=80,
    contam_prob=0.1,
    innov="gaussian",
    do_fit=TRUE,
    fit_args=list(
      r=40, K1=12, K2=12,
      cv=FALSE,
      lambda1 = 2,
      lambda2 = 0.6,
      n_iter=30,
      dia=0.4, gamma=0.2,
      plotA=FALSE
    )
  )
  
  I = log(res2$fit$I)
  B = res2$fit$B
  
  I_mt = matrix(0,dim(I)[1]+1,dim(I)[2])
  I_ke = matrix(0,dim(I)[1]+1,dim(I)[2])
  for(j in 1:dim(I)[2]){
    II =  sdfhat_vec1599_to_mat(res2$fit$I[,j],n1=80,n2=80)
    II[1] = 0
    II = fftshift_mat(II)
    
    IIMT = mt_smooth_periodogram_2d(II)
    IIMT = as.vector(IIMT)
    I_mt[,j] = IIMT
    
    IIKE = image.smooth(II)$z
    IIKE = as.vector(IIKE)
    I_ke[,j] = IIKE
  }
  ISP = B %*% solve(t(B)%*%B) %*% t(B) %*% I
  ISP = exp(ISP)
  SRCSE = res2$fit$sdfhat
  SE = res2$fit0$sdfhat
  
  tru = res2$true_block
  
  d = dist(t(SRCSE))
  h = hclust(d, method = 'ward.D')
  res = cutree(h,4)
  ari[1] = ari[1] + cluster_similarity(res,tru)/nrep #ari
  sim[1] = sim[1] + cluster.evaluation(res,tru)/nrep #sim
  
  d = dist(t(SE))
  h = hclust(d, method = 'ward.D')
  res = cutree(h,4)
  ari[2] = ari[2] + cluster_similarity(res,tru)/nrep #ari
  sim[2] = sim[2] + cluster.evaluation(res,tru)/nrep #sim
  
  
  d = dist(t(ISP))
  h = hclust(d, method = 'ward.D')
  res = cutree(h,4)
  ari[3] = ari[3] + cluster_similarity(res,tru)/nrep #ari
  sim[3] = sim[3] + cluster.evaluation(res,tru)/nrep #sim
  
  d = dist(t(I_mt))
  h = hclust(d, method = 'ward.D')
  res = cutree(h,4)
  ari[4] = ari[4] + cluster_similarity(res,tru)/nrep #ari
  sim[4] = sim[4] + cluster.evaluation(res,tru)/nrep #sim
  
  d = dist(t(I_ke))
  h = hclust(d, method = 'ward.D')
  res = cutree(h,4)
  ari[5] = ari[5] + cluster_similarity(res,tru)/nrep #ari
  sim[5] = sim[5] + cluster.evaluation(res,tru)/nrep #sim
  print(i)
  }
  

sdfhat = res2$fit$sdfhat
d = dist(t(sdfhat))
h = hclust(d, method = 'ward.D')
res = cutree(h,4)
plot(res2$location,col=res,pch=20,cex=1.5
     ,xlab='',ylab='',xaxt='n',yaxt='n', main='GGCV')
title(xlab=expression(s[1]), ylab=expression(s[2]), 
      line=2, cex.lab=1.2)
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)


plot(res2$location,col=res2$used_block,pch=20,cex=1.5
     ,xlab='',ylab='',xaxt='n',yaxt='n', main='GGCV')
title(xlab=expression(s[1]), ylab=expression(s[2]), 
      line=2, cex.lab=1.2)
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)

#############################plot sim2
res2 <- sim2_2DAR(
  grid_nr= 400, grid_nc=1600,
  sample_nr= 40, sample_nc=40,
  contam_prob=0.15,
  innov="gaussian",
  do_fit=TRUE,
  fit_args=list(
    r=10, K1=6, K2=6,
    cv=FALSE,
    lambda1 = 2,
    lambda2 = 0.4,
    n_iter=30,
    dia=0.4, gamma=0.2,
    plotA=FALSE
  )
)

I = log(res2$fit$I)
B = res2$fit$B

I_mt = matrix(0,dim(I)[1]+1,dim(I)[2])
I_ke = matrix(0,dim(I)[1]+1,dim(I)[2])
for(j in 1:dim(I)[2]){
  II =  sdfhat_vec1599_to_mat(res2$fit$I[,j])
  II[1] = 0
  IIMT = mt_smooth_periodogram_2d(II)
  IIMT = as.vector(IIMT)
  I_mt[,j] = IIMT
  
  IIKE = image.smooth(II)$z
  IIKE = as.vector(IIKE)
  I_ke[,j] = IIKE
}
ISP = B %*% solve(t(B)%*%B) %*% t(B) %*% I
ISP = exp(ISP)
I_mt = exp(I_mt)
I_ke = exp(I_ke)
SRCSE = res2$fit$sdfhat
SE = res2$fit0$sdfhat


set.panel(4,1)
par(mar=c(3.5,3.5,1.5,1.5),mgp=c(2,0.5,0))
d = dist(t(SRCSE))
h = hclust(d, method = 'ward.D')
res = cutree(h,4)
plot(res2$location,col=res,pch=15,cex=2
     ,xlab='',ylab='',xaxt='n',yaxt='n', main='SRCSE',xaxs="i",xlim=c(-0.01,4.01))
title(xlab=expression(s[1]), ylab=expression(s[2]), 
      line=2, cex.lab=1.2)
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)


d = dist(t(SE))
h = hclust(d, method = 'ward.D')
res = cutree(h,4)
plot(res2$location,col=res,pch=15,cex=2
     ,xlab='',ylab='',xaxt='n',yaxt='n', main='SE',xaxs="i",xlim=c(-0.01,4.01))
title(xlab=expression(s[1]), ylab=expression(s[2]), 
      line=2, cex.lab=1.2)
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)


d = dist(t(ISP))
h = hclust(d, method = 'ward.D')
res = cutree(h,4)
res[which(res==1)] = 5
res[which(res==2)] = 1
res[which(res==5)] = 2
plot(res2$location,col=res,pch=15,cex=2
     ,xlab='',ylab='',xaxt='n',yaxt='n', main='SMO1',xaxs="i",xlim=c(-0.01,4.01))
title(xlab=expression(s[1]), ylab=expression(s[2]), 
      line=2, cex.lab=1.2)
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)

d = dist(t(I_mt))
h = hclust(d, method = 'ward.D')
res = cutree(h,4)
plot(res2$location,col=res,pch=15,cex=2
     ,xlab='',ylab='',xaxt='n',yaxt='n', main='MT',xaxs="i",xlim=c(-0.01,4.01))
title(xlab=expression(s[1]), ylab=expression(s[2]), 
      line=2, cex.lab=1.2)
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)



#######################################case 
sim2_2DAR <-
  function(
    grid_nr = 400, grid_nc = 1600,      # overall grid size (rows, cols)
    sample_nr = 40, sample_nc = 40,     # each replicate size
    n_block_col = 4,                    # 4 blocks left->right
    ar_list = NULL,                     # length n_block_col; each list(ax, ay, axy, sd)
    contam_prob = 0.15,                 # contamination prob
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
        list(ax = 0.5,  ay = 0.10, axy = 0.05, sd = 1.00),
        list(ax = 0.25, ay = 0.20, axy = 0.05, sd = 1.00),
        list(ax = 0.18, ay = 0.30, axy = 0.05, sd = 1.00),
        list(ax = 0.1,  ay = 0.40, axy = 0.05, sd = 1.00)
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
set.seed(1)
res2 <- sim2_2DAR(
  grid_nr= 800, grid_nc=3200,
  sample_nr= 80, sample_nc=80,
  contam_prob=0.15,
  innov="gaussian",
  do_fit=TRUE,
  fit_args=list(
    r=40, K1=12, K2=12,
    cv=FALSE,
    lambda1 = 2,
    lambda2 = 0.6,
    n_iter=30,
    dia=0.4, gamma=0.2,
    plotA=FALSE
  )
)
ar_list <- res2$ar_list
true_spec_list <- lapply(ar_list, function(par) {
  true_spectrum_2dar(par$ax, par$ay, par$axy, sd = par$sd, n1 = 80, n2 = 80, shifted =F)
})

set.panel(1,4)
par(mar=c(3.8,3.8,1.5,1),mgp=c(2,0.5,0))
omega1 = seq(from=-pi,to=pi,length.out=80)
omega2 = seq(from=-pi,to=pi,length.out=80)
image(omega1/pi/2, omega2/pi/2 , true_spec_list[[1]]$f,xaxt="n",yaxt="n",xlab='',ylab='',col=tim.colors(),
      main = 'True SDF: block 1',xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab=expression(paste(omega[1], " (×", " 2", pi, ")")),ylab=expression(paste(omega[2], " (×", " 2", pi, ")")), line=2, cex.lab=1.2)

image(omega1/pi/2, omega2/pi/2 , true_spec_list[[2]]$f,xaxt="n",yaxt="n",xlab='',ylab='',col=tim.colors(),
      main = 'True SDF: block 2',xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab=expression(paste(omega[1], " (×", " 2", pi, ")")),ylab=expression(paste(omega[2], " (×", " 2", pi, ")")), line=2, cex.lab=1.2)

contam = 7
srcse = sdfhat_vec1599_to_mat(res2$fit$sdfhat[,contam],n1=80,n2=80)
srcse[1] = 0
srcse = fftshift_mat(srcse)
image(omega1/pi/2, omega2/pi/2 , srcse, xaxt="n",yaxt="n",xlab='',ylab='',col=tim.colors(),
      main = 'SRCSE (contaminated)',xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab=expression(paste(omega[1], " (×", " 2", pi, ")")),ylab=expression(paste(omega[2], " (×", " 2", pi, ")")), line=2, cex.lab=1.2)

srcse = sdfhat_vec1599_to_mat(res2$fit0$sdfhat[,contam],n1=80,n2=80)
srcse[1] = 0
srcse = fftshift_mat(srcse)
image(omega1/pi/2, omega2/pi/2 , srcse, xaxt="n",yaxt="n",xlab='',ylab='',col=tim.colors(),
      main = 'SE (contaminated)',xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab=expression(paste(omega[1], " (×", " 2", pi, ")")),ylab=expression(paste(omega[2], " (×", " 2", pi, ")")), line=2, cex.lab=1.2)






