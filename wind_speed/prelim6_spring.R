library(tidyverse)
library(lubridate)
library(viridis)
library(fields)
library(Hmisc)
library(orthopolynom)
library(geoR)
library(sp)
library(TruncatedNormal)
library(emulator)
library(MASS)
rm(list = ls())


dat = read.csv("wind_speed_hourly.csv")



dat$date_time = ymd_hms(dat$date_time)
dat$use_time = as.POSIXct(dat$date_time)
dat$windd_bin = cut(dat$windd,breaks = c(0,45,90,135,180,225,270,315,360))
dat$time_numeric = (as.numeric(dat$use_time) - min(as.numeric(dat$use_time)))/(24*3600)
dat$log_winds = log(dat$winds + 1)

month_use = 3:5

dat_use = dat[which(complete.cases(dat[,c("windd","airt_max","rh_max","solar")]) & 
                      month(dat$date_time) %in% month_use),]

n = nrow(dat_use)

dd_angle = acos(cos(rdist(dat_use$windd * 2 * pi / 360)))
dd_season = acos(cos(rdist(dat_use$time_numeric * 2 * pi)))
dd = rdist(dat_use$time_numeric)

poly_terms = 6
n_poly = poly_terms^2


gegenbauer_list = gegenbauer.polynomials(poly_terms - 1,0, normalized=TRUE)


###### empirical variogram fitting

dat_geo1 = as.geodata(cbind(0,jitter(dat_use$windd)*2*pi/360,dat_use$log_winds),
                      coords.col=1:2, data.col=3)
v1 = variog(dat_geo1,estimator.type="classical",max.dist = 2,
            breaks = seq(0,2,length = 100))

#plot(v1)

dat_geo2 = as.geodata(cbind(0,jitter(dat_use$time_numeric),dat_use$log_winds),
                      coords.col=1:2, data.col=3)
v2 = variog(dat_geo2,estimator.type="classical",max.dist = 7,
            breaks = seq(0,7,length = 200))

#plot(v2)


cov_empir = outer(max(v1$v) - v1$v, max(v2$v) - v2$v)

dd_angle_empir = simplify2array(lapply(gegenbauer_list,predict,cos(v1$uvec)))
dd_season_empir = simplify2array(lapply(gegenbauer_list,predict,cos(acos(cos(v2$uvec * 2 * pi)))))

gegenbauer_prod_empir = array(0,c(length(v1$v),length(v2$v),n_poly))

count = 1

for(i in 1:poly_terms){
  for(j in 1:poly_terms){   #### seasonal terms first
    
    gegenbauer_prod_empir[,,count] = outer(dd_angle_empir[,i],dd_season_empir[,j] * exp(-v2$uvec / 24))
    count = count +1
    
  }
}

X_gegen = matrix(0, length(v1$v)*length(v2$v), n_poly)

Y_cov = c(cov_empir)   #### angle vectors

for(i in 1:length(v2$v)){
  X_gegen[((i-1)*length(v1$v) + 1):(i*length(v1$v) ),] = gegenbauer_prod_empir[,i,]
}


lm_semivario = lm(Y_cov ~ X_gegen)
b_now = coef(lm_semivario)[-1]
covar_empir_pos = vcov(lm_semivario)[-1,-1]

log_b_now_sim = log(rtmvnorm(1e4, b_now, covar_empir_pos, rep(0,n_poly),rep(Inf,n_poly)))



############### data setup       

gegenbauer_dist_angle = simplify2array(lapply(gegenbauer_list,predict,cos(dd_angle)))
gegenbauer_dist_season = simplify2array(lapply(gegenbauer_list,predict,cos(dd_season)))

gegenbauer_prod = array(0,c(n_poly,n,n))

count = 1

for(i in 1:poly_terms){
  for(j in 1:poly_terms){  ####### seasonal terms -- summed first
    
    gegenbauer_prod[count,,] = gegenbauer_dist_angle[,,i] * gegenbauer_dist_season[,,j]
    count = count +1
    
  }
}

gegenbauer_distm = lapply(seq(dim(gegenbauer_prod)[1]), function(x) gegenbauer_prod[x , , ])

Y = dat_use$log_winds
X = cbind(1,scale(cbind(dat_use$airt_max,dat_use$rh_max,dat_use$solar)))
p = ncol(X)
lm_xb = lm(Y ~ X[,-1])

############################################################# 
############################################################# 
################## model 
############################################################# 
############################################################# 

my_dmvnorm = function(Y,sigma2,log_det,prec){
  n_local = length(Y) 
  -n_local/2 * log(2 * pi) - 
    1/2 * (log_det + n_local * log(sigma2)) - 
    quad.form(prec,Y) / (2 * sigma2)
}

my_dmvlnorm = function(logY,logmean,sigma2,log_det,prec){
  n_local = length(logY) 
  -n_local/2 * log(2 * pi) - 
    1/2 * (log_det + n_local * log(sigma2)) - 
    quad.form(prec,logY - logmean) / (2 * sigma2) - sum(logY)
}


likelihood = function(Y,log_det,prec){
  n_local = length(Y) 
  -n_local/2 * log(2 * pi) - 
    1/2 * (log_det) - 
    quad.form(prec,Y) / 2
}

reps = 3e4
burn = 3e4
tune = 100

b = array(0,c(reps ,n_poly)); b_now = exp(apply(log_b_now_sim,2,mean))
phi = rep(1/24,reps ); phi_now = 1/24
tau2 = rep(0.01,reps); tau2_now = 0.01
#sig2 = rep(0.5,reps ) ; sig2_now = 0.5
beta_reg = matrix(0,reps,p) ; beta_reg_now = coef(lm_xb)
like_save = numeric(reps)

pars_now = c(b_now,phi_now,tau2_now)
pars_save = matrix(0,reps + burn,n_poly + 2)
pars_save[1,] = pars_now

xb = X %*% beta_reg_now

Sig_now_terms = lapply(1:n_poly,function(qq){
  pars_now[qq] * gegenbauer_distm[[qq]]
})

Sig_now = exp( -pars_now[n_poly+1] * dd) * Reduce('+', Sig_now_terms) + pars_now[n_poly+2] * diag(n)
chol_sig_now = chol(Sig_now)
log_sig_det = 2 * sum(log(diag(chol_sig_now)))
inv_sig = chol2inv(chol_sig_now)

XsigX = quad.form(inv_sig,X)

like_now = likelihood(Y - xb,log_sig_det,inv_sig)


scale_fac1 = 2.38^2 / (n_poly+2)
cand_keep =  0.0001 * diag(n_poly+2)
cand_keep[1:n_poly,1:n_poly] = cov(log_b_now_sim)/10
cand_keep[n_poly+1,n_poly+1] = 0.01
cand_keep[n_poly+2,n_poly+2] = 0.1

cand_var1 = scale_fac1 *cand_keep


count1 = 0 ; chol1 = chol(cand_var1); inv1 = chol2inv(chol1); log_det1 = 2 * sum(log(diag(chol1)))

st = proc.time()

for(i in 2:(reps + burn)){
  
  ######################### update Beta Regression
  
  v_bet = solve(XsigX + diag(p)/100)
  m_bet = t(X) %*% c(inv_sig %*% Y)
  
  beta_reg_now = mvrnorm(1,v_bet %*% m_bet, v_bet) 
  xb = X %*% beta_reg_now
  
  like_now = likelihood(Y - xb,log_sig_det,inv_sig)
  
  
  ######################### update all pars
  
  log_now = log(pars_now)
  log_cand = log_now + t(chol1) %*% rnorm(n_poly + 2)
  cand = exp(log_cand)
  
  Sig_cand_terms = lapply(1:n_poly,function(qq){
    cand[qq] * gegenbauer_distm[[qq]]
  })
  
  Sig_cand = exp( -cand[n_poly + 1] * dd) * Reduce('+', Sig_cand_terms) + 
    cand[n_poly+2] * diag(n)
  
  chol_sig_cand= chol(Sig_cand)
  log_sig_det_cand = 2 * sum(log(diag(chol_sig_cand)))
  inv_sig_cand = chol2inv(chol_sig_cand)
  
  like_cand = likelihood(Y - xb,log_sig_det_cand,inv_sig_cand)
  
  prior_dif = dgamma(cand[n_poly + 1],2,72, log = TRUE ) +  
    sum(dcauchy(cand[-(n_poly + 1)],0,1, log = TRUE )) - 
    dgamma(pars_now[n_poly + 1],2,72, log = TRUE ) - 
    sum(dcauchy(pars_now[-(n_poly + 1)],0,1, log = TRUE ))
  
  MH_dif = my_dmvlnorm(log_now,log_cand,1,log_det1,inv1) -
    my_dmvlnorm(log_cand,log_now,1,log_det1,inv1)
  
  if(like_cand - like_now + prior_dif + MH_dif > log(runif(1))){
    
    pars_now = cand 
    
    Sig_now = Sig_cand
    Sig_now_terms = Sig_cand_terms
    log_sig_det= log_sig_det_cand
    inv_sig = inv_sig_cand
    
    like_now = like_cand
    
    XsigX = quad.form(inv_sig,X)
    
    count1 = count1 + 1
  }
  
  
  
  pars_save[i,] = pars_now
  
  
  
  if(i > burn){
    
    b[i - burn,] = pars_now[1:n_poly]
    phi[i - burn] = pars_now[n_poly+1]
    tau2[i - burn] = pars_now[n_poly+2]
    beta_reg[i - burn,] = beta_reg_now
    like_save[i - burn] = like_now 
    
  }  
  
  
  if(i %% tune == 0){
    
    if(i < burn){
      
      acc1 = count1 / tune; count1 = 0
      
      if(acc1 > 0.05 & acc1 < 0.9){

        if(acc1 < 0.1){
          scale_fac1 = scale_fac1 / 1.25
        }

        if(acc1 > 0.6){
          scale_fac1 = scale_fac1 * 1.1

        }

      } else if(acc1 < 0.05){

        scale_fac1 = scale_fac1 / 1.5

      } else{

        scale_fac1 = scale_fac1 * 1.25

      }
      
      
      if(i > burn/10){

        cand_var1 = scale_fac1 * cov(log(unique(pars_save[1:i,])))
        cand_var1 = cand_var1 + diag(diag(cand_var1))/10

      } else{

        cand_var1 = scale_fac1 * cand_keep

      }
      
      
    } else{
      
      cand_var1 = scale_fac1 * cov( unique(log(pars_save[(burn/2):i,]) ))
      cand_var1 = cand_var1 + diag(diag(cand_var1))/10
    }
    
    chol1 = chol(cand_var1); inv1 = chol2inv(chol1); log_det1 = 2 * sum(log(diag(chol1)))
    
    # time_its <- (proc.time() - st)[3] / (i)
    # time_used <- round((proc.time() - st)[3]/(60),digits=4)
    # time_left <- round(time_its * (reps + burn - i )/(60),digits=4)
    # cat("\r", i, " of ", reps + burn,"||| Time left: ",floor(time_left/60),
    #     " hours",time_left%%60," minutes") 
    # flush.console()
    # cat("\n,   acceptance rate:", acc1,"scale=",scale_fac1)
    
  }
  
  chol1 = chol(cand_var1); inv1 = chol2inv(chol1); log_det1 = 2 * sum(log(diag(chol1)))
  
  time_its <- (proc.time() - st)[3] / (i)
  time_used <- round((proc.time() - st)[3]/(60),digits=4)
  time_left <- round(time_its * (reps + burn - i )/(60),digits=4)
  cat("\r", i, " of ", reps + burn,"||| Time left: ",floor(time_left/60),
      " hours",time_left%%60," minutes")
  flush.console()
  cat("\n,   count1:", c(count1), "log-like = ",like_now,"scale_fac= ",scale_fac1 )
  
  
}

# saveRDS(pars_save,"par_start6.rds")

save.image("mod6_spring..RData")

