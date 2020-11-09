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
library(globe)
library(geosphere)
library(sf)
library(spData)
library(splines)
library(MCMCpack)
library(Matrix)

rm(list = ls())

dat_use = read.csv("cloud_use.csv")[,-1]

dat_use$y = scale(dat_use$cloud_cover)

library(ggplot2)




n_t = max(dat_use$month)
n = nrow(dat_use)

dd_sphere = distm(cbind(dat_use$lon,dat_use$lat)) / 6.371e6
dd_season = acos(cos(rdist(dat_use$month * 2 * pi/n_t)))
dd = rdist(dat_use$month)

poly_terms1 = 20
poly_terms2 = 4

n_poly = poly_terms1 * poly_terms2

gegenbauer_list1 = gegenbauer.polynomials(poly_terms1 - 1,1/2,normalized = TRUE)
gegenbauer_list2 = gegenbauer.polynomials(poly_terms2 - 1,0,normalized = TRUE)


Y = c(dat_use$y)

lm_temp = lm(y ~ is_land*scale(abs(lat)) + bs(lat,knots = c(-30,30))*is_land,data = dat_use)

lm_temp_month = lm(y ~ factor(dat_use$month)*(is_land*scale(abs(lat)) + bs(lat,knots = c(-30,30))*is_land),data = dat_use)
anova(lm_temp,lm_temp_month)

X = cbind(1,model.matrix(lm_temp)[,-1])
p = ncol(X)

X_t = lapply(1:n_t,function(qq){X[which(dat_use$month == qq),]})
X_block = as.matrix(bdiag(X_t))

dat_use$resid = resid(lm_temp_month)

###### empirical variogram fitting

cor_use = cor(matrix(dat$resid,ncol = n_t))


temp_cor = numeric(n_t)
temp_cor[1] = 1
for(i in 1:11){
  temp_cor[i+1] = mean(cor_use[row(cor_use) == (col(cor_use) - i)])
  
}

plot(0:11,temp_cor,type = "o",xlab = "Difference in Months",
     ylab = "Average Autocorrelation",pch = 16,cex = 2,lwd =2)

ggplot(dat_use,aes(x = proj_x,y=proj_y)) + 
  geom_point(aes(col = resid),size = 3) + 
  scale_color_viridis()

fit_v = vector(length = n_t,mode = "list")
vario_temp = vector(length = n_t,mode = "list")

for(i in 1:n_t){
  dat_geo1 = as.geodata(cbind(dat_use$proj_x,
                              dat_use$proj_y,
                              scale(dat_use$resid))[which(dat_use$month == i),],
                        coords.col=1:2, data.col=3)
  
  vario_temp[[i]] = variog(dat_geo1,estimator.type="classical",max.dist = pi/4,
                           breaks = seq(0.05,pi/4,length = 100))
}


vario_temp_all = vario_temp[[1]]
vario_temp_all$v = apply(sapply(1:n_t,function(x){vario_temp[[x]]$v}),1,mean)

# fit_all = variofit(vario_temp_all, c(1,.2), cov.model = "matern",
#                    nugget = .1,fix.kappa = TRUE,kappa = 1/2)

plot(vario_temp_all)

# lines.variomodel(fit_all)

spat_cov = max(vario_temp_all$v) - vario_temp_all$v

cov_empir = outer(spat_cov, c(temp_cor))

dd_sphere_empir = simplify2array(lapply(gegenbauer_list1,predict,cos(vario_temp_all$u)))
dd_season_empir = simplify2array(lapply(gegenbauer_list2,predict,cos(acos(cos(0:11 * 2 * pi/n_t)))))


# b_season = coef(lm(temp_cor ~ 0 + dd_season_empir))
# b_sphere = coef(lm( spat_cov ~ 0 + dd_sphere_empir))

gegenbauer_prod_empir = array(0,c(length(vario_temp_all$v),length(temp_cor),n_poly))

count = 1

for(i in 1:poly_terms1){
  for(j in 1:poly_terms2){   #### seasonal terms first
    
    gegenbauer_prod_empir[,,count] = outer(dd_sphere_empir[,i],dd_season_empir[,j] * exp(-(0:11)/(6)))
    count = count +1
    
  }
}


X_gegen = matrix(0, length(vario_temp_all$v)*(length(temp_cor)), n_poly)

Y_cov = c(cov_empir)   #### angle vectors

for(i in 1:length(temp_cor)){
  X_gegen[((i-1)*length(vario_temp_all$v) + 1):(i*length(vario_temp_all$v)),] = gegenbauer_prod_empir[,i,]
}

library(glmnet)
ridge_mod = glmnet(X_gegen, Y_cov, alpha = 0, lambda = 1)
b_now = predict(ridge_mod, type = "coefficients")[-1,]

lm_semivario = lm(Y_cov ~ 0 + X_gegen)
covar_empir_pos = summary(lm_semivario)$sigma^2 * solve(t(X_gegen) %*% X_gegen + diag(ncol(X_gegen)))

log_b_now_sim = log(rtmvnorm(1e4, b_now, covar_empir_pos, rep(0,n_poly),rep(Inf,n_poly)))


############### data setup       

gegenbauer_dist_angle = simplify2array(lapply(gegenbauer_list1,predict,cos(dd_sphere)))
gegenbauer_dist_season = simplify2array(lapply(gegenbauer_list2,predict,cos(dd_season)))

gegenbauer_prod = array(0,c(n_poly,n,n))

count = 1

for(i in 1:poly_terms1){
  for(j in 1:poly_terms2){  ####### seasonal terms -- summed first
    
    gegenbauer_prod[count,,] = gegenbauer_dist_angle[,,i] * gegenbauer_dist_season[,,j]
    count = count +1
    
  }
}

gegenbauer_distm = lapply(seq(dim(gegenbauer_prod)[1]), function(x) gegenbauer_prod[x , , ])


############################################################# 
############################################################# 
################## model 
############################################################# 
############################################################# 

prec_ou = function(difs, phi, sig2, n){
  
  prec = matrix(0,ncol = n,nrow = n)
  diags = numeric(n)
  
  if(length(unique(difs))==1){
    
    dif_u = unique(difs)
    diags[c(1,n)] = 1 / (1 - exp(-2 * phi *dif_u ))
    diags[2:(n-1)] = 1 / (1 - exp(-2 * phi *dif_u)) + exp(-2 * phi *dif_u) / 
      (1 - exp(-2 * phi * dif_u))
    
    diag(prec) = diags
    
    cross = - exp(-phi * dif_u) / (1 - exp(-2 * phi * dif_u))
    
    for(i in 1:(n-1)){
      prec[i,i+1] <- prec[i+1,i] <- cross
    }
    
  } else{
    
    diags[1] = 1 / (1 - exp(-2 * phi * difs[1] ))
    diags[2:(n-1)] = 1 / (1 - exp(-2 * phi * difs[1:(n-2)])) +
      exp(-2 * phi * difs[2:(n-1)]) / (1 - exp(-2 * phi * difs[2:(n-1)]))
    diags[n] = 1 / (1 - exp(-2 * phi * difs[n-1] ))
    
    diag(prec) = diags
    
    for(i in 1:(n-1)){
      prec[i+1,i] <- prec[i,i+1] <- - exp(-phi * difs[i]) / (1 - exp(-2 * phi * difs[i]))
    }
    
  }
  
  
  return(prec/sig2)
}


dmatnorm = function(X,mu,U_inv,V_inv,log_det_U,log_det_V){
  n = nrow(X)
  p = ncol(X)
  X_scale = sweep(X,2,mu,"-")
  left = V_inv %*% t(X_scale)
  right = U_inv %*% X_scale
  -n * p * log(2 * pi) / 2  - n/2 * log_det_V - p/2 * log_det_U - 1/2 * sum(diag(left %*% right))
}


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

reps = 1e4
burn = 3e4
tune = 100

b = array(0,c(reps ,n_poly)); b_now = exp(apply(log_b_now_sim,2,mean))
phi = rep(1/6,reps ); phi_now = 1/6
phi_beta = rep(1/6,reps ); phi_beta_now = 1/6
tau2 = rep(0.1,reps); tau2_now = 0.1
#sig2 = rep(0.5,reps ) ; sig2_now = 0.5
beta_reg = array(0,c(reps,p,n_t)) ; beta_reg_now = matrix(0,p,n_t)
beta_tot = matrix(0,reps,p) ; beta_tot_now = rep(0,p)
V_reg = array(0,c(reps,p,p)) ; V_reg_now = matrix(0,p,p)
like_save = numeric(reps)

for(i in 1:n_t){
  lm_temp = lm(Y[dat_use$month == i] ~ 0 + X[dat_use$month == i,])
  beta_reg_now[,i] = coef(lm_temp)
  V_reg_now = vcov(lm_temp)
}


log_det_V = determinant(V_reg_now)$modulus[1]
V_reg_inv_now = solve(V_reg_now)


pars_now = c(b_now,phi_now,tau2_now)
pars_save = matrix(0,reps + burn,n_poly + 2)
pars_save[1,] = pars_now

xb = unlist(lapply(1:n_t,function(qq){ X_t[[qq]] %*% beta_reg_now[,qq]}))

Sig_now_terms = lapply(1:n_poly,function(qq){
  pars_now[qq] * gegenbauer_distm[[qq]]
})


inv_R = prec_ou(rep(1,n_t),phi_beta_now,1,n_t)
R = solve(inv_R)
log_det_R = 2 * sum(log(diag(chol(R))))
R_row_sum = apply(inv_R,1,sum)
R_quad_form = sum(inv_R)

Sig_now = exp(-pars_now[n_poly+1] * dd) * Reduce('+', Sig_now_terms) + pars_now[n_poly+2] * diag(n)
chol_sig_now = chol(Sig_now)
log_sig_det = 2 * sum(log(diag(chol_sig_now)))
inv_sig = chol2inv(chol_sig_now)

XsigX = quad.form(inv_sig,X_block)

like_now = likelihood(Y - xb,log_sig_det,inv_sig)


scale_fac1 = 2.38^2 / (n_poly+2)
cand_keep =  0.0001 * diag(n_poly+2)
cand_keep[1:n_poly,1:n_poly] = cov(log_b_now_sim)/10
cand_keep[n_poly+1,n_poly+1] = 0.01
cand_keep[n_poly+2,n_poly+2] = 0.1

cand_var1 = scale_fac1 *cand_keep


count1 = 0 ; chol1 = chol(cand_var1); inv1 = chol2inv(chol1); log_det1 = 2 * sum(log(diag(chol1)))

cand_var_phi = 0.1
count_phi = 0

st = proc.time()

for(i in 2:(reps + burn)){
  
  ######################### update Beta
  inv_prior = inv_R %x% V_reg_inv_now
  
  v_bet = solve(XsigX + inv_prior)
  m_bet = t(X_block) %*% c(inv_sig %*% Y) + inv_prior %*% rep(beta_tot_now,times = n_t)
  
  beta_reg_now = matrix(mvrnorm(1,v_bet %*% m_bet, v_bet),ncol = n_t)
  xb = unlist(lapply(1:n_t,function(qq){ X_t[[qq]] %*% beta_reg_now[,qq]}))
  
  like_now = likelihood(Y - xb,log_sig_det,inv_sig)
  
  ##### update hierarchical mean and variance  -- V and R
  
  v_bet = solve(V_reg_inv_now * R_quad_form + diag(p)/100)
  m_bet = (V_reg_inv_now %*% (beta_reg_now %*% R_row_sum))
  beta_tot_now = mvrnorm(1,v_bet %*% m_bet, v_bet)
  
  
  psi = quad.tform(inv_R,sweep(beta_reg_now,1,beta_tot_now,"-")) + 1 * diag(p)
  nu_new = n_t + p + 1
  
  V_reg_now = riwish(nu_new,psi)
  log_det_V = determinant(V_reg_now)$modulus[1]
  V_reg_inv_now = solve(V_reg_now)
  
  phi_cand = phi_beta_now
  phi_cand = rlnorm(1,log(phi_beta_now),cand_var_phi)
  
  if(phi_cand < 10){
    
    inv_R_cand = prec_ou(rep(1,n_t),phi_cand,1,n_t)
    R_cand = solve(inv_R_cand)
    log_det_R_cand = 2 * sum(log(diag(chol(R_cand))))
    R_row_sum_cand = apply(inv_R_cand,1,sum)
    R_quad_form_cand = sum(inv_R_cand)
    
    MH_dif = dlnorm(phi_beta_now,log(phi_cand),cand_var_phi,log = TRUE) - 
      dlnorm(phi_cand,log(phi_beta_now),cand_var_phi,log = TRUE)  
    
    prior_dif = dmatnorm(t(beta_reg_now),beta_tot_now,inv_R_cand,V_reg_inv_now,log_det_R_cand,log_det_V) - 
      dmatnorm(t(beta_reg_now),beta_tot_now,inv_R,V_reg_inv_now,log_det_R,log_det_V) + 
      dgamma(phi_cand,1,2,log = TRUE) - dgamma(phi_beta_now,1,2,log = TRUE)
    
    if(prior_dif + MH_dif > log(runif(1))){
      
      phi_beta_now = phi_cand
      
      inv_R = inv_R_cand
      R = R_cand
      log_det_R = log_det_R_cand
      R_row_sum = R_row_sum_cand
      R_quad_form = R_quad_form_cand
      
      count_phi = count_phi + 1
    }
    
  }
  
  
  ######################### update all pars
  
  log_now = log(pars_now)
  log_cand = log_now + t(chol1) %*% rnorm(n_poly + 2)
  cand = exp(log_cand)
  
  Sig_cand_terms = lapply(1:n_poly,function(qq){
    cand[qq] * gegenbauer_distm[[qq]]
  })
  
  Sig_cand = exp( -cand[n_poly + 1] * dd) * Reduce('+', Sig_cand_terms) + 
    (cand[n_poly+2] + 1e-3) * diag(n)
  
  chol_sig_cand= try(chol(Sig_cand),silent = TRUE)
  
  if(is.matrix(chol_sig_cand)){
    
    log_sig_det_cand = 2 * sum(log(diag(chol_sig_cand)))
    inv_sig_cand = chol2inv(chol_sig_cand)
    
    like_cand = likelihood(Y - xb,log_sig_det_cand,inv_sig_cand)
    
    prior_dif = dgamma(cand[n_poly + 1],2,24, log = TRUE ) +  
      sum(dcauchy(cand[-(n_poly + 1)],0,1, log = TRUE )) - 
      dgamma(pars_now[n_poly + 1],2,24, log = TRUE ) - 
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
      
      XsigX = quad.form(inv_sig,X_block)
      
      count1 = count1 + 1
    }
    
  }
  
  pars_save[i,] = pars_now
  
  
  
  if(i > burn){
    
    b[i - burn,] = pars_now[1:n_poly]
    phi[i - burn] = pars_now[n_poly+1]
    tau2[i - burn] = pars_now[n_poly+2]
    
    beta_reg[i - burn,,] = beta_reg_now
    beta_tot[i - burn,] = beta_tot_now
    V_reg[i - burn,,] = V_reg_now
    phi_beta[i - burn] = phi_beta_now
    
    like_save[i - burn] = like_now 
    
  }  
  
  
  if(i %% tune == 0){
    
    if(i < burn){
      
      acc1 = count1 / tune; count1 = 0
      
      acc_phi = count_phi / tune; count_phi = 0
      
      cand_var_phi = ifelse(acc_phi > .6, 2 * cand_var_phi,
                            ifelse(acc_phi < 0.15, cand_var_phi/3 ,cand_var_phi))
      
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

rm(list = setdiff(ls(),c("b","phi","tau2","beta_reg","beta_tot","V_reg","phi_beta","like_save")))


