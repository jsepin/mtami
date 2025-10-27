
rm(list = ls());gc()
library(tidyverse)
library(mice)
library(multcomp)
library(parallel)

mygetcor <- function(Vbetahat,C){
  Sigma <- diag(1 / sqrt(diag(C %*% Vbetahat %*% t(C))))
  Cor <- Sigma %*% (C %*% Vbetahat %*% t(C)) %*% t(Sigma)
  Cor
}



exp_fact <- expand.grid(n = seq(from =400, to = 800, by = 50)
                        ,rho = seq(from = 0, to = 0.8, length.out =5) # (0:9)/10
                        ,beta_trt = c(0,0.25) #0.25
                        ,beta_prog = 0.25
                        ,p = c(2:5)
                        ,id = 1:1000)%>%
  mutate(id = 1:nrow(.))
print(nrow(exp_fact))

# set.seed(21321)
cores <- parallel::detectCores()
res <- mclapply(mc.cores = cores, 1:nrow(exp_fact),function(i){

  # i = 1
  n <- exp_fact$n[i]
  rho <- exp_fact$rho[i]
  beta_trt <- exp_fact$beta_trt[i]
  beta_prog <- exp_fact$beta_prog[i]
  p <- exp_fact$p[i]
  
  trt <- rep(c(0,1),each = n/2)#rbinom(n, 1, 0.5)
  age <- rep(c(0,1),n/2)
  Y <- replicate(p, trt)%*%diag(rep(beta_trt,p))
  # age_mat <- cbind(age - mean(age), age - mean(age))%*%diag(c(beta_prog,beta_prog))
  age_mat <- replicate(p, age - mean(age))%*%diag(rep(beta_prog,p))
  # sigma=matrix(c(1,rho,rho,1),2,2)
  sigma <- diag(p)
  sigma[upper.tri(sigma)] <- sigma[lower.tri(sigma)] <- rho
  noise_mat <- mvtnorm::rmvnorm(n,sigma=sigma)
  noise_mat <- noise_mat - colMeans(noise_mat)
  Y <- Y + age_mat + noise_mat
  colnames(Y) <- paste0("Outcome", LETTERS[1:p])
  dd <- amp <- data.frame(Y,trt,age)
  
  # summary(lm(data = dd, OutcomeA ~ trt*age))$coef
  # summary(lm(data = dd, OutcomeB ~ trt))$coef

  
  mi <- rbinom(n=nrow(dd),size=1,prob=dd$age*0.5)==1
  mi[trt==1] <- FALSE
  amp[mi,"OutcomeA"] <- NA # ~30% missing
  
  # mi <- rbinom(n=n,size=1
  #              ,prob=plogis((bmi-mean(bmi))-2))==1
  # #mi[trt==1] <- FALSE
  # amp[mi,2] <- NA # ~30% missing

  # fit analysis models
  f <- paste0("Outcome", LETTERS[1:p], collapse = ", ")
  f <- paste0("cbind(", f, ") ~ 1 + trt")
  mfit_full <- lm(formula = formula(f),data = dd)
  mfit_cc <- lm(formula = formula(f),data = amp)
  
  # Contrast matrix
  C_index <- str_detect(colnames(vcov(mfit_full)), "trt") & !str_detect(colnames(vcov(mfit_full)), "trt:")
  C <- lapply(1:length(C_index),function(j){
    out <- rep(0,length(C_index))
    out[j] <- as.integer(C_index[j])
    out
  })
  C <- do.call("rbind", C)
  C <- C[rowSums(C)!=0,]
  rownames(C) <- LETTERS[1:p]
  
  
  # correct for multiplicity
  get_adj_pval <- function(fit,C,types = c("none","free","bonferroni", "holm", "Westfall")){
    # fit = mfit_full; types = c("none", "free")
    name <-  deparse(substitute(fit))
    if(any(class(fit)=="mira")){
      test <- mypool_glht(fit, linfct = C)
    }else{
      test <- glht(fit, linfct = C)
      test$coef <- c(test$coef) # quick fix...
    }
    res <- lapply(types,function(type){
      # type = "none"
      s <- summary(test,test = adjusted(type = type))
      names(s$test$pval) <- names(s$test$sigma)# not always defined...
      
      data.frame("type" = s$test$type
                 ,"fit" = name
                 ,"est" = t(s$test$coefficients) #t(C%*%test$coef)
                 ,"se" = t(s$test$sigma) #t(C%*%sqrt(diag(test$vcov)))
                 ,"pval" = t(s$test$pval)
                 )
    })
    res <- do.call("rbind",res)
    res
  }

  # return outcome
  out <- bind_rows(
     get_adj_pval(mfit_full,C)
    ,get_adj_pval(mfit_cc,C)
  ) %>%
    mutate(id = i, p_missing =   mean(rowSums(is.na(amp))>0)
           , rho_actual = mygetcor(vcov(mfit_full),C)[1,2]
           , .before = 1)
  # gc()
  out
})
gc()
res <- do.call("bind_rows",res)
res <- left_join(res, exp_fact,by = "id")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Visualization ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k <- 1
alpha <- 0.05

# fwer_k_of_p <- function(p, k, alpha=0.05){
#   s <- 0
#   for(i in k:p) s <- s + choose(p,i) * alpha^i * (1-alpha)^(p-i)
#   return(s)
# }
# fwer_k_of_p(p=2, k = 1)
# fwer_k_of_p(p=3, k = 1)
# fwer_k_of_p(p=3, k = 3)
# # 1-( 1-alpha)^p

res$sig <- as.integer(rowSums(res[,str_detect(colnames(res), "pval")] < alpha,na.rm = TRUE) >= k )
my_cov <- function(x){
  out <- prop.test(x = sum(x,na.rm=TRUE), n = length(x))
  out <- data.frame(out$estimate, t(out$conf.int), row.names = NULL)
  names(out) <- c("cov_estimate", "cov_lower", "cov_upper")
  return(out)
}


res_sum <- res %>%
  group_by(type, fit, n, rho, beta_trt,p)%>%
  summarise(my_cov(sig),Nsim=n())

dodge <- 15#15#50
# making a type 1 error assuming correlation is 0
line_df <- data.frame(beta_trt = 0,p = unique(res$p), yintercept = c( 1-( 1-alpha)^unique(res$p) ))
line_df <- bind_rows(line_df,data.frame(beta_trt = 0,p = unique(res$p), yintercept = alpha ))
line_df <- bind_rows(line_df
                     ,data.frame(beta_trt = setdiff(unique(res$beta_trt),0), p = unique(res$p), yintercept = 0.8 ))

ggplot(data = res_sum %>% filter(fit == "mfit_full")
       , aes(x = n,y = cov_estimate,col = type,shape = fit, linetype = fit))+
  geom_point(position=position_dodge(width=dodge))+
  geom_errorbar(aes(ymin = cov_lower,ymax = cov_upper),position=position_dodge(width=dodge))+
  geom_line()+
  #geom_smooth(se = F)+
  geom_hline(data = line_df, aes(yintercept = yintercept ))+
  facet_grid(p+beta_trt~rho, labeller = label_both,scale = "free")+
  theme_bw()+
  labs(y = "Rejection of Null Hypothesis")+
  scale_x_continuous(breaks = seq(min(res_sum$n), max(res_sum$n), by = 50))+
  scale_y_continuous(breaks = round(seq(0, max(res_sum$cov_estimate), by = 0.05),2))


# ggplot(data = res_sum
#        , aes(x = n,y= rho,z = cov_estimate, linetype = type))+
#   # geom_contour(bins = 20)+
#   # geom_contour(aes(colour = after_stat(level)))+
#   geom_contour_filled()+
#   facet_grid(beta_trt~fit, labeller = label_both,scale = "free")+
#   theme_bw()
# 
# 
# ggplot(data = res %>% filter(type == "none")%>%pivot_longer(cols = c(est.A, est.B))
#        , aes(x = as.factor(n),y = value,fill = paste0(name,"+", fit)))+
#   geom_boxplot()+
#   facet_grid(beta_trt~rho, labeller = label_both,scale = "free")+
#   geom_hline(aes(yintercept = beta_trt))+
#   theme_bw()
# 
# ggplot(data = res %>% filter(type == "none",fit == "mfit_full")
#        , aes(x = n,y = rho_actual,col = type,shape = fit))+
#   geom_point()+
#   geom_line()+
#   facet_grid(beta_trt~rho, labeller = label_both,scale = "free")+
#   theme_bw()
# 



