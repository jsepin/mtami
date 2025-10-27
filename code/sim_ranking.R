# Not really a big difference between Holm and "free"

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


exp_fact <- expand.grid(n = seq(from = 200, to = 800, by = 50)
                        ,rho = seq(from = 0, to = 0.8, length.out =5) # (0:9)/10
                        ,beta_mindiff = c(0.25,0) #0.25
                        ,beta_prog = 0.25
                        ,p = 3 # c(2,3)
                        ,id = 1:1000)%>%
  mutate(id = 1:nrow(.))

# set.seed(21321)
cores <- parallel::detectCores()
res <- mclapply(mc.cores = cores, 1:nrow(exp_fact),function(i){
  
  # i = 1
  n <- exp_fact$n[i]
  rho <- exp_fact$rho[i]
  beta_mindiff <- exp_fact$beta_mindiff[i]
  beta_prog <- exp_fact$beta_prog[i]
  p <- exp_fact$p[i]
  
  
  
  trt <- rep(LETTERS[1:p],each = n)#rbinom(n, 1, 0.5)
  age <- rbinom(n = n*p, 1, p = 0.5)
  dd <- data.frame(trt = factor(trt), age, age_scaled = age - mean(age))
  Xmat <- model.matrix(data = dd, ~ trt + trt*age_scaled)
  dd$y <- Xmat%*%c(cumsum(rep(beta_mindiff,p)), beta_prog, rep(beta_prog,p-1 ))+rnorm(n=n*p)
  amp <- dd
  # summary(lm(data = dd, y ~ trt*age))$coef
  # summary(lm(data = dd, y ~ trt))$coef

  mi <- rbinom(n=nrow(dd),size=1,prob=dd$age*0.5)==1
  amp[mi,"y"] <- NA # ~30% missing
  
  summary(lm(data = dd, y ~ trt))$coef
  summary(lm(data = amp, y ~ trt))$coef
  
  
  # mi <- rbinom(n=n,size=1
  #              ,prob=plogis((bmi-mean(bmi))-2))==1
  # #mi[trt==1] <- FALSE
  # amp[mi,2] <- NA # ~30% missing
  
  # fit analysis models
  mfit_full <- lm(y ~ trt,data = dd)
  mfit_cc <- lm(y ~ trt,data = amp)
  

  # correct for multiplicity
  get_adj_pval <- function(fit,types = c("none","free", "Westfall","bonferroni", "holm")){
    # fit = mfit_full; types = c("none", "free","Westfall")
    name <-  deparse(substitute(fit))
    if(any(class(fit)=="mira")){
      test <- mypool_glht(fit, linfct = C)
    }else{
      test <- glht(fit, linfct = mcp(trt = "Tukey"))
    }
    res <- lapply(types,function(type){
      # type = "free"
      # type = "Westfall"
      s <- summary(test,test = adjusted(type = type))
      names(s$test$pval) <- names(s$test$sigma)# not defined...
      
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
     get_adj_pval(mfit_full)
    ,get_adj_pval(mfit_cc)
  ) %>%
    mutate(id = i, p_missing =   mean(rowSums(is.na(amp))>0)
           #, rho_actual = mygetcor(vcov(mfit_full),C)[1,2]
           , .before = 1)
  # gc()
  out
})

res <- do.call("bind_rows",res)
res <- left_join(res, exp_fact,by = "id")



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Visualization ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha <- 0.05
res$sig <- as.integer(rowSums(res[,str_detect(colnames(res), "pval")] < alpha,na.rm = TRUE) >= ifelse(res$beta_mindiff==0, 1, res$p) )
my_cov <- function(x){
  out <- prop.test(x = sum(x,na.rm=TRUE), n = length(x))
  out <- data.frame(out$estimate, t(out$conf.int), row.names = NULL)
  names(out) <- c("cov_estimate", "cov_lower", "cov_upper")
  return(out)
}


res_sum <- res %>%
  group_by(type, fit, n, rho, beta_mindiff,p)%>%
  summarise(my_cov(sig),Nsim=n())

dodge <- 0 # 50
# making a type 1 error assuming correlation is 0
line_df <- data.frame(beta_mindiff = 0,p = unique(res$p), yintercept = c( 1-( 1-alpha)^unique(res$p) ))
line_df <- bind_rows(line_df,data.frame(beta_mindiff = 0,p = unique(res$p), yintercept = alpha ))
line_df <- bind_rows(line_df
                     ,data.frame(beta_mindiff = setdiff(unique(res$beta_mindiff),0), p = unique(res$p), yintercept = 0.8 ))

ggplot(data = res_sum #%>% filter(fit == "mfit_full")
       , aes(x = n,y = cov_estimate,col = type,shape = fit, linetype = fit))+
  geom_point(position=position_dodge(width=dodge))+
  geom_errorbar(aes(ymin = cov_lower,ymax = cov_upper),position=position_dodge(width=dodge))+
  geom_line()+
  #geom_smooth(se = F)+
  geom_hline(data = line_df, aes(yintercept = yintercept ))+
  facet_grid(beta_mindiff+p~rho, labeller = label_both,scale = "free")+
  theme_bw()+
  labs(y = "Rejection of Null Hypothesis")+
  scale_x_continuous(breaks = seq(min(res_sum$n), max(res_sum$n), by = 50))+
  scale_y_continuous(breaks = round(seq(0, max(res_sum$cov_estimate), by = 0.05),2))
