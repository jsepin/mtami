rm(list = ls());gc()
library(tidyverse)
library(mice)
library(multcomp)
library(parallel)
library(mmrm)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Helper function ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# adjusted autoregressive correlation matrix
corar1 <- function(rho,p){
  times <- 1:p
  sigma <- 1
  H <- abs(outer(times, times, "-"))
  V <- sigma * rho^H
  p <- nrow(V)
  V[cbind(1:p, 1:p)] <- V[cbind(1:p, 1:p)] * sigma
  return(V)
}

# correct for multiplicity
get_adj_pval <- function(fit,C,types = c("none","free","bonferroni", "holm", "Westfall")){
  # fit = mfit_full; types = c("none", "free")
  name <-  deparse(substitute(fit))
  if(any(class(fit)=="mira")){
    test <- mypool_glht(fit, linfct = C)
  }else{
    test <- glht(fit, linfct = C)
    test$coef <- c(test$coef) # quick fix for some...
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


mygetcor <- function(Vbetahat,C){
  Sigma <- diag(1 / sqrt(diag(C %*% Vbetahat %*% t(C))))
  Cor <- Sigma %*% (C %*% Vbetahat %*% t(C)) %*% t(Sigma)
  Cor
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Simulation ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exp_fact <- expand.grid(n = seq(from = 200, to = 600, by = 20)
                        ,rho = seq(from = 0, to = 0.8, length.out = 3) # (0:9)/10
                        ,beta_trt  = c(1,0) #0.25
                        ,beta_time = 1
                        ,beta_prog = 0.1
                        ,p = c(3:4)
                        )%>%
  mutate(exp_fact_id = 1:nrow(.))
print(nrow(exp_fact))


# set.seed(21321)
cores <- parallel::detectCores()
Nsim <- 5000
res <- lapply(1:nrow(exp_fact),function(e){
  # e = 1
  print(paste0(e, " out of ", nrow(exp_fact)))
  n <- exp_fact$n[e]
  rho <- exp_fact$rho[e]
  beta_trt <- exp_fact$beta_trt[e]
  beta_time <- exp_fact$beta_time[e]
  beta_prog <- exp_fact$beta_prog[e]
  p <- exp_fact$p[e]
  

  out <- mclapply(mc.cores = cores, 1:Nsim,function(i){
    gc()
    # i = 1
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## DGP
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # n = 2000; rho = 0; beta_trt = 1; beta_time = 1; beta_prog = 1
    # age <- rnorm(n=n,mean = 50,sd = 3)
    # trt <- rbinom(n=n,size = 1, prob = 0.5)
    trt <- rep(c(0,1),each = n/2)#rbinom(n, 1, 0.5)
    age <- rep(c(0,1),n/2)
    age <- rnorm(n = n, mean = 50, sd = 2)#2
    sigma <- corar1(rho = rho, p = p)*16
    
    mod_mat <- replicate(p,rep(1, n))%*%diag(cumsum(rep(beta_time,p)))+replicate(p,trt)*beta_trt
    noise_mat <- mvtnorm::rmvnorm(n = n, sigma = sigma)
    prog_mat  <- mvtnorm::rmvnorm(n = n, sigma = sigma)
    age <- prog_mat[,1]
    total_unexp <- prog_mat%*%diag((rep(beta_prog,p)))+noise_mat
    Y <- mod_mat + total_unexp
    dd_long <- dd_long_amp <- data.frame(Y,USUBJID = as.factor(1:n)
                                         , "TRT" = as.factor(trt) 
                                         , "AGE" = age
    )
    dd_long <- dd_long %>%
      pivot_longer(cols = -c(USUBJID,TRT,AGE), names_to = "AVISIT", values_to = "Y")  %>%
      mutate(AVISIT = as.factor(gsub("[^0-9.-]", "", AVISIT)))
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## MGP (probably quite high!)
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u <- matrix(runif(n = n*p), ncol = p)
    mi <- u > plogis(scale(prog_mat, scale =F))
    mi[trt==0,] <- FALSE
    colMeans(mi)
    dd_long_amp[which(mi, arr.ind = TRUE)] <- NA
    
    # mi <- rbinom(n=n,size=1,prob=plogis(2*(age - mean(age))))==1
    # mi[trt==0] <- FALSE # max 12.5% missing
    # mean(mi)
    # dd_long_amp[which(mi),1:p] <- NA # ~30% missing
    dd_long_amp <- dd_long_amp %>%
      pivot_longer(cols = -c(USUBJID,TRT,AGE), names_to = "AVISIT", values_to = "Y")  %>%
      mutate(AVISIT = as.factor(gsub("[^0-9.-]", "", AVISIT)))
    
  
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## MMRM
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    # ggplot(dd_long, aes(x =AGE, y= Y, col =TRT))+
    #   geom_point()+
    #   geom_smooth(method ="lm")+
    #   facet_wrap(~AVISIT)
    

    fit_full <- mmrm(
      formula = Y ~ TRT*AVISIT + us(AVISIT | USUBJID),
      data = dd_long
    )
    fit_cc <- mmrm(
      formula = Y ~ TRT*AVISIT + us(AVISIT | USUBJID),
      data = dd_long_amp
    )
    # summary(fit_full)$coef
    # summary(fit_cc)$coef
    # summary(lm(data = dd_long,Y~TRT*AVISIT))$coef
    # summary(lm(data = dd_long%>%filter(AVISIT==1),Y~TRT))$coef; summary(lm(data = dd_long_amp%>%filter(AVISIT==1),Y~TRT))$coef
    # summary(lm(data = dd_long%>%filter(AVISIT==2),Y~TRT))$coef; summary(lm(data = dd_long_amp%>%filter(AVISIT==2),Y~TRT))$coef
    # summary(lm(data = dd_long%>%filter(AVISIT==3),Y~TRT))$coef; summary(lm(data = dd_long_amp%>%filter(AVISIT==3),Y~TRT))$coef
    
    
    # Contrast matrix
    C_index <- str_detect(colnames(vcov(fit_full)), "TRT")
    C <- lapply(1:length(C_index),function(j){
      # j = 1
      out <- rep(0,length(C_index))
      out[j] <- as.integer(C_index[j])
      out
    })
    C <- do.call("rbind", C)
    C <- C[rowSums(C)!=0,]
    C[,C[1,]==1]<- 1
    rownames(C) <- paste0("TRT", 1:p)
    C
    # rownames(C) <- colnames(vcov(fit_full))[str_detect(colnames(vcov(fit_full)), "TRT")]
    rho_actual <- mygetcor(vcov(fit_full),C)
    rho_actual <- rho_actual[upper.tri(rho_actual)]
    names(rho_actual) <- apply(t(combn(1:p,2)),1,function(x) paste(x,collapse = ""))
    
    out <- bind_rows(
      get_adj_pval(fit_full,C)
      ,get_adj_pval(fit_cc,C)
    ) %>%
      mutate(id = i
             #, p_missing =   mean(rowSums(is.na(amp))>0)
             #, rho_actual = t(rho_actual)
             , .before = 1)
    # gc()
    out
  })
  out <- do.call("rbind",out)%>%
    mutate(exp_fact_id = e)%>%
    left_join(.,exp_fact, by = "exp_fact_id")
  out

})
gc()
res <- do.call("bind_rows",res)
# res <- left_join(res, exp_fact,by = "id")
# saveRDS(res, "../data/res.RDS")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Visualization ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k <- 1
alpha <- 0.05

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
# line_df <- bind_rows(line_df
#                      ,data.frame(beta_trt = setdiff(unique(res$beta_trt),0), p = unique(res$p), yintercept = 0.8 ))

ggplot(data = res_sum #%>% filter(fit == "fit_full")
       , aes(x = n,y = cov_estimate,col = type,shape = fit, linetype = fit))+
  # geom_point(position=position_dodge(width=dodge))+
  # geom_errorbar(aes(ymin = cov_lower,ymax = cov_upper),position=position_dodge(width=dodge))+
  geom_line()+
  #geom_smooth(se = F)+
  geom_hline(data = line_df, aes(yintercept = yintercept ))+
  facet_grid(p+beta_trt~rho, labeller = label_both,scale = "free")+
  theme_bw()+
  labs(y = "P(Rejecting at least one Null Hypothesis)")+
  # scale_x_continuous(breaks = seq(min(res_sum$n), max(res_sum$n), by = 50))+
  scale_y_continuous(breaks = round(seq(0, max(res_sum$cov_estimate), by = 0.05),2))+
  geom_hline(yintercept = c(0.8),col = "gray")




ggplot(data = res %>% filter(type == "none")%>%
         dplyr::select(-starts_with(c("se","pval"))) %>%
         pivot_longer(cols = starts_with("est"))
       , aes(x = as.factor(n),y = (value-beta_trt),fill = paste0(name,"+", fit)))+
  geom_boxplot()+
  #coord_cartesian(ylim = c(0.5,-0.5))+
  facet_grid(p+beta_trt~rho, labeller = label_both,scale = "free")+
  geom_hline(aes(yintercept = 0))+
  theme_bw()

