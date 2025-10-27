rm(list = ls())
library(tidyverse)
library(mice)
library(multcomp)
library(parallel)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Helper functions ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mypool <- function(fit){
  vcov_list <- lapply(fit$analyses,function(x)vcov(x))
  coef_list <- lapply(fit$analyses,function(x)c(coef(x)))
  m <- length(coef_list)
  U_bar <- Reduce("+", vcov_list) / m
  coef_mat <- do.call("rbind",coef_list)
  colnames(coef_mat) <- colnames(U_bar)
  Qbar <- colMeans(coef_mat)
  B <- var(coef_mat)
  TotVar <- U_bar + (1 + 1/m)*B
  list("coef" = Qbar, "vcov" = TotVar)
}


mypool_glht <- function(fit,...){
  out <- glht(fit$analyses[[1]], ...)# copy structure
  pooled_est <- mypool(fit)
  out$coef <- pooled_est$coef
  out$vcov <- pooled_est$vcov
  out$df <- Inf #FIXME! normal approximation or really using the pooled df?
  out$model <- NULL
  return(out)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Very simple co-primary endpoint (one enough to show efficacy) ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exp_fact <- expand.grid(n = (3:10)*100
                        ,rho = 0 # c(0, 0.3, 0.5, 0.7)
                        ,beta_trt = 0#c(0,0.5) #0.25
                        ,p = 3
                        ,m = 5 # 20 
                        ,id = 1:10000)%>%
  mutate(id = 1:nrow(.))

set.seed(21321)
cores <- parallel::detectCores()
res <- mclapply(mc.cores = cores, 1:nrow(exp_fact),function(i){
  # i = 1
  print(i)
  n <- exp_fact$n[i]
  rho <- exp_fact$rho[i]
  beta_trt <- exp_fact$beta_trt[i]
  p <- exp_fact$p[i]
  m <- exp_fact$m[i]
  
  # RCT with missingness depending on age
  beta <- rep(beta_trt,p)
  sigma <- diag(p)
  sigma[upper.tri(sigma)] <- sigma[lower.tri(sigma)] <- rho
  Y <- mvtnorm::rmvnorm(n=n,sigma = sigma)
  colnames(Y) <- paste0("Outcome",LETTERS[1:p])
  trt <- rbinom(n = n, size = 1, prob = 0.5)
  age <- rnorm(n=n,mean = 50,sd = 3)
  Y <- Y + trt%*%t(beta) + 0.5*age + 0.5*trt*(age-mean(age))
  dd <- amp <- data.frame(Y, trt, age)
  mi <- rbinom(n=n,size=1
               ,prob=plogis((age-mean(age))-2))==1
  amp[mi,1:p] <- NA # ~30% missing
  
  # Contrast matrix
  C <- lapply(1:p,function(j){
    out <- rep(0,2*p)
    out[2*(j-1) + 2] <- 1
    out
  })
  C <- do.call("rbind", C)
  rownames(C) <- LETTERS[1:p]
  
  # multiple imputation with treatment interaction
  imp <- mice(data = amp %>% mutate(int = trt*age)
              , method = "norm", m = 1, maxit = 0, printFlag = F,eps = 0)
  pm <- imp$predictorMatrix
  pm[1:p,1:p] <- 0
  imp <- mice(data = amp %>% mutate(int = trt*age)
              , method = "norm", m = m, maxit = 1, printFlag = F,eps = 0
              ,predictorMatrix = pm)
  
  # if(sum(is.na(complete(imp, include = F)))>0){
  #   stop("Still missing!")
  # }
  
  # fit analysis models
  f <- paste0("Outcome", LETTERS[1:p], collapse = ", ")
  f <- paste0("cbind(", f, ") ~ trt")
  mfit_mi <- with(data = imp, exp = lm(formula = formula(f)) )
  mfit_full <- lm(formula = formula(f),data = dd)
  mfit_cc  <- lm(formula = formula(f),data = amp)
  
  # correct for multiplicity
  get_adj_pval <- function(fit,C,types = c("none", "bonferroni", "holm","free")){
    # fit = mfit_full; types = c("none", "bonferroni", "holm","free")
    name <-  deparse(substitute(fit))
    if(any(class(fit)=="mira")){
      test <- mypool_glht(fit, linfct = C)
    }else{
      test <- glht(fit, linfct = C)
      test$coef <- c(test$coef) # quick fix..
    }
    res <- lapply(types,function(type){
      # type = "none"
      s <- summary(test,test = adjusted(type = type))
      
      data.frame("type" = s$test$type
                 , "fit" = name
                 ,"est" = s$test$coefficients #t(C%*%test$coef)
                 ,"se" = s$test$sigma #t(C%*%sqrt(diag(test$vcov)))
                 ,"pval" = t(s$test$pval) )
    })
    res <- do.call("rbind",res)
    res
  }
  
  # return outcome
  out <- bind_rows(
    get_adj_pval(mfit_mi,C)
    ,get_adj_pval(mfit_full,C)
    ,get_adj_pval(mfit_cc,C)
  )%>%
    mutate(id = i, p_missing = mean(mi), .before = 1)
  out
})
res <- do.call("rbind",res)
res <- left_join(res, exp_fact,by = "id")


# saveRDS(res, "data/res.RDS")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Visualization ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha <- 0.05
k <- 1 # number of successes
# res <- res %>% mutate(sig = as.integer(pval.A < alpha | pval.B<alpha))
res$sig <- as.integer(rowSums(res[,str_detect(colnames(res), "pval")] < alpha) >= k)

# all.equal(with(res, as.integer(pval.A < alpha | pval.B<alpha))
# ,as.integer(rowSums(res[,str_detect(colnames(res), "pval")] < alpha) >= k)
# )

# dd_check <- expand.grid(pval.A = c(0.1, 0.04), pval.B = c(0.1, 0.04))%>%
#   mutate(sig = as.integer(pval.A < alpha | pval.B<alpha))

my_cov <- function(x){
  out <- prop.test(x = sum(x,na.rm=TRUE), n = length(x))
  out <- data.frame(out$estimate, t(out$conf.int), row.names = NULL)
  names(out) <- c("cov_estimate", "cov_lower", "cov_upper")
  return(out)
}


res_sum <- res %>%
  group_by(type, fit, n, rho, beta_trt)%>%
  summarise(my_cov(sig))

dodge <- 0 #15#50
line_df <- data.frame(beta_trt = 0, yintercept = c(alpha, 1-( 1-alpha)^unique(res$p) ))
ggplot(data = res_sum %>% filter(fit == "mfit_full")
       , aes(x = n,y = cov_estimate,col = type,shape = fit))+
  geom_point(position=position_dodge(width=dodge))+
  geom_errorbar(aes(ymin = cov_lower,ymax = cov_upper),position=position_dodge(width=dodge))+
  geom_line()+
  geom_hline(data = line_df, aes(yintercept = yintercept ))+
  facet_grid(beta_trt~rho, labeller = label_both,scale = "free")+
  theme_bw()+
  labs(y = "Rejection of Null Hypothesis")+
  scale_x_continuous(breaks = unique(res_sum$n))+
  scale_y_continuous(breaks = round(seq(0, max(res_sum$cov_estimate), by = 0.05),2))
