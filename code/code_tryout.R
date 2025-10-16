# Jäger2022-Obesity in the context of migration and socio-economic risk factors – a multivariate epidemiologic analysis
# Schmidt-Kraepelin2022-Amisulpride and olanzapine combination treatment versus each monotherapy in acutely ill patients with schizophrenia in Germany (COMBINE) a double-blind randomised controlled trial

# FIMD Chapter 5.3 Multi-parameter inference
# https://stefvanbuuren.name/fimd/sec-multiparameter.html


rm(list = ls())
library(tidyverse)
library(mice)
library(multcomp)

n <- 50
sigma <- c( 1.0,0.5,0.5
           ,0.5,1.0,0.5
           ,0.5,0.5,1.0)
m <- 5
est <- lapply(1:m,function(...){
  sigma <- matrix(sigma,ncol =3)
  beta <- mvtnorm::rmvnorm(n = 1, mean = c(1,2,3),sigma = sigma)
  X <- mvtnorm::rmvnorm(n, mean = c(0,0,0))
  y <- c(X%*%c(beta) + rnorm(n))
  fit <- lm(y ~ X-1)
  list( "coef" = coef(fit),"vcov" = vcov(fit)) 
  })

coef_list <- lapply(est,function(x)x$coef)
vcov_list <- lapply(est,function(x)x$vcov)

Qbar <- Reduce("+", coef_list) /length(coef_list)
Qbar

U_bar <- Reduce("+", vcov_list) / length(vcov_list)
coef_mat <- do.call("rbind",coef_list)
B <- var(coef_mat)
(T <- U_bar + (1 + 1/m)*B)


library(multcomp)
data("warpbreaks")
m <- lm(breaks ~ tension-1, data = warpbreaks)
coef(m)
vcov(m)

Xmat <- model.matrix(m)
y <- Xmat%*%coef(m)+rnorm(n=nrow(Xmat), sd = sigma(m))
dd <- amp <- data.frame("tension" = warpbreaks$tension, "breaks" = y)
fit_full <- lm(breaks ~ tension, data = dd)
amp$breaks[rbinom(n=nrow(amp), size=1, p = 0.3)==1] <- NA
fit_drop <- lm(breaks ~ tension, data = amp)

imp <- mice(data = amp, method = "norm", m = 50, maxit = 1, printFlag = F)
fit <- with(data = imp, exp = lm(breaks ~ tension))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions ---- 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# matrix(1:6, ncol= 2); c(matrix(1:6, ncol= 2))

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
# Simulation ---- 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigma <- c(1.0,0.9,0.1
          ,0.9,1.0,0.5
          ,0.1,0.5,1.0)
sigma <- matrix(sigma,ncol =3)


p <- 3
sigma <- diag(p)
d <- 0.1
n <- 300

mean <- cumsum(rep(d,ncol(sigma)))
Y <- mvtnorm::rmvnorm(n = n, mean = mean, sigma = sigma)
# plot(Y[,1:2])
# dd_long <- data.frame(Y,ID = 1:n) %>%
#   pivot_longer(cols = -ID)  %>%
#   mutate(name = as.factor(name))
# fit_long <- lm(data = dd_long, value ~ name -1 ) #+ as.factor(ID)
# C <- diag(length(coef(fit_long)))
# diag(C) <- 0
# C[1,1] <- C[2,2] <- C[3,3] <- 1
# test_long <- glht(fit_long, linfct = mcp(name = "Dunnett"))
# test_long <- glht(fit_long, linfct = C )
# test_long$linfct
# Vbetahat <- vcov(fit_long)
# betahat <- coef(fit_long)
# betahat
# Sigma <- diag(1 / sqrt(diag(C %*% Vbetahat %*% t(C))))
# t <- Sigma %*% C %*% betahat
# Cor <- Sigma %*% (C %*% Vbetahat %*% t(C)) %*% t(Sigma)
# Cor

trt <- sample(1:length(mean), n, replace = TRUE)
dd <- amp <-  data.frame("y"=Y[cbind(1:length(trt),trt)], "trt" = factor(LETTERS[trt]))
amp$y[rbinom(n=nrow(amp), size=1, p = 0.3)==1] <- NA

fit_full <- with(data = dd, exp = lm(formula( f) ))
test_full <- glht(fit_full, linfct = mcp(trt = "Tukey"))

C <- (test_full$linfct)
Vbetahat <- vcov(fit_full)
betahat <- coef(fit_full)
Sigma <- diag(1 / sqrt(diag(C %*% Vbetahat %*% t(C))))
# t <- Sigma %*% C %*% betahat
Cor <- Sigma %*% (C %*% Vbetahat %*% t(C)) %*% t(Sigma)
Cor

f <- "y ~ trt"
imp <- mice(data = amp, method = "norm", m = 50, maxit = 1, printFlag = F,eps = 0)
imp$loggedEvents
if(sum(is.na(complete(imp, include = F)))>0){
  stop("Still missing!")
}
fit <- with(data = imp, exp = lm(formula( f) ))



test <- mypool_glht(fit, linfct = mcp(trt = "Tukey"))
types <- c("none", "holm", "free", "Westfall")
out <- lapply(types,function(type){
  # type = "holm"
  summary(test,test = adjusted(type = type))$test$pval
})
out <- data.frame(type = types, do.call("rbind", out))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Very simple co-primary endpoint (one enough to show efficacy) ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exp_fact <- expand.grid(n = (3:10)*100
                        ,rho = c(0.3, 0.5, 0.7)
                        ,beta_trt = c(0,0.5) #0.25
                        ,p = 2
                        ,id = 1:1000
                        )%>%
  mutate(id = 1:nrow(.))

set.seed(21321)
res <- lapply(1:nrow(exp_fact),function(i){
  # i = 1
  print(i)
  n <- exp_fact$n[i]
  rho <- exp_fact$rho[i]
  beta_trt <- exp_fact$beta_trt[i]
  p <- exp_fact$p[i]
  

  beta <- rep(beta_trt,p)
  sigma <- diag(p)
  sigma[upper.tri(sigma)] <- sigma[lower.tri(sigma)] <- rho
  Y <- mvtnorm::rmvnorm(n=n,sigma = sigma)
  colnames(Y) <- paste0("Outcome",LETTERS[1:p])
  trt <- rbinom(n = n, size = 1, prob = 0.5)
  age <- rnorm(n=n,mean = 50,sd = 3)
  Y <- Y + trt%*%t(beta) + 0.5*age + 0.5*trt*(age-mean(age))
  dd <- amp <- data.frame(Y, trt, age)
  # RCT with missingness depending on age
  mi <- rbinom(n=n,size=1
               ,prob=plogis((age-mean(age))-2))==1
  amp[mi,1:p] <- NA # ~30% missing
  
  # Contrast matrix
  # C <- rbind("A" = c(0,1,0,0)
  #            ,"B" = c(0,0,0,1)
  # )
  # C <- matrix(0,nrow = p, ncol = 2*p)
  # rownames(C) <- LETTERS[1:p]
  C <- lapply(1:p,function(j){
    out <- rep(0,2*p)
    out[2*(j-1) + 2] <- 1
    out
  })
  C <- do.call("rbind", C)
  rownames(C) <- LETTERS[1:p]
  
  
  imp <- mice(data = amp %>% mutate(int = trt*age)
              , method = "norm", m = 1, maxit = 0, printFlag = F,eps = 0)
  pm <- imp$predictorMatrix
  pm[1:p,1:p] <- 0
  imp <- mice(data = amp %>% mutate(int = trt*age)
              , method = "norm", m = 50, maxit = 1, printFlag = F,eps = 0
              ,predictorMatrix = pm)
  
  imp$loggedEvents
  if(sum(is.na(complete(imp, include = F)))>0){
    stop("Still missing!")
  }
  f <- paste0("Outcome", LETTERS[1:p], collapse = ", ")
  #f <- as.formula(paste0("cbind(", f, ") ~ trt"))
  f <- paste0("cbind(", f, ") ~ trt")
  
  # mfit_mi <- with(data = imp, exp = lm(cbind(OutcomeA,OutcomeB) ~ trt ))
  # mfit_full  <- lm(cbind(OutcomeA,OutcomeB) ~ trt,data = dd)
  # mfit_cc  <- lm(cbind(OutcomeA,OutcomeB) ~ trt,data = amp)
  mfit_mi <- with(data = imp, exp = lm(formula = formula(f)) )
  mfit_full <- lm(formula = formula(f),data = dd)
  mfit_cc  <- lm(formula = formula(f),data = amp)
  
  
  get_adj_pval <- function(fit,C,types = c("none", "bonferroni", "holm","free")){
    # fit = mfit_full
    name <-  deparse(substitute(fit))
    if(any(class(fit)=="mira")){
      test <- mypool_glht(fit, linfct = C)
    }else{
      test <- glht(fit, linfct = C)
      test$coef <- c(test$coef) # quick fix..
    }
    res <- lapply(types,function(type){
      # type = "none"
      data.frame("type" = type, "fit" =name
                 ,"est" = t(C%*%test$coef)
                 ,"se" = t(C%*%sqrt(diag(test$vcov)))
                 ,"pval" = t(summary(test,test = adjusted(type = type))$test$pval))
    })
    res <- do.call("rbind",res)
    res
  }
  
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
saveRDS(res, "../data/res.RDS")
alpha <- 0.05
k <- 1 # number of successes
# res <- res %>% mutate(sig = as.integer(pval.A < alpha | pval.B<alpha))
res$sig <- as.integer(rowSums(res[,str_detect(colnames(res), "pval")] < alpha) >= k)

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

dodge <- 15#50
ggplot(data = res_sum %>% filter(fit == "mfit_full")
       , aes(x = n,y = cov_estimate,col = type,shape = fit))+
  geom_point(position=position_dodge(width=dodge))+
  geom_errorbar(aes(ymin = cov_lower,ymax = cov_upper),position=position_dodge(width=dodge))+
  geom_line()+
  geom_hline(yintercept =alpha )+
  facet_grid(beta_trt~rho, labeller = label_both,scale = "free")+
  theme_bw()+
  labs(y = "Rejection of Null Hypothesis (one of both)")+
  scale_x_continuous(breaks = unique(res_sum$n))+
  scale_y_continuous(breaks = round(seq(0, max(res_sum$cov_estimate), by = 0.05),2))



# # Get back correlation (to check)
# C <- diag(p)
# Vbetahat <- vcov(mfit)[-c(1,3),-c(1,3)]
# betahat <- coef(mfit)["trt",]
# Sigma <- diag(1 / sqrt(diag(C %*% Vbetahat %*% t(C))))
# # t <- Sigma %*% C %*% betahat
# Cor <- Sigma %*% (C %*% Vbetahat %*% t(C)) %*% t(Sigma)
# Cor




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# MMRM ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(mmrm)
require(emmeans)
fit_mmrm <- mmrm(
  formula = FEV1 ~ ARMCD*AVISIT + us(AVISIT | USUBJID),
  data = fev_data
)
coef(fit_mmrm)
emmeans(fit_mmrm, ~ ARMCD | AVISIT)
pairs(emmeans(fit_mmrm, ~ ARMCD | AVISIT), reverse = TRUE)


fit_mmrm <- mmrm(
  formula = FEV1 ~ ARMCD:AVISIT -1 + us(AVISIT | USUBJID),
  data = fev_data
)
coef(fit_mmrm)


C <- matrix(0,nrow = 4, ncol = length(coef(fit_mmrm)))
C[1,1] <- 1; C[1,2] <- -1
C[2,3] <- 1; C[2,4] <- -1
C[3,5] <- 1; C[3,6] <- -1
C[4,6] <- 1; C[4,8] <- -1
C%*%coef(fit_mmrm)
# C[2,2] <- 1
# C <- diag(length(coef(fit_mmrm)))

test_mmrm <- glht(fit_mmrm, linfct = C)
dim(test_mmrm$linfct)
summary(test_mmrm,test = adjusted(type = "none"))$test$pval
summary(test_mmrm,test = adjusted(type = "bonferroni"))$test$pval
summary(test_mmrm,test = adjusted(type = "free"))$test$pval



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
# set.seed(213)
sigma <- corar1(rho = 0.7, p = 3)
beta_trt <- 0.5
beta_time <- 0.2
n <- 300

res <-lapply(1:5000, function(i){
  age <- rnorm(n=n,mean = 50,sd = 3)
  trt <- rbinom(n=n,size = 1, prob = 0.5)
  Y <- mvtnorm::rmvnorm(n = n, mean = cumsum(rep(beta_time,ncol(sigma))), sigma = sigma)
  Y[,1] <- Y[,1]+trt*0.2 + 0.2*trt*(age-mean(age))
  Y[,2] <- Y[,2]+trt*0.4 + 0.4*trt*(age-mean(age))
  Y[,3] <- Y[,3]+trt*0.6 + 0.6*trt*(age-mean(age))
  dd_long <- dd_long_amp <- data.frame(Y,USUBJID = as.factor(1:n)
                        , "TRT" = as.factor(trt) 
                        , "AGE" = age
                        )
  dd_long <- dd_long %>%
    pivot_longer(cols = -c(USUBJID,TRT,AGE), names_to = "AVISIT", values_to = "Y")  %>%
    mutate(AVISIT = as.factor(gsub("[^0-9.-]", "", AVISIT)))
  
  mi <- rbinom(n=n,size=1,prob=plogis((age-mean(age))-2))==1
  dd_long_amp[which(mi),1:3] <- NA # ~30% missing
  dd_long_amp <- dd_long_amp %>%
    pivot_longer(cols = -c(USUBJID,TRT,AGE), names_to = "AVISIT", values_to = "Y")  %>%
    mutate(AVISIT = as.factor(gsub("[^0-9.-]", "", AVISIT)))
  
  
  fi <- mmrm(
    formula = Y ~ TRT*AVISIT + us(AVISIT | USUBJID),
    data = dd_long
  )
  
  fi_amp <- mmrm(
    formula = Y ~ TRT*AVISIT + us(AVISIT | USUBJID),
    data = dd_long_amp
  )
  
  out <- data.frame(rbind( summary(fi)$coef["TRT1",], summary(fi_amp)$coef["TRT1",]))
  out$type <- c("without", "dropping")
  out
})
res <- do.call("rbind",res)
ggplot(data = res, aes(x = Estimate ,fill = type))+
  geom_histogram(position = "identity", col = "black",alpha = 0.6)+
  geom_vline(xintercept = res %>% group_by(type) %>% summarise(mean_est = mean(Estimate)) %>% pull(mean_est),col = "blue")



ggplot(data = dd_long, aes(x = AVISIT,y = Y,col = TRT))+
  geom_boxplot()


library(mmrm)
fit_mmrm <- mmrm(
  formula = Y ~ TRT:AVISIT-1 + us(AVISIT | USUBJID),
  data = dd_long
)
coef(fit_mmrm)
summary(fit_mmrm)

C <- matrix(0,nrow = 3, ncol = length(coef(fit_mmrm)))
C[1,1] <- -1; C[1,2] <- 1
C[2,3] <- -1; C[2,4] <- 1
C[3,5] <- -1; C[3,6] <- 1
C%*%coef(fit_mmrm)

test_mmrm <- glht(fit_mmrm, linfct = C)
dim(test_mmrm$linfct)
summary(test_mmrm,test = adjusted(type = "none"))$test$pval
summary(test_mmrm,test = adjusted(type = "bonferroni"))$test$pval
summary(test_mmrm,test = adjusted(type = "free"))$test$pval


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# OS-Surv via time-series ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(survival)

km <- survfit(Surv(time, status) ~ 1, data = lung)
plot(km)
dd <- data.frame("time" = km$time, "surv" = km$surv)%>%arrange(time)
ggplot(data = dd %>% mutate(diff = c(diff(surv),NA))
         , aes(x = time,y = (diff)))+
  geom_point()+
  geom_smooth()

time_new <- base::setdiff(seq(min(dd$time),max(dd$time)), unique(dd$time))
# linear interpolation
interp <- approxfun( km$time, km$surv, method = "linear", rule = 2)
surv_new <- interp(time_new)
dd_extended <- bind_rows(dd
                         ,data.frame("time" = time_new , "surv" = surv_new)
                         ,data.frame("time" = 0 , "surv" = 1)
                         )%>%
  arrange(time)


dd_extended$diff <- c(diff(dd_extended$surv),NA)
plot(y = log(-dd_extended$diff),x = dd_extended$time,type = "o")

plot(dd_extended$diff)
plot(dd_extended$surv)
lines(cumsum(dd_extended$diff)+1,col = "red")



ggplot(data = dd_extended, aes(x = time,y = surv))+
  geom_point()+
  geom_smooth()+
  geom_hline(yintercept = 0)
ggplot(data = dd_extended, aes(x = time,y = diff))+
  geom_point()+
  geom_smooth()
ggplot(data = dd_extended, aes(x = time,y = cumsum(diff)))+
  geom_point()+
  geom_smooth()

minmax <- function(x){
  min_x <- min(x,na.rm = T)
  max_x <- max(x,na.rm = T)
  out <- (x - min_x)/(max_x - min_x)
  attr(out,"min") <- min_x
  attr(out,"max") <- max_x
  out
}
dd_extended$scaled_diff <- scaled_diff <- minmax(dd_extended$diff)

back_trans <- function(x, min, max){
  x*(max - min) + min
}

all.equal(
dd_extended$diff
  ,back_trans(scaled_diff,attributes(scaled_diff)$min,attributes(scaled_diff)$max)
  ,check.attributes = F
)

cumsum(back_trans(scaled_diff,attributes(scaled_diff)$min,attributes(scaled_diff)$max))

ggplot(data = dd_extended, aes(x = time,y = scaled_diff))+
  geom_point()+
  geom_smooth()

# scaled_diff <- scaled_diff[1:10]
lagged <- lapply(1:length(scaled_diff), function(i){
  # i = 1
  out <- rep(1, length(scaled_diff))
  out[(length(scaled_diff) - i + 1):length(scaled_diff)] <- scaled_diff[1:i]
  out
})
lagged <- data.frame(do.call("rbind",lagged))
colnames(lagged)[ncol(lagged)] <- "y"

library(ranger)
lagged$y <- factor(lagged$y )
rf <- ranger(y ~ ., data = na.omit(lagged))
summary(rf)
rf
my_num <- function(x){ as.numeric(as.character(x))}
ypred <- predict(rf,data = lagged)
plot(x = my_num(lagged$y)
     , y = my_num(ypred$predictions) );abline(a = 0,b = 1)


extend <- lagged[nrow(lagged),,drop = F] %>% mutate(y=NULL)
y_new <- unname(unlist(extend[1,,drop =F]))
steps <- 500
for(i in 1:steps){
  # i = 1
  pred <- c(predict(rf,data = extend, predict.all = TRUE)$predictions)
  pred <- rf$forest$levels[pred]
  pred <- table(pred)/sum(table(pred))
  new_pred <- as.numeric(sample(rownames(pred),1,p = pred))
  extend$y <- new_pred
  y_new <- c(y_new, new_pred)
  if(min(cumsum(back_trans(y_new,attributes(scaled_diff)$min,attributes(scaled_diff)$max))+1)< 0){
    break
  }
  extend <- extend[1,2:ncol(extend)]
  # print(pred)
  colnames(extend) <- paste0("X", 1:ncol(extend))
}
y_recover <- cumsum(back_trans(y_new,attributes(scaled_diff)$min,attributes(scaled_diff)$max))
y_recover <- y_recover + 1 - y_recover[1]
plot(y_recover,col = "red",type = "l");abline(h=0)
lines(dd_extended$surv,type = "l")







