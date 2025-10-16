# Multiplicity Adjustment for Multiple Imputation 

Covariance-Based Multiplicity Adjustment for Multiply Imputed Data

Pooling Variance-Covariance matrizes for Holm-based multiplicity correction (with multcomp package).

Motivation: e.g WASH-trial (https://tlverse.org/tmlcimx2021-workshop/data.html) has 6 intervention groups. But missing values need to be introduced.




library(multcomp)
amod <- aov(breaks ~ tension, data = warpbreaks)
test <- glht(amod, linfct = mcp(tension = "Tukey"))

test$coef # can be adapted (take the pooled estimates)
test$vcov # can be adapted (take the pooled estimates)
test$focus

summary(glht(amod, linfct = mcp(tension = "Tukey")),test = adjusted(type = "bonferroni"))
summary(glht(amod, linfct = mcp(tension = "Tukey")),test = adjusted(type = "holm"))
summary(glht(amod, linfct = mcp(tension = "Tukey")),test = adjusted(type = "none"))
summary(glht(amod, linfct = mcp(tension = "Tukey")),test = adjusted(type = "free"))


Does Bonferroni-Holm use the covariance? Most likely not.