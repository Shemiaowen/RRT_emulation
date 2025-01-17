# Impute the data
# Let's explore some data before that
colSums(is.na(patients.proof.dup))

# MICE
set.seed(123)
target_imp<-subset(patients.proof.dup[,c("gcs", "HR", "RR", "SBP", "DBP", "MBP", "temp", "SpO2",
                                         "pH", "bicarbonate", "potassium", "lactate", "neutrophils", "hemoglobin", "platelet", 
                                         "bilirubin", "urea_nitrogen", "creatinine")])
imp_data <- mice(target_imp,m=5)
Ind.Na11.imput = apply(is.na(target_imp),2,which)
median.imputed.data = target_imp
for(i in 1:length(median.imputed.data)){
  aa2.mimic = matrix(unlist(imp_data$imp[i]),ncol=5,byrow = F)
  median.imputed.data[unlist(Ind.Na11.imput[i]),i] = apply(aa2.mimic,1,median)
}
target.unimputed.to.merge<-patients.proof.dup%>%select(-c("gcs", "HR", "RR", "SBP", "DBP", "MBP", "temp", "SpO2",
                                                          "pH", "bicarbonate", "potassium", "lactate", "neutrophils", "hemoglobin", "platelet", 
                                                          "bilirubin", "urea_nitrogen", "creatinine"))
imputed.final.data<-cbind(target.unimputed.to.merge,median.imputed.data)
colSums(is.na(imputed.final.data))

# #let's study data before and after imputation
apply(patients.proof.dup[,c("gcs", "pH", "bicarbonate", "creatinine", "urea_nitrogen", "lactate")],2,mean,na.rm=T)
apply(patients.proof.dup[,c("gcs", "pH", "bicarbonate", "creatinine", "urea_nitrogen", "lactate")],2,sd,na.rm=T)
apply(patients.proof.dup[,c("gcs", "pH", "bicarbonate", "creatinine", "urea_nitrogen", "lactate")],2,quantile,na.rm=T)

apply(imputed.final.data[,c("gcs", "pH", "bicarbonate", "creatinine", "urea_nitrogen", "lactate")],2,mean,na.rm=T)
apply(imputed.final.data[,c("gcs", "pH", "bicarbonate", "creatinine", "urea_nitrogen", "lactate")],2,sd,na.rm=T)
apply(imputed.final.data[,c("gcs", "pH", "bicarbonate", "creatinine", "urea_nitrogen", "lactate")],2,quantile,na.rm=T)

# we'll filter out observations beyond the 72 repetitions
imputed.final.data <- imputed.final.data %>%
  group_by(id) %>%
  mutate(reset.hour = row_number()) %>%
  filter(reset.hour < 73) %>%
  ungroup()
length(unique(imputed.final.data$id))

initiation.proof<-imputed.final.data%>%filter(crrt==1)
length(unique(initiation.proof$id))
dim(initiation.proof) # 7607 patients initiation out of 825 patients
min(imputed.final.data$hr) 
quantile(imputed.final.data$hr)

# 5TH STEP:IPTW
names(imputed.final.data)
imputed.final.data$gender <- as.factor(imputed.final.data$gender)
imputed.final.data$race <- as.factor(imputed.final.data$race)
imputed.final.data$crrt <- as.factor(imputed.final.data$crrt)
imputed.final.data$ever.initiation <- as.factor(imputed.final.data$ever.initiation)
imputed.final.data$sepsis <- as.factor(imputed.final.data$sepsis)
imputed.final.data$mechanical_vent <- as.factor(imputed.final.data$mechanical_vent)
imputed.final.data$vasopressor <- as.factor(imputed.final.data$vasopressor)
imputed.final.data$mortality_90 <- as.factor(imputed.final.data$mortality_90)
imputed.final.data$mortality_30 <- as.factor(imputed.final.data$mortality_30)
table(imputed.final.data$reset.hour,imputed.final.data$crrt)

cumulative_crrt <- imputed.final.data %>%
  group_by(reset.hour) %>%
  summarise(cumulative_count = cumsum(sum(as.numeric(crrt)))) %>%
  ungroup()

# lagging to allow modelling the deltas
imputed.final.data <- imputed.final.data %>%
  arrange(id, reset.hour) %>%  
  group_by(id) %>%  
  mutate(
    prev_HR = lag(HR, 1, order_by = reset.hour),
    prev_RR = lag(RR, 1, order_by = reset.hour),
    prev_SBP = lag(SBP, 1, order_by = reset.hour),
    prev_DBP = lag(DBP, 1, order_by = reset.hour),
    prev_MBP = lag(MBP, 1, order_by = reset.hour),
    prev_SpO2 = lag(SpO2, 1, order_by = reset.hour),
    prev_temp = lag(temp, 1, order_by = reset.hour),
    prev_gcs = lag(gcs, 1, order_by = reset.hour) #'1' is the lag interval; you can change it to create lagged variables over different intervals
  ) %>%
  ungroup()

imputed.final.data <- imputed.final.data %>%
  arrange(id, reset.hour) %>%  
  group_by(id) %>%              
  fill(                                
    prev_HR, prev_RR, prev_SBP, prev_DBP, prev_MBP,prev_SpO2, prev_temp, prev_gcs,
    .direction = "up"  
  ) %>%
  ungroup()

set.seed(123)
imputed.final.data<-as.data.frame(imputed.final.data)

# calculate weights
# for hour 1
imputed.final.data.hour1<-imputed.final.data%>%filter(reset.hour==1)
ps_model1<-glm(crrt~age+gender+sofa+charlson+HR+RR+DBP+SpO2+sepsis+vasopressor,family="binomial",data=imputed.final.data.hour1)
probability.of.initiation<-sum(imputed.final.data.hour1$crrt==1)/dim(imputed.final.data.hour1)[1]
imputed.final.data.hour1$ps<-ps_model1$fitted.values
imputed.final.data.hour1$weights<-ifelse(imputed.final.data.hour1$crrt==1,probability.of.initiation/imputed.final.data.hour1$ps,(1-probability.of.initiation)/(1-imputed.final.data.hour1$ps))

# for the remaining hours
imputed.final.data.hours<-imputed.final.data%>%filter(reset.hour>1)
ps_model2<-glm(crrt~reset.hour+age+gender+sofa+charlson+HR+RR+DBP+SpO2+sepsis+vasopressor+prev_HR+prev_RR+prev_DBP+prev_SpO2,family="binomial",data=imputed.final.data.hours)
imputed.final.data.hours$ps<-ps_model2$fitted.values
marginal_probs <- imputed.final.data.hours %>%
  group_by(reset.hour) %>%
  summarize(marginal_treatment_prob = mean(as.numeric(ever.initiation)-1),
            marginal_no_treatment_prob = 1 - marginal_treatment_prob) %>%
  ungroup()

# Join the marginal probabilities back to the original dataframe
imputed.final.data.hours <- imputed.final.data.hours %>%
  left_join(marginal_probs%>% select(marginal_treatment_prob,marginal_no_treatment_prob,reset.hour), by = "reset.hour")

# calculate the stabilized weights
imputed.final.data.hours <- imputed.final.data.hours %>%
  mutate(
    weights = ifelse(ever.initiation == 1, 
                     marginal_treatment_prob / ps,
                     marginal_no_treatment_prob / (1 - ps))
  )

# merge the cohorts
imputed.final.data.hours<-imputed.final.data.hours%>%select(-c("marginal_treatment_prob","marginal_no_treatment_prob"))
imputed.final.data<-rbind(imputed.final.data.hour1,imputed.final.data.hours)

# weights accross time
weight.plot<-ggplot(imputed.final.data,aes(x=as.factor(reset.hour),y=weights))+geom_boxplot()
weight.plot

# truncation
# perform some weight truncation at 99%
truncate<-imputed.final.data
quantile(imputed.final.data$weights)
truncate$weights[truncate$weights>quantile(truncate$weights,0.99)]<-quantile(truncate$weights,0.99)
quantile(truncate$weights)

# weights accross time
weight.plot.truncated<-ggplot(truncate,aes(x=as.factor(reset.hour),y=weights))+geom_boxplot()
weight.plot.truncated

truncate<-truncate%>% arrange(desc(id), hr)
truncate<-truncate%>%group_by(id)%>%mutate(initiation= ifelse(any(crrt == 1), 1, 0))%>%ungroup()
initiation.proof<-truncate%>%filter(initiation==1)
length(unique(initiation.proof$id))

# 90d main analysis
# survival object
set.seed(123)
surv_object.truncated <- with(truncate, Surv(original_survive_days, mortality_90 == 1))

# fit a weighted Cox proportional hazards model
cox_model.truncated <- coxph(surv_object.truncated ~ initiation+SBP+MBP+temp+pH+bicarbonate+lactate+hemoglobin+urea_nitrogen+creatinine,weights=weights,cluster=id,data = truncate,robust=T)
summary(cox_model.truncated)

# fit a stratified object to use it for a survival curve
cox_model.truncated.strata <- coxph(surv_object.truncated ~ strata(initiation)+SBP+MBP+temp+pH+bicarbonate+lactate+hemoglobin+urea_nitrogen+creatinine,cluster=id,data = truncate,weights = weights,robust=T)
summary(cox_model.truncated.strata)

# proportional hazards (PH) test based on Schoenfeld Residuals
ph_test <- cox.zph(cox_model.truncated)
print(ph_test)
ggcoxzph(ph_test, 
         var = c("initiation"), 
         font.main = 8,
         ggtheme=theme(axis.title.y = element_text(size=8)))

# survival curves
surv_fit_truncated <- survfit(cox_model.truncated.strata, data = truncate)
end_time <- 90
summary(surv_fit_truncated ,times=end_time)

cox.plot.truncated <- ggsurvplot(
  surv_fit_truncated, 
  data = truncate,
  pval = "Hazard ratio: 0.653, 95% CI 0.512-0.834, p < 0.001",
  censor = FALSE,
  risk.table = FALSE,
  legend.labs = c("Initiation = delay", "Initiation = early"),
  title = "90-day Survival Curves by Initiation Group",
  xlab = "Time (days)",
  ylab = "Survival Probability"
)
cox.plot.truncated$plot+geom_vline(aes(xintercept = 30), linetype = "dashed", color = "blue") +scale_x_continuous(breaks = c(0,30,90),expand = c(0,0))

# 30d
# survival object
set.seed(123)
surv_object.truncated.30 <- with(truncate, Surv(survive_days, mortality_30 == 1))

# fit a weighted Cox proportional hazards model
cox_model.truncated <- coxph(surv_object.truncated.30 ~ initiation+SBP+MBP+temp+pH+bicarbonate+lactate+hemoglobin+urea_nitrogen+creatinine,weights=weights,cluster=id,data = truncate,robust=T)
summary(cox_model.truncated)

# fit a stratified object to use it for a survival curve
cox_model.truncated.strata.30 <- coxph(surv_object.truncated.30 ~ strata(initiation)+SBP+MBP+temp+pH+bicarbonate+lactate+hemoglobin+urea_nitrogen+creatinine,cluster=id,data = truncate, weights = weights,robust=T)
summary(cox_model.truncated.strata.30)

# plot a survival curve
surv_fit_truncated.30d <- survfit(cox_model.truncated.strata.30, data = truncate)
cox.plot.truncated.30d <- ggsurvplot(
  surv_fit_truncated.30d, 
  data = truncate,
  pval = "Hazard ratio: 0.649, 95% CI 0.504-0.835, p < 0.001",
  censor = FALSE,
  risk.table = FALSE,
  legend.labs = c("Initiation = delay", "Initiation = early"),
  title = "30-day Survival Curves by Initiation Group",
  xlab = "Time (days)",
  ylab = "Survival Probability"
)
print(cox.plot.truncated.30d)

# AFT
# 90d
surv_object <- Surv(truncate$original_survive_days, truncate$mortality_90)
aft_model <- survreg(surv_object ~ initiation, data = truncate, dist = "weibull")
summary(aft_model)
coef(aft_model)
exp(coef(aft_model))
exp(confint(aft_model))[-1,]

# 30d
surv_object <- Surv(truncate$original_survive_days, truncate$mortality_30)
aft_model <- survreg(surv_object ~ initiation,data = truncate,dist = "weibull")
summary(aft_model)
coef(aft_model)
exp(coef(aft_model))
exp(confint(aft_model))[-1,]


# AIPW
cov<-c("age","gender","sofa","charlson","HR","RR","DBP","SpO2","sepsis","vasopressor","SBP","MBP","pH","bicarbonate","lactate","hemoglobin", "urea_nitrogen","creatinine")
cov2<-c("age","gender","sofa","charlson","HR","RR","DBP","SpO2","sepsis","vasopressor")

# 90d
aipw.mortality.lr <- AIPW$new(
  Y = as.numeric(truncate$mortality_90)-1, 
  A = as.numeric(truncate$initiation) ,
  W.Q = subset(truncate, select = cov),
  W.g = subset(truncate, select = cov2),
  Q.SL.library = "SL.glm",
  g.SL.library = "SL.glm",
  verbose = TRUE,
  k_split = 10
)
aipw.mortality.lr$fit()
aipw.mortality.lr$plot.p_score()
suppressWarnings({
  aipw.mortality.lr$stratified_fit()$summary()
})

# 30d
aipw.mortality.lr.30 <- AIPW$new(
  Y = as.numeric(truncate$mortality_30)-1, 
  A = as.numeric(truncate$initiation),
  W.Q = subset(truncate, select = cov),
  W.g = subset(truncate, select = cov2),
  Q.SL.library = "SL.glm",
  g.SL.library = "SL.glm",
  verbose = TRUE,
  k_split = 10
)
aipw.mortality.lr.30$fit()
aipw.mortality.lr$plot.p_score()
suppressWarnings({
  aipw.mortality.lr.30$stratified_fit()$summary()
})


