## SOFA subgroup analysis
data <- truncate[truncate$reset.hour == 1, ]
ROC<-survivalROC(Stime=data$original_survive_days,  
                 status=data$mortality_90,          
                 marker=data$sofa,                  
                 predict.time=90,                   
                 method="KM")                       
cutoff<-ROC$cut.values[which.max(ROC$TP-ROC$FP)]
cutoff #The optimal cutoff value for SOFA score is 8

# Plotting the ROC Curve.
plot(ROC$FP,    
     ROC$TP,    
     type="l",  
     xlim=c(0,1),
     ylim=c(0,1),
     xlab=paste("FP","\n","AUC=",round(ROC$AUC,3)),
     ylab="TP",  
     main="90-day survival ROC",col="red")+
  abline(0,1)+
  legend("bottomright",c("sofa"),col="red",lty=c(1,1))

# Categorizing SOFA Scores
data$new_sofa<-ifelse(data$sofa>8,"high","low")
value_distribution <- table(data$new_sofa)
print(value_distribution)

truncate <- truncate %>%
  group_by(id) %>%
  mutate(new_sofa = ifelse(any(new_sofa == "high"), "high", new_sofa)) %>%
  ungroup() 

# Survival Analysis
set.seed(123)
surv_object <- with(data, Surv(original_survive_days, mortality_90 == 1))
cutoff_cox<-coxph(surv_object ~ new_sofa,data=data,weights = weights,robust=T)
summary(cutoff_cox)

cutoff_cox.strata <- coxph(surv_object ~ strata(new_sofa),data = data,weights = weights,robust=T)
summary(cutoff_cox.strata)

cox_model <- survfit(cutoff_cox.strata, data = data)
end_time <- 90
summary(cox_model,times=end_time)

ggsurvplot(cox_model, 
           data = data,
           pval = "Hazard ratio: 0.244, 95% CI 0.22-0.27, p < 0.001",
           censor = FALSE,
           risk.table = FALSE,
           title = "Survival Curves by sofa Group",
           xlab = "Time",
           ylab = "Survival Probability")



## lactate subgroup analysis
ROC<-survivalROC(Stime=data$original_survive_days,  
                 status=data$mortality_90,          
                 marker=data$lactate,                  
                 predict.time=90,                   
                 method="KM")                       
cutoff<-ROC$cut.values[which.max(ROC$TP-ROC$FP)]
cutoff #The optimal cutoff value for lactate score is 2.1

plot(ROC$FP,    
     ROC$TP,    
     type="l",  
     xlim=c(0,1),
     ylim=c(0,1),
     xlab=paste("FP","\n","AUC=",round(ROC$AUC,3)),
     ylab="TP",  
     main="90-day survival ROC",col="red")+
  abline(0,1)+
  legend("bottomright",c("lactate"),col="red",lty=c(1,1))

data$new_lactate<-ifelse(data$lactate>2.1,"high","low")
value_distribution <- table(data$new_lactate)
print(value_distribution)

set.seed(123)
surv_object <- with(data, Surv(original_survive_days, mortality_90 == 1))
cutoff_cox<-coxph(surv_object ~ new_lactate,data=data,weights = weights,robust=T)
summary(cutoff_cox)

cutoff_cox.strata <- coxph(surv_object ~ strata(new_lactate),data = data,weights = weights,robust=T)
summary(cutoff_cox.strata)

cox_model <- survfit(cutoff_cox.strata, data = data)
end_time <- 90
summary(cox_model,times=end_time)

ggsurvplot(cox_model, 
           data = data,
           pval = "Hazard ratio: 0.89, 95% CI 0.80-0.98, p = 0.021",
           censor = FALSE,
           risk.table = FALSE,
           title = "Survival Curves by lactate Group",
           xlab = "Time",
           ylab = "Survival Probability")



## urea nitrogen subgroup analysis
ROC<-survivalROC(Stime=data$original_survive_days,  
                 status=data$mortality_90,          
                 marker=data$urea_nitrogen,                  
                 predict.time=90,                   
                 method="KM")                       
cutoff<-ROC$cut.values[which.max(ROC$TP-ROC$FP)]
cutoff #The optimal cutoff value for urea nitrogen score is 30

plot(ROC$FP,    
     ROC$TP,    
     type="l",  
     xlim=c(0,1),
     ylim=c(0,1),
     xlab=paste("FP","\n","AUC=",round(ROC$AUC,3)),
     ylab="TP",  
     main="90-day survival ROC",col="red")+
  abline(0,1)+
  legend("bottomright",c("urea_nitrogen"),col="red",lty=c(1,1))

data$new_urea_nitrogen<-ifelse(data$urea_nitrogen>30,"high","low")
value_distribution <- table(data$new_urea_nitrogen)
print(value_distribution)

surv_object <- with(data, Surv(original_survive_days, mortality_90 == 1))
cutoff_cox<-coxph(surv_object ~ new_urea_nitrogen,data=data,weights = weights,robust=T)
summary(cutoff_cox)

cutoff_cox.strata <- coxph(surv_object ~ strata(new_urea_nitrogen),data = data,weights = weights,robust=T)
summary(cutoff_cox.strata)

cox_model <- survfit(cutoff_cox.strata, data = data)
end_time <- 90
summary(cox_model,times=end_time)

ggsurvplot(cox_model, 
           data = data,
           pval = "Hazard ratio: 0.897, 95% CI 0.808-0.996, p = 0.042",
           censor = FALSE,
           risk.table = FALSE,
           title = "Survival Curves by urea nitrogen Group",
           xlab = "Time",
           ylab = "Survival Probability")



## CCI subgroup analysis
ROC<-survivalROC(Stime=data$original_survive_days,  
                 status=data$mortality_90,          
                 marker=data$charlson,                  
                 predict.time=90,                   
                 method="KM")                       
cutoff<-ROC$cut.values[which.max(ROC$TP-ROC$FP)]
cutoff #The optimal cutoff value for charlson score is 8

plot(ROC$FP,    
     ROC$TP,    
     type="l",  
     xlim=c(0,1),
     ylim=c(0,1),
     xlab=paste("FP","\n","AUC=",round(ROC$AUC,3)),
     ylab="TP",  
     main="90-day survival ROC",col="red")+
  abline(0,1)+
  legend("bottomright",c("charlson"),col="red",lty=c(1,1))

data$new_charlson<-ifelse(data$charlson>8,"high","low")
value_distribution <- table(data$new_charlson)
print(value_distribution)

surv_object <- with(data, Surv(original_survive_days, mortality_90 == 1))
cutoff_cox<-coxph(surv_object ~ new_charlson,data=data,weights = weights,robust=T)
summary(cutoff_cox)

cutoff_cox.strata <- coxph(surv_object ~ strata(new_charlson),data = data,weights = weights,robust=T)
summary(cutoff_cox.strata)

cox_model <- survfit(cutoff_cox.strata, data = data)
end_time <- 90
summary(cox_model,times=end_time)

ggsurvplot(cox_model, 
           data = data,
           pval = "Hazard ratio: 0.613, 95% CI 0.555-0.677, p < 0.001",
           censor = FALSE,
           risk.table = FALSE,
           title = "Survival Curves by charlson Group",
           xlab = "Time",
           ylab = "Survival Probability")

