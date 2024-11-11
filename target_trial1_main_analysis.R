setwd("D:/code")

#load data
patients <- read.csv("patients_aki.csv")

#make sure the data is ordered by id and then by time
patients<-patients%>% arrange(desc(stay_id), hr)

#rename vars
patients <- patients %>%
  rename(
    id = stay_id,
    age = admission_age,
    race = ethnicity,
    HR = heartrate,
    RR = respiratoryrate,
    SBP = nibp_systolic,
    DBP = nibp_diastolic,
    MBP = nibp_mean,
    temp = temperature,
    SpO2 = spo2,
    vasopressor = vaso_used,
    sepsis = sepsis3,
    crrt = crrt_start_time,
    charlson = charlson_comorbidity_index,
    survive_days = death_days,
    los = los_icu,
    pH = do_50820,
    bicarbonate = do_50882,
    potassium = do_50971,
    lactate = do_50813,
    neutrophils = do_51256,
    hemoglobin = do_51222,
    platelet = do_51265,  
    bilirubin = do_50885,
    ALT = do_50861,
    AST = do_50878,
    albumin = do_50862,
    urea_nitrogen = do_51006,  
    creatinine = do_50912
  )

#remove unnecessary variables
patients <- 
  patients %>%
  select(id,
    hr,
    admission_type,
    sofa,
    sapsii,
    gcs,
    charlson,
    gender,
    age,
    race,
    HR,
    RR,
    SBP,
    DBP,
    MBP,
    temp,
    SpO2,
    pH,
    bicarbonate,
    potassium,
    lactate,
    neutrophils,
    hemoglobin,
    platelet,  
    bilirubin,
    urea_nitrogen,  
    creatinine,
    crrt,
    sepsis,
    mechanical_vent,
    vasopressor,
    los,
    survive_days
  )

#####define outcome variable
patients$original_survive_days <- patients$survive_days
patients$survive_days[which(patients$survive_days>90)]=90
patients$survive_days[which(is.na(patients$survive_days))]=90
patients$mortality_90=ifelse(patients$survive_days>=90,0,1)
patients$survive_days[which(patients$survive_days>30)]=30
patients$mortality_30=ifelse(patients$survive_days>=30,0,1)

#gender
patients$gender[patients$gender == "F"] <- 1
patients$gender[patients$gender == "M"] <- 0
unique(patients$gender)

#sepsis
unique(patients$sepsis)
patients$sepsis[patients$sepsis == ""] <- NA
patients <- patients %>%
  mutate(sepsis = case_when(
    sepsis == "t" ~ 1,  
    is.na(sepsis) ~ 0   
  ))

#mechanical_vent
patients$mechanical_vent[is.na(patients$mechanical_vent)] <- 0
unique(patients$mechanical_vent)

#create categories for race
unique(patients$race)
patients <- patients %>%
  mutate(race = case_when(
    race %in% c('OTHER','UNKNOWN','UNABLE TO OBTAIN') ~ 'OTHER',
    race %in% c('WHITE') ~ 'WHITE',
    race %in% c('BLACK/AFRICAN AMERICAN') ~ 'BLACK',
    race %in% c('ASIAN') ~ 'ASIAN',
    race %in% c('AMERICAN INDIAN/ALASKA NATIVE') ~ 'AMERICAN INDIAN/ALASKA NATIVE',
    race %in% c('HISPANIC/LATINO') ~ 'HISPANIC/LATINO'
  ))

#lag time-varying covariates since modelling will depend on past covariate values
patients.proof<-patients%>%group_by(id) %>%
  mutate(HR = lag(HR),
         RR = lag(RR),
         SBP = lag(SBP),
         DBP = lag(DBP),
         MBP = lag(MBP),
         temp = lag(temp),
         SpO2 = lag(SpO2),
         gcs = lag(gcs))

patients.proof<-patients.proof%>%group_by(id)%>%mutate(ever.initiation= ifelse(any(crrt == 1), 1, 0))%>%ungroup()
patients.proof<-patients.proof%>%group_by(id)%>%mutate(cumulative.crrt=cumsum(crrt))%>%ungroup()
patients.proof$RR[patients.proof$RR<10]<-NA
patients.proof <- patients.proof %>%
  group_by(id) %>%
  fill(gcs, HR, RR, SBP, DBP, MBP, temp, SpO2) %>%
  ungroup()
patients.proof<-patients.proof%>%filter(hr>-1)

initiation.proof<-patients.proof%>%filter(crrt==1)
dim(initiation.proof)  # 1101 patients initiation out of 7607 patients
min(patients.proof$hr) 
quantile(patients.proof$hr)
patients.proof.dup <- patients.proof
patients.proof.dup$SBP<-as.integer(patients.proof$SBP)
patients.proof.dup$DBP<-as.integer(patients.proof$DBP)
patients.proof.dup$MBP<-as.integer(patients.proof$MBP)
patients.proof.dup$HR<-as.integer(patients.proof$HR)
patients.proof.dup$RR<-as.integer(patients.proof$RR)
patients.proof.dup$SpO2<-as.integer(patients.proof$SpO2)

#We'll change the flag to reset it to hr=1 when the conditions are met
patients.proof.dup <- patients.proof.dup %>%
  group_by(id) %>% 
  mutate(reset.hour = row_number())%>%
  ungroup()

###Impute the data
##Let's explore some data before that
#MICE
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

missing_percentage <- sapply(imputed.final.data, function(x) mean(is.na(x)) * 100)
print(missing_percentage)

write.csv(imputed.final.data, "D:/RRT_emulation/aki/imputed.final.data.csv")

imputed.final.data <- read.csv("imputed.final.data.csv")

#let's study data before and after imputation
apply(patients.proof.dup[,c("gcs", "pH", "bicarbonate", "creatinine", "urea_nitrogen", "lactate")],2,mean,na.rm=T)
apply(patients.proof.dup[,c("gcs", "pH", "bicarbonate", "creatinine", "urea_nitrogen", "lactate")],2,sd,na.rm=T)
apply(patients.proof.dup[,c("gcs", "pH", "bicarbonate", "creatinine", "urea_nitrogen", "lactate")],2,quantile,na.rm=T)

apply(imputed.final.data[,c("gcs", "pH", "bicarbonate", "creatinine", "urea_nitrogen", "lactate")],2,mean,na.rm=T)
apply(imputed.final.data[,c("gcs", "pH", "bicarbonate", "creatinine", "urea_nitrogen", "lactate")],2,sd,na.rm=T)
apply(imputed.final.data[,c("gcs", "pH", "bicarbonate", "creatinine", "urea_nitrogen", "lactate")],2,quantile,na.rm=T)

#we'll filter out observations beyond the 72 repetitions
imputed.final.data <- imputed.final.data %>%
  group_by(id) %>%
  mutate(reset.hour = row_number()) %>%
  filter(reset.hour < 73) %>%
  ungroup()
length(unique(imputed.final.data$id))

initiation.proof<-imputed.final.data%>%filter(crrt==1)
length(unique(initiation.proof$id))
dim(initiation.proof) # 825 patients initiation out of 7607 patients
min(imputed.final.data$hr) 
quantile(imputed.final.data$hr)

#5TH STEP:IPTW
colSums(is.na(imputed.final.data))
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

#lagging to allow modelling the deltas
imputed.final.data <- imputed.final.data %>%
  arrange(id, reset.hour) %>% #Make sure the data is ordered by id and then by time 
  group_by(id) %>% #Group the data by the 'id' variable 
  mutate(
    prev_HR = lag(HR, 1, order_by = reset.hour),
    prev_RR = lag(RR, 1, order_by = reset.hour),
    prev_SBP = lag(SBP, 1, order_by = reset.hour),
    prev_DBP = lag(DBP, 1, order_by = reset.hour),
    prev_MBP = lag(MBP, 1, order_by = reset.hour),
    prev_SpO2 = lag(SpO2, 1, order_by = reset.hour),
    prev_gcs = lag(gcs, 1, order_by = reset.hour) #'1' is the lag interval; you can change it to create lagged variables over different intervals
  ) %>%
  ungroup()

imputed.final.data <- imputed.final.data %>%
  arrange(id, reset.hour) %>%  
  group_by(id) %>%              
  fill(                                
    prev_HR, prev_RR, prev_SBP, prev_DBP, prev_MBP, prev_SpO2, prev_gcs,
    .direction = "up"  
  ) %>%
  ungroup()
colSums(is.na(imputed.final.data))
set.seed(123)
imputed.final.data<-as.data.frame(imputed.final.data)

#calculate weights
#for hour 1
imputed.final.data.hour1<-imputed.final.data%>%filter(reset.hour==1)
ps_model1<-glm(crrt~age+gender+gcs+sofa+charlson+HR+RR+SBP+DBP+MBP,family="binomial",data=imputed.final.data.hour1)
probability.of.initiation<-sum(imputed.final.data.hour1$crrt==1)/dim(imputed.final.data.hour1)[1]
imputed.final.data.hour1$ps<-ps_model1$fitted.values
imputed.final.data.hour1$weights<-ifelse(imputed.final.data.hour1$crrt==1,probability.of.initiation/imputed.final.data.hour1$ps,(1-probability.of.initiation)/(1-imputed.final.data.hour1$ps))

#for the remaining hours
imputed.final.data.hours<-imputed.final.data%>%filter(reset.hour>1)
ps_model2<-glm(crrt~reset.hour+age+gender+gcs+sofa+charlson+HR+RR+SBP+DBP+MBP+prev_HR+prev_RR+prev_SBP+prev_DBP+prev_MBP+prev_gcs,family="binomial",data=imputed.final.data.hours)
imputed.final.data.hours$ps<-ps_model2$fitted.values
marginal_probs <- imputed.final.data.hours %>%
  group_by(reset.hour) %>%
  summarize(marginal_treatment_prob = mean(as.numeric(ever.initiation)-1),
            marginal_no_treatment_prob = 1 - marginal_treatment_prob) %>%
  ungroup()

#Join the marginal probabilities back to the original dataframe
imputed.final.data.hours <- imputed.final.data.hours %>%
  left_join(marginal_probs%>% select(marginal_treatment_prob,marginal_no_treatment_prob,reset.hour), by = "reset.hour")

#Then calculate the stabilized weights
imputed.final.data.hours <- imputed.final.data.hours %>%
  mutate(
    weights = ifelse(ever.initiation == 1, 
                     marginal_treatment_prob / ps,
                     marginal_no_treatment_prob / (1 - ps))
  )
#merge the cohorts
imputed.final.data.hours<-imputed.final.data.hours%>%select(-c("marginal_treatment_prob","marginal_no_treatment_prob"))
imputed.final.data<-rbind(imputed.final.data.hour1,imputed.final.data.hours)

weight.plot<-ggplot(imputed.final.data,aes(x=as.factor(reset.hour),y=weights))+geom_boxplot()
weight.plot

#truncation
#Perform some weight truncation at 99%
truncate<-imputed.final.data
quantile(imputed.final.data$weights)
truncate$weights[truncate$weights>quantile(truncate$weights,0.99)]<-quantile(truncate$weights,0.99)
quantile(truncate$weights)

weight.plot.truncated<-ggplot(truncate,aes(x=as.factor(reset.hour),y=weights))+geom_boxplot()
weight.plot.truncated

truncate<-truncate%>% arrange(desc(id), hr)
truncate<-truncate%>%group_by(id)%>%mutate(initiation= ifelse(any(crrt == 1), 1, 0))%>%ungroup()
initiation.proof<-truncate%>%filter(initiation==1)
length(unique(initiation.proof$id))

##MAIN ANALYSIS
#SURVIVAL OBJECT
set.seed(123)
surv_object.truncated <- with(truncate, Surv(original_survive_days, mortality_90 == 1))

#Fit a weighted Cox proportional hazards model
cox_model.truncated <- coxph(surv_object.truncated ~ initiation+gcs+sofa+charlson+HR+RR+SBP+DBP+MBP+pH+bicarbonate+potassium+lactate+hemoglobin+platelet+urea_nitrogen+creatinine,weights=weights,cluster=id,data = truncate,robust=T)
summary(cox_model.truncated)

#We fit a stratified object to use it for a survival curve
cox_model.truncated.strata <- coxph(surv_object.truncated ~ strata(initiation)+gcs+sofa+charlson+HR+RR+SBP+DBP+MBP+pH+bicarbonate+potassium+lactate+hemoglobin+platelet+urea_nitrogen+creatinine,cluster=id,data = truncate,weights = weights,robust=T)
summary(cox_model.truncated.strata)

#Proportional Hazards (PH) test based on Schoenfeld Residuals
ph_test <- cox.zph(cox_model.truncated)
print(ph_test)
ggcoxzph(ph_test, 
         var = c("initiation", "gcs", "sofa", "SBP","pH","lactate"), 
         font.main = 8,
         ggtheme=theme(axis.title.y = element_text(size=8)))

#survival curve
surv_fit_truncated <- survfit(cox_model.truncated.strata, data = truncate)
end_time <- 90
summary(surv_fit_truncated ,times=end_time)

#Plot
cox.plot.truncated <- ggsurvplot(
  surv_fit_truncated, 
  data = truncate,
  pval = "Hazard ratio: 0.706, 95% CI 0.56-0.89, p= 0.003",
  censor = FALSE,
  risk.table = FALSE,
  legend.labs = c("Initiation = No", "Initiation = Yes"),
  title = "Survival Curves by Initiation Group",
  xlab = "Time",
  ylab = "Survival Probability"
)

cox.plot.truncated$plot+geom_vline(aes(xintercept = 30), linetype = "dashed", color = "blue") +scale_x_continuous(breaks = c(0,30,90),expand = c(0,0))

# 30d
set.seed(123)
surv_object.truncated.30 <- with(truncate, Surv(survive_days, mortality_30 == 1))
cox_model.truncated <- coxph(surv_object.truncated.30 ~ initiation+gcs+sofa+charlson+HR+RR+SBP+DBP+MBP+pH+bicarbonate+potassium+lactate+hemoglobin+platelet+urea_nitrogen+creatinine,weights=weights,cluster=id,data = truncate,robust=T)
summary(cox_model.truncated)
cox_model.truncated.strata.30 <- coxph(surv_object.truncated.30 ~ strata(initiation)+gcs+sofa+charlson+HR+RR+SBP+DBP+MBP+pH+bicarbonate+potassium+lactate+hemoglobin+platelet+urea_nitrogen+creatinine,cluster=id,data = truncate, weights = weights,robust=T)
summary(cox_model.truncated.strata.30)
surv_fit_truncated.30d <- survfit(cox_model.truncated.strata.30, data = truncate)

cox.plot.truncated.30d <- ggsurvplot(
  surv_fit_truncated.30d, 
  data = truncate,
  pval = "Hazard ratio: 0.71, 95% CI 0.56-0.90, p= 0.005",
  censor = FALSE,
  risk.table = FALSE,
  legend.labs = c("Initiation = No", "Initiation = Yes"),
  title = "Survival Curves by Initiation Group",
  xlab = "Time (days)",
  ylab = "Survival Probability"
)

print(cox.plot.truncated.30d)

# AFT model
truncate$original_survive_days[truncate$original_survive_days == 0] <- 1e-10
surv_object <- Surv(truncate$original_survive_days, truncate$mortality_90)
aft_model <- survreg(surv_object ~ initiation+gcs+sofa+charlson+HR+RR+SBP+DBP+MBP+pH+bicarbonate+potassium+
                       lactate+hemoglobin+platelet+urea_nitrogen+creatinine,data = truncate, dist = "weibull")
aft_model <- survreg(surv_object ~ initiation,data = truncate, dist = "weibull")

summary(aft_model)
coef(aft_model)
exp(coef(aft_model))
exp(confint(aft_model))[-1,]

# Survival Function Plot
d<- seq(0, 90, length.out=10000)
h<-1-pweibull(d,shape=1/2.49,scale=exp(8.184))
df<-data.frame(t=d,s=h)
ggplot(df, aes(x=t, y=s)) + 
  geom_line(colour="green") + 
  ggtitle("Survival Function s(t)") +
  xlab("Time(day)") + 
  ylab("Survival Probability")

######AIPW#######
mortality.lr<-geeglm(as.numeric(mortality_90) ~ initiation + gcs + sofa + charlson + HR + RR + SBP + DBP + MBP + pH + bicarbonate + potassium + lactate + hemoglobin + platelet + urea_nitrogen + creatinine,id=id,data = truncate,family = binomial,corstr = "exchangeable")
summary(mortality.lr)
exp(mortality.lr$coefficients) #OR
truncate$predicted_outcome <- predict(mortality.lr, type = "response")
cov<-c("age","gcs","sofa","charlson","HR","RR","SBP","pH","bicarbonate","lactate", "urea_nitrogen")
cov2<-c("age","gcs","sofa","charlson","HR","RR","SBP")

aipw.mortality.lr <- AIPW$new(
  Y = as.numeric(truncate$mortality_90) , 
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
  Y = as.numeric(truncate$mortality_30), 
  A = as.numeric(truncate$initiation),
  W.Q = subset(truncate, select = cov),
  W.g = subset(truncate, select = cov2),
  Q.SL.library = "SL.glm",
  g.SL.library = "SL.glm",
  verbose = TRUE,
  k_split = 10
)
aipw.mortality.lr.30$fit()
suppressWarnings({
  aipw.mortality.lr.30$stratified_fit()$summary()
})