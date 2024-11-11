##### hour < 9#####
prova<-truncate%>%filter(reset.hour<9)
table(prova$crrt)
prova <- prova %>% 
  mutate(
    mortality.30d.censored = ifelse(original_survive_days <= 30, 1, 0),
    time.fu = ifelse(original_survive_days <= 30, original_survive_days, 30)
  )

### 90-day mortality
surv_object.truncated.prova <- with(prova, Surv(original_survive_days, mortality_90 == 1))

#Fit a weighted Cox proportional hazards model
set.seed(123)
cox_model.truncated.prova <- coxph(surv_object.truncated.prova ~ crrt,weights=weights,cluster=id,data = prova,robust=T)
summary(cox_model.truncated.prova)

### 30-day mortality
surv_object.truncated.prova.30 <- with(prova, Surv(time.fu, mortality.30d.censored == 1))

# Fit a weighted Cox proportional hazards model
cox_model.truncated.prova.30 <- coxph(surv_object.truncated.prova.30 ~ crrt,weights=weights,cluster=id,data = prova,robust=T)
summary(cox_model.truncated.prova.30)
#######



##### hour 9 to 24#######
prova<-truncate%>%filter(reset.hour>8&reset.hour<25)
table(prova$crrt)
prova <- prova %>% 
  mutate(
    mortality.30d.censored = ifelse(original_survive_days <= 30, 1, 0),
    time.fu = ifelse(original_survive_days <= 30, original_survive_days, 30)
  )

### 90-day mortality
set.seed(123)
surv_object.truncated.prova <- with(prova, Surv(original_survive_days, mortality_90 == 1))

#Fit a weighted Cox proportional hazards model
cox_model.truncated.prova <- coxph(surv_object.truncated.prova ~ crrt,weights=weights,cluster=id,data = prova,robust=T)
summary(cox_model.truncated.prova)

#### 30-day mortality
surv_object.truncated.prova.30 <- with(prova, Surv(time.fu, mortality.30d.censored == 1))

# Fit a weighted Cox proportional hazards model
cox_model.truncated.prova.30 <- coxph(surv_object.truncated.prova.30 ~ crrt,weights=weights,cluster=id,data = prova,robust=T)
summary(cox_model.truncated.prova.30)
###########



#####hour 25 to 48#######
prova<-truncate%>%filter(reset.hour>24&reset.hour<49)
table(prova$crrt)
prova <- prova %>% 
  mutate(
    mortality.30d.censored = ifelse(original_survive_days <= 30, 1, 0),
    time.fu = ifelse(original_survive_days <= 30, original_survive_days, 30)
  )
### 90-day mortality
set.seed(123)
surv_object.truncated.prova <- with(prova, Surv(original_survive_days, mortality_90 == 1))

#Fit a weighted Cox proportional hazards model
cox_model.truncated.prova <- coxph(surv_object.truncated.prova ~ crrt,weights=weights,cluster=id,data = prova,robust=T)
summary(cox_model.truncated.prova)

#### 30-day mortality
surv_object.truncated.prova.30 <- with(prova, Surv(time.fu, mortality.30d.censored == 1))

# Fit a weighted Cox proportional hazards model
cox_model.truncated.prova.30 <- coxph(surv_object.truncated.prova.30 ~ crrt,weights=weights,cluster=id,data = prova,robust=T)
summary(cox_model.truncated.prova.30)
###########



#####hour 49 to 72#######
prova<-truncate%>%filter(reset.hour>48&reset.hour<73)
table(prova$crrt)
prova <- prova %>% 
  mutate(
    mortality.30d.censored = ifelse(original_survive_days <= 30, 1, 0),
    time.fu = ifelse(original_survive_days <= 30, original_survive_days, 30)
  )
### 90-day mortality
set.seed(123)
surv_object.truncated.prova <- with(prova, Surv(original_survive_days, mortality_90 == 1))

#Fit a weighted Cox proportional hazards model
cox_model.truncated.prova <- coxph(surv_object.truncated.prova ~ crrt,weights=weights,cluster=id,data = prova,robust=T)
summary(cox_model.truncated.prova)

### 30-day mortality
surv_object.truncated.prova.30 <- with(prova, Surv(time.fu, mortality.30d.censored == 1))

#Fit a weighted Cox proportional hazards model
cox_model.truncated.prova.30 <- coxph(surv_object.truncated.prova.30 ~ crrt,weights=weights,cluster=id,data = prova,robust=T)
summary(cox_model.truncated.prova.30)
###########

