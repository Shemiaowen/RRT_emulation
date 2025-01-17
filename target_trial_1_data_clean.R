# load data
patients <- read.csv("D:/code/patients_aki.csv")

# organize the dataframe in descending order
patients<-patients%>% arrange(desc(stay_id), hr)

# rename vars
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

# select vars
patients <- 
  patients %>%
  select(
    id,
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

# define outcome variable
patients$original_survive_days <- patients$survive_days
patients$survive_days[which(patients$survive_days>90)]=90
patients$survive_days[which(is.na(patients$survive_days))]=90
patients$mortality_90=ifelse(patients$survive_days>=90,0,1)
patients$survive_days[which(patients$survive_days>30)]=30
patients$mortality_30=ifelse(patients$survive_days>=30,0,1)

# gender
patients$gender[patients$gender == "F"] <- 1
patients$gender[patients$gender == "M"] <- 0
unique(patients$gender)

# sepsis
patients$sepsis[patients$sepsis == ""] <- NA
patients <- patients %>%
  mutate(sepsis = case_when(
    sepsis == "t" ~ 1,  
    is.na(sepsis) ~ 0   
  ))
unique(patients$sepsis)

# mechanical_vent
patients$mechanical_vent[is.na(patients$mechanical_vent)] <- 0
unique(patients$mechanical_vent)

# race
patients <- patients %>%
  mutate(race = case_when(
    race %in% c('OTHER','UNKNOWN','UNABLE TO OBTAIN') ~ 'OTHER',
    race %in% c('WHITE') ~ 'WHITE',
    race %in% c('BLACK/AFRICAN AMERICAN') ~ 'BLACK',
    race %in% c('ASIAN') ~ 'ASIAN',
    race %in% c('AMERICAN INDIAN/ALASKA NATIVE') ~ 'AMERICAN INDIAN/ALASKA NATIVE',
    race %in% c('HISPANIC/LATINO') ~ 'HISPANIC/LATINO'
  ))
unique(patients$race)

# lag all time-varying covariates since modelling will depend on past covariate values
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
patients.proof<-patients.proof%>%filter(hr>-1)
length(unique(patients.proof$id))
patients.initiated<-patients.proof%>%filter(crrt==1)
length(unique(patients.initiated$id))


patients.proof.dup <- patients.proof
patients.proof.dup$SBP<-as.integer(patients.proof$SBP)
patients.proof.dup$DBP<-as.integer(patients.proof$DBP)
patients.proof.dup$MBP<-as.integer(patients.proof$MBP)
patients.proof.dup$HR<-as.integer(patients.proof$HR)
patients.proof.dup$RR<-as.integer(patients.proof$RR)
patients.proof.dup$SpO2<-as.integer(patients.proof$SpO2)

# We'll change the flag to reset it to hr=1 when the conditions are met
patients.proof.dup <- patients.proof.dup %>%
  group_by(id) %>% 
  mutate(reset.hour = row_number())%>%
  ungroup()
