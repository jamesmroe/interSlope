rm(list=ls())
version = "v2"

# Q - How do multiple variables (early and late life)
# associate with variance in intercept versus change of brain and cognition across the lifespan? 


#===cohort inputs=========================================================================

#--- load FreeSurfer data with variables linked
# DF = fread("path/to/FSdata.csv")
#---

#--- var of interest
VAR = "GCA" # name of general cog ability var
long_only = TRUE # use multi timepoint factor with multi timepoint brain
saveres = FALSE
#---

cohort = "LCBC"

# set necessary inputs to varnames in cohort data
varObsID = "uid" # unique obs ID var
varSubID = "subject_id" # subID var
varBrainL = "left_hippocampus" # hippocampus var L
varBrainR = "right_hippocampus" # hippocampus var R
varSex = "subject_sex" # sex var (coded "Female" / "Male")
varAge = "visit_age" # age var
varEdu = "edu_years" # edu var
varICV = "estimatedtotalintracranialvol" # ICV var

# optional inputs
correct_for_scanner = TRUE
varScanner = "mri_site_name" # scanner var

#===cohort inputs=========================================================================


#---load packages
loadPackages = function() {
  packages = c("dplyr","tidyverse", "magrittr","here","metagam","mgcv","gamm4", "tictoc", "data.table")
  new.packages = packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) {
    install.packages(new.packages)
  }
  sapply(packages, require, character.only = T)
}
loadPackages()
setwd("/Users/jamesroe/LCBC/Users/jamesroe/interSlope/interSlope")
here()


#---testing
# runAll = 1; if (runAll == 1) {
cohort = "LCBC"
# cohort = "BACS"
if (cohort == "LCBC") {
  load(here("data/LCBC/DF_LCBC_noas.Rda"))
  load(here("data/LCBC/test_LCBC.Rda"))
  varObsID = "mri_info_folder"
  varSubID = "subject_id"
  varBrainL = "mri_aseg_volume_left_hippocampus"
  varBrainR = "mri_aseg_volume_right_hippocampus"
  varSex = "subject_sex"
  varAge = "visit_age"
  varEdu = "edu_years"
  varICV = "mri_aseg_volume_estimatedtotalintracranialvol"
  varScanner = "mri_site_name"
} else if (cohort == "BACS") {
  VAR = "em_cfa"
  load(here("data/BACS/test_BACS.Rda"))
  DF$sex[DF$sex == "M"] = "Male"
  varObsID = "imageLink"
  varSubID = "BACSID"
  varBrainL = "vol.Left.Hippocampus"
  varBrainR = "vol.Right.Hippocampus"
  varSex = "sex"
  varAge = "age_at_visit"
  varEdu = "age_at_visit"
  varICV = "vol.EstimatedTotalIntraCranialVol"
  varScanner = "scanner_var"
}


analysis = "X-GCA_Y-Hippo"


# harmonize vars -----
DF %<>% rename(uid = varObsID,
               id = varSubID)
DF$id = as.character(DF$id)
DF$brainvar = ( DF[[varBrainL]] + DF[[varBrainR]] ) / 2
DF$sex = factor(DF[[varSex]])
DF$age = DF[[varAge]]
DF$edu = DF[[varEdu]]
DF$scanner = factor(DF[[varScanner]])
DF$icv = ( DF[[varICV]] - mean(DF[[varICV]], na.rm = T) ) / sd(DF[[varICV]], na.rm = T)


#--- sets var of interest
DF$VAR = DF[[VAR]]
#--- sets var of interest


# derive time vars
if (long_only == T) {
  #use multi timepoint factor with multi timepoint brain
  DF %<>% arrange(id, age) %>% filter(
                                        !is.na(VAR) & !is.na(brainvar)) %>% group_by(id) %>% mutate(
                                          age_bsl = min(age),
                                          age_mean = mean(age),
                                          time_bsl = age - min(age),
                                          nTimepoints = length(unique(age)),
                                          nTimepointsVAR = length(unique(age[!is.na(VAR)])),
                                          nTime = row_number(),
                                          intervalFirstLast = max(age) - min(age)
                                        ) %>% ungroup() %>% filter(nTimepoints >= 2)
  
  #standardize VAR by first timepoint
  mean_tp1 = mean(DF$VAR[DF$nTime == 1])
  sd_tp1 = sd(DF$VAR[DF$nTime == 1])
  DF$VAR_Z = (DF$VAR - mean_tp1) /sd_tp1
  
} else {
  #use one timepoint factor with multi timepoint brain
  DF %<>% arrange(id, age) %>% filter(
                                        !is.na(brainvar)) %>% group_by(id) %>% mutate(
                                          age_bsl = min(age),
                                          age_mean = mean(age),
                                          time_bsl = age - min(age),
                                          nTimepoints = length(unique(age)),
                                          nTimepointsVAR = length(unique(age[!is.na(VAR)])),
                                          nTime = row_number(),
                                          intervalFirstLast = max(age) - min(age)
                                        ) %>% ungroup() %>% filter(nTimepoints >= 2) %>% group_by(id) %>% 
    #remove subjects with no observations on VAR
    mutate(anyVAR = as.integer(any(!is.na(VAR)))) %>% ungroup() %>% filter(anyVAR == 1) %>% select(-anyVAR)
  
  #standardize VAR
  mean_all = mean(DF$VAR, na.rm = TRUE)
  sd_all = sd(DF$VAR, na.rm = TRUE)
  DF$VAR_Z = (DF$VAR - mean_all) /sd_all
}
# DF %>% filter(nTime == 1) %>% select(id) %>% unlist() %>% unname() %>% unique() %>% length()
# length(unique(DF$id))


# check derived vars
core_cols = c("uid", "id", "brainvar", "sex", "age", "edu")
derived_cols = c("age_bsl", "age_mean", "time_bsl", "nTimepoints", "nTime", "intervalFirstLast")
head(DF) %>% select(all_of(c(core_cols, derived_cols)))



# split VAR to categorical factors corrected for age and sex
VAR_Z_mean =
  DF %>% 
  filter(!is.na(VAR_Z)) %>% 
  group_by(id) %>% 
  summarize(VAR_Z_mean = mean(VAR_Z, na.rm = TRUE))


# (U)nique subject baselines
U = DF %>% filter(nTime == 1)


if (!"VAR_Z_mean" %in% names(U)) {
  U = left_join(U, VAR_Z_mean)
}
sum(is.na(U$VAR_Z_mean))


# remove if no observations on VAR 
U %<>% filter(!is.na(VAR_Z_mean))


#split VAR to categorical factors corrected for age and sex
gmod1 = gamm4::gamm4(VAR_Z_mean ~
                       s(age_mean) + sex, data = U)


addPred = function(dat, gam_mod, modstring) {
  # add GAMM predictions to data
  predictions = dat %>%
    mutate(sex = "Female") %>%
    select(age, sex, age_mean) %>%
    predict(gam_mod$gam, newdata = ., se.fit = T)
  dat$residuals = residuals(gam_mod$mer)
  dat$partial_residuals = predictions$fit + dat$residuals
  dat$fit = predictions$fit
  dat$sefit = predictions$se.fit
  dat$cifit = predictions$se.fit * 1.96
  names(dat)[names(dat) == "residuals"] = paste(modstring, "residuals", sep = "_")
  names(dat)[names(dat) == "partial_residuals"] = paste(modstring, "partial_residuals", sep = "_")
  names(dat)[names(dat) == "fit"] = paste(modstring, "fit", sep = "_")
  names(dat)[names(dat) == "sefit"] = paste(modstring, "sefit", sep = "_")
  names(dat)[names(dat) == "cifit"] = paste(modstring, "cifit", sep = "_")
  return(list(data = dat, model = gam_mod))
}


U = addPred(U, gmod1, "gmod1")[[1]]
# plot(U$age_mean, U$gmod1_partial_residuals)
# points(U$age_mean, U$gmod1_fit,col="blue")
# plot(U$age_mean, U$VAR_Z_mean)


# smooth function of age for categorical split -----
U %<>% mutate(VAR_Z_mean_meansplit = ifelse(VAR_Z_mean > U$gmod1_fit, 1, 0))
U %<>% mutate(VAR_Z_mean_tertile = ntile(U$gmod1_residuals, 3))
U %<>% mutate(VAR_Z_mean_quartile = ntile(U$gmod1_residuals, 4))


# check splits
ggplot() +
  geom_point(data = U, aes(x = age_mean,VAR_Z_mean, color = factor(VAR_Z_mean_meansplit)), stat = "identity", alpha = 1, size = 0.5) +
  geom_line(data = U, aes(x = age_mean, y = gmod1_fit), col = "black") +
  theme_minimal()
ggplot() +
  geom_point(data = U, aes(x = age_mean, VAR_Z_mean, color = factor(VAR_Z_mean_quartile)), stat = "identity", alpha = 1, size = 0.5) +
  geom_line(data = U, aes(x = age_mean, y = gmod1_fit), col = "black") +
  theme_minimal()
ggplot() +
  geom_point(data = U, aes(x = age_mean, VAR_Z_mean, color = factor(VAR_Z_mean_tertile)), stat="identity", alpha = 1, size = 0.5) +
  geom_smooth(data = U, aes(x = age_mean, VAR_Z_mean, color = factor(VAR_Z_mean_tertile))) +
  theme_minimal()


#add in categorical splits to DF
sum(U$id %in% DF$id) == nrow(U); sum(DF$id %in% U$id) == nrow(DF)
if (! "VAR_Z_mean_quartile" %in% names(DF)) {
  
  link = U %>% select(id, VAR_Z_mean_meansplit, contains("VAR_Z_mean")) %>% distinct()
  for (i in 1:length(link$id)) {
    if (i == 1) {
      DF$VAR_Z_mean_meansplit = NA
      DF$VAR_Z_mean_tertile = NA
      DF$VAR_Z_mean_quartile = NA
    }
    sub = link$id[i]
    DF$VAR_Z_mean_meansplit[DF$id == sub] = link$VAR_Z_mean_meansplit[link$id == sub]
    DF$VAR_Z_mean_tertile[DF$id == sub] = link$VAR_Z_mean_tertile[link$id == sub]
    DF$VAR_Z_mean_quartile[DF$id == sub] = link$VAR_Z_mean_quartile[link$id == sub]
  }
  
  # all should be there
  if (sum(is.na(DF$VAR_Z_mean_meansplit))) { stop("not all linked. Check this.")}
  if (sum(is.na(DF$VAR_Z_mean_quartile))) { stop("not all linked. Check this.")}
  
  DF %<>% mutate(
    fVAR_Z_mean_meansplit = factor(VAR_Z_mean_meansplit),
    fVAR_Z_mean_tertile = factor(DF$VAR_Z_mean_tertile),
    fVAR_Z_mean_quartile = factor(DF$VAR_Z_mean_quartile),
    oVAR_Z_mean_meansplit = factor(DF$VAR_Z_mean_meansplit, levels = c(0, 1), ordered = T),
    oVAR_Z_mean_tertile = factor(DF$VAR_Z_mean_tertile, levels = c(1, 2, 3), ordered = T),
    oVAR_Z_mean_quartile = factor(DF$VAR_Z_mean_quartile, levels = c(1, 2, 3, 4), ordered = T))
}


if (correct_for_scanner == TRUE) {
  covariates = c("sex" , "scanner")
  covariatesICV = c("sex" , "scanner", "icv")
} else {
  covariates = c("sex")
  covariatesICV = c("sex", "icv")
}



# main categorical models -------
# (including one timepoint VAR observations)
gammf_mod2_mean = as.formula(paste("brainvar ~ s(age, bs = \"cr\", by = fVAR_Z_mean_meansplit) + fVAR_Z_mean_meansplit","+",paste(covariates, collapse = " + ")))
gammf_mod3_tert = as.formula(paste("brainvar ~ s(age, bs = \"cr\", by = fVAR_Z_mean_tertile) + fVAR_Z_mean_tertile","+",paste(covariates, collapse = " + ")))
gammf_mod4_quart = as.formula(paste("brainvar ~ s(age, bs = \"cr\", by = fVAR_Z_mean_quartile) + fVAR_Z_mean_quartile","+",paste(covariates, collapse = " + ")))
gammf_mod2_omean = as.formula(paste("brainvar ~ fVAR_Z_mean_meansplit + s(age, bs = \"cr\") + s(age, bs = \"cr\", by = oVAR_Z_mean_meansplit)","+",paste(covariates, collapse = " + ")))
gammf_mod3_otert = as.formula(paste("brainvar ~ fVAR_Z_mean_tertile + s(age, bs = \"cr\") + s(age, bs = \"cr\", by = oVAR_Z_mean_tertile)","+",paste(covariates, collapse = " + ")))
gammf_mod4_oquart = as.formula(paste("brainvar ~ fVAR_Z_mean_quartile + s(age, bs = \"cr\") + s(age, bs = \"cr\", by = oVAR_Z_mean_quartile)","+",paste(covariates, collapse = " + ")))

gammf_mod2_mean_icv = as.formula(paste("brainvar ~ s(age, bs = \"cr\", by = fVAR_Z_mean_meansplit) + fVAR_Z_mean_meansplit","+",paste(covariatesICV, collapse = " + ")))
gammf_mod3_tert_icv = as.formula(paste("brainvar ~ s(age, bs = \"cr\", by = fVAR_Z_mean_tertile) + fVAR_Z_mean_tertile","+",paste(covariatesICV, collapse = " + ")))
gammf_mod4_quart_icv = as.formula(paste("brainvar ~ s(age, bs = \"cr\", by = fVAR_Z_mean_quartile) + fVAR_Z_mean_quartile","+",paste(covariatesICV, collapse = " + ")))
gammf_mod2_omean_icv = as.formula(paste("brainvar ~ fVAR_Z_mean_meansplit + s(age, bs = \"cr\") + s(age, bs = \"cr\", by = oVAR_Z_mean_meansplit)","+",paste(covariatesICV, collapse = " + ")))
gammf_mod3_otert_icv = as.formula(paste("brainvar ~ fVAR_Z_mean_tertile + s(age, bs = \"cr\") + s(age, bs = \"cr\", by = oVAR_Z_mean_tertile)","+",paste(covariatesICV, collapse = " + ")))
gammf_mod4_oquart_icv = as.formula(paste("brainvar ~ fVAR_Z_mean_quartile + s(age, bs = \"cr\") + s(age, bs = \"cr\", by = oVAR_Z_mean_quartile)","+",paste(covariatesICV, collapse = " + ")))


simple_mod = gamm4::gamm4(
  brainvar ~ s(age, bs = "cr"),# + sex + scanner,
  # random = list(id = ~1),
  random = ~ (1 |id),
  data = DF
)

gamm_mod2_mean = gamm4::gamm4(
  gammf_mod2_mean,
  random = ~ (1 |id),
  # random = list(id = ~1),
  data = DF
)
# plot.gam(gamm_mod2_mean$gam) #plots age across factor0 then factor1
gamm_mod2_omean = gamm4::gamm4(
  gammf_mod2_omean,
  random = ~ (1 |id),
  data = DF
)
# plot.gam(gamm_mod2_omean$gam) #plots average age across factors then difference


# summary of main categorical models ------
summary(gamm_mod2_mean$gam)
summary(gamm_mod2_omean$gam)


gamm_mod2_mean_icv = gamm4::gamm4(
  gammf_mod2_mean_icv,
  random = ~ (1 |id),
  # random = list(id = ~1),
  data = DF
)
gamm_mod2_omean_icv = gamm4::gamm4(
  gammf_mod2_omean_icv,
  random = ~ (1 |id),
  data = DF
)


# summary of main categorical models corrected for ICV ------
summary(gamm_mod2_mean_icv$gam)
summary(gamm_mod2_omean_icv$gam)


# alt models (just in case) -------
gamm_mod3_tert = gamm4::gamm4(
  gammf_mod3_tert,
  random = ~ (1 |id),
  data = DF
)
gamm_mod3_otert = gamm4::gamm4(
  gammf_mod3_otert,
  random = ~ (1 |id),
  data = DF
)
gamm_mod3_tert_icv = gamm4::gamm4(
  gammf_mod3_tert_icv,
  random = ~ (1 |id),
  data = DF
)
gamm_mod3_otert_icv = gamm4::gamm4(
  gammf_mod3_otert_icv,
  random = ~ (1 |id),
  data = DF
)
# gamm_mod4_quart = gamm4::gamm4(gammf_mod4_quart, random = ~ (1 |id), data = DF)
# gamm_mod4_oquart = gamm4::gamm4(gammf_mod4_oquart, random = ~ (1 |id), data = DF)


# summary of alt categorical models ------
summary(gamm_mod3_tert$gam)
summary(gamm_mod3_otert$gam)
summary(gamm_mod3_tert_icv$gam)
summary(gamm_mod3_otert_icv$gam)
# summary(gamm_mod4_quart$gam); summary(gamm_mod4_oquart$gam)


# main longitudinal model -------
# (including only multi timepoint VAR observations)
DF_long = DF %>% filter(!is.na(VAR_Z)) %>% group_by(id) %>% mutate(
  age_bsl = min(age),
  age_mean = mean(age),
  time_bsl = age - min(age),
  nTimepoints = length(unique(age)),
  nTime = row_number(),
  intervalFirstLast = max(age) - min(age)) %>% ungroup() %>% filter(nTimepoints >= 2)


#NEWSECTION ----
gammf_mod_long = as.formula(paste("brainvar ~ VAR_Z * time_bsl + s(age_bsl, bs = \"cr\") + ",paste(covariates, collapse = " + ")))
lmmf_mod_long = as.formula(paste("brainvar ~ VAR_Z * time_bsl + age_bsl + ",paste(covariates, collapse = " + ")))
gammf_mod_groupmean = as.formula(paste("brainvar ~ VAR_Z_mean_meansplit * time_bsl + s(age_bsl, bs = \"cr\") + ",paste(covariates, collapse = " + ")))
lmmf_mod_groupmean = as.formula(paste("brainvar ~ VAR_Z_mean_meansplit * time_bsl + age_bsl + ",paste(covariates, collapse = " + ")))

gammf_mod_long_icv = as.formula(paste("brainvar ~ VAR_Z * time_bsl + s(age_bsl, bs = \"cr\") + ",paste(covariatesICV, collapse = " + ")))
lmmf_mod_long_icv = as.formula(paste("brainvar ~ VAR_Z * time_bsl + age_bsl + ",paste(covariatesICV, collapse = " + ")))
gammf_mod_groupmean_icv = as.formula(paste("brainvar ~ VAR_Z_mean_meansplit * time_bsl + s(age_bsl, bs = \"cr\") + ",paste(covariatesICV, collapse = " + ")))
lmmf_mod_groupmean_icv = as.formula(paste("brainvar ~ VAR_Z_mean_meansplit * time_bsl + age_bsl + ",paste(covariatesICV, collapse = " + ")))
 

gamm_mod_long = gamm4::gamm4(
  gammf_mod_long,
  random = ~ (1 |id), #gamm4
  REML = T, #gamm4
  # random = list(id = ~1), #mgcv
  data = DF_long
)


lmm_mod_long = gamm4::gamm4(
  lmmf_mod_long,
  random = ~ (1 |id),
  REML = T,
  data = DF_long
)


gamm_mod_long_icv = gamm4::gamm4(
  gammf_mod_long_icv,
  random = ~ (1 |id),
  REML = T,
  data = DF_long
)


lmm_mod_long_icv = gamm4::gamm4(
  lmmf_mod_long_icv,
  random = ~ (1 |id),
  REML = T,
  data = DF_long
)


# alt models (just in case) -------
gamm_mod_groupmean = gamm4::gamm4(
  gammf_mod_groupmean,
  random = ~ (1 |id),
  REML = T,
  data = DF
)

lmm_mod_groupmean = gamm4::gamm4(
  lmmf_mod_groupmean,
  random = ~ (1 |id),
  REML = T,
  data = DF
)

gamm_mod_groupmean_icv = gamm4::gamm4(
  gammf_mod_groupmean_icv,
  random = ~ (1 |id),
  REML = T,
  data = DF
)

lmm_mod_groupmean_icv = gamm4::gamm4(
  lmmf_mod_groupmean_icv,
  random = ~ (1 |id),
  REML = T,
  data = DF
)


# NEWSECTION - precorrect for scanner ----
if (correct_for_scanner == TRUE) {

  # accounting for scan distribution across age/sex
  DF$scanner = droplevels(DF$scanner)
  most_scans = names(table(DF$scanner)[which.max(as.data.frame(table(DF$scanner))$Freq)])
  scanner_levels = levels(DF$scanner)
  
  # change scanner baseline coding to common var
  DF$scanner = as.character(DF$scanner)
  DF$scanner[DF$scanner == most_scans] = "base" 
  DF$scanner = factor(DF$scanner)
  DF$scanner = relevel(DF$scanner, ref = "base")
  scanner_correction = gamm4::gamm4(
    brainvar ~ s(age, bs = "cr") + sex + scanner,
    random = ~ (1 |id),
    data = DF
  )
  
  
  fixed_effects = summary(scanner_correction$gam)$p.table
  scanner_names = sub("scanner", "", rownames(fixed_effects)[grepl("scanner", rownames(fixed_effects))])
  scanner_effects = fixed_effects[grepl("scanner", rownames(fixed_effects)), "Estimate"]
  for (n in 1:length(scanner_names)) {
    if (n == 1) {
      DF$brainvarcor = DF$brainvar
    }
    scanner = scanner_names[n]
    DF$brainvarcor = ifelse(DF$scanner == scanner,
                             DF$brainvarcor - scanner_effects[n],
                             DF$brainvarcor)
  }

} else {
  # models then duplicated (easier for saving)
  DF$brainvarcor = DF$brainvar
  scanner_levels = NA
}


# models with brain precorrected for scanner ----
gammf_mod2_mean_cor = as.formula(paste("brainvarcor ~ s(age, bs = \"cr\", by = fVAR_Z_mean_meansplit) + fVAR_Z_mean_meansplit + sex"))
gammf_mod3_tert_cor = as.formula(paste("brainvarcor ~ s(age, bs = \"cr\", by = fVAR_Z_mean_tertile) + fVAR_Z_mean_tertile + sex"))
gammf_mod4_quart_cor = as.formula(paste("brainvarcor ~ s(age, bs = \"cr\", by = fVAR_Z_mean_quartile) + fVAR_Z_mean_quartile + sex"))
gammf_mod2_omean_cor = as.formula(paste("brainvarcor ~ fVAR_Z_mean_meansplit + s(age, bs = \"cr\") + s(age, bs = \"cr\", by = oVAR_Z_mean_meansplit) + sex"))
gammf_mod3_otert_cor = as.formula(paste("brainvarcor ~ fVAR_Z_mean_tertile + s(age, bs = \"cr\") + s(age, bs = \"cr\", by = oVAR_Z_mean_tertile) + sex"))
gammf_mod4_oquart_cor = as.formula(paste("brainvarcor ~ fVAR_Z_mean_quartile + s(age, bs = \"cr\") + s(age, bs = \"cr\", by = oVAR_Z_mean_quartile) + sex"))

gammf_mod2_mean_icv_cor = as.formula(paste("brainvarcor ~ s(age, bs = \"cr\", by = fVAR_Z_mean_meansplit) + fVAR_Z_mean_meansplit + sex + icv"))
gammf_mod3_tert_icv_cor = as.formula(paste("brainvarcor ~ s(age, bs = \"cr\", by = fVAR_Z_mean_tertile) + fVAR_Z_mean_tertile + sex + icv"))
gammf_mod4_quart_icv_cor = as.formula(paste("brainvarcor ~ s(age, bs = \"cr\", by = fVAR_Z_mean_quartile) + fVAR_Z_mean_quartile + sex + icv"))
gammf_mod2_omean_icv_cor = as.formula(paste("brainvarcor ~ fVAR_Z_mean_meansplit + s(age, bs = \"cr\") + s(age, bs = \"cr\", by = oVAR_Z_mean_meansplit) + sex + icv"))
gammf_mod3_otert_icv_cor = as.formula(paste("brainvarcor ~ fVAR_Z_mean_tertile + s(age, bs = \"cr\") + s(age, bs = \"cr\", by = oVAR_Z_mean_tertile) + sex + icv"))
gammf_mod4_oquart_icv_cor = as.formula(paste("brainvarcor ~ fVAR_Z_mean_quartile + s(age, bs = \"cr\") + s(age, bs = \"cr\", by = oVAR_Z_mean_quartile) + sex + icv"))

gammf_mod_long_cor = as.formula(paste("brainvarcor ~ VAR_Z * time_bsl + s(age_bsl, bs = \"cr\") + sex"))
lmmf_mod_long_cor = as.formula(paste("brainvarcor ~ VAR_Z * time_bsl + age_bsl + sex"))
gammf_mod_groupmean_cor = as.formula(paste("brainvarcor ~ VAR_Z_mean_meansplit * time_bsl + s(age_bsl, bs = \"cr\") + sex"))
lmmf_mod_groupmean_cor = as.formula(paste("brainvarcor ~ VAR_Z_mean_meansplit * time_bsl + age_bsl + sex"))

gammf_mod_long_icv_cor = as.formula(paste("brainvarcor ~ VAR_Z * time_bsl + s(age_bsl, bs = \"cr\") + sex + icv"))
lmmf_mod_long_icv_cor = as.formula(paste("brainvarcor ~ VAR_Z * time_bsl + age_bsl + sex + icv"))
gammf_mod_groupmean_icv_cor = as.formula(paste("brainvarcor ~ VAR_Z_mean_meansplit * time_bsl + s(age_bsl, bs = \"cr\") + sex + icv"))
lmmf_mod_groupmean_icv_cor = as.formula(paste("brainvarcor ~ VAR_Z_mean_meansplit * time_bsl + age_bsl + sex + icv"))


DF_long = DF %>% filter(!is.na(VAR_Z)) %>% group_by(id) %>% mutate(
  age_bsl = min(age),
  age_mean = mean(age),
  time_bsl = age - min(age),
  nTimepoints = length(unique(age)),
  nTime = row_number(),
  intervalFirstLast = max(age) - min(age)) %>% ungroup() %>% filter(nTimepoints >= 2)


simple_mod_cor = gamm4::gamm4(
  brainvarcor ~ s(age, bs = "cr"),
  # random = list(id = ~1),
  random = ~ (1 |id),
  data = DF
)

gamm_mod2_mean_cor = gamm4::gamm4(
  gammf_mod2_mean_cor,
  random = ~ (1 |id),
  # random = list(id = ~1),
  data = DF
)
# plot.gam(gamm_mod2_mean_cor$gam) #plots age across factor0 then factor1
gamm_mod2_omean_cor = gamm4::gamm4(
  gammf_mod2_omean_cor,
  random = ~ (1 |id),
  data = DF
)
# plot.gam(gamm_mod2_omean_cor$gam) #plots average age across factors then difference


summary(gamm_mod2_mean_cor$gam)
summary(gamm_mod2_omean_cor$gam)


gamm_mod2_mean_icv_cor = gamm4::gamm4(
  gammf_mod2_mean_icv_cor,
  random = ~ (1 |id),
  # random = list(id = ~1),
  data = DF
)
gamm_mod2_omean_icv_cor = gamm4::gamm4(
  gammf_mod2_omean_icv_cor,
  random = ~ (1 |id),
  data = DF
)


summary(gamm_mod2_mean_icv_cor$gam)
summary(gamm_mod2_omean_icv_cor$gam)


gamm_mod3_tert_cor = gamm4::gamm4(
  gammf_mod3_tert_cor,
  random = ~ (1 |id),
  data = DF
)
gamm_mod3_otert_cor = gamm4::gamm4(
  gammf_mod3_otert_cor,
  random = ~ (1 |id),
  data = DF
)
gamm_mod3_tert_icv_cor = gamm4::gamm4(
  gammf_mod3_tert_icv_cor,
  random = ~ (1 |id),
  data = DF
)
gamm_mod3_otert_icv_cor = gamm4::gamm4(
  gammf_mod3_otert_icv_cor,
  random = ~ (1 |id),
  data = DF
)
# gamm_mod4_quart_cor = gamm4::gamm4(gammf_mod4_quart_cor, random = ~ (1 |id), data = DF)
# gamm_mod4_oquart_cor = gamm4::gamm4(gammf_mod4_oquart_cor, random = ~ (1 |id), data = DF)


summary(gamm_mod3_tert_cor$gam)
summary(gamm_mod3_otert_cor$gam)
summary(gamm_mod3_tert_icv_cor$gam)
summary(gamm_mod3_otert_icv_cor$gam)
# summary(gamm_mod4_quart_cor$gam); summary(gamm_mod4_oquart_cor$gam)


gamm_mod_long_cor = gamm4::gamm4(
  gammf_mod_long_cor,
  random = ~ (1 |id), #gamm4
  REML = T, #gamm4
  # random = list(id = ~1), #mgcv
  data = DF_long
)


lmm_mod_long_cor = gamm4::gamm4(
  lmmf_mod_long_cor,
  random = ~ (1 |id),
  REML = T,
  data = DF_long
)


gamm_mod_long_icv_cor = gamm4::gamm4(
  gammf_mod_long_icv_cor,
  random = ~ (1 |id),
  REML = T,
  data = DF_long
)


lmm_mod_long_icv_cor = gamm4::gamm4(
  lmmf_mod_long_icv_cor,
  random = ~ (1 |id),
  REML = T,
  data = DF_long
)


gamm_mod_groupmean_cor = gamm4::gamm4(
  gammf_mod_groupmean_cor,
  random = ~ (1 |id),
  REML = T,
  data = DF
)

lmm_mod_groupmean_cor = gamm4::gamm4(
  lmmf_mod_groupmean_cor,
  random = ~ (1 |id),
  REML = T,
  data = DF
)

gamm_mod_groupmean_icv_cor = gamm4::gamm4(
  gammf_mod_groupmean_icv_cor,
  random = ~ (1 |id),
  REML = T,
  data = DF
)

lmm_mod_groupmean_icv_cor = gamm4::gamm4(
  lmmf_mod_groupmean_icv_cor,
  random = ~ (1 |id),
  REML = T,
  data = DF
)


# describe sample
descVar = function(dat, var, rm_NA = F) {
  if (rm_NA == F) {
    desc = data.frame("M" = mean(dat[[var]]), "SD" = sd(dat[[var]]), "min" = min(dat[[var]]), "max" = max(dat[[var]]))
  } else if (rm_NA == T) {
    desc = data.frame("M" = mean(dat[[var]], na.rm = T), 
             "SD" = sd(dat[[var]], na.rm = T), 
             "min" = min(dat[[var]], na.rm = T),
             "max" = max(dat[[var]], na.rm = T))
  }
  desc$nObs = sum(!is.na(dat[[var]]))
  desc$var = var
  return(desc)
}


(descU = rbind(
  descVar(U, "age_bsl"), #baseline age
  descVar(U, "age_mean"), #mean age
  descVar(U, "intervalFirstLast"), #scan interval
  # descVar(U, "edu", rm_NA = T), #edu
  descVar(U, "VAR_Z_mean"), #VAR at baseline
  descVar(U, "nTimepoints"), #n timepoints
  descVar(U, "nTimepointsVAR") #n timepoints VAR
) %>% mutate(df = "U", 
             n = length(unique(U$id)),
             nF = length(unique(U$id[U$sex == "Female"])),
             nM = length(unique(U$id[U$sex == "Male"]))))


#categorical split models
(descDF = rbind(
  descVar(DF, "time_bsl"),
  descVar(DF, "age")
) %>% mutate(df = "DF",
             n = length(unique(DF$id)),
             nF = length(unique(DF$id[DF$sex == "Female"])),
             nM = length(unique(DF$id[DF$sex == "Male"])),
             obsVAR = sum(!is.na(DF$VAR_Z))))


#full longitudinal
(descDF_long = rbind(
  descVar(DF_long, "time_bsl"),
  descVar(DF_long, "age")
) %>% mutate(df = "DF_long",
             n = length(unique(DF_long$id)),
             nF = length(unique(DF_long$id[DF$sex == "Female"])),
             nM = length(unique(DF_long$id[DF$sex == "Male"]))))


nTimeY = table(U$nTimepoints)
nTimeVAR = table(U$nTimepointsVAR)


# strip data and save models for sharing -----
cohort_fitsimple = metagam::strip_rawdata(simple_mod)
cohort_fit2_mean = metagam::strip_rawdata(gamm_mod2_mean)
cohort_fit2_omean = metagam::strip_rawdata(gamm_mod2_omean)
cohort_fit3_tert = metagam::strip_rawdata(gamm_mod3_tert)
cohort_fit3_otert = metagam::strip_rawdata(gamm_mod3_otert)

cohort_fit2_mean_icv = metagam::strip_rawdata(gamm_mod2_mean_icv)
cohort_fit2_omean_icv = metagam::strip_rawdata(gamm_mod2_omean_icv)
cohort_fit3_tert_icv = metagam::strip_rawdata(gamm_mod3_tert_icv)
cohort_fit3_otert_icv = metagam::strip_rawdata(gamm_mod3_otert_icv)

cohort_gamm_mod_long = metagam::strip_rawdata(gamm_mod_long)
cohort_lmm_mod_long = metagam::strip_rawdata(lmm_mod_long)
cohort_gamm_mod_long_icv = metagam::strip_rawdata(gamm_mod_long_icv)
cohort_lmm_mod_long_icv = metagam::strip_rawdata(lmm_mod_long_icv)

cohort_gamm_mod_groupmean = metagam::strip_rawdata(gamm_mod_groupmean)
cohort_lmm_mod_groupmean = metagam::strip_rawdata(lmm_mod_groupmean)
cohort_gamm_mod_groupmean_icv = metagam::strip_rawdata(gamm_mod_groupmean_icv)
cohort_lmm_mod_groupmean_icv = metagam::strip_rawdata(lmm_mod_groupmean_icv)
#________

#NEWSECTION ----
cohort_fitsimple_cor = metagam::strip_rawdata(simple_mod_cor)
cohort_fit2_mean_cor = metagam::strip_rawdata(gamm_mod2_mean_cor)
cohort_fit2_omean_cor = metagam::strip_rawdata(gamm_mod2_omean_cor)
cohort_fit3_tert_cor = metagam::strip_rawdata(gamm_mod3_tert_cor)
cohort_fit3_otert_cor = metagam::strip_rawdata(gamm_mod3_otert_cor)

cohort_fit2_mean_icv_cor = metagam::strip_rawdata(gamm_mod2_mean_icv_cor)
cohort_fit2_omean_icv_cor = metagam::strip_rawdata(gamm_mod2_omean_icv_cor)
cohort_fit3_tert_icv_cor = metagam::strip_rawdata(gamm_mod3_tert_icv_cor)
cohort_fit3_otert_icv_cor = metagam::strip_rawdata(gamm_mod3_otert_icv_cor)

cohort_gamm_mod_long_cor = metagam::strip_rawdata(gamm_mod_long_cor)
cohort_lmm_mod_long_cor = metagam::strip_rawdata(lmm_mod_long_cor)
cohort_gamm_mod_long_icv_cor = metagam::strip_rawdata(gamm_mod_long_icv_cor)
cohort_lmm_mod_long_icv_cor = metagam::strip_rawdata(lmm_mod_long_icv_cor)

cohort_gamm_mod_groupmean_cor = metagam::strip_rawdata(gamm_mod_groupmean_cor)
cohort_lmm_mod_groupmean_cor = metagam::strip_rawdata(lmm_mod_groupmean_cor)
cohort_gamm_mod_groupmean_icv_cor = metagam::strip_rawdata(gamm_mod_groupmean_icv_cor)
cohort_lmm_mod_groupmean_icv_cor = metagam::strip_rawdata(lmm_mod_groupmean_icv_cor)


timestamp = as.character(Sys.time())
timestamp_str = gsub(" ", "_", timestamp)
timestamp_str = gsub(":", "-", timestamp_str)


version = "v2"

if (saveres == TRUE) {
  resdir = here("results", cohort)
  if (! dir.exists(resdir)) {
    dir.create(resdir, recursive = T)
  }
  save('cohort',
       'timestamp',
       
       'long_only',
       'correct_for_scanner',
       
       'descU', 
       'descDF', 
       'descDF_long', 
       'nTimeY',
       'nTimeVAR',
       
       'cohort_fitsimple',
       'cohort_fit2_mean',
       'cohort_fit2_omean',
       'cohort_fit3_tert',
       'cohort_fit3_otert',
      
       'cohort_fit2_mean_icv',
       'cohort_fit2_omean_icv',
       'cohort_fit3_tert_icv',
       'cohort_fit3_otert_icv',
       
       'cohort_gamm_mod_long',
       'cohort_lmm_mod_long',
       'cohort_gamm_mod_long_icv',
       'cohort_lmm_mod_long_icv',
       
       'cohort_gamm_mod_groupmean',
       'cohort_lmm_mod_groupmean',
       'cohort_gamm_mod_groupmean_icv',
       'cohort_lmm_mod_groupmean_icv',
       
       #NEWSECTION ----
       'cohort_fitsimple_cor',
       'cohort_fit2_mean_cor',
       'cohort_fit2_omean_cor',
       'cohort_fit3_tert_cor',
       'cohort_fit3_otert_cor',
       
       'cohort_fit2_mean_icv_cor',
       'cohort_fit2_omean_icv_cor',
       'cohort_fit3_tert_icv_cor',
       'cohort_fit3_otert_icv_cor',
       
       'cohort_gamm_mod_long_cor',
       'cohort_lmm_mod_long_cor',
       'cohort_gamm_mod_long_icv_cor',
       'cohort_lmm_mod_long_icv_cor',
       
       'cohort_gamm_mod_groupmean_cor',
       'cohort_lmm_mod_groupmean_cor',
       'cohort_gamm_mod_groupmean_icv_cor',
       'cohort_lmm_mod_groupmean_icv_cor',
       'scanner_levels',
       
       file = paste0(resdir, "/results.", 
                     cohort, ".", analysis, ".", 
                     "optlong", as.numeric(long_only), ".", 
                     "optscan", as.numeric(correct_for_scanner), 
                     ".sharesafe_", version, ".Rda")
       )
}

# }
#END #######















##########################################################################################################-
##########################################################################################################-
#################################### CODE BELOW HERE IS FOR TESTING ######################################
##########################################################################################################-
##########################################################################################################-


# splitSample = function(DF, U, seed, nSplits = 3) {
#   set.seed(seed)
#   ids = unique(DF$id)
#   n = length(ids)
#   indices <- sample(1:n)
# 
#   # Calculate the size of each part
#   if (nSplits == 3) {
#     part_size <- ceiling(n / 3)
#   } else if (nSplits == 2) {
#     part_size <- ceiling(n / 2)
#   } else {
#     stopifnot("invalid split size")
#   }
#   split1 <- indices[1:part_size]
#   split1_ids <- ids[split1]
#   split2 <- indices[(part_size + 1):(part_size * 2)]
#   split2_ids <- ids[split2]
# 
#   if (nSplits == 3) {
#     split3 <- indices[(part_size * 2 + 1):n]
#     split3_ids <- ids[split3]
#   }
# 
#   # Create the three subsets
#   DF1 <- DF[DF$id %in% split1_ids, ]
#   DF2 <- DF[DF$id %in% split2_ids, ]
#   if (nSplits == 3) {
#     DF3 <- DF[DF$id %in% split3_ids, ]
#   }
# 
#   DF1 %<>% select(-contains("mod"), contains("meansplit"), contains("quintile"))
#   DF2 %<>% select(-contains("mod"), contains("meansplit"), contains("quintile"))
#   if (nSplits == 3) {
#     DF3 %<>% select(-contains("mod"), contains("meansplit"), contains("quintile"))
# 
#     return(list(DF1 = DF1, DF2 = DF2, DF3 = DF3))
# 
#   } else if (nSplits == 2) {
# 
#     return(list(DF1 = DF1, DF2 = DF2))
#   }
# 
# }
# 
# 
# testFirst = function(DF_num) {
# 
#   # DF_num = DF1
# 
#   if (long_only == T) {
#     #use multi timepoint factor with multi timepoint brain
#     DF_num %<>% arrange(id, age) %>% filter(
#       !is.na(VAR) & !is.na(brainvar)) %>% group_by(id) %>% mutate(
#         age_bsl = min(age),
#         age_mean = mean(age),
#         time_bsl = age - min(age),
#         nTimepoints = length(unique(age)),
#         nTimepointsVAR = length(unique(age[!is.na(VAR)])),
#         nTime = row_number(),
#         intervalFirstLast = max(age) - min(age)
#       ) %>% ungroup() %>% filter(nTimepoints >= 2)
# 
#     #standardize VAR by first timepoint
#     mean_tp1 = mean(DF_num$VAR[DF_num$nTime == 1])
#     sd_tp1 = sd(DF_num$VAR[DF_num$nTime == 1])
#     DF_num$VAR_Z = (DF_num$VAR - mean_tp1) /sd_tp1
# 
#   } else {
#     #use one timepoint factor with multi timepoint brain
#     DF_num %<>% arrange(id, age) %>% filter(
#       !is.na(brainvar)) %>% group_by(id) %>% mutate(
#         age_bsl = min(age),
#         age_mean = mean(age),
#         time_bsl = age - min(age),
#         nTimepoints = length(unique(age)),
#         nTimepointsVAR = length(unique(age[!is.na(VAR)])),
#         nTime = row_number(),
#         intervalFirstLast = max(age) - min(age)
#       ) %>% ungroup() %>% filter(nTimepoints >= 2) %>% group_by(id) %>%
#       #remove subjects with no observations on VAR
#       mutate(anyVAR = as.integer(any(!is.na(VAR)))) %>% ungroup() %>% filter(anyVAR == 1) %>% select(-anyVAR)
# 
#     #standardize VAR
#     mean_all = mean(DF_num$VAR, na.rm = TRUE)
#     sd_all = sd(DF_num$VAR, na.rm = TRUE)
#     DF_num$VAR_Z = (DF_num$VAR - mean_all) /sd_all
#   }
#   # DF_num %>% filter(nTime == 1) %>% select(id) %>% unlist() %>% unname() %>% unique() %>% length()
#   # length(unique(DF_num$id))
# 
# 
#   # check derived vars
#   core_cols = c("uid", "id", "brainvar", "sex", "age", "edu")
#   derived_cols = c("age_bsl", "age_mean", "time_bsl", "nTimepoints", "nTime", "intervalFirstLast")
#   head(DF_num) %>% select(all_of(c(core_cols, derived_cols)))
# 
# 
# 
#   # split VAR to categorical factors corrected for age and sex
#   VAR_Z_mean =
#     DF_num %>%
#     filter(!is.na(VAR_Z)) %>%
#     group_by(id) %>%
#     summarize(VAR_Z_mean = mean(VAR_Z, na.rm = TRUE))
# 
# 
#   # (U)nique subject baselines
#   U = DF_num %>% filter(nTime == 1)
# 
# 
#   if (!"VAR_Z_mean" %in% names(U)) {
#     U = left_join(U, VAR_Z_mean)
#   }
#   sum(is.na(U$VAR_Z_mean))
# 
# 
#   # remove if no observations on VAR
#   U %<>% filter(!is.na(VAR_Z_mean))
# 
# 
#   #split VAR to categorical factors corrected for age and sex
#   gmod1 = gamm4::gamm4(VAR_Z_mean ~
#                          s(age_mean) + sex, data = U)
# 
# 
# 
#   U = addPred(U, gmod1, "gmod1")[[1]]
#   # plot(U$age_mean, U$gmod1_partial_residuals)
#   # points(U$age_mean, U$gmod1_fit,col="blue")
#   # plot(U$age_mean, U$VAR_Z_mean)
# 
# 
#   # smooth function of age for categorical split -----#
#   U %<>% mutate(VAR_Z_mean_meansplit = ifelse(VAR_Z_mean > U$gmod1_fit, 1, 0))
#   U %<>% mutate(VAR_Z_mean_tertile = ntile(U$gmod1_residuals, 3))
#   U %<>% mutate(VAR_Z_mean_quartile = ntile(U$gmod1_residuals, 4))
#   U %<>% mutate(VAR_Z_mean_quintile = ntile(U$gmod1_residuals, 5))
# 
# 
#   # check splits
#   ggplot() +
#     geom_point(data = U, aes(x = age_mean,VAR_Z_mean, color = factor(VAR_Z_mean_meansplit)), stat = "identity", alpha = 1, size = 0.5) +
#     geom_line(data = U, aes(x = age_mean, y = gmod1_fit), col = "black") +
#     theme_minimal()
#   ggplot() +
#     geom_point(data = U, aes(x = age_mean, VAR_Z_mean, color = factor(VAR_Z_mean_quintile)), stat = "identity", alpha = 1, size = 0.5) +
#     geom_line(data = U, aes(x = age_mean, y = gmod1_fit), col = "black") +
#     theme_minimal()
#   ggplot() +
#     geom_point(data = U, aes(x = age_mean, VAR_Z_mean, color = factor(VAR_Z_mean_tertile)), stat="identity", alpha = 1, size = 0.5) +
#     geom_smooth(data = U, aes(x = age_mean, VAR_Z_mean, color = factor(VAR_Z_mean_tertile))) +
#     theme_minimal()
# 
# 
#   #add in categorical splits to DF_num
#   sum(U$id %in% DF_num$id) == nrow(U); sum(DF_num$id %in% U$id) == nrow(DF_num)
#   if (! "VAR_Z_mean_quintile" %in% names(DF_num)) {
# 
#     link = U %>% select(id, VAR_Z_mean_meansplit, contains("VAR_Z_mean")) %>% distinct()
#     for (i in 1:length(link$id)) {
#       if (i == 1) {
#         DF_num$VAR_Z_mean_meansplit = NA
#         DF_num$VAR_Z_mean_tertile = NA
#         DF_num$VAR_Z_mean_quartile = NA
#         DF_num$VAR_Z_mean_quintile = NA
#       }
#       sub = link$id[i]
#       DF_num$VAR_Z_mean_meansplit[DF_num$id == sub] = link$VAR_Z_mean_meansplit[link$id == sub]
#       DF_num$VAR_Z_mean_tertile[DF_num$id == sub] = link$VAR_Z_mean_tertile[link$id == sub]
#       DF_num$VAR_Z_mean_quartile[DF_num$id == sub] = link$VAR_Z_mean_quartile[link$id == sub]
#       DF_num$VAR_Z_mean_quintile[DF_num$id == sub] = link$VAR_Z_mean_quintile[link$id == sub]
#     }
# 
#     # all should be there
#     if (sum(is.na(DF_num$VAR_Z_mean_meansplit))) { stop("not all linked. Check this.")}
#     if (sum(is.na(DF_num$VAR_Z_mean_quintile))) { stop("not all linked. Check this.")}
# 
#     DF_num %<>% mutate(
#       fVAR_Z_mean_meansplit = factor(VAR_Z_mean_meansplit),
#       fVAR_Z_mean_tertile = factor(DF_num$VAR_Z_mean_tertile),
#       fVAR_Z_mean_quartile = factor(DF_num$VAR_Z_mean_quartile),
#       fVAR_Z_mean_quintile = factor(DF_num$VAR_Z_mean_quintile),
#       oVAR_Z_mean_meansplit = factor(DF_num$VAR_Z_mean_meansplit, levels = c(0, 1), ordered = T),
#       oVAR_Z_mean_tertile = factor(DF_num$VAR_Z_mean_tertile, levels = c(1, 2, 3), ordered = T),
#       oVAR_Z_mean_quartile = factor(DF_num$VAR_Z_mean_quartile, levels = c(1, 2, 3, 4), ordered = T),
#       oVAR_Z_mean_quintile = factor(DF_num$VAR_Z_mean_quintile, levels = c(1, 2, 3, 4, 5), ordered = T))
#   }
# 
#   return(DF_num)
# }
# 
# testSecond = function(DF_num, iteration = 1) {
# 
#   
#   
#   #test changing scanner coding
#   if (iteration == 2) {
#     # DF_num %<>% mutate(scanner = case_when(
#     #   scanner == "ousAvanto" ~ "ascanner1",
#     #   scanner == "ousSkyra" ~ "scanner2",
#     #   scanner == "ousPrisma" ~ "scanner3",
#     #   scanner == "ntnuAvanto" ~ "scanner3",
#     # )) %>% mutate(scanner = factor(scanner))
#   }
#   if (iteration == 1) {
#     correct_for_scanner = F
#     # DF_num$scanner = as.character(DF_num$scanner)
#     # DF_num$scanner[DF_num$scanner == "ousAvanto"] = "ascanner1"
#     # DF_num$scanner = factor(DF_num$scanner)
#   }
#   if (iteration == 3) {
#     correct_for_scanner = F
#   }
#   
#   print("running models. May take a sec")
# 
#   # DF_num = DF1_num
#   if (correct_for_scanner == T) {
#     covariates = c("sex" , "scanner")
#     covariatesICV = c("sex" , "scanner", "icv")
#   } else {
#     covariates = c("sex", "icv")
#     covariatesICV = c("sex", "icv")
#   }
# 
# 
#   # main categorical models -------#
#   # (including one timepoint VAR observations)
#   gammf_mod2_mean = as.formula(paste("brainvar ~ s(age, bs = \"cr\", by = fVAR_Z_mean_meansplit) + fVAR_Z_mean_meansplit","+",paste(covariates, collapse = " + ")))
#   gammf_mod3_tert = as.formula(paste("brainvar ~ s(age, bs = \"cr\", by = fVAR_Z_mean_tertile) + fVAR_Z_mean_tertile","+",paste(covariates, collapse = " + ")))
#   gammf_mod4_quart = as.formula(paste("brainvar ~ s(age, bs = \"cr\", by = fVAR_Z_mean_quartile) + fVAR_Z_mean_quartile","+",paste(covariates, collapse = " + ")))
#   gammf_mod5_quint = as.formula(paste("brainvar ~ s(age, bs = \"cr\", by = fVAR_Z_mean_quintile) + fVAR_Z_mean_quintile","+",paste(covariates, collapse = " + ")))
#   gammf_mod2_omean = as.formula(paste("brainvar ~ fVAR_Z_mean_meansplit + s(age, bs = \"cr\") + s(age, bs = \"cr\", by = oVAR_Z_mean_meansplit)","+",paste(covariates, collapse = " + ")))
#   gammf_mod3_otert = as.formula(paste("brainvar ~ fVAR_Z_mean_tertile + s(age, bs = \"cr\") + s(age, bs = \"cr\", by = oVAR_Z_mean_tertile)","+",paste(covariates, collapse = " + ")))
#   gammf_mod4_oquart = as.formula(paste("brainvar ~ fVAR_Z_mean_quartile + s(age, bs = \"cr\") + s(age, bs = \"cr\", by = oVAR_Z_mean_quartile)","+",paste(covariates, collapse = " + ")))
#   gammf_mod5_oquint = as.formula(paste("brainvar ~ fVAR_Z_mean_quintile + s(age, bs = \"cr\") + s(age, bs = \"cr\", by = oVAR_Z_mean_quintile)","+",paste(covariates, collapse = " + ")))
# 
#   gammf_mod2_mean_icv = as.formula(paste("brainvar ~ s(age, bs = \"cr\", by = fVAR_Z_mean_meansplit) + fVAR_Z_mean_meansplit","+",paste(covariatesICV, collapse = " + ")))
#   gammf_mod3_tert_icv = as.formula(paste("brainvar ~ s(age, bs = \"cr\", by = fVAR_Z_mean_tertile) + fVAR_Z_mean_tertile","+",paste(covariatesICV, collapse = " + ")))
#   gammf_mod4_quart_icv = as.formula(paste("brainvar ~ s(age, bs = \"cr\", by = fVAR_Z_mean_quartile) + fVAR_Z_mean_quartile","+",paste(covariatesICV, collapse = " + ")))
#   gammf_mod5_quint_icv = as.formula(paste("brainvar ~ s(age, bs = \"cr\", by = fVAR_Z_mean_quintile) + fVAR_Z_mean_quintile","+",paste(covariatesICV, collapse = " + ")))
#   gammf_mod2_omean_icv = as.formula(paste("brainvar ~ fVAR_Z_mean_meansplit + s(age, bs = \"cr\") + s(age, bs = \"cr\", by = oVAR_Z_mean_meansplit)","+",paste(covariatesICV, collapse = " + ")))
#   gammf_mod3_otert_icv = as.formula(paste("brainvar ~ fVAR_Z_mean_tertile + s(age, bs = \"cr\") + s(age, bs = \"cr\", by = oVAR_Z_mean_tertile)","+",paste(covariatesICV, collapse = " + ")))
#   gammf_mod4_oquart_icv = as.formula(paste("brainvar ~ fVAR_Z_mean_quartile + s(age, bs = \"cr\") + s(age, bs = \"cr\", by = oVAR_Z_mean_quartile)","+",paste(covariatesICV, collapse = " + ")))
#   gammf_mod5_oquint_icv = as.formula(paste("brainvar ~ fVAR_Z_mean_quintile + s(age, bs = \"cr\") + s(age, bs = \"cr\", by = oVAR_Z_mean_quintile)","+",paste(covariatesICV, collapse = " + ")))
# 
# 
#   simple_mod = gamm4::gamm4(
#     brainvar ~ s(age, bs = "cr"),# + sex + scanner,
#     # random = list(id = ~1),
#     random = ~ (1 |id),
#     data = DF_num
#   )
# 
#   gamm_mod2_mean = gamm4::gamm4(
#     gammf_mod2_mean,
#     random = ~ (1 |id),
#     # random = list(id = ~1),
#     data = DF_num
#   )
#   # plot.gam(gamm_mod2_mean$gam) #plots age across factor0 then factor1
#   gamm_mod2_omean = gamm4::gamm4(
#     gammf_mod2_omean,
#     random = ~ (1 |id),
#     data = DF_num
#   )
#   # plot.gam(gamm_mod2_omean$gam) #plots average age across factors then difference
# 
# 
#   # summary of main categorical models ------#
#   summary(gamm_mod2_mean$gam)
#   summary(gamm_mod2_omean$gam)
# 
# 
#   gamm_mod2_mean_icv = gamm4::gamm4(
#     gammf_mod2_mean_icv,
#     random = ~ (1 |id),
#     # random = list(id = ~1),
#     data = DF_num
#   )
#   gamm_mod2_omean_icv = gamm4::gamm4(
#     gammf_mod2_omean_icv,
#     random = ~ (1 |id),
#     data = DF_num
#   )
# 
# 
#   # summary of main categorical models corrected for ICV ------#
#   summary(gamm_mod2_mean_icv$gam)
#   summary(gamm_mod2_omean_icv$gam)
# 
# 
#   # alt models (just in case) -------#
#   gamm_mod3_tert = gamm4::gamm4(
#     gammf_mod3_tert,
#     random = ~ (1 |id),
#     data = DF_num
#   )
#   gamm_mod3_otert = gamm4::gamm4(
#     gammf_mod3_otert,
#     random = ~ (1 |id),
#     data = DF_num
#   )
#   gamm_mod3_tert_icv = gamm4::gamm4(
#     gammf_mod3_tert_icv,
#     random = ~ (1 |id),
#     data = DF_num
#   )
#   gamm_mod3_otert_icv = gamm4::gamm4(
#     gammf_mod3_otert_icv,
#     random = ~ (1 |id),
#     data = DF_num
#   )
#   gamm_mod4_quart = gamm4::gamm4(gammf_mod4_quart, random = ~ (1 |id), data = DF_num)
#   gamm_mod4_oquart = gamm4::gamm4(gammf_mod4_oquart, random = ~ (1 |id), data = DF_num)
#   gamm_mod5_quint = gamm4::gamm4(gammf_mod5_quint, random = ~ (1 |id), data = DF_num)
#   gamm_mod5_oquint = gamm4::gamm4(gammf_mod5_oquint, random = ~ (1 |id), data = DF_num)
# 
# 
#   # summary of alt categorical models ------#
#   summary(gamm_mod3_tert$gam)
#   summary(gamm_mod3_otert$gam)
#   summary(gamm_mod3_tert_icv$gam)
#   summary(gamm_mod3_otert_icv$gam)
#   # summary(gamm_mod4_quart$gam); summary(gamm_mod4_oquart$gam); summary(gamm_mod5_quint$gam); summary(gamm_mod5_oquint$gam)
# 
# 
# 
#   # main longitudinal model -------#
#   # (including only multi timepoint VAR observations)
#   DF_long = DF_num %>% filter(!is.na(VAR_Z)) %>% group_by(id) %>% mutate(
#     age_bsl = min(age),
#     age_mean = mean(age),
#     time_bsl = age - min(age),
#     nTimepoints = length(unique(age)),
#     nTime = row_number(),
#     intervalFirstLast = max(age) - min(age)) %>% ungroup() %>% filter(nTimepoints >= 2)
# 
# 
#   gamm_mod_long = gamm4::gamm4(
#     brainvar ~ VAR_Z * time_bsl + s(age_bsl, bs = "cr") + sex + scanner,
#     random = ~ (1 |id), #gamm4
#     REML = T, #gamm4
#     # random = list(id = ~1), #mgcv
#     data = DF_long
#   )
# 
# 
#   lmm_mod_long = gamm4::gamm4(
#     brainvar ~ VAR_Z * time_bsl + age_bsl + sex + scanner,
#     random = ~ (1 |id), #gamm4
#     REML = T, #gamm4
#     data = DF_long
#   )
# 
#   
#   # NEWSECTION - precorrect for scanner ----
#   if (correct_for_scanner == TRUE) {
#     
#     # accounting for scan distribution across age/sex
#     DF_num$scanner = droplevels(DF_num$scanner)
#     most_scans = names(table(DF_num$scanner)[which.max(as.data.frame(table(DF_num$scanner))$Freq)])
#     scanner_levels = levels(DF_num$scanner)
#     
#     # change scanner baseline coding to common var
#     DF_num$scanner = as.character(DF_num$scanner)
#     DF_num$scanner[DF_num$scanner == most_scans] = "base" 
#     DF_num$scanner = factor(DF_num$scanner)
#     DF_num$scanner = relevel(DF_num$scanner, ref = "base")
#     scanner_correction = gamm4::gamm4(
#       brainvar ~ s(age, bs = "cr") + sex + scanner,
#       random = ~ (1 |id),
#       data = DF_num
#     )
#     
#     
#     fixed_effects = summary(scanner_correction$gam)$p.table
#     scanner_names = sub("scanner", "", rownames(fixed_effects)[grepl("scanner", rownames(fixed_effects))])
#     scanner_effects = fixed_effects[grepl("scanner", rownames(fixed_effects)), "Estimate"]
#     for (n in 1:length(scanner_names)) {
#       if (n == 1) {
#         DF$brainvarcor = DF$brainvar
#       }
#       scanner = scanner_names[n]
#       DF$brainvarcor = ifelse(DF$scanner == scanner,
#                               DF$brainvarcor - scanner_effects[n],
#                               DF$brainvarcor)
#     }
#     
#   } else {
#     # models then duplicated (easier for saving)
#     DF_num$brainvarcor = DF_num$brainvar
#   }
#   
#   
#   # models with brain precorrected for scanner ----
#   gammf_mod2_mean_cor = as.formula(paste("brainvarcor ~ s(age, bs = \"cr\", by = fVAR_Z_mean_meansplit) + fVAR_Z_mean_meansplit + sex"))
#   gammf_mod3_tert_cor = as.formula(paste("brainvarcor ~ s(age, bs = \"cr\", by = fVAR_Z_mean_tertile) + fVAR_Z_mean_tertile + sex"))
#   gammf_mod4_quart_cor = as.formula(paste("brainvarcor ~ s(age, bs = \"cr\", by = fVAR_Z_mean_quartile) + fVAR_Z_mean_quartile + sex"))
#   gammf_mod2_omean_cor = as.formula(paste("brainvarcor ~ fVAR_Z_mean_meansplit + s(age, bs = \"cr\") + s(age, bs = \"cr\", by = oVAR_Z_mean_meansplit) + sex"))
#   gammf_mod3_otert_cor = as.formula(paste("brainvarcor ~ fVAR_Z_mean_tertile + s(age, bs = \"cr\") + s(age, bs = \"cr\", by = oVAR_Z_mean_tertile) + sex"))
#   gammf_mod4_oquart_cor = as.formula(paste("brainvarcor ~ fVAR_Z_mean_quartile + s(age, bs = \"cr\") + s(age, bs = \"cr\", by = oVAR_Z_mean_quartile) + sex"))
#   
#   gammf_mod2_mean_icv_cor = as.formula(paste("brainvarcor ~ s(age, bs = \"cr\", by = fVAR_Z_mean_meansplit) + fVAR_Z_mean_meansplit + sex + icv"))
#   gammf_mod3_tert_icv_cor = as.formula(paste("brainvarcor ~ s(age, bs = \"cr\", by = fVAR_Z_mean_tertile) + fVAR_Z_mean_tertile + sex + icv"))
#   gammf_mod4_quart_icv_cor = as.formula(paste("brainvarcor ~ s(age, bs = \"cr\", by = fVAR_Z_mean_quartile) + fVAR_Z_mean_quartile + sex + icv"))
#   gammf_mod2_omean_icv_cor = as.formula(paste("brainvarcor ~ fVAR_Z_mean_meansplit + s(age, bs = \"cr\") + s(age, bs = \"cr\", by = oVAR_Z_mean_meansplit) + sex + icv"))
#   gammf_mod3_otert_icv_cor = as.formula(paste("brainvarcor ~ fVAR_Z_mean_tertile + s(age, bs = \"cr\") + s(age, bs = \"cr\", by = oVAR_Z_mean_tertile) + sex + icv"))
#   gammf_mod4_oquart_icv_cor = as.formula(paste("brainvarcor ~ fVAR_Z_mean_quartile + s(age, bs = \"cr\") + s(age, bs = \"cr\", by = oVAR_Z_mean_quartile) + sex + icv"))
#   
#   gammf_mod_long_cor = as.formula(paste("brainvarcor ~ VAR_Z * time_bsl + s(age_bsl, bs = \"cr\") + sex"))
#   lmmf_mod_long_cor = as.formula(paste("brainvarcor ~ VAR_Z * time_bsl + age_bsl + sex"))
#   gammf_mod_groupmean_cor = as.formula(paste("brainvarcor ~ VAR_Z_mean_meansplit * time_bsl + s(age_bsl, bs = \"cr\") + sex"))
#   lmmf_mod_groupmean_cor = as.formula(paste("brainvarcor ~ VAR_Z_mean_meansplit * time_bsl + age_bsl + sex"))
#   
#   gammf_mod_long_icv_cor = as.formula(paste("brainvarcor ~ VAR_Z * time_bsl + s(age_bsl, bs = \"cr\") + sex + icv"))
#   lmmf_mod_long_icv_cor = as.formula(paste("brainvarcor ~ VAR_Z * time_bsl + age_bsl + sex + icv"))
#   gammf_mod_groupmean_icv_cor = as.formula(paste("brainvarcor ~ VAR_Z_mean_meansplit * time_bsl + s(age_bsl, bs = \"cr\") + sex + icv"))
#   lmmf_mod_groupmean_icv_cor = as.formula(paste("brainvarcor ~ VAR_Z_mean_meansplit * time_bsl + age_bsl + sex + icv"))
#   
#   
#   DF_long = DF_num %>% filter(!is.na(VAR_Z)) %>% group_by(id) %>% mutate(
#     age_bsl = min(age),
#     age_mean = mean(age),
#     time_bsl = age - min(age),
#     nTimepoints = length(unique(age)),
#     nTime = row_number(),
#     intervalFirstLast = max(age) - min(age)) %>% ungroup() %>% filter(nTimepoints >= 2)
#   
#   
#   simple_mod_cor = gamm4::gamm4(
#     brainvarcor ~ s(age, bs = "cr"),
#     # random = list(id = ~1),
#     random = ~ (1 |id),
#     data = DF_num
#   )
#   
#   gamm_mod2_mean_cor = gamm4::gamm4(
#     gammf_mod2_mean_cor,
#     random = ~ (1 |id),
#     # random = list(id = ~1),
#     data = DF_num
#   )
#   # plot.gam(gamm_mod2_mean_cor$gam) #plots age across factor0 then factor1
#   gamm_mod2_omean_cor = gamm4::gamm4(
#     gammf_mod2_omean_cor,
#     random = ~ (1 |id),
#     data = DF_num
#   )
#   # plot.gam(gamm_mod2_omean_cor$gam) #plots average age across factors then difference
#   
#   
#   summary(gamm_mod2_mean_cor$gam)
#   summary(gamm_mod2_omean_cor$gam)
#   
#   
#   gamm_mod2_mean_icv_cor = gamm4::gamm4(
#     gammf_mod2_mean_icv_cor,
#     random = ~ (1 |id),
#     # random = list(id = ~1),
#     data = DF_num
#   )
#   gamm_mod2_omean_icv_cor = gamm4::gamm4(
#     gammf_mod2_omean_icv_cor,
#     random = ~ (1 |id),
#     data = DF_num
#   )
#   
#   
#   summary(gamm_mod2_mean_icv_cor$gam)
#   summary(gamm_mod2_omean_icv_cor$gam)
#   
#   
#   gamm_mod3_tert_cor = gamm4::gamm4(
#     gammf_mod3_tert_cor,
#     random = ~ (1 |id),
#     data = DF_num
#   )
#   gamm_mod3_otert_cor = gamm4::gamm4(
#     gammf_mod3_otert_cor,
#     random = ~ (1 |id),
#     data = DF_num
#   )
#   gamm_mod3_tert_icv_cor = gamm4::gamm4(
#     gammf_mod3_tert_icv_cor,
#     random = ~ (1 |id),
#     data = DF_num
#   )
#   gamm_mod3_otert_icv_cor = gamm4::gamm4(
#     gammf_mod3_otert_icv_cor,
#     random = ~ (1 |id),
#     data = DF_num
#   )
#   # gamm_mod4_quart_cor = gamm4::gamm4(gammf_mod4_quart_cor, random = ~ (1 |id), data = DF_num)
#   # gamm_mod4_oquart_cor = gamm4::gamm4(gammf_mod4_oquart_cor, random = ~ (1 |id), data = DF_num)
#   
#   
#   summary(gamm_mod3_tert_cor$gam)
#   summary(gamm_mod3_otert_cor$gam)
#   summary(gamm_mod3_tert_icv_cor$gam)
#   summary(gamm_mod3_otert_icv_cor$gam)
#   # summary(gamm_mod4_quart_cor$gam); summary(gamm_mod4_oquart_cor$gam)
#   
#   
#   gamm_mod_long_cor = gamm4::gamm4(
#     gammf_mod_long_cor,
#     random = ~ (1 |id), #gamm4
#     REML = T, #gamm4
#     # random = list(id = ~1), #mgcv
#     data = DF_long
#   )
#   
#   
#   lmm_mod_long_cor = gamm4::gamm4(
#     lmmf_mod_long_cor,
#     random = ~ (1 |id),
#     REML = T,
#     data = DF_long
#   )
#   
#   
#   gamm_mod_long_icv_cor = gamm4::gamm4(
#     gammf_mod_long_icv_cor,
#     random = ~ (1 |id),
#     REML = T,
#     data = DF_long
#   )
#   
#   
#   lmm_mod_long_icv_cor = gamm4::gamm4(
#     lmmf_mod_long_icv_cor,
#     random = ~ (1 |id),
#     REML = T,
#     data = DF_long
#   )
#   
#   
#   gamm_mod_groupmean_cor = gamm4::gamm4(
#     gammf_mod_groupmean_cor,
#     random = ~ (1 |id),
#     REML = T,
#     data = DF_num
#   )
#   
#   lmm_mod_groupmean_cor = gamm4::gamm4(
#     lmmf_mod_groupmean_cor,
#     random = ~ (1 |id),
#     REML = T,
#     data = DF_num
#   )
#   
#   gamm_mod_groupmean_icv_cor = gamm4::gamm4(
#     gammf_mod_groupmean_icv_cor,
#     random = ~ (1 |id),
#     REML = T,
#     data = DF_num
#   )
#   
#   lmm_mod_groupmean_icv_cor = gamm4::gamm4(
#     lmmf_mod_groupmean_icv_cor,
#     random = ~ (1 |id),
#     REML = T,
#     data = DF_num
#   )
# 
#   #__________
#   cohort_fitsimple = metagam::strip_rawdata(simple_mod)
#   cohort_fit2_mean = metagam::strip_rawdata(gamm_mod2_mean)
#   cohort_fit2_omean = metagam::strip_rawdata(gamm_mod2_omean)
#   cohort_fit3_tert = metagam::strip_rawdata(gamm_mod3_tert)
#   cohort_fit3_otert = metagam::strip_rawdata(gamm_mod3_otert)
#   cohort_fit4_quart = metagam::strip_rawdata(gamm_mod4_quart)
#   cohort_fit4_oquart = metagam::strip_rawdata(gamm_mod4_oquart)
#   cohort_fit5_quint = metagam::strip_rawdata(gamm_mod5_quint)
#   cohort_fit5_oquint = metagam::strip_rawdata(gamm_mod5_oquint)
# 
#   cohort_fit2_mean_icv = metagam::strip_rawdata(gamm_mod2_mean_icv)
#   cohort_fit2_omean_icv = metagam::strip_rawdata(gamm_mod2_omean_icv)
#   cohort_fit3_tert_icv = metagam::strip_rawdata(gamm_mod3_tert_icv)
#   cohort_fit3_otert_icv = metagam::strip_rawdata(gamm_mod3_otert_icv)
# 
#   cohort_gamm_mod_long = metagam::strip_rawdata(gamm_mod_long)
#   cohort_lmm_mod_long = metagam::strip_rawdata(lmm_mod_long)
#   #________
# 
#   #NEWSECTION ----
#   cohort_fitsimple_cor = metagam::strip_rawdata(simple_mod_cor)
#   cohort_fit2_mean_cor = metagam::strip_rawdata(gamm_mod2_mean_cor)
#   cohort_fit2_omean_cor = metagam::strip_rawdata(gamm_mod2_omean_cor)
#   cohort_fit3_tert_cor = metagam::strip_rawdata(gamm_mod3_tert_cor)
#   cohort_fit3_otert_cor = metagam::strip_rawdata(gamm_mod3_otert_cor)
#   
#   cohort_fit2_mean_icv_cor = metagam::strip_rawdata(gamm_mod2_mean_icv_cor)
#   cohort_fit2_omean_icv_cor = metagam::strip_rawdata(gamm_mod2_omean_icv_cor)
#   cohort_fit3_tert_icv_cor = metagam::strip_rawdata(gamm_mod3_tert_icv_cor)
#   cohort_fit3_otert_icv_cor = metagam::strip_rawdata(gamm_mod3_otert_icv_cor)
#   
#   cohort_gamm_mod_long_cor = metagam::strip_rawdata(gamm_mod_long_cor)
#   cohort_lmm_mod_long_cor = metagam::strip_rawdata(lmm_mod_long_cor)
#   cohort_gamm_mod_long_icv_cor = metagam::strip_rawdata(gamm_mod_long_icv_cor)
#   cohort_lmm_mod_long_icv_cor = metagam::strip_rawdata(lmm_mod_long_icv_cor)
#   
#   cohort_gamm_mod_groupmean_cor = metagam::strip_rawdata(gamm_mod_groupmean_cor)
#   cohort_lmm_mod_groupmean_cor = metagam::strip_rawdata(lmm_mod_groupmean_cor)
#   cohort_gamm_mod_groupmean_icv_cor = metagam::strip_rawdata(gamm_mod_groupmean_icv_cor)
#   cohort_lmm_mod_groupmean_icv_cor = metagam::strip_rawdata(lmm_mod_groupmean_icv_cor)
# 
#   return(list(
#     cohort_fitsimple = cohort_fitsimple,
#     cohort_fit2_mean = cohort_fit2_mean,
#     cohort_fit2_omean = cohort_fit2_omean,
#     cohort_fit3_tert = cohort_fit3_tert,
#     cohort_fit3_otert = cohort_fit3_otert,
#     cohort_fit4_quart = cohort_fit4_quart,
#     cohort_fit4_oquart = cohort_fit4_oquart,
#     cohort_fit5_quint = cohort_fit5_quint,
#     cohort_fit5_oquint = cohort_fit5_oquint,
# 
#     cohort_fit2_mean_icv = cohort_fit2_mean_icv,
#     cohort_fit2_omean_icv = cohort_fit2_omean_icv,
#     cohort_fit3_tert_icv = cohort_fit3_tert_icv,
#     cohort_fit3_otert_icv = cohort_fit3_otert_icv,
# 
#     cohort_gamm_mod_long = cohort_gamm_mod_long,
#     cohort_lmm_mod_long = cohort_lmm_mod_long,
#     
#     cohort_fitsimple_cor = cohort_fitsimple_cor,
#     cohort_fit2_mean_cor = cohort_fit2_mean_cor,
#     cohort_fit2_omean_cor = cohort_fit2_omean_cor,
#     cohort_fit3_tert_cor = cohort_fit3_tert_cor,
#     cohort_fit3_otert_cor = cohort_fit3_otert_cor,
#     
#     cohort_fit2_mean_icv_cor = cohort_fit2_mean_icv_cor,
#     cohort_fit2_omean_icv_cor = cohort_fit2_omean_icv_cor,
#     cohort_fit3_tert_icv_cor = cohort_fit3_tert_icv_cor,
#     cohort_fit3_otert_icv_cor = cohort_fit3_otert_icv_cor,
#     
#     cohort_gamm_mod_long_cor = cohort_gamm_mod_long_cor,
#     cohort_lmm_mod_long_cor = cohort_lmm_mod_long_cor
#   ))
#   print("finished")
# }
# 
# 
# rerunMeta = 1
# seed = 17
# nSplits = 3
# if (rerunMeta == 1) {
# 
#   print(paste("splitting DF by", nSplits))
#   if (nSplits == 2) {
#     DF1_num = splitSample(DF, U, seed, nSplits = 2)$DF1
#     DF2_num = splitSample(DF, U, seed, nSplits = 2)$DF2
#   } else {
#     DF1_num = splitSample(DF, U, seed)$DF1
#     DF2_num = splitSample(DF, U, seed)$DF2
#     DF3_num = splitSample(DF, U, seed)$DF3
#   }
# 
#   print("prepping")
#   DF1_num = testFirst(DF1_num)
#   DF2_num = testFirst(DF2_num)
#   if (nSplits == 3) {
#     DF3_num = testFirst(DF3_num)
#   }
# 
#   #models
#   meta1 = testSecond(DF1_num)
#   meta2 = testSecond(DF2_num, iteration = 2)
#   if (nSplits == 3) {
#     meta3 = testSecond(DF3_num, iteration = 3)
#   }
# 
#   # DF1_num$id %in% DF2_num$id
#   # DF1_num$id %in% DF3_num$id
# 
#   if (nSplits == 2) {
#     simple_mod = list(meta1$cohort_fitsimple,
#                       meta2$cohort_fitsimple)
# 
#     cohort_fit2_mean = list(meta1$cohort_fit2_mean,
#                             meta2$cohort_fit2_mean)
# 
#     cohort_fit2_omean = list(meta1$cohort_fit2_omean,
#                              meta2$cohort_fit2_omean)
# 
#     cohort_fit3_tert = list(meta1$cohort_fit3_tert,
#                             meta2$cohort_fit3_tert)
# 
#     cohort_fit3_otert = list(meta1$cohort_fit3_otert,
#                              meta2$cohort_fit3_otert)
# 
#     cohort_fit4_quart = list(meta1$cohort_fit4_quart,
#                              meta2$cohort_fit4_quart)
# 
#     cohort_fit5_quint = list(meta1$cohort_fit5_quint,
#                              meta2$cohort_fit5_quint)
# 
#     cohort_fit5_oquint = list(meta1$cohort_fit5_oquint,
#                               meta2$cohort_fit5_oquint)
# 
#     cohort_fit2_mean_icv = list(
#       meta1$cohort_fit2_mean_icv,
#       meta2$cohort_fit2_mean_icv
#     )
#     cohort_fit2_omean_icv = list(
#       meta1$cohort_fit2_omean_icv,
#       meta2$cohort_fit2_omean_icv
#     )
#     cohort_fit3_tert_icv = list(
#       meta1$cohort_fit3_tert_icv,
#       meta2$cohort_fit3_tert_icv
#     )
#     cohort_fit3_otert_icv = list(
#       meta1$cohort_fit3_otert_icv,
#       meta2$cohort_fit3_otert_icv
#     )
#   } else if (nSplits == 3) {
#     simple_mod = list(meta1$cohort_fitsimple,
#                    meta2$cohort_fitsimple,
#                    meta3$cohort_fitsimple)
#     cohort_fit2_mean = list(meta1$cohort_fit2_mean,
#                             meta2$cohort_fit2_mean,
#                             meta3$cohort_fit2_mean)
#     cohort_fit2_omean = list(meta1$cohort_fit2_omean,
#                              meta2$cohort_fit2_omean,
#                              meta3$cohort_fit2_omean)
#     cohort_fit3_tert = list(meta1$cohort_fit3_tert,
#                             meta2$cohort_fit3_tert,
#                             meta3$cohort_fit3_tert)
#     cohort_fit3_otert = list(meta1$cohort_fit3_otert,
#                              meta2$cohort_fit3_otert,
#                              meta3$cohort_fit3_otert)
#     cohort_fit4_quart = list(meta1$cohort_fit4_quart,
#                              meta2$cohort_fit4_quart,
#                              meta3$cohort_fit4_quart)
#     cohort_fit5_quint = list(meta1$cohort_fit5_quint,
#                              meta2$cohort_fit5_quint,
#                              meta3$cohort_fit5_quint)
#     cohort_fit5_oquint = list(meta1$cohort_fit5_oquint,
#                               meta2$cohort_fit5_oquint,
#                               meta3$cohort_fit5_oquint)
#     cohort_fit2_mean_icv = list(
#       meta1$cohort_fit2_mean_icv,
#       meta2$cohort_fit2_mean_icv,
#       meta3$cohort_fit2_mean_icv
#     )
#     cohort_fit2_omean_icv = list(
#       meta1$cohort_fit2_omean_icv,
#       meta2$cohort_fit2_omean_icv,
#       meta3$cohort_fit2_omean_icv
#     )
#     cohort_fit3_tert_icv = list(
#       meta1$cohort_fit3_tert_icv,
#       meta2$cohort_fit3_tert_icv,
#       meta3$cohort_fit3_tert_icv
#     )
#     cohort_fit3_otert_icv = list(
#       meta1$cohort_fit3_otert_icv,
#       meta2$cohort_fit3_otert_icv,
#       meta3$cohort_fit3_otert_icv
#     )
#     
#     simple_mod_cor = list(meta1$cohort_fitsimple_cor,
#                       meta2$cohort_fitsimple_cor,
#                       meta3$cohort_fitsimple_cor)
#     
#     cohort_fit2_mean_cor = list(meta1$cohort_fit2_mean_cor,
#                             meta2$cohort_fit2_mean_cor,
#                             meta3$cohort_fit2_mean_cor)
#     cohort_fit2_omean_cor = list(meta1$cohort_fit2_omean_cor,
#                              meta2$cohort_fit2_omean_cor,
#                              meta3$cohort_fit2_omean_cor)
#     cohort_fit3_tert_cor = list(meta1$cohort_fit3_tert_cor,
#                             meta2$cohort_fit3_tert_cor,
#                             meta3$cohort_fit3_tert_cor)
#     cohort_fit3_otert_cor = list(meta1$cohort_fit3_otert_cor,
#                              meta2$cohort_fit3_otert_cor,
#                              meta3$cohort_fit3_otert_cor)
#     
#     cohort_fit2_mean_icv_cor = list(
#       meta1$cohort_fit2_mean_icv_cor,
#       meta2$cohort_fit2_mean_icv_cor,
#       meta3$cohort_fit2_mean_icv_cor
#     )
#     cohort_fit2_omean_icv_cor = list(
#       meta1$cohort_fit2_omean_icv_cor,
#       meta2$cohort_fit2_omean_icv_cor,
#       meta3$cohort_fit2_omean_icv_cor
#     )
#     cohort_fit3_tert_icv_cor = list(
#       meta1$cohort_fit3_tert_icv_cor,
#       meta2$cohort_fit3_tert_icv_cor,
#       meta3$cohort_fit3_tert_icv_cor
#     )
#     cohort_fit3_otert_icv_cor = list(
#       meta1$cohort_fit3_otert_icv_cor,
#       meta2$cohort_fit3_otert_icv_cor,
#       meta3$cohort_fit3_otert_icv_cor
#     )
#     
#   }
#   grid0 = data.frame(age = seq(4, 90),
#                      sex = factor("Female", levels = c("Female", "Male")),
#                      fVAR_Z_mean_meansplit = factor(0, levels = c(0, 1)),
#                      VAR_Z_mean_meansplit = 0,
#                      oVAR_Z_mean_meansplit = ordered(0, levels = c(0, 1)),
#                      # scanner = factor("ousSkyra", levels = c("ousAvanto",  "ousSkyra",   "ousPrisma",  "ntnuAvanto")),
#                      # scanner = factor("ascanner1", levels = c("ascanner1",  "ousSkyra",   "ousPrisma",  "ntnuAvanto")),
#                      scanner = factor("base", levels = c("base",  "another")),
#                      icv = 0
#   )
#   grid1 = data.frame(age = seq(4, 90),
#                      sex = factor("Female", levels = c("Female", "Male")),
#                      fVAR_Z_mean_meansplit = factor(1, levels = c(0, 1)),
#                      VAR_Z_mean_meansplit = 1,
#                      oVAR_Z_mean_meansplit = ordered(1, levels = c(0, 1)),
#                      # scanner = factor("ascanner1", levels = c("ascanner1",  "ousSkyra",   "ousPrisma",  "ntnuAvanto")),
#                      scanner = factor("base", levels = c("base",  "another")),
#                      icv = 0
#   )
#   grid = rbind(grid0, grid1)
# }
# 
# test = metagam::metagam(simple_mod, grid = grid, terms = "s(age)")
# plot(test, ci = "pointwise")
# test = metagam::metagam(simple_mod_cor, grid = grid, terms = "s(age)")
# plot(test, ci = "pointwise")
# 
# meta_fit2_mean0 = metagam::metagam(cohort_fit2_mean, grid0, terms = "s(age):fVAR_Z_mean_meansplit0", nsim=1000)
# meta_fit2_mean1 = metagam::metagam(cohort_fit2_mean_cor, grid1, terms = "s(age):fVAR_Z_mean_meansplit1", nsim=1000)
# plot(meta_fit2_mean0, ci = "pointwise")
# plot(meta_fit2_mean1, ci = "pointwise")
# 
# 
# meta_fit2_omean = metagam::metagam(cohort_fit2_omean_cor, grid1, terms = "s(age):oVAR_Z_mean_meansplit1",  nsim=1000)
# plot(meta_fit2_omean, ci = "pointwise")
# 
# 
# meta_fit3_tert1 = metagam::metagam(cohort_fit3_tert, grid0 %>% mutate(fVAR_Z_mean_tertile = factor(1, levels = c(1, 2, 3))), terms = "s(age):fVAR_Z_mean_tertile1", nsim = 1000)
# meta_fit3_tert2 = metagam::metagam(cohort_fit3_tert, grid0 %>% mutate(fVAR_Z_mean_tertile = factor(2, levels = c(1, 2, 3))), terms = "s(age):fVAR_Z_mean_tertile2", nsim = 1000)
# meta_fit3_tert3 = metagam::metagam(cohort_fit3_tert, grid0 %>% mutate(fVAR_Z_mean_tertile = factor(3, levels = c(1, 2, 3))), terms = "s(age):fVAR_Z_mean_tertile3", nsim = 1000)
# plot(meta_fit3_tert1, ci = "pointwise")
# plot(meta_fit3_tert2, ci = "pointwise")
# plot(meta_fit3_tert3, ci = "pointwise")
# # 
# # 
# # meta_fit3_otert2 = metagam::metagam(cohort_fit3_otert, grid0 %>% mutate(oVAR_Z_mean_tertile = factor(2, levels = c(1, 2, 3)),
# #                                                                         fVAR_Z_mean_tertile = factor(2, levels = c(1, 2, 3))), terms = "s(age):oVAR_Z_mean_tertile2", nsim = 1000)
# # meta_fit3_otert3 = metagam::metagam(cohort_fit3_otert, grid0 %>% mutate(oVAR_Z_mean_tertile = factor(3, levels = c(1, 2, 3)),
# #                                                                         fVAR_Z_mean_tertile = factor(3, levels = c(1, 2, 3))), terms = "s(age):oVAR_Z_mean_tertile3", nsim = 1000)
# # plot(meta_fit3_otert2, ci = "pointwise") + geom_hline(yintercept = 0, linetype = 3)
# # plot(meta_fit3_otert3, ci = "pointwise") + geom_hline(yintercept = 0, linetype = 3)
# # 
# # # plot_dominance(simple_mod)
# # # plot_heterogeneity(simple_mod)
# # 
# # 
# # #-----goal
# # #relative effects on intercept and slope of various lifespan factors
# # #-----goal
