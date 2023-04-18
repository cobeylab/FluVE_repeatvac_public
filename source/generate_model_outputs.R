# load utility file
source(here::here("source/utils.R"))
reload_source()


#### functions
##'
get_repvac_OR_wrt01 <- function(model_output){
  coef_index11 <- which(names(coef(model_output))=="vac_group11")
  coef_index01 <- which(names(coef(model_output))=="vac_group01")
  coef <- coef(model_output)[coef_index11] - coef(model_output)[coef_index01]
  se <- deltamethod(as.formula(paste("~x",coef_index11,"-x",coef_index01, sep="")), 
                    coef(model_output), vcov(model_output))
  c(exp(c(coef, coef-1.96*se, coef+1.96*se)))
}
##'
get_lastinf_glm <- function(model_output){
  coef_index01 <- which(substr(names(coef(model_output)),1,18)=="last_infection_cat")
  data.frame(matrix(exp(cbind(coef(model_output),confint.default(model_output))[coef_index01,]),ncol=3),
             infyr=substr(names(coef(model_output)),24,30)[coef_index01])
}
##'
get_lastinf_glm_neg <- function(model_output){
  coef_index01 <- which(substr(names(coef(model_output)),1,18)=="last_infection_cat")
  data.frame(matrix(exp(cbind(coef(model_output),confint.default(model_output))[coef_index01,]),ncol=3),
             infyr=substr(names(coef(model_output)),24+4,30+4)[coef_index01])
}
##'
get_repvac_glm_0117 <- function(model_output){
  coef_index <- which(names(coef(model_output))=="vac_group11")
  exp(cbind(coef(model_output),confint.default(model_output))[coef_index,])
}
##'
get_waning_glm_0117 <- function(model_output){
  coef_index <- which(substr(names(coef(model_output)),1,21)=="time_VAX_to_ONSET_cat")
  exp(cbind(coef(model_output),confint.default(model_output))[coef_index,])
}

#filter_type_by_season <- unique(dat1_nofil[,c("SEASON","H3N2_dom","B_year")]) %>% 
#  mutate(SEASON = as.character(SEASON))

##'
get_se_coef2 <- function(model_output){
  coef_index1 <- which(names(coef(model_output))=="RES_any")
  coef_index2 <- which(substr(names(coef(model_output)),1,8)=="RES_any:")
  coef <- c(coef(model_output)[coef_index1], coef(model_output)[coef_index1] + coef(model_output)[coef_index2])
  se <- c(sqrt(vcov(model_output)[coef_index1,coef_index1]), 
          unlist(lapply(1:length(coef_index2),
                        function(i)deltamethod(as.formula(paste("~x",coef_index1,"+x",coef_index2[i], sep="")),
                                               coef(model_output), vcov(model_output)))))
  vac <- c("no","yes")
  return(list(coef,se,vac))
}
##'
get_se_coef <- function(model_output){
  coef_index1 <- which(names(coef(model_output))=="RES_any")
  coef_index2 <- which(substr(names(coef(model_output)),1,8)=="RES_any:")
  coef <- c(coef(model_output)[coef_index1], coef(model_output)[coef_index1] + coef(model_output)[coef_index2])
  se <- c(sqrt(vcov(model_output)[coef_index1,coef_index1]), 
          unlist(lapply(1:length(coef_index2),
                        function(i)deltamethod(as.formula(paste("~x",coef_index1,"+x",coef_index2[i], sep="")),
                                               coef(model_output), vcov(model_output)))))
  age <- c("0-4",substr(names(coef),18,22)[-1])
  return(list(coef,se,age))
}


##' this generates part of the data used for figure 1
##'
##' @param dat_glm cleaned FLU VE Network data
##'
##' @return model outputs before adjusting for waning
##' estimates are available by season and site
##' 
##'
get_output_glm <- function(dat_glm){
  all_seasons <- as.character(unique(dat_glm$SEASON))
  all_sites <- as.character(unique(dat_glm$SITE))
  # create output
  out <- list()
  # 
  list_names <- NULL
  list_i <- 0
  # all sites and seasons
  for (var in c("RES_B", "RES_H1", "RES_H3")){
    var_index <- which(names(dat_glm) == var)
    dat_index <- dat_glm[,var_index] <=1
    model <- NULL
    model <- glm(as.formula(paste(var, "~ vac_group + age_cat_5 + P_SEX + ANY_HR_PY + 
                                    factor(ONSET_MONTH_shift) + as.factor(SEASON) +
                                    as.factor(SITE)")),
                 family="binomial",
                 data = dat_glm[dat_index,])
    list_names <- c(list_names, paste(substr(var,5,6),"all", "all", sep="_"))
    list_i <- list_i + 1
    out[[list_i]] <- model
  }

  # by site
  for (var in c("RES_B", "RES_H1", "RES_H3")){
    for (site_i in all_sites){
      var_index <- which(names(dat_glm) == var)
      var_index2 <- which(names(dat_glm) == "SITE")
      dat_index <- (dat_glm[,var_index] <=1) & (dat_glm[,var_index2] == site_i)
      model <- NULL
      model <- glm(as.formula(paste(var, "~ vac_group + age_cat_5 + P_SEX + ANY_HR_PY + 
                                    as.factor(ONSET_MONTH_shift) + as.factor(SEASON)")),
                   family="binomial",
                   data = dat_glm[dat_index,])
      list_names <- c(list_names, paste(substr(var,5,6),site_i, "all", sep="_"))
      list_i <- list_i + 1
      out[[list_i]] <- model
    }
  }
  
  # by season
  for (var in c("RES_B", "RES_H1", "RES_H3")){
    for (season_i in all_seasons){
      var_index <- which(names(dat_glm) == var)
      var_index2 <- which(names(dat_glm) == "SEASON")
      dat_index <- (dat_glm[,var_index] <=1) & (dat_glm[,var_index2] == season_i)
      model <- NULL
      model <- glm(as.formula(paste(var, "~ vac_group + age_cat_5 + P_SEX + ANY_HR_PY + 
                                    as.factor(ONSET_MONTH_shift) + as.factor(SITE)")),
                   family="binomial",
                   data = dat_glm[dat_index,])
      list_names <- c(list_names, paste(substr(var,5,6),"all", season_i, sep="_"))
      list_i <- list_i + 1
      out[[list_i]] <- model
    }
  }
  
  # by season and site
  for (var in c("RES_B", "RES_H1", "RES_H3")){
    for (season_i in all_seasons){
      for (site_i in all_sites){
        var_index <- which(names(dat_glm) == var)
        var_index2 <- which(names(dat_glm) == "SEASON")
        var_index3 <- which(names(dat_glm) == "SITE")
        dat_index <- (dat_glm[,var_index] <=1) & (dat_glm[,var_index2] == season_i)& (dat_glm[,var_index3] == site_i)
        model <- NULL
        model <- glm(as.formula(paste(var, "~ vac_group + age_cat_5 + P_SEX + ANY_HR_PY + 
                                      as.factor(ONSET_MONTH_shift)")),
                     family="binomial",
                     data = dat_glm[dat_index,])
        list_names <- c(list_names, paste(substr(var,5,6),site_i, season_i, sep="_"))
        list_i <- list_i + 1
        out[[list_i]] <- model
      }
    }
  }
  names(out) <- list_names
  return(out)
}

waning_not_adj_main_output <- get_output_glm(dat_glm)

waning_not_adj_main_output_ray <- get_output_glm(
  dat_glm %>% filter( (CONFIRMED_VAX1_WEEK_shift < first_case_week)|
                        is.na(CONFIRMED_VAX1_WEEK_shift))
)



##' this generates part of the data used for figure 1
##'
##' @param dat_glm cleaned FLU VE Network data
##'
##' @return model outpus after adjusting for waning
##' estimates are available by season and site
##' 
##' 
get_waning_adj_main_output <- function(dat_glm){
  all_seasons <- unique(dat_glm$SEASON)
  all_sites <- unique(dat_glm$SITE)
  out <- list()
  list_names <- NULL
  list_i <- 0
  # all season site combined
  for (var in c("RES_B", "RES_H1", "RES_H3")){
    var_index <- which(names(dat_glm) == var)
    dat_index <- dat_glm[,var_index] <= 1
    model <- NULL
    model <- glm(as.formula(paste(var, "~ P_SEX + ANY_HR_PY + age_cat_5 + 
                                      vac_group11 + vac_group10 + 
                                      time_VAX_to_ONSET_cat_2_9 + 
                                      time_VAX_to_ONSET_cat_10_13 + 
                                      time_VAX_to_ONSET_cat_14_17 + 
                                      time_VAX_to_ONSET_cat_18_21 +
                                      time_VAX_to_ONSET_cat_22 +
                                      as.factor(ONSET_MONTH_shift) + as.factor(SEASON) +
                                      as.factor(SITE)")),
                 family="binomial",
                 data = dat_glm[dat_index,])
    list_names <- c(list_names, paste(substr(var,5,6),"all","all", sep="_"))
    list_i <- list_i + 1
    out[[list_i]] <- model
  }
  
  # by season
  for (var in c("RES_B", "RES_H1", "RES_H3")){
    for (season_i in all_seasons){
      var_index <- which(names(dat_glm) == var)
      var_index2 <- which(names(dat_glm) == "SEASON")
      dat_index <- (dat_glm[,var_index] <=1) & (dat_glm[,var_index2] == season_i)
      model <- NULL
      model <- glm(as.formula(paste(var, "~ P_SEX + ANY_HR_PY + age_cat_5 + 
                                      vac_group11 + vac_group10 + 
                                      time_VAX_to_ONSET_cat_2_9 + time_VAX_to_ONSET_cat_10_13 + 
                                      time_VAX_to_ONSET_cat_14_17 + time_VAX_to_ONSET_cat_18_21 +
                                      time_VAX_to_ONSET_cat_22 +
                                      as.factor(ONSET_MONTH_shift) +
                                    as.factor(SITE)")),
                   family="binomial",
                   data = dat_glm[dat_index,])
      list_names <- c(list_names, paste(substr(var,5,6),"all", season_i, sep="_"))
      list_i <- list_i + 1
      out[[list_i]] <- model
    }
  }
  
  # by site
  for (var in c("RES_B", "RES_H1", "RES_H3")){
    for (site_i in all_sites){
      var_index <- which(names(dat_glm) == var)
      var_index2 <- which(names(dat_glm) == "SITE")
      dat_index <- (dat_glm[,var_index] <=1) & (dat_glm[,var_index2] == site_i)
      model <- NULL
      model <- glm(as.formula(paste(var, "~ P_SEX + ANY_HR_PY + age_cat_5 + vac_group11 + vac_group10 + 
                                      time_VAX_to_ONSET_cat_2_9 + time_VAX_to_ONSET_cat_10_13 + 
                                      time_VAX_to_ONSET_cat_14_17 + time_VAX_to_ONSET_cat_18_21 +
                                      time_VAX_to_ONSET_cat_22 +
                                      as.factor(ONSET_MONTH_shift) + as.factor(SEASON)")),
                   family="binomial",
                   data = dat_glm[dat_index,])
      list_names <- c(list_names, paste(substr(var,5,6),  site_i, "all",sep="_"))
      list_i <- list_i + 1
      out[[list_i]] <- model
    }
  }
  # by season and by site
  for (var in c("RES_B", "RES_H1", "RES_H3")){
    for (site_i in all_sites){
      for (season_i in all_seasons){
        var_index <- which(names(dat_glm) == var)
        var_index2 <- which(names(dat_glm) == "SITE")
        var_index3 <- which(names(dat_glm) == "SEASON")
        dat_index <- (dat_glm[,var_index] <=1) & (dat_glm[,var_index2] == site_i)& (dat_glm[,var_index3] == season_i)
        model <- NULL
        model <- glm(as.formula(paste(var, "~ P_SEX + ANY_HR_PY + age_cat_5 + vac_group11 + vac_group10 + 
                                      time_VAX_to_ONSET_cat_2_9 + time_VAX_to_ONSET_cat_10_13 + 
                                      time_VAX_to_ONSET_cat_14_17 + time_VAX_to_ONSET_cat_18_21 +
                                      time_VAX_to_ONSET_cat_22 +
                                      as.factor(ONSET_MONTH_shift)")),
                     family="binomial",
                     data = dat_glm[dat_index,])
        list_names <- c(list_names, paste(substr(var,5,6), site_i, season_i, sep="_"))
        list_i <- list_i + 1
        out[[list_i]] <- model
      }
    }
  }
  names(out) <- list_names
  return(out)
}

waning_adj_main_output <- get_waning_adj_main_output(dat_glm)

waning_adj_main_output_ray <- get_waning_adj_main_output(
  dat_glm %>% 
    filter( (CONFIRMED_VAX1_WEEK_shift < first_case_week)|
              is.na(CONFIRMED_VAX1_WEEK_shift))
)


##' (base model in figure 2c)
##'
##' @param dat_t1 current season data
##' @param all_seasons_mf seasons analyzed
##'
##' @return model outputs;
##' no adjustment for waning; no weights;
##' data were pooled across all sites and seasons
##' 
##'
get_output_glm_mf <- function(dat_t1, all_seasons_mf){
  out <- list()
  # base model, not adjusting for seasons since the last clinical infection
  for (var in c("RES_B", "RES_H1", "RES_H3")){
    var_index <- which(names(dat_t1) == var)
    dat_index <- dat_t1[,var_index] <=1
    
    model <- NULL
    model <- glm(as.formula(paste(var, "~ vac_group + age_cat_5 + gender + HR + 
                                  as.factor(ONSET_MONTH) + as.factor(year)")),
                 family="binomial",
                 data = dat_t1[dat_index,] )
    if (var=="RES_H1") {out$H1 <- model}
    if (var=="RES_H3") {out$H3 <- model}
    if (var=="RES_B") {out$B <- model}
  }
  return(out)
}

# main analyses
output_glm_mf <- get_output_glm_mf(
  msm_long %>% filter(time==1) %>% filter(!year %in% c("2008","2009Pan","2010")), 
  unique(msm_long  %>% filter(!year %in% c("2008","2009Pan","2010")) %>% 
           pull(year))
  )

# sensitivity analyses: exclude data from those who refused enrollment in the previous season 
output_glm_mf_refused_incl <- get_output_glm_mf(
  msm_long %>% filter(time==1) %>% 
    #filter(refused_lastyear!=1|is.na(refused_lastyear)) %>% 
    mutate(refused_lastyear = replace_na(refused_lastyear, 4)) %>%
    filter(!(maari_year_prior == 1 & refused_lastyear == 1)) %>%
    filter(!year %in% c("2008","2009Pan","2010")), 
  unique(msm_long  %>% filter(!year %in% c("2008","2009Pan","2010")) %>%
           pull(year))
  )

##' (figure 2a)
##' @param dat_t1 current season data
##' @param all_seasons_mf seasons analyzed
##'
##' @return model output
##' models adjusted for years since the last 
##' clinical infection of the homologous (sub)type 
##' data were pooled across all sites and seasons
##' 
##'
get_output_glm_mf_lastinf_sametype <- function(dat_t1, all_seasons_mf){
  #
  out <- list()
  # glm base
  for (var in c("RES_B", "RES_H1", "RES_H3")){
    var_index <- which(names(dat_t1) == var)
    dat_index <- dat_t1[,var_index] <=1
    
    model <- NULL
    if (var=="RES_B"){
      model <- glm(as.formula(paste(var, "~ vac_group + age_cat_5 + gender + HR + last_infection_cat_BB_2 +
                                  as.factor(ONSET_MONTH) + as.factor(year)")),
                   family="binomial",
                   data = dat_t1[dat_index,])
      out$B_base <- model
    }
    if (var=="RES_H1"){
      model <- glm(as.formula(paste(var, "~ vac_group + age_cat_5 + gender + HR + last_infection_cat_H1_2 +
                                  as.factor(ONSET_MONTH) + as.factor(year)")),
                   family="binomial",
                   data = dat_t1[dat_index,])
      out$H1_base <- model
    }
    if (var=="RES_H3"){
      model <- glm(as.formula(paste(var, "~ vac_group + age_cat_5 + gender + HR + last_infection_cat_H3_2 +
                                  as.factor(ONSET_MONTH) + as.factor(year)")),
                   family="binomial",
                   data = dat_t1[dat_index,])
      out$H3_base <- model
    }
  }
  return(out)
}

output_glm_mf_lastinf_sametype <- get_output_glm_mf_lastinf_sametype(
  msm_long %>% filter(time==1) %>% filter(!year %in% c("2008")),
  unique(msm_long %>% filter(year!=2008) %>% pull(year))
)


# data for the sensitivity analyses: exclude data from those who refused enrollment in the previous season
output_glm_mf_lastinf_sametype_refused_incl <- get_output_glm_mf_lastinf_sametype(
  msm_long %>% filter(time==1) %>% filter(!year %in% c("2008","2009Pan","2010")) %>% 
    #filter(refused_lastyear!=1|is.na(refused_lastyear)) %>% 
    mutate(refused_lastyear = replace_na(refused_lastyear, 4)) %>%
    filter(!(maari_year_prior == 1 & refused_lastyear == 1)) ,
  unique(msm_long %>% filter(!year %in% c("2008","2009Pan","2010")) %>% pull(year))
)


##' MSM model output
##' no waning in weighted outcome model
get_output_msm_nowaning <- function(dat, all_seasons_mf){
  out <- list()
  for (var in c("RES_B", "RES_H1", "RES_H3")){
    var_index <- which(names(dat) == var)
    dat_index <- dat[,var_index] <= 1
    dat_short <- dat[dat_index,]
    model <- NULL
    cond_stabilized <- ipwtm(
      exposure = vac,
      family = "binomial",
      link = "logit",
      numerator = ~ lag_vac + time + age_cat_5 + HR + gender,
      denominator = ~ lag_vac*inf_pastseason + time + age_cat_5 + HR + gender + year,
      id = studyid,
      timevar = time,
      type = "all",
      data = as.data.frame(dat_short)
    )
    dat_short$cond_stab_weight <- cond_stabilized$ipw.weights
    model <- glm(as.formula(paste(var,"~ age_cat_5 + gender + HR + vac_group")), 
                 weights = cond_stab_weight,
                 family="binomial",
                 data = dat_short %>% filter(time == 1))
    if(var=="RES_B"){
      out$B_MSM <- model
      out$B_MSM_weights <-  cond_stabilized$ipw.weights
    }
    if(var=="RES_H3"){
      out$H3_MSM <- model
      out$H3_MSM_weights <-  cond_stabilized$ipw.weights
    }
    if(var=="RES_H1"){
      out$H1_MSM <- model
      out$H1_MSM_weights <-  cond_stabilized$ipw.weights
    }
  }
  return(out)
}

# # data for the main analyses
# output_msm_nowanings <- get_output_msm_nowaning(
#   msm_long %>% filter(!year %in% c("2008","2009Pan","2010")) %>%
#     filter( (vd1_WEEK_shift < first_case_week)|is.na(vd1_WEEK_shift)), 
#   unique(msm_long %>% filter(!year %in% c("2008","2009Pan","2010")) %>% pull(year))
# )


output_msm_nowanings_refused_incl <- get_output_msm_nowaning(
  msm_long %>% 
    #filter(refused_lastyear!=1|is.na(refused_lastyear)) %>% 
    mutate(refused_lastyear = replace_na(refused_lastyear, 4)) %>%
    filter(!(maari_year_prior == 1 & refused_lastyear == 1)) %>%
    filter(!year %in% c("2008","2009Pan","2010")) %>%
    filter( (vd1_WEEK_shift < first_5pct_case_week)|is.na(vd1_WEEK_shift)), 
  unique(msm_long %>% filter(!year %in% c("2008","2009Pan","2010")) %>% 
           pull(year))
)



##' MSM model output
##' adjusted for waning in weighted outcome model
get_output_msm_waning <- function(dat, all_seasons_mf){
  out <- list()
  list_i <- 0
  list_names <- NULL
  # results by season
  for (var in c("RES_B", "RES_H1", "RES_H3")){
    for (season_i in all_seasons_mf){
        var_index <- which(names(dat) == var)
        var_index2 <- which(names(dat) == "year")
        dat_index <- (dat[,var_index] <=1) & (dat[,var_index2] == season_i)
        dat_short <- dat[dat_index,]
        model <- NULL
        cond_stabilized <- ipwtm(
          exposure = vac,
          family = "binomial",
          link = "logit",
          numerator = ~ lag_vac + time + age_cat_5 + HR + gender,
          denominator = ~ lag_vac*inf_pastseason + time + age_cat_5 + HR + gender,
          id = studyid,
          timevar = time,
          type = "all",
          data = as.data.frame(dat_short)
        )
        dat_short$cond_stab_weight <- cond_stabilized$ipw.weights
        model <- glm(as.formula(paste(var,"~ age_cat_5 + gender + HR + 
                                  vac_group11 + vac_group10 +  
                                  time_VAX_to_ONSET_cat_2_9 + time_VAX_to_ONSET_cat_10_13 + 
                                  time_VAX_to_ONSET_cat_14_17 + time_VAX_to_ONSET_cat_18_21 +
                                  time_VAX_to_ONSET_cat_22
                                  ")), 
                     family="binomial",
                     data = dat_short %>% filter(time == 1))
        
        list_names <- c(list_names, paste(substr(var,5,6), season_i, sep="_"))
        list_i <- list_i + 1
        out[[list_i]] <- model
    }
  }
  names(out) <- list_names
  
  for (var in c("RES_B", "RES_H1", "RES_H3")){
    var_index <- which(names(dat) == var)
    dat_index <- dat[,var_index] <=1
    dat_short <- dat[dat_index,]
    model <- NULL
    cond_stabilized <- ipwtm(
      exposure = vac,
      family = "binomial",
      link = "logit",
      numerator = ~ lag_vac + time + age_cat_5 + HR + gender,
      denominator = ~ lag_vac*inf_pastseason + time + age_cat_5 + HR + gender + year,
      id = studyid,
      timevar = time,
      type = "all",
      data = as.data.frame(dat_short)
    )
    dat_short$cond_stab_weight <- cond_stabilized$ipw.weights
    model <- glm(as.formula(paste(var,"~ age_cat_5 + gender + HR + 
                                  vac_group11 + vac_group10 +  
                                  time_VAX_to_ONSET_cat_2_9 + time_VAX_to_ONSET_cat_10_13 + 
                                  time_VAX_to_ONSET_cat_14_17 + time_VAX_to_ONSET_cat_18_21 +
                                  time_VAX_to_ONSET_cat_22
                                  ")), 
                 family="binomial",
                 data = dat_short %>% filter(time == 1))
    if(var=="RES_B"){
      out$B_MSM <- model
      out$B_MSM_weights <-  cond_stabilized$ipw.weights
    }
    if(var=="RES_H3"){
      out$H3_MSM <- model
      out$H3_MSM_weights <-  cond_stabilized$ipw.weights
    }
    if(var=="RES_H1"){
      out$H1_MSM <- model
      out$H1_MSM_weights <-  cond_stabilized$ipw.weights
    }
  }
  return(out)
}


output_msm_wanings_5pct <- get_output_msm_waning(
  msm_long %>% 
    filter(!year %in% c("2008","2009Pan","2010")) %>%
    filter( (vd1_WEEK_shift < first_5pct_case_week)|is.na(vd1_WEEK_shift)), 
  unique(msm_long %>% filter(!year %in% c("2008","2009Pan","2010")) %>% pull(year))
)



output_msm_nowanings_5pct <- get_output_msm_nowaning(
  msm_long %>% filter(!year %in% c("2008","2009Pan","2010")) %>%
    filter( (vd1_WEEK_shift < first_5pct_case_week)|is.na(vd1_WEEK_shift)), 
  unique(msm_long %>% filter(!year %in% c("2008","2009Pan","2010")) %>% pull(year))
)


# we want to delete those who have a maari visit but refused enrollment
output_msm_wanings_refused_incl <- get_output_msm_waning(
  msm_long %>% 
    #filter(refused_lastyear!=1|is.na(refused_lastyear)) %>% 
    mutate(refused_lastyear = replace_na(refused_lastyear, 4)) %>%
    filter(!(maari_year_prior == 1 & refused_lastyear == 1)) %>%
    filter(!year %in% c("2008","2009Pan","2010")) %>%
    filter( (vd1_WEEK_shift < first_5pct_case_week)|is.na(vd1_WEEK_shift)), 
  unique(msm_long %>% filter(!year %in% c("2008","2009Pan","2010")) %>% 
           pull(year))
)

##' (figure Sfig3.2)
##'
##' @param dat_t1 data from the current season only
##' @param all_seasons_mf seasons analyzed
##'
##' @return model output
##' models adjusted for years since the last 
##' clinical infection of the HETEROLOGOUS (sub)type 
##' data were pooled across all sites and seasons
get_output_glm_mf_lastinf_heterotype <- function(dat_t1, all_seasons_mf){
  out <- list()
  # glm base
  for (var in c("RES_B", "RES_H1", "RES_H3")){
    var_index <- which(names(dat_t1) == var)
    dat_index <- dat_t1[,var_index] <=1

    model <- NULL
    if (var=="RES_B"){
      model <- glm(as.formula(paste(var, "~ vac_group + age_cat_5 + gender + HR + last_infection_cat_BB_2_neg +
                                  as.factor(ONSET_MONTH) + as.factor(year)")),
                   family="binomial",
                   data = dat_t1[dat_index,])
      out$B_base <- model
    }
    if (var=="RES_H1"){
      model <- glm(as.formula(paste(var, "~ vac_group + age_cat_5 + gender + HR + last_infection_cat_H1_2_neg +
                                  as.factor(ONSET_MONTH) + as.factor(year)")),
                   family="binomial",
                   data = dat_t1[dat_index,])
      out$H1_base <- model
    }
    if (var=="RES_H3"){
      model <- glm(as.formula(paste(var, "~ vac_group + age_cat_5 + gender + HR + last_infection_cat_H3_2_neg +
                                  as.factor(ONSET_MONTH) + as.factor(year)")),
                   family="binomial",
                   data = dat_t1[dat_index,])
      out$H3_base <- model
    }
  }
  return(out)
}

output_glm_mf_lastinf_heterotype <- get_output_glm_mf_lastinf_heterotype(
  msm_long %>% filter(time==1) %>% filter(!year %in% c("2008")),
  unique(msm_long %>% filter(!year %in% c("2008")) %>% 
           pull(year))
  )
