library(ggplot2)
library(data.table)

# check the PAM data
setwd("~/Documents/Bioinformatics/Projects/TURNT/PAM_data/")

# read in day 3 data (first good set, I think)
sample_key = fread("TURNT_sample_sheet - PAM.csv")
raw_pam  = fread("raw pam data - all_raw_pam.csv") 
reps = fread("TURNT_sample_sheet - replicates.csv")

# drop 24hr time point (incorrect calibration)
sample_key = sample_key[sample_key$dat != 2022-04-02,]

# drop leading 1s from colnames
names(raw_pam) <- gsub("[1:()/ ]", "", names(raw_pam))

# get curve summary information from raw PAM data
curve_list = lapply(1:nrow(sample_key),  function(i){
  
  pam_id = sample_key[i, pam_id]
  pam_curve_no = sample_key[i, curve_number] 
  pam_date = sample_key[i, date]

  
  # identify the start of a rapid light curve in the raw PAM data
  first_point = which(raw_pam$No. == pam_curve_no & raw_pam$Date == pam_date)
  
  # pull out Fv.fm at the first light curve data point
  # should be same for all points in the curve
  FvFm = raw_pam[first_point, FvFm]
  
  # pull out parameters from first rapid light curve regression generate PAM software
  curve_params = unlist(strsplit( raw_pam[first_point - 1, `F`], split = ":|,"))
  
  a_curve <- data.table(
    pam_id = pam_id[[1]],
    alpha = as.numeric(curve_params[3]),
    ETRm = as.numeric(curve_params[5]),
    Ik = as.numeric(curve_params[7]),
    FvFm = as.numeric(FvFm)
  )
  
  return(a_curve)
})

curve_dt = rbindlist(curve_list)
# merge
pam_dat = merge.data.table(curve_dt, sample_key, all.y = T, by = "pam_id")

pam_dat = merge(pam_dat, reps, all.x = T, by = "rep")


# plot all data
good_dat = pam_dat[quality %in% c("good", "okay") 
                   & date %in% as.Date(c("2022-04-04", "2022-04-05"))]

ggplot(data = good_dat, aes(x = treatment, y = ETRm, color = treatment))+
  geom_boxplot()+
  facet_wrap(~date)

ggplot(data = good_dat, aes(x = treatment, y = alpha, color = treatment))+
  geom_boxplot()+
  facet_wrap(~date)


ggplot(data = good_dat, aes(x = treatment, y = ETRm, color = treatment))+
  geom_point()+
  facet_wrap(~date)

ggplot(data = good_dat, aes(x = treatment, y = FvFm, color = treatment))+
  geom_boxplot()+
  facet_wrap(~date)


ggplot(data = good_dat, aes(x = date,
                            y = FvFm,
                            color = treatment,
                            group = individual,
                            label = individual))+
  geom_line()+
  geom_label()+
  facet_wrap(~treatment)

ggplot(data = good_dat, aes(x = date,
                            y = ETRm,
                            color = treatment,
                            group = individual,
                            label = individual))+
  geom_line()+
  geom_label()+
  facet_wrap(~treatment)

ggplot(data = good_dat, aes(x = date,
                            y = alpha,
                            color = treatment,
                            group = individual,
                            label = individual))+
  geom_line()+
  geom_label()+
  facet_wrap(~treatment)



mod <- lm(FvFm ~ treatment , good_dat)
anova(mod)

mod2<- lm(FvFm ~ treatment*date, good_dat)
anova(mod2)

mod3 <- lm(FvFm ~ treatment, good_dat[date == "2022-04-05"])  
anova(mod3)

mod4 <- lm(FvFm ~ treatment * date, good_dat[date == "2022-04-04"])
anova(mod4)
