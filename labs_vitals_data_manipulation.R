
hrs = c('01', '03', '06', '12', '24')

## Want to retrieve data with multiple hours, but it is apparently too big to do in 1 shot and in R. 
## Have written a sql script (in dataManip dir) to do this in netezza. Trying to get this to write to a txt file.
# labs_vitals_all_hrs = fread(input = '../data/labs_vitals_all_hrs_output_from_netezza.csv'
#                             , sep = '|'
#                             , header = TRUE
#                             , na.strings = c('NA', ",,", "")
#                             , stringsAsFactors = FALSE
#                             , drop = exclude_these_strings
#                             , encoding = "UTF-8"
#                             , blank.lines.skip = TRUE                       # May not need this. There shouldn't be any blank lines.
#                             , key = toupper(c('encounterid', 'hour'))
#                             , showProgress = TRUE
#                             , logical01 = FALSE)
## Read in the csvs that were output from netezza.
## Then write each one as an R object.
fname_beginning = '../data/labs_vitals'
fname_end = paste(substr(date_lb, 2, 11), "_", substr(date_ub, 2, 11), sep = '')

output_from_netezza = '_output_from_netezza_' 
r_obj = '_netezza_to_R_'

 read_names = paste(fname_beginning, '_hr_', hrs, output_from_netezza, fname_end, '.csv', sep = "")
r_obj_names = paste(fname_beginning, '_hr_', hrs, r_obj, fname_end, '.Rdata', sep = "")


     labs_vitals_hr_01 =  fread(read_names[1])
save(labs_vitals_hr_01, file = r_obj_names[1])
     labs_vitals_hr_03 =  fread(read_names[2])  
save(labs_vitals_hr_03, file = r_obj_names[2])
     labs_vitals_hr_06 =  fread(read_names[3]) 
save(labs_vitals_hr_06, file = r_obj_names[3])
     labs_vitals_hr_12 =  fread(read_names[4]) 
save(labs_vitals_hr_12, file = r_obj_names[4])
     labs_vitals_hr_24 =  fread(read_names[5])
save(labs_vitals_hr_24, file = r_obj_names[5])
gc()


# system.time(load(paste(fname_beginning, '_hr_01', r_obj, fname_end, '.Rdata', sep = ""), verbose = TRUE), gcFirst = FALSE)

nms_I_want = names(labs_vitals_hr_01)
setnames(labs_vitals_hr_03, nms_I_want)
setnames(labs_vitals_hr_06, nms_I_want)
setnames(labs_vitals_hr_24, nms_I_want)



labs_vitals_all_hrs = rbindlist(list(labs_vitals_hr_01, labs_vitals_hr_06), use.names = TRUE)
#labs_vitals_all_hrs = rbindlist(list(labs_vitals_hr_01, labs_vitals_hr_03, labs_vitals_hr_06, labs_vitals_hr_12, labs_vitals_hr_24), use.names = TRUE)


dim(labs_vitals_all_hrs)
setnames(labs_vitals_all_hrs, names(labs_vitals_all_hrs), tolower(names(labs_vitals_all_hrs)))
setnames(labs_vitals_all_hrs, 'encounterid', 'encounter_id')

# Should include age_at_admit? Make sure to exclude when combining with other data which already includes age at admit.
exclude_these_strings = Cs(
    encntr_id
  , hour_start_ts
  , hour_end_ts
  , admit_ts
  , discharge_ts
  , birth_dt_tm
  , encounterid_daac
  , sepsis_narrow
  , sepsis_icd
  , sepsis_icd2
  , v7chfc
  , fact_clinical_encounter_unique_id)
# Remove all the time variables.
time_vars = c(grep('(_ts)|(dt_tm)$', names(labs_vitals_all_hrs), value = TRUE)
  , grep('time[12]$', names(labs_vitals_all_hrs), value = TRUE))
#grep('(_ts)|(dt_tm)$', names(sepsis_dev_dat), value = TRUE)
#grep('time[12]$', names(sepsis_dev_dat), value = TRUE)
labs_vitals_vars =
  setdiff(
    names(labs_vitals_all_hrs), 
    c(exclude_these_strings, time_vars))
labs_vitals_all_hrs = labs_vitals_all_hrs[ , ..labs_vitals_vars]



# Check to see if these are already numeric!!!
# Convert all labs and vitals to numeric
#labs_vitals_dat[ , (labs_vitals_vars) := lapply(.SD, as.numeric), .SDcols = labs_vitals_vars]

## "Impute" missing numeric values with an extreme value
labs_vitals_all_hrs[is.na(labs_vitals_all_hrs)] <- -9999
#repMiss =  function(x) ifelse(is.na(x), -9999, x)
#labs_vitals_dat[ , (labs_vitals_vars) := lapply(.SD, repMiss), .SDcols = labs_vitals_vars]


save(labs_vitals_all_hrs, file = paste(fname_beginning, '_all_hrs_', r_obj, fname_end, '.Rdata', sep = ""))
#rm(labs_vitals_hr_01, labs_vitals_hr_03, labs_vitals_hr_06, labs_vitals_hr_12, labs_vitals_hr_24)




























#######################################################################################
## Cut
## This takes quite a long time to pull down from netezza. Going to save and load.

uatConnection = odbcConnect("UAT", uid="joann_alvarez", pwd="ChangeMe2017", believeNRows=FALSE)
# This data should go back to Jan 2016. Oct 2017.


labs_vitals_query = paste(
"select 
  distinct
  *
from 
  clinical_ops..dataframe_wide
where 
  --hour = 1
  age_at_admit >= 18
  and discharge_ts >= duration_subtract(",
  date_lb, 
  ", 00000002::numeric(8,0))
  and discharge_ts <= duration_add(",
  date_ub,
  ", 00000002::numeric(8,0))",  
  sep = "")

ss = Sys.time()
#labs_vitals_dat <- data.table(sqlQuery(uatConnection, 
#  labs_vitals_query, 
#  as.is = TRUE))   
print(Sys.time() - ss)

#####################################################################
#####################################################################








#??results3 <- results3[, lapply(.SD, as.numeric), by=ENCOUNTERID]

#mainDat[!is.na(hour) & ip_encounters == FALSE, .N]

# About 13% of people in main dat did not have a row matching from dataframe_wide.
# Now only 7.7%!
# ('MOD:28359357', 'MOD:28360002', 'MOD:28360037', 'MOD:28360355', 'MOD:28360606')
# 'AHD:300007028''MOD:28360720' These did not have a match from dataframe_wide because their discharge date was after 2017-10-30 00000, but before Nov. 1st. Ryan's table includes everyone discharged during the day of the 31st. I will adjust my dates in the future. 
# The dates in my source are different from the dates in their source.
# so, the dates are different, and the dates in the dataframe_wide is after oct 31 0000.
# still don't know what the issue is for these other five.

#main_w_prev_enc[is.na(hour), .N]
# 15129
#main_w_prev_enc[is.na(hour) & discharge_ts > '2017-10-31', .N]
# 4170
# Of the people in my data that weren't in dataframe_wide, about 27% had dates after 2017-10-31.
# I still don't know why the others are being excluded.

# 1) discharge dates are slightly different between idm tables and dataframe_wide _Because the times in Ryan's data are in different time zones._
# 2) Lots of missing (13%) row in dataframe_wide _um, what_
# 3) Some of the missing rows (27%) are caused by the combination of difference in the discharge dates between the two tables combined with the fact that my date range ended at 10-31 time zero, but I don't know the reason the others were not in the dataframe_wide.
# Another reason is a difference in the time that dataframe wide includes. Jan 2016 - Oct 2017

#main_w_prev_enc[!is.na(hour) & discharge_ts == discharge_ts.labs_vitals_dat, .N]
#main_w_prev_enc[!is.na(hour), .N]
#main_w_prev_enc[!is.na(hour) & discharge_ts < discharge_ts.labs_vitals_dat, .N]
#[1] 118173
# None of them are equal. Almost all of the ones in Ryan's table are bigger. Maybe a few hours earlier.
