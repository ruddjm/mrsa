options(stringsAsFactors = FALSE)
library(RODBC)
library(Hmisc)
library(data.table)
library(recipes)
library(xgboost)
library(caret)
library(pROC)
library(plyr)
library(dplyr)
library(ggplot2)
library(ROCR)

setwd("E:/MRSA/dataManipulation")

########################################################################################################
########################################################################################################
##
##
##  Prepare data to examine factors related to systemic MRSA infections. Will hopefully then make a predictive model.
##  
##
########################################################################################################
########################################################################################################
# Population: ED and inpatients
#########################################################################################################

makoConnection = odbcConnect("MAKO", uid="joann_alvarez", pwd="ChangeMe2017", believeNRows=FALSE)
uatConnection = odbcConnect("UAT", uid="joann_alvarez", pwd="ChangeMe2017", believeNRows=FALSE)

# For *training,* need to only get patients who have already been discharged and had a chance to get all dxs coded, so they have the outcome.
#
date_lb = "'2016-01-01 00:00:00'"
date_ub = "'2017-12-01 00:00:00'"
whereClause = paste(" (emergency_room_flg = 1 or ip_encounters = 1)
  and age >= 18
  and discharge_dt >= ",
  date_lb, 
  " and discharge_dt <= ",
  date_ub, 
  " and fact_encounter.source_system = 'DAAC'
  --and dmor..emrgncy_dept.outcome_location = 'Admit'
  -- and fact_diagnosis.source_system = 'DAAC'
  -- and admission_dt >= '2018-02-01 00:00:00'
  --and admission_dt <= '2018-02-29 00:00:00'
  --and diag_typ_cd = 'ICD10'
  --limit 100
  ", 
  sep = ' ')

# (powerchart_flg = 1
#   or powerchart_flg is NULL
#   or (powerchart_flg = 0 and contributor_system <> 'PowerChart'))

facility_filter = 
   " facility_type_cd_desc = 'Acute'
    and group_nm is not NULL
    and group_nm not in ('DIVESTED GROUP')
    and dim_facility.facility_cd not in ('DCR', 'DRA')"   # DCR is a children's hosp. DRA is a rehab hosp. Should I exclude RHB??



# mrsa_results_table_q = 
#   "select lower(result_name), lower(result_value), count(1) 
#   from idm..fact_clinical_result
#   where 
#   (powerchart_flg = 1
#     or powerchart_flg is NULL
#     or (powerchart_flg = 0 and contributor_system <> 'PowerChart'))
#   and result_ts > '2017-11-01 00:00:00'
#   group by lower(result_name), lower(result_value)
#   having lower(result_name) like ('%mrsa%')
#   order by lower(result_name), lower(result_value)
#   limit 100;"
#  
# 
# 
# mrsa_results_table <- data.table(sqlQuery(makoConnection, 
#   mrsa_results_table_q, 
#   as.is = TRUE))
# dim(mrsa_results_table)
# setnames(mrsa_results_table, names(mrsa_results_table), tolower(names(mrsa_results_table)))


mainQuery = paste(
"select 
  distinct
  ---------------------------------------------------------------
  -- Data to select patient population
  fact_encounter.facility_cd
  , fact_encounter.FACILITY_CD||':'||fact_encounter.PATIENT_ACCOUNT_NBR as ENCOUNTER_ID
  , fact_encounter.FACILITY_CD||':'||dim_patient.medical_record_nbr as patient_id -- Can use this to join this to previous encounter queries.
  , admission_dt + admission_tm as admit_ts
  , ed_arrival_ts
  , lower(ed_arrival_method) as ed_arrival_method
  --, case when emergency_room_flg = 0 and ip_encounters = 1  then  admit_ts else ed_arrival_ts   end as presentation_time
  , case when ed_arrival_ts is not null then min(admit_ts, ed_arrival_ts) else admit_ts end as presentation_time  

  , discharge_dt + discharge_tm as discharge_ts -- some 0001-01-01 00:00:00
  , case when discharge_ts = '0001-01-01 00:00:00' then null else discharge_ts end as discharge_ts_corrected
  -- los so far:
  , case when discharge_ts = '0001-01-01 00:00:00' then current_date - admission_dt else discharge_dt - admission_dt end as los_so_far 
  , emergency_room_flg
  , encounter_cd_desc -- From fact_encounter
  , encntr_class   -- from fact_clinical_encounter
  , encntr_type
  , ip_encounters
  , icu_flg   -- This variable not available real time, but can use it to look at trends
  --, trauma_flg 
  --, acute_flg  
  --, ip_days  
  --------------------------------------------------------------
  -- Features
  , round(months_between (admit_ts, date_of_birth)/12, 1) as age --(used as both filter var and feature)

  , ltrim(lower(derived_primary_visit_reason)) as derived_primary_visit_reason  -- fact_clinical_encounter
  , ltrim(lower(reason_for_visit)) as reason_for_visit -- fact_clinical_encounter  

  , sex_cd
  , rtrim(lower(race_desc)) as race_desc
  , rtrim(lower(ethnicity_desc)) as ethnicity_desc
  --, lower(marital_status_cd) as marital_status_cd
  
  , case when lower(marital_status_desc) = 'single' then 'single' else 
    case when lower(marital_status_desc) in ('married', 'life partner', 'significant other') then 'married' else 
    case when lower(marital_status_desc) = 'divorced' then 'divorced' else  
    case when lower(marital_status_desc) like '%widow%' then 'widow' else 
    case when lower(marital_status_desc) = 'legally separated' then 'separated' else 
      NULL  
      end end end end end
      as marital_status
  
  , lower(financial_class) as financial_class -- from fact_clinical_encounter. Needs to be grouped. 89 categories
  , case when vip is NULL then 0 else 1 end as vip
  , srg9_desc
  , do_not_resuscitate_ind -- missing for half
  , admit_type
  , longitude
  , latitude
  , regn_nm
  , licensed_bed_count
  , sub_mrkt_nm
  , case 
      when lower(ub_admit_source_desc) in ('transfer from hospital', 'transfer from a hospital') 
      then 'hospital transfer'
      else 
        case 
          when lower(ub_admit_source_desc) in ('information not available') 
          then NULL
          else lower(ub_admit_source_desc) 
          end
      end  
      as ub_admit_source_desc --fact_encounter  
      
  , case 
      when lower(ub_admit_source_desc) in 
        ('transfer from a snf inpatient', 'hospital transfer', 'transfer from hospital', 'transfer from a hospital', 'transfer from hospice', 'transfer from another health care facility') 
      then 'healthcare contacts'
      else 
        case 
          when lower(ub_admit_source_desc) in ('information not available') 
          then NULL
          else lower(ub_admit_source_desc) 
          end
      end  
      as ub_admit_source_desc_cat --fact_encounter        
      
      
  , case 
      when lower(ub_admit_source_desc) = 'transfer from a snf inpatient'
      then 1
      else 0
      end 
      as snf
  --, admit_src -- fact_clinical_encounter
  --, ub_discharge_status_desc as discharge_status_fact_encounter
  , fact_clinical_encounter.med_service as med_service_fce -- What is this capturing? I think we want to capture where is the patient now
  ----------------------------------------------------------------
  -- Outcomes/labels  
  , case when mrsa_sepsis = 'A4102' then 1 else 0 end as mrsa_sepsis
  , case when mrsa_pneum = 'J15212' then 1 else 0 end as mrsa_pneum
  , case when mrsa_unspec = 'A4902' then 1 else 0 end as mrsa_unspec
  , case when mrsa_cause_other_disease = 'B9562' then 1 else 0 end as mrsa_cause_other_disease
  , case when mrsa_carrier = 'Z22322' then 1 else 0 end as mrsa_carrier
  , mrsa_sepsis_poa
  , mrsa_pneum_poa
  , mrsa_unspec_poa
  , mrsa_cod_poa
  , mrsa_carrier_poa    

  , diag_cd as primary_dx_cd
  , diag_desc as primary_dx_desc

  ----------------------------------------------------------------

from idm..fact_encounter  
  --left join dmor..emrgncy_dept
  --  on fact_encounter.facility_cd = emrgncy_dept.hosp_cd
  --  and fact_encounter.patient_account_nbr = emrgncy_dept.pat_acct_id
  left join idm..fact_clinical_encounter
    on fact_clinical_encounter.facility_cd = fact_encounter.facility_cd
    and fact_clinical_encounter.patient_account_nbr = fact_encounter.patient_account_nbr

  left join idm..dim_diag
    on fact_encounter.dim_primary_icd10_diag_sk = dim_diag.dim_diag_sk

  left join idm..dim_encounter_type
    on dim_encounter_type.dim_encounter_type_sk = fact_encounter.dim_encounter_type_sk 
  left join idm..dim_patient
    on dim_patient.dim_patient_sk = fact_encounter.dim_patient_sk
  left join idm..dim_facility
    on dim_facility.facility_cd = fact_encounter.facility_cd
  left join idm..ref_admit_source 
    on fact_encounter.admit_source_cd = ref_admit_source.ub_admit_source_cd
  left join idm..dim_std_rev_grp
    on dim_std_rev_grp.dim_std_rev_grp_sk = fact_encounter.dim_std_rev_grp_sk

  left join (
    select * from (
      select lvl1_hosp_cd
        , patient_account_nbr
        , arr_dt as ed_arrival_ts
        , arr_method as ed_arrival_method
        , row_number() over (partition by lvl1_hosp_cd, patient_account_nbr order by arr_dt asc nulls last) as arr_dt_rank
      from dmor..emrgncy_dept
        join dmor..dim_fac using (dim_fac_key)
      where 
        ed_arrival_ts is not NULL
      ) f where arr_dt_rank = 1
    ) emergency_dept
  on fact_encounter.facility_cd = emergency_dept.lvl1_hosp_cd
    and fact_encounter.patient_account_nbr = emergency_dept.patient_account_nbr

  left join 
     (select
        diag_cd as mrsa_sepsis
        , present_on_admit as mrsa_sepsis_poa
        , facility_cd
        , patient_account_nbr
      from idm..fact_diagnosis
        left join idm..dim_diag
        on fact_diagnosis.dim_diag_sk = dim_diag.dim_diag_sk
      where diag_cd = 'A4102'
        and diag_typ_cd = 'ICD10') fact_diag_mrsa_sepsis
    on fact_diag_mrsa_sepsis.facility_cd = fact_encounter.facility_cd
      and fact_diag_mrsa_sepsis.patient_account_nbr = fact_encounter.patient_account_nbr
      
  left join 
     (select
        diag_cd as mrsa_pneum
        , present_on_admit as mrsa_pneum_poa
        , facility_cd
        , patient_account_nbr
      from idm..fact_diagnosis
        left join idm..dim_diag
        on fact_diagnosis.dim_diag_sk = dim_diag.dim_diag_sk
      where diag_cd = 'J15212'
        and diag_typ_cd = 'ICD10') fact_diag_mrsa_pneum
    on fact_diag_mrsa_pneum.facility_cd = fact_encounter.facility_cd
      and fact_diag_mrsa_pneum.patient_account_nbr = fact_encounter.patient_account_nbr      
    
  left join 
     (select
        diag_cd as mrsa_unspec
        , present_on_admit as mrsa_unspec_poa        
        , facility_cd
        , patient_account_nbr
      from idm..fact_diagnosis
        left join idm..dim_diag
        on fact_diagnosis.dim_diag_sk = dim_diag.dim_diag_sk
      where diag_cd = 'A4902'
        and diag_typ_cd = 'ICD10') fact_diag_mrsa_unspec
    on fact_diag_mrsa_unspec.facility_cd = fact_encounter.facility_cd
      and fact_diag_mrsa_unspec.patient_account_nbr = fact_encounter.patient_account_nbr      
    
  left join 
     (select
        diag_cd as mrsa_cause_other_disease
        , present_on_admit as mrsa_cod_poa        
        , facility_cd
        , patient_account_nbr
      from idm..fact_diagnosis
        left join idm..dim_diag
        on fact_diagnosis.dim_diag_sk = dim_diag.dim_diag_sk
      where diag_cd = 'B9562'
        and diag_typ_cd = 'ICD10') fact_diag_mrsa_cause_other_disease
    on fact_diag_mrsa_cause_other_disease.facility_cd = fact_encounter.facility_cd
      and fact_diag_mrsa_cause_other_disease.patient_account_nbr = fact_encounter.patient_account_nbr     
    
  left join 
     (select
        diag_cd as mrsa_carrier
        , present_on_admit as mrsa_carrier_poa        
        , facility_cd
        , patient_account_nbr
      from idm..fact_diagnosis
        left join idm..dim_diag
        on fact_diagnosis.dim_diag_sk = dim_diag.dim_diag_sk
      where diag_cd = 'Z22322'
        and diag_typ_cd = 'ICD10') fact_diag_mrsa_carrier
    on fact_diag_mrsa_carrier.facility_cd = fact_encounter.facility_cd
      and fact_diag_mrsa_carrier.patient_account_nbr = fact_encounter.patient_account_nbr  
      
where",
  whereClause, 
   " and ",
   facility_filter,
   ";")

ss = Sys.time()
mainDatInit <- data.table(sqlQuery(makoConnection, 
  mainQuery, 
  as.is = TRUE),
  #na.strings = c(""),
  key = 'ENCOUNTER_ID')
Sys.time() - ss
mainDat = copy(mainDatInit)
dim(mainDat)
setnames(mainDat, names(mainDat), tolower(names(mainDat)))
mainDat = mainDat[!duplicated(encounter_id)]
# Maindat has distinct encounter_ids but duplicated patients (probably)


mainDat[ , (names(mainDat)) := lapply(.SD, function(x) {x[x %in% c('', 'declined', 'declined to specify', 'unknown', 'Information Not Available')] <- NA; x})]
dim(mainDat)
# 148,000 for jan 2018
# 548,000 for 4 months.
mainDat[ , age := as.numeric(age)]
mainDat[ , licensed_bed_count := as.numeric(licensed_bed_count)]
mainDat[tolower(sex_cd) == "u", sex_cd := NA]

# Make an interaction term between ed flag and ip flag
mainDat[ip_encounters == 1 & emergency_room_flg == 0, encounter_ed_ip := 'ip_not_through_ed']
mainDat[ip_encounters == 1 & emergency_room_flg == 1, encounter_ed_ip := 'ed_admitted']
mainDat[ip_encounters == 0 & emergency_room_flg == 1, encounter_ed_ip := 'ed_not_admitted']

# Clean financial_class
# There were ~100 different values of financial_class, so I grouped them in the following way.
# I got a table combinations in netezza of financial_class by srg9_desc using dim_std_rev_grp, which has 9 categories.
# I also made a few different categories.
# I output the netezza data to a csv and manually changed some of the categories. Then I saved the new file with a different name to prevent overwriting. These files are in data_documentation.

# Are there major different goups in self pay/uninsured
# For people who had fin class = NA, but were then resolved to medicare managed care: Depends on why they started with a missing value. Is there more info in the fact that they started with a missing value, or that they ended up medicaid managed care? ## Using missing category, since, if you look at ALL patients who started with fin class = NULL they end up with an array of srg9_desc. 
fin_class_mapping = data.table(read.csv('../data_documentation/fin_class_grouping_annotated.csv',
                                        na.strings = c(NA, '', 'NULL')))
setnames(fin_class_mapping, 'count', 'countt')
fin_class_mapping[ , fin_class_cat := tolower(srg9_desc)]
fin_class_mapping[!is.na(grouping), fin_class_cat := grouping]
#fin_class_mapping[ , sum(countt), by = fin_class_cat][order(-V1)]

mainDat[ , financial_class_cat := fin_class_mapping$fin_class_cat[match(financial_class, fin_class_mapping$fin_class)]]
mainDat[!is.na(financial_class) & is.na(financial_class_cat), financial_class_cat := 'other']



ed_arrival_mapping = data.table(read.csv('../data_documentation/ed_arrival_grouping.csv',
                                        na.strings = c(NA, '', 'NULL')))
# This is just so data.table doesnt think I want the count() *function*, but rather the column name.
setnames(ed_arrival_mapping, 'count', 'countt')


#mainDat[ , ed_arrival_cat := ed_arrival_mapping$ed_arrival_cat[match(ed_arrival_method, ed_arrival_mapping$ed_arrival_cat)]]
mainDat[!is.na(ed_arrival_method), ed_arrival_ambulatory := ed_arrival_mapping$ambulatory[match(ed_arrival_method, ed_arrival_mapping$method)]]
mainDat[!is.na(ed_arrival_method), ed_arrival_fire := ed_arrival_mapping$fire[match(ed_arrival_method, ed_arrival_mapping$method)]]
mainDat[!is.na(ed_arrival_method), ed_arrival_ambulance := ed_arrival_mapping$ambulance[match(ed_arrival_method, ed_arrival_mapping$method)]]
mainDat[!is.na(ed_arrival_method), ed_arrival_assisted  := ed_arrival_mapping$assisted[match(ed_arrival_method, ed_arrival_mapping$method)]]
mainDat[!is.na(ed_arrival_method), ed_arrival_air       := ed_arrival_mapping$air[match(ed_arrival_method, ed_arrival_mapping$method)]]
mainDat[!is.na(ed_arrival_method), ed_arrival_stretcher := ed_arrival_mapping$stretcher[match(ed_arrival_method, ed_arrival_mapping$method)]]
#describe(mainDat[ed_arrival_stretcher == 1 , .(ed_arrival_method, ed_arrival_ambulatory, ed_arrival_ambulance, ed_arrival_assisted, ed_arrival_fire, ed_arrival_air, ambulance_or_fire_not_walking)])
mainDat[ , ambulance_or_fire_not_walking := as.numeric(ed_arrival_stretcher | (ed_arrival_ambulatory == 0 & (ed_arrival_fire == 1 | ed_arrival_ambulance == 1)))]
mainDat[is.na(ambulance_or_fire_not_walking), ambulance_or_fire_not_walking := -9999]



mainDat[ , any_systemic_mrsa := as.numeric(mrsa_sepsis == 1 | mrsa_pneum == 1 | mrsa_cause_other_disease == 1 | mrsa_unspec == 1)]
mainDat[ , mrsa_fac := factor(any_systemic_mrsa, levels = c(1, 0), labels = c('yes', 'no'))]

mainDat[race_desc %in% c('native american', 'pacific islander'), race_desc := 'other']


mainDat[ , race := race_desc]
mainDat[is.na(race_desc) & ethnicity_desc == 'asian/oriental', race := 'asian']
mainDat[is.na(race_desc) & ethnicity_desc == 'african american/black', race := 'black']
mainDat[is.na(race_desc) & ethnicity_desc == 'white/caucasian', race := 'white']
mainDat[(is.na(race_desc) | race %in% c('other', 'unknown', 'white')) & ethnicity_desc == 'hispanic', race := 'hispanic']
mainDat[ , nonwhite := factor(race != 'white', levels = c(FALSE, TRUE), labels = c('White', 'Non-white'))]


#mainDat[ , race := interaction(race_desc, ethnicity_desc)]

# See 2018-06-12 email and dataDocumentation notes 
mainDat[ , do_not_resuscitate := factor(do_not_resuscitate_ind, levels = c(1, 2), labels = c('Yes', 'No'))]

#######################
# Strings for different kinds of catheters/ports
#   central line: central line and insert  
#    picc: picc and insert
#    cvc insert and cvc, insert and cv cath
#    urinary cath: urin and cath and insert, uret and cath and insert
#    ? foley cath foley + cath + insert or 'foley catheter'
#    port + cath + insert

# These are orders from the *same* encounter as the encounter in mainDat
# There should be one row for every encounter, and encounters with no orders for central line, etc, have a missing value which must be populated with 0.
ordersQuery = paste(
"
  select 
  distinct
    fact_encounter.facility_cd||':'||fact_encounter.patient_account_nbr as encounter_id
    , fact_encounter.facility_cd
    , admission_dt + admission_tm as admit_ts
    , discharge_dt + discharge_tm as discharge_ts -- some 0001-01-01 00:00:00
    , case when discharge_ts = '0001-01-01 00:00:00' then null else discharge_ts end as discharge_ts_corrected
    , round(months_between (admit_ts, date_of_birth)/12, 1) as age
    , urinary_cath
    , urinary_cath_ts
    , central_line
    , central_line_ts
  from idm..fact_encounter  
    left join idm..dim_patient
      on dim_patient.dim_patient_sk = fact_encounter.dim_patient_sk
    left join idm..dim_facility
      on dim_facility.facility_cd = fact_encounter.facility_cd

        -- urinary cath orders
        left join 
          (select
            distinct
              facility_cd||':'||patient_account_nbr as encounter_id
            , 'urinary_cath' as urinary_cath
            , min(order_written_ts) as urinary_cath_ts
          from idm..fact_clinical_order
          where (
          (((lower(order_name) like '%urin%' or lower(order_name) like '%uret%' or lower(order_name) like '%foley%')
          and lower(order_name) like '%cath%'
          and lower(order_name) like '%insert%')
          or lower(order_name) = 'foley catheterization')
          and lower(status) in ('completed', 'comp')
          )
          and 
            (powerchart_flg = 1
             or powerchart_flg is NULL)
          group by encounter_id
          ) fact_urinary_cath
          on fact_urinary_cath.encounter_id = fact_encounter.facility_cd||':'||fact_encounter.patient_account_nbr

        -- catheter orders
        left join 
          (select
            distinct
              facility_cd||':'||patient_account_nbr as encounter_id
            , 'central_line' as central_line
            , min(order_written_ts) as central_line_ts
          from idm..fact_clinical_order
          where           
            (
              (lower(order_name) like '%central line%' 
            or lower(order_name) like '%picc%' 
            or lower(order_name) like '%cvc%' 
            or (lower(order_name) like '%cv%' and lower(order_name) like '%cath%')
            or (lower(order_name) like '%port%' and lower(order_name) like '%cath%') 
            )
            and lower(order_name) like '%insert%')
            and lower(status) in ('completed', 'comp')
            and (powerchart_flg = 1
               or powerchart_flg is NULL)
          group by encounter_id
          ) fact_central_line
          on fact_central_line.encounter_id = fact_encounter.facility_cd||':'||fact_encounter.patient_account_nbr

  where",
    whereClause, 
     " and ",
     facility_filter,
    ";") 

ss = Sys.time()
ordersDat <- data.table(sqlQuery(makoConnection, 
  ordersQuery, 
  as.is = TRUE),
  #na.strings = c(""),
  key = 'ENCOUNTER_ID')
Sys.time() - ss
dim(ordersDat)
setnames(ordersDat, names(ordersDat), tolower(names(ordersDat)))
# These are orders from the *same* encounter
# I think this one should have one row for every encounter, and encounters with no orders for central line, etc, would have a missing value which must be populated with 0

# Merge!
main_w_orders = merge(
  mainDat,
  ordersDat[ , c("encounter_id", setdiff(names(ordersDat), names(mainDat))), with = FALSE],
  by.x = 'encounter_id', 
  by.y = 'encounter_id',
  all.x = TRUE,
  all.y = FALSE
  #, suffixes = c('', '.labs_vitals_dat')
  )

main_w_orders[is.na(central_line), central_line := '0']
main_w_orders[is.na(urinary_cath), urinary_cath := '0']
###########################################################
##########################################################

# The index_encounter table should have a row for every encounter in mainDat. 
# Then, the prev_encounter table should have a row for many encounters, including index encounters, but specified that discharge time is not 01:01:00 (excluding some inhouse patients.)
# Finanlly, these two tables are joined by *patient*, and only those where the prev encounter presentation time is strictly less than the index encounter presentation_time are retained.
# So prev_encounters_dat has one row per prev encounter. Not all pats in mainDat had a prev encounter, and some patients had more than 1. 
# prev_hiv_dx hais either 1 or 0 and no missing.

prev_encounters_query =
paste("
with index_encounter as (
  ---------------------------------------------------------------
  ---------------------------------------------------------------
  select 
  distinct
    fact_encounter.facility_cd||':'||fact_encounter.patient_account_nbr as index_encounter_id
    , fact_encounter.facility_cd
    , medical_record_nbr as mrn
    , admission_dt + admission_tm as admit_ts
    , case when ed_arrival_ts is not null then min(admit_ts, ed_arrival_ts) else admit_ts end as index_presentation_time  
  --, case when     emergency_room_flg = 0 and ip_encounters = 1  then  admit_ts else ed_arrival_ts   end as index_presentation_time
  --, discharge_dt + discharge_tm as index_encounter_discharge_ts -- some 0001-01-01 00:00:00
    , round(months_between (admit_ts, date_of_birth)/12, 1) as age
  from idm..fact_encounter  
    left join idm..dim_patient
      on dim_patient.dim_patient_sk = fact_encounter.dim_patient_sk
    left join idm..dim_facility
      on dim_facility.facility_cd = fact_encounter.facility_cd
    left join (
      select * from (
        select lvl1_hosp_cd
          , patient_account_nbr
          , arr_dt as ed_arrival_ts
          , arr_method as ed_arrival_method
          , row_number() over (partition by lvl1_hosp_cd, patient_account_nbr order by arr_dt asc nulls last) as arr_dt_rank
        from dmor..emrgncy_dept
          join dmor..dim_fac using (dim_fac_key)
        where 
          ed_arrival_ts is not NULL
        ) f where arr_dt_rank = 1
      ) emergency_dept
      on fact_encounter.facility_cd = emergency_dept.lvl1_hosp_cd
      and fact_encounter.patient_account_nbr = emergency_dept.patient_account_nbr
  where",
    whereClause, 
     " and ",
     facility_filter,
    "), 
  ---------------------------------------------------------------
  ---------------------------------------------------------------
prev_encounter as (
      select 
      distinct
        fact_encounter.facility_cd||':'||fact_encounter.patient_account_nbr as prev_encounter_id
        , fact_encounter.facility_cd
        , medical_record_nbr as mrn
        , discharge_dt + discharge_tm as prev_encounter_discharge_ts
        --, row_number() over (partition by prev_encounter_id order by prev_encounter_discharge_ts asc nulls last) as prev_encounter_discharge_ts -- Can use this to select the most recent discharge
        , length_of_stay as prev_encounter_los
        , emergency_room_flg as prev_emergency_room_flg
        , ip_encounters as prev_ip_encounters
        , icu_flg as prev_icu_flg
        , case when prev_hiv_dx = 'B20' then 1 else 0 end as prev_hiv_dx
        , case when regexp_like(prev_transplant_dx, '^Z94') then 1 else 0 end as prev_transplant_dx  
        , case when regexp_like(prev_dep_on_ren_dialysis_dx, '^Z992') then 1 else 0 end as prev_dep_on_ren_dialysis_dx        
      
      from idm..fact_encounter
        join idm..dim_patient using (dim_patient_sk)
        
        -- hiv
        left join 
          (select
            distinct
            diag_cd as prev_hiv_dx
            , facility_cd
            , patient_account_nbr
          from idm..fact_diagnosis
            left join idm..dim_diag
              on fact_diagnosis.dim_diag_sk = dim_diag.dim_diag_sk
          where (diag_cd = 'B20' or diag_cd = 'Z21')
            and diag_typ_cd = 'ICD10') fact_hiv
          on fact_hiv.facility_cd = fact_encounter.facility_cd
          and fact_hiv.patient_account_nbr = fact_encounter.patient_account_nbr

        -- transplant
        left join 
          (select
            distinct
            diag_cd as prev_transplant_dx
            , facility_cd
            , patient_account_nbr
          from idm..fact_diagnosis
            left join idm..dim_diag
              on fact_diagnosis.dim_diag_sk = dim_diag.dim_diag_sk
          where staging..regexp_like(diag_cd,  '^Z94')
            and diag_typ_cd = 'ICD10') fact_transplant
          on fact_transplant.facility_cd = fact_encounter.facility_cd
          and fact_transplant.patient_account_nbr = fact_encounter.patient_account_nbr   
          
        -- dialysis
        left join 
          (select
            distinct
            diag_cd as prev_dep_on_ren_dialysis_dx
            , facility_cd
            , patient_account_nbr
          from idm..fact_diagnosis
            left join idm..dim_diag
              on fact_diagnosis.dim_diag_sk = dim_diag.dim_diag_sk
          where staging..regexp_like(diag_cd,  '^Z992')
            and diag_typ_cd = 'ICD10') fact_dialysis
          on fact_dialysis.facility_cd = fact_encounter.facility_cd
          and fact_dialysis.patient_account_nbr = fact_encounter.patient_account_nbr  
        
      where 
        discharge_dt <> '0001-01-01'  -- If the discharge date equals 0001-01-01, that means theyre still in the hospital, and then it wouldnt be a previous encounter!
        and (ip_encounters = 1 or emergency_room_flg = 1)
        ) 
  ---------------------------------------------------------------
  ---------------------------------------------------------------
select 
  distinct
  index_encounter_id
  , index_encounter.facility_cd||':'||index_encounter.mrn as patient_id
  , index_encounter.facility_cd
  , index_presentation_time
  , index_encounter.mrn

  , prev_encounter_id
  , prev_encounter_discharge_ts
  , prev_encounter_los
  , prev_ip_encounters
  , prev_emergency_room_flg
  , prev_icu_flg
  , prev_hiv_dx
  , prev_transplant_dx
  , prev_dep_on_ren_dialysis_dx
  --, index_presentation_time - prev_encounter_discharge_ts as time_prev_encounter_to_index_encounter
  , days_between(prev_encounter_discharge_ts , index_presentation_time) as time_prev_encounter_to_index_encounter
  
from index_encounter
  left join prev_encounter 
    on index_encounter.facility_cd = prev_encounter.facility_cd
    and index_encounter.mrn = prev_encounter.mrn
where
  prev_encounter_discharge_ts < index_presentation_time
;")
ss = Sys.time()
prev_encounters_dat <- data.table(sqlQuery(makoConnection, 
  prev_encounters_query, 
  as.is = TRUE),
  #na.strings = c(""),
  key = 'INDEX_ENCOUNTER_ID')
Sys.time() - ss
dim(prev_encounters_dat)
# 13,364,829 (For 2016-01--2017-11-01)

setnames(prev_encounters_dat, names(prev_encounters_dat), tolower(names(prev_encounters_dat)))
prev_encounters_dat[ , (names(prev_encounters_dat)) := lapply(.SD, function(x) {x[x == ''] <- NA; x})]
archive_prev = copy(prev_encounters_dat)
gc()
#describe(prev_encounters_dat$index_encounter_id)
#prev_encounters_dat$index_encounter_id 
#       n  missing distinct 
# 8088293        0  1019916 

# Of about 1.6 million encounters in mainDat, about 1.0 million had at least one previous encounter. However, there were only 500,000 distinct *patients* made up these 1 million prev encounters
# WAIT>>> There are ~8 million rows in prev_encounters_dat, but only ~ 2 million distinct prev encounter ids in the table!!
# Some patients have multiple index encounters
#describe(prev_encounters_dat[patient_id == 'AHD:547140'])
# This pat:(AHD:547140) has five encounters in mainDat. Each one of those has multiple previous encounters.

#setorder(prev_encounters_dat, patient_id, -index_encounter_id, time_prev_encounter_to_index_encounter)
# any prev encounters with hiv
prev_encounters_dat[ , any_prev_hiv_dx := as.integer(any(prev_hiv_dx == 1)), by = index_encounter_id]
prev_encounters_dat[ , any_prev_dep_on_ren_dialysis_dx := as.integer(any(prev_dep_on_ren_dialysis_dx == 1)), by = index_encounter_id]
prev_encounters_dat[ , any_prev_transplant :=             as.integer(any(prev_transplant_dx == 1)), by = index_encounter_id]
prev_encounters_dat[ , prev_encounter_los := as.numeric(prev_encounter_los)]
#prev_encounters_dat[  , .(prev_encounter_discharge_ts, index_presentation_time, prev_encounter_los, time_prev_encounter_to_index_encounter)]
prev_encounters_dat[ , prev_encounter_within_year := as.numeric(time_prev_encounter_to_index_encounter <= 365)]

# Aggregations to get key features.
prev_encounters_dat[ , num_encounters_past_year := sum(prev_encounter_within_year ==  1), by = index_encounter_id]  # patient_id is facility, mrn.
suppressWarnings(
prev_encounters_dat[ , time_since_last_ed := as.integer(min(time_prev_encounter_to_index_encounter[prev_emergency_room_flg == 1])), by = index_encounter_id])
suppressWarnings(
prev_encounters_dat[ , time_since_last_ip := as.integer(min(time_prev_encounter_to_index_encounter[prev_ip_encounters == 1])), by = index_encounter_id])
suppressWarnings(
prev_encounters_dat[ , time_since_last_ic := as.integer(min(time_prev_encounter_to_index_encounter[prev_icu_flg == 1])), by = index_encounter_id])
#describe(prev_encounters_dat)
#"AHD:302879192"
# los of most recent IP visit.

setorder(prev_encounters_dat, patient_id, index_encounter_id, time_prev_encounter_to_index_encounter)
## This must be sorted to work!!
# For a given encounter in the denominator data set:
#   Take all the rows of prev encounters that were IP
#   Take the first one, and get it's length_of_stay. (They are already sorted by date above.)
prev_encounters_dat[ , last_ip_los := prev_encounter_los[prev_ip_encounters == 1][1], by = index_encounter_id]
#prev_encounters_dat[ , .(index_encounter_id, patient_id, time_prev_encounter_to_index_encounter, prev_ip_encounters,
#                         prev_encounter_los, last_ip_los)][300:330]

# Before there was one row for every prev encounter. Now there is one per index encounter
encounter_lev_prev_encounters = unique(prev_encounters_dat[ , .(index_encounter_id, num_encounters_past_year, time_since_last_ed, time_since_last_ip, time_since_last_ic, last_ip_los, any_prev_hiv_dx, any_prev_transplant, any_prev_dep_on_ren_dialysis_dx)])

main_w_prev_enc = merge(
  main_w_orders, 
  encounter_lev_prev_encounters,
  by.x = "encounter_id",
  by.y = "index_encounter_id", 
  all.x = TRUE) # This should work with both alls = TRUE (all = TRUE).

# For people who don't have any records for previous ED encounters, how to handle time_since_last_ed? The more recent the ED visit, the worse. Here, I'm assigning a very LARGE number if they did not have any previous ED visits.
suppressWarnings(
main_w_prev_enc[is.na(time_since_last_ed), time_since_last_ed := runif(.N, min = 2000, max = 4000)])
suppressWarnings(
main_w_prev_enc[is.na(time_since_last_ip), time_since_last_ip := runif(.N, min = 2000, max = 4000)])
suppressWarnings(
main_w_prev_enc[is.na(time_since_last_ic), time_since_last_ic := runif(.N, min = 2000, max = 4000)])
main_w_prev_enc[is.na(num_encounters_past_year), num_encounters_past_year := 0]
main_w_prev_enc[is.na(any_prev_hiv_dx), any_prev_hiv_dx := 0]
main_w_prev_enc[is.na(any_prev_transplant), any_prev_transplant := 0]
main_w_prev_enc[is.na(any_prev_dep_on_ren_dialysis_dx), any_prev_dep_on_ren_dialysis_dx := 0]
main_w_prev_enc[is.na(last_ip_los), last_ip_los := 0]



















































# In terms of determining whether someone is in trauma vs medical vs surgical, a good way to handle it may be current location. 


# take info from locations within the 1st 24 hours, and see (1) if they were in a trauma ward.
# For "current location," in the model development (train/test data), we can use current location at 24 hours after presentation. For production, we would use current time.
##!! Still need to take the last observation.
location_query = 
paste(
"select 
  distinct
  fact_encounter.facility_cd||':'||fact_encounter.patient_account_nbr as encounter_id
  , admission_dt + admission_tm as admit_ts
  , ed_arrival_ts
  , case when ed_arrival_ts is not null then min(admit_ts, ed_arrival_ts) else admit_ts end as presentation_time  
  , round(months_between (admit_ts, date_of_birth)/12, 1) as age
  , lower(three_hour_loc.medical_service) as location_at_3_hours      
  , three_hour_loc.critical_care_flg as in_critical_care_at_3_hours 
  
  , lower(first_location) as first_location
  , first_location_critical_care
  , first_location_start_ts
  
  from idm..fact_encounter
    left join idm..dim_facility
      on dim_facility.facility_cd = fact_encounter.facility_cd
    left join idm..dim_patient
      on dim_patient.dim_patient_sk = fact_encounter.dim_patient_sk      
    left join (
      select * from ( 
        select
          facility_cd
          , patient_account_nbr
          , start_ts as first_location_start_ts
          , lower(medical_service) as first_location
          , critical_care_flg as first_location_critical_care
          , row_number() over (partition by facility_cd, patient_account_nbr order by start_ts asc nulls last) as start_ts_rank
        from idm..fact_clinical_location
        
        ) fff where start_ts_rank = 1
        ) first_loc
      on fact_encounter.facility_cd = first_loc.facility_cd
      and fact_encounter.patient_account_nbr = first_loc.patient_account_nbr      
      
    left join idm..fact_clinical_location three_hour_loc
      on fact_encounter.facility_cd = three_hour_loc.facility_cd
      and fact_encounter.patient_account_nbr = three_hour_loc.patient_account_nbr  

    -- Get hospital, pat number, and ed arrival time
    left join (
      select * from (
        select 
            lvl1_hosp_cd
          , patient_account_nbr
          , arr_dt as ed_arrival_ts
          , row_number() over (partition by lvl1_hosp_cd, patient_account_nbr order by arr_dt asc nulls last) as arr_dt_rank
        from dmor..emrgncy_dept
          join dmor..dim_fac using (dim_fac_key)
        where 
          ed_arrival_ts is not NULL
        ) f where arr_dt_rank = 1
      ) emergency_dept
      on fact_encounter.facility_cd = emergency_dept.lvl1_hosp_cd
      and fact_encounter.patient_account_nbr = emergency_dept.patient_account_nbr      
      
  where
    
    ",
  whereClause, 
  " and ",
  facility_filter,
  " and (powerchart_flg = 1
  --or powerchart_flg is NULL
  --or (powerchart_flg = 0 and contributor_system <> 'PowerChart')
  )
 
  and three_hour_loc.start_ts < presentation_time + cast(3||' Hour' as interval)
  and three_hour_loc.end_ts > presentation_time + cast(3||' Hour' as interval)              
  ;")


ss = Sys.time()
location_dat <- data.table(sqlQuery(makoConnection, 
  location_query, 
  as.is = TRUE))
Sys.time() - ss
dim(location_dat)

setnames(location_dat, names(location_dat), tolower(names(location_dat)))

location_dat[ , trauma_ward_by_3_hours := as.numeric(grepl('trauma', first_location) | grepl('trauma', location_at_3_hours))]
location_dat[ , ob_location_by_3_hours := as.numeric(grepl('obstetrics|labor', first_location) | grepl('obstetrics|labor', location_at_3_hours))]
location_dat[ , surg_location_by_3_hours := as.numeric(grepl('surgical|cardiac surg|surg/inv diag|surg/diag obs', first_location) | grepl('surgical|cardiac surg|surg/inv diag|surg/diag obs', location_at_3_hours))]
location_dat[ , rehab_by_3_hours := as.numeric(grepl('rehab', first_location) | grepl('rehab', location_at_3_hours))]



main_w_prev_enc = merge(main_w_prev_enc
      , location_dat[!duplicated(encounter_id), .(encounter_id, first_location, location_at_3_hours, in_critical_care_at_3_hours, trauma_ward_by_3_hours, rehab_by_3_hours, surg_location_by_3_hours, ob_location_by_3_hours)]
      , by = "encounter_id"
      , all.x = TRUE
      , all.y = FALSE)

main_w_prev_enc[is.na(in_critical_care_at_3_hours), in_critical_care_at_3_hours := -9999]
main_w_prev_enc[is.na(trauma_ward_by_3_hours), trauma_ward_by_3_hours := -9999]
main_w_prev_enc[is.na(rehab_by_3_hours), rehab_by_3_hours := -9999]
main_w_prev_enc[is.na(surg_location_by_3_hours), surg_location_by_3_hours := -9999]
main_w_prev_enc[is.na(ob_location_by_3_hours), ob_location_by_3_hours := -9999]


# There are about 700 extra rows that I wasn't expecting. I expected unique rows.
# Could we use oncology as a location?
# 95 - SNF/Extended Care n = 
# What is
# non-obsv op in bed






fnamee = paste('../data/saved_R_workspaces_or_objects/key_datatables_end_of_dataManipulation_', substr(date_lb, 2, 11), "_", substr(date_ub, 2, 11), '.Rdata', sep = "")

save(mainDat, date_lb, date_ub, encounter_lev_prev_encounters, main_w_orders, main_w_prev_enc, ordersDat, location_dat, file = fnamee)






































# Medication info
query_string <-"select DISTINCT A.FACILITY_CD, A.PATIENT_ACCOUNT_NBR, A.ORDERING_PHYSICIAN, A.MEDICATION_NAME, A.DETAILS,

IDM..DIM_DRG.DRG_CD,IDM..DIM_APR_DRG.DRG_CD, SEVERITY_CD, DIAG_CD, ORDER_ID,HEALTH_SYSTEM_SOURCE_ID
FROM

(select idm..FACT_DISCHARGE_MEDICATION.*, DISCHARGE_DT, CERNER_ORDERS_HIST.ORDER_ID,  
CERNER_ORDERS_HIST.HEALTH_SYSTEM_SOURCE_ID, DIM_PRIMARY_ICD10_DIAG_SK, DIM_DRG_SK, DIM_APR_DRG_SK
from idm..FACT_DISCHARGE_MEDICATION
inner join staging..CERNER_ORDERS_HIST
on CERNER_ORDERS_HIST.HEALTH_SYSTEM_SOURCE_ID||'~'|| CERNER_ORDERS_HIST.ORDER_ID=FACT_DISCHARGE_MEDICATION.unique_id

inner join IDM..DIM_MEDICATION
on CERNER_ORDERS_HIST.HEALTH_SYSTEM_SOURCE_ID||'~'||CERNER_ORDERS_HIST.CATALOG_CD = DIM_MEDICATION.MEDICATION_CD

inner join IDM..REF_MEDICATION_THEREPAUTIC_CLASS 
on REF_MEDICATION_THEREPAUTIC_CLASS.MEDICATION_CD = DIM_MEDICATION.MEDICATION_CD

where upper(TOP_LEVEL_THERAPEUTIC_CLASS) =  'ANALGESICS'
and DISCHARGE_DT BETWEEN '20171001' AND '20171231'
and STATUS IN ('Completed','Ordered')) A,

IDM..DIM_DIAG,"


# Example of geting 24 hours worth of data.
results1 <- sqlQuery(channel, "SELECT DISTINCT A.FACILITY_DISP, A.FIN_NBR, A.FACILITY_DISP||':'||A.FIN_NBR as ENCOUNTERID
                    FROM RSTUDIES..DATA_SEPSIS_LAB_RESULTS2017 A, 
                    RSTUDIES.DATA_SEPSIS2017 B
                    where A.FACILITY_DISP=B.FACILITY_DISP AND
                    A.FIN_NBR=B.FIN_NBR AND
                    B.ADMIT_AGE>=18 AND
                    A.COLLECTION_DT_TM < B.ADMIT_DT_TM + cast(24||' Hour' as interval) AND 
                    A.COLLECTION_DT_TM > B.ADMIT_DT_TM - cast(24||' Hour' as interval);" , as.is=TRUE)



mainDat[ , age := as.numeric(age)]
#mainDat[ , length_of_stay := as.numeric(length_of_stay)]
#mainDat[ , ip_days := as.numeric(ip_days)]

mainDat[ , ageCat := factor(cut(age, breaks = seq(0, 110, 10), right = FALSE), ordered = TRUE)]
label(mainDat$age) = 'Age'
# calculate who has each of the five codes




mrsa_conditions = Cs(mrsa_carrier, mrsa_sepsis, mrsa_pneum, mrsa_unspec, mrsa_cause_other_disease, any_systemic_mrsa)

label(mainDat$mrsa_sepsis) = "A41.02: Sepsis due to MRSA"
label(mainDat$mrsa_pneum) = "J15.212: Pneumonia due to MRSA"
label(mainDat$mrsa_unspec) = "A49.02: MRSA infection of unspecified site"
label(mainDat$mrsa_cause_other_disease) = "B95.62: MRSA infection causing disease classified elsewhere"
label(mainDat$mrsa_carrier) = "Z22.322: Carrier or suspected MRSA carrier"
label(mainDat$any_systemic_mrsa) = "Any systemic MRSA infection"


#A41.02 MRSA sepsis
#J15.212 MRSA pneumonia
#A49.02 MRSA infection of an unspecified site
#B95.62 Methicillin resistant S. aureus (MRSA) infection as the cause of diseases classified elsewhere
#Z22.322, Carrier or suspected carrier of MRSA



mrsa_rates_by_age = mainDat[ , .(.N, mrsa_rate = mean(any_systemic_mrsa == 1)), 
                    by = ageCat]
mrsa_rates_by_age[ , lower := mrsa_rate - qnorm(0.975)*sqrt(mrsa_rate*(1 - mrsa_rate)/N)]
mrsa_rates_by_age[ , upper := mrsa_rate + qnorm(0.975)*sqrt(mrsa_rate*(1 - mrsa_rate)/N)]


# Find out how mrsa many dxs people have
systemic_mrsa_icds = Cs(mrsa_sepsis, mrsa_pneum, mrsa_unspec, mrsa_cause_other_disease)
mainDat[ , num_sys_mrsa_icds := rowSums(.SD), .SDcols = systemic_mrsa_icds]











if(FALSE){
# POA codes: 
# y = yes
# E = yes
# n = no
# u = missing due to insufficient documentation
# w = missing due to clinically not determinable
# 1 = missing because 'Unreported/Not used.  Exempt from POA reporting.  This code is equivalent to a blank on the UB-04, however; it was determined that blanks are undesirable when submitting this data via the 4010A'

# Among the codes, determine if they were present on admit.
mainDat[ , mrsa_sepsis_poa_rc := ifelse(mrsa_sepsis_poa == 'Y' | mrsa_sepsis_poa == 'E', 1, ifelse(mrsa_sepsis_poa == "N", 0, NA))]
mainDat[ , mrsa_pneum_poa_rc := ifelse(mrsa_pneum_poa == 'Y' | mrsa_sepsis_poa == 'E', 1, ifelse(mrsa_pneum_poa == "N", 0, NA))]
mainDat[ , mrsa_cod_poa_rc := ifelse(mrsa_cod_poa == 'Y' | mrsa_sepsis_poa == 'E', 1, ifelse(mrsa_cod_poa == "N", 0, NA))]
mainDat[ , mrsa_unspec_poa_rc := ifelse(mrsa_unspec_poa == 'Y' | mrsa_sepsis_poa == 'E', 1, ifelse(mrsa_unspec_poa == "N", 0, NA))]
mainDat[ , mrsa_carrier_poa_rc := ifelse(mrsa_carrier_poa == 'Y' | mrsa_sepsis_poa == 'E', 1, ifelse(mrsa_carrier_poa == "N", 0, NA))]

mainDat[mrsa_sepsis == 1, 
        round(100*mean(mrsa_sepsis_poa_rc == 1, na.rm = TRUE), 0)]
mainDat[mrsa_pneum == 1, 
        round(100*mean(mrsa_pneum_poa_rc == 1, na.rm = TRUE), 0)]
mainDat[mrsa_unspec == 1, 
        round(100*mean(mrsa_unspec_poa_rc == 1, na.rm = TRUE), 0)]
mainDat[mrsa_cause_other_disease == 1, 
        round(100*mean(mrsa_cod_poa_rc == 1, na.rm = TRUE), 0)]
mainDat[mrsa_carrier == 1, 
        round(100*mean(mrsa_carrier_poa_rc == 1, na.rm = TRUE), 0)]
}

" --, present_on_admit  --fact_diagnosis
  --, diagnosis_seq_nbr  --fact_diagnosis    
  --, diagnosis_entered_ts --fact_clinical_diagnosis
  --, clinical_dx   --fact_clinical_diagnosis
  --, condition_name   --fact_clinical_diagnosis
  --, vocabulary --fact_clinical_diagnosis
  --, dx_type  --fact_clinical_diagnosis
  --, clinical_service  --fact_clinical_diagnosis
  --, fact_clinical_diagnosis.contributor_system"
# RD's workflow:
# Begin with eta at 0.1-0.5 to get features over a few runs (e.g. 200 features  100  50  25 )
# Max depth 3-6 – 3 works just fine. More offers negligible improvement. 3 does fine.
# Subsample 0.3 – 0.7 depending on trainset size (e.g. 2M rows use 0.3, 10k-500k rows use 0.7)
# Once you have final feature list, turn eta down to  0.001- 0.1 to improve model auc. 0.1 works fine and is fast enough.

if(FALSE){


set.seed(1234)
xgb_cv_1 = xgb.cv(params = xgb_params_1,
                  data = sparse_matrix,
                  label = trainset$OUTCOME,
                  nrounds = 12000, 
                  nfold = 5,                                                   # number of folds in K-fold
                  prediction = TRUE,                                           # return the prediction using the final model 
                  showsd = TRUE,                                               # standard deviation of loss across folds
                  stratified = TRUE,                                           # sample is unbalanced; use stratified sampling
                  verbose = TRUE,
                  print.every.n = 1, 
                  early.stop.round = 20,
                  missing="NAN" 
)

max_auc = max(xgb_cv_1$dt$test.auc.mean)
max_auc_index = which.max(xgb_cv_1$dt$test.auc.mean==max_auc)

# fit the model with the arbitrary parameters specified above
xgb_1 = xgboost(data = sparse_matrix,
                label = trainset$OUTCOME,
                params = xgb_params_1,
                nrounds = max_auc_index,                                                 # max number of trees to build                verbose = TRUE,                                         
                print.every.n = 1,
                early.stop.round = 20                                # stop if no improvement within 20 trees
                
)

pred <- predict(xgb_1, sparse_matrix)
somers2(pred, trainset$OUTCOME)     #AUC = 0.986 with difftime and COMFORT_MEASURES
}

