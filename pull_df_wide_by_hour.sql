


select 
  distinct
  *
from 
  clinical_ops..dataframe_wide
where 
  hour = 1
  and age_at_admit >= 18
  and discharge_ts >= duration_subtract('2016-01-01 00:00:00', 00000002::numeric(8,0))
  and discharge_ts <= duration_add('2017-11-01 00:00:00', 00000002::numeric(8,0))
 ;
  
  


select 
  distinct
  *
from 
  clinical_ops..dataframe_wide
where 
  hour = 3
  and age_at_admit >= 18
  and discharge_ts >= duration_subtract('2016-01-01 00:00:00', 00000002::numeric(8,0))
  and discharge_ts <= duration_add('2017-11-01 00:00:00', 00000002::numeric(8,0))
 ;
  
 
select 
  distinct
  *
from 
  clinical_ops..dataframe_wide
where 
  hour = 6
  and age_at_admit >= 18
  and discharge_ts >= duration_subtract('2016-01-01 00:00:00', 00000002::numeric(8,0))
  and discharge_ts <= duration_add('2017-11-01 00:00:00', 00000002::numeric(8,0))
 ;
   

select 
  distinct
  *
from 
  clinical_ops..dataframe_wide
where 
  hour = 12
  and age_at_admit >= 18
  and discharge_ts >= duration_subtract('2016-01-01 00:00:00', 00000002::numeric(8,0))
  and discharge_ts <= duration_add('2017-11-01 00:00:00', 00000002::numeric(8,0))
 ;
 
 
select 
  distinct
  *
from 
  clinical_ops..dataframe_wide
where 
  hour = 24
  and age_at_admit >= 18
  and discharge_ts >= duration_subtract('2016-01-01 00:00:00', 00000002::numeric(8,0))
  and discharge_ts <= duration_add('2017-11-01 00:00:00', 00000002::numeric(8,0))
 ;
  
  