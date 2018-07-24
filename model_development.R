

#######################################################################################
##
##
##  Train model to predict any systemic MRSA ICD10
##
##
#######################################################################################
setwd("E:/MRSA/model_development")
options(stringsAsFactors = FALSE)
library(Matrix)
library(RODBC)
library(Hmisc)
library(data.table)
library(recipes)
library(xgboost)
library(caret)
library(pROC)
library(ggplot2)
library(ROCR)
######################################################################################################
## Load data  
######################################################################################################
#fnamee = paste('../data/saved_R_workspaces_or_objects/key_datatables_end_of_dataManipulation.Rdata', sep = "")
#load(file = fnamee, verbose = TRUE)

## Load labs and vitals data
# Individual hour files
#for(i in seq_along(r_obj_names)){
#  load(r_obj_names[i], verbose = TRUE)}
# Combined file
#load(paste(fname_beginning, '_all_hrs_', r_obj, fname_end, '.Rdata', sep = ""), verbose = TRUE)

# Make sure there is only one age variable. There are age columns in both dataframe wide and mainDat.
## Features to use from main data
fewer_vars = Cs(
                    #  health_system_source_id
                    # facility_cd
                     ambulance_or_fire_not_walking
                    , emergency_room_flg    # maybe want to change this to something more timely later? It should reflect whether they entered through the ed. Would change it to something like encntr_class, but 
                       # encounter_cd_desc will not work for past data because after patients are discharged, the become ip-deceased or ip_discharged. Not sure how to handle this.
                    #, age
                    # los_so_far ## Need to think more about how to deal with this. 
                    # time of day of presentation?
                    , trauma_ward_by_3_hours
                    , ob_location_by_3_hours
                    , surg_location_by_3_hours
                    , rehab_by_3_hours
                    #, first_location
                    #, first_location_critical_care
                    #, location_at_3_hours
                    , in_critical_care_at_3_hours 
                    , regn_nm
                    , licensed_bed_count
                    , do_not_resuscitate
                    , race # interaction between race and ethnicity 
                    , sex_cd
                    , marital_status
                    , financial_class_cat
                    , ub_admit_source_desc_cat      # This is available early in encounter. We were originally using ub_admit_source_desc, without collapsing groups.
                    #snf,       # The factorized admit source will provide an identical variable.
                    , num_encounters_past_year
                    , last_ip_los
                    , time_since_last_ed
                    , time_since_last_ic
                    , time_since_last_ip
                    
                    , any_prev_hiv_dx
                    , any_prev_transplant
                    , any_prev_dep_on_ren_dialysis_dx
                    , urinary_cath
                    , central_line)
#fewer_vars_subset = main_w_prev_enc[ , ..fewer_vars]

ryans_vars = tolower(
   c('HEIGHT_M'
  ,'CELLULITIS'
  ,'LACTIC_ACID_CHANGE'
  ,'WEIGHT_KG'
  ,'ABSCESS'
  ,'LACTIC_ACID_MIN'
  ,'BMI'
  ,'RSLT_GLUCOSE_CHANGE'
  ,'ALBUMIN_LVL_MAX'
  ,'BUN_LAST'
  ,'TEMP_MAX'
  ,'PLATELET_COUNTC_FIRST'
  ,'AGE_AT_ADMIT'
  ,'WBC_FIRST'
  ,'RSLT_MAP_FIRST'
  ,'SIA_END'))            # shock index by age.
labs_vitals_all_hrs = labs_vitals_all_hrs[ , c('encounter_id', 'hour', ryans_vars), with = FALSE]
# Add outcome
labs_vitals_all_hrs[ , any_systemic_mrsa := main_w_prev_enc$any_systemic_mrsa[match(encounter_id, main_w_prev_enc$encounter_id)]]
labs_vitals_all_hrs = labs_vitals_all_hrs[!is.na(any_systemic_mrsa)]
labs_vitals_all_hrs[ , mrsa_fac := factor(any_systemic_mrsa, levels = c(1, 0), labels = c('yes', 'no'))]


ggplot(labs_vitals_all_hrs, aes(age_at_admit, any_systemic_mrsa)) +
geom_smooth() +
coord_cartesian(xlim = c(18, 100), ylim = c(0, 0.10))

ggplot(labs_vitals_all_hrs[weight_kg != -9999], aes(weight_kg*2.2046, any_systemic_mrsa)) +
geom_smooth()  +
coord_cartesian(xlim = c(75, 400), ylim = c(0, 0.10)) +
geom_line(y = 0.003)
# Why are their probs elevated? Because their weight was measured in the first place.


ggplot(labs_vitals_all_hrs[weight_kg != -9999], aes(weight_kg*2.2046)) +
geom_density()  



# less interesting
ggplot(labs_vitals_all_hrs[!is.na(bmi)], aes(bmi, any_systemic_mrsa)) +
geom_smooth() +
coord_cartesian(xlim = c(15, 40), ylim = c(0, 0.04)) +
geom_line(y = 0.003)

######################################################################################################
## Merge lab/vit and other IDM data.
######################################################################################################
all_all_hrs = merge(
   main_w_prev_enc[ , c('encounter_id', fewer_vars), with = FALSE],
   labs_vitals_all_hrs,
   by = 'encounter_id', 
   
   # Records must be in both files.
   all.x = FALSE,
   all.y = FALSE,
   
   # There should not be any of these.
   suffixes = c('', '.labs_vitals_dat'))
######################################################################################################
##  Separate data into training and test sets. I think we should use the same index for every different model we try.
######################################################################################################
# trainIndex = createDataPartition(main_w_prev_enc$any_systemic_mrsa,  # stratifies by mrsa
#                                  p = 0.70,
#                                  list = FALSE,
#                                  times = 1) 
# trainEncounters = main_w_prev_enc[trainIndex, encounter_id]
# testEncounters = main_w_prev_enc[-trainIndex, encounter_id]
set.seed(124800)
trainIndex = createDataPartition(all_all_hrs[hour == 1, any_systemic_mrsa],  # stratifies by mrsa
                                 p = 0.70,
                                 list = FALSE,
                                 times = 1) 
trainEncounters = all_all_hrs[hour == 1][trainIndex, encounter_id]
testEncounters = all_all_hrs[hour == 1][-trainIndex, encounter_id]

## For the **LAB DATA**, choose one hour (one row) for each encounter. Then filter to only train patients. This will be used to train on.
# (We don't need to choose an hour for each *test* patient, because we will probably assess performance separately on each hour.)
setkey(labs_vitals_all_hrs, encounter_id, hour)
setkey(all_all_hrs, encounter_id, hour)
# Randomly sort the hours:
hours = labs_vitals_all_hrs[ , sort(unique(hour))]    #   c(1, 3, 6, 12, 24)
hrz = sample(hours, length(unique(labs_vitals_all_hrs$encounter_id)), replace = TRUE)

labs_vitals_random_hour = labs_vitals_all_hrs[.(unique(encounter_id), hrz)]
all_random_hour = all_all_hrs[.(unique(encounter_id), hrz)]

labs_vitals_random_hour_train = labs_vitals_random_hour[trainEncounters]
all_random_hour_train                 = all_random_hour[trainEncounters]
# I think we don't need a test set from these data, because it only has one hour per encounter. 
# For our test data, I think we need the all hours data. 

###################################################################################################################
###################################################################################################################
## Prepare data that includes factors, etc.
all_vars_RHS = paste(c('hour', ryans_vars, fewer_vars), collapse = ' + ')
all_vars_recipe_obj = recipe(
  as.formula(paste("any_systemic_mrsa ~ ", all_vars_RHS)),
             data = all_random_hour_train) %>% 
  step_modeimpute(all_predictors(), -all_numeric()) %>%
 # step_bin2factor(all_outcomes()) %>%           # If want to use *classification* with caret, I think you need to use a factor for the outcome. 
  # I could exclude this step if I use the xgboost training function. ## What??
  step_novel(all_nominal(), -all_outcomes()) %>%   # If there are any new levels in the new data, this will make a new level for it.
  step_dummy(all_nominal(), -all_outcomes())
trained_recipe_all = prep(all_vars_recipe_obj, 
                    training = all_random_hour_train, 
                    stringsAsFactors = FALSE,      # You need this...
                    retain = TRUE)
# Apply the recipe to the original data
prepared_training_dat_all = data.table(juice(trained_recipe_all))
# Prepare test data
prepared_test_dat_all = data.table(bake(trained_recipe_all, newdata = all_all_hrs[testEncounters]))
setcolorder(prepared_test_dat_all, names(prepared_training_dat_all))
###################################################################################################################
###################################################################################################################
## 
##  Preliminary models with Ryan's pre-selected lab features.
##
###################################################################################################################
###################################################################################################################
## Make sparse matrices for xgboost input
# The -1 in the formula below means to not include an intercept.
sparse_mat_train = sparse.model.matrix(any_systemic_mrsa ~ . - 1, 
                                 data = labs_vitals_random_hour_train[ , c('any_systemic_mrsa', 'hour', ryans_vars), with = FALSE])
# This test set includes all hours (among test encounters)
sparse_mat_test  = sparse.model.matrix(any_systemic_mrsa ~ . - 1, 
                                 data = labs_vitals_all_hrs[testEncounters, c('any_systemic_mrsa', 'hour', ryans_vars), with = FALSE])  ## Cannot use -trainIndex because that is only applicable to data[hour = 1].

#sparse_mat_all_train = sparse.model.matrix(any_systemic_mrsa ~ . - 1, 
#                                 data = all_random_hour_train[ , c('any_systemic_mrsa', 'hour', ryans_vars, fewer_vars), with = #FALSE])
sparse_mat_all_train = sparse.model.matrix(any_systemic_mrsa ~ . - 1, 
                                 data = prepared_training_dat_all)
# This test set includes all hours (among test encounters)
sparse_mat_all_test  = sparse.model.matrix(any_systemic_mrsa ~ . - 1, 
                                 data = prepared_test_dat_all)
###################################################################################################################
## Model tuning/hyperparameters
ryan_orig_parameters = list(
   objective = "binary:logistic", 
   eval_metric = "auc", 
   eta = 0.1,                     # default = 0.3               
   max.depth = 4,                 # default = 6         
   subsample=0.5                  # default = 1
   # Defaults
   # min_child_weight  default = 1
   # colsample_bytree  default = 1
)
# This set of params gave these oos stats: 
# For both 1 hr and 6 hour together, for a > 0.01 cut off:
      ROC      Sens      Spec 
0.9112544 0.6512089 0.9344103 But very bad ppv.
a_la_verga_parameters = list(
   objective = "binary:logistic", 
   eval_metric = "auc", 
   
   eta = 0.03,                     # default = 0.3               
   max.depth = 5,                 # default = 6         
   subsample = 0.7,                  # default = 1
   
   # Defaults
   # min_child_weight  default = 1
   colsample_bytree = 0.8 # default = 1
)
# I wrote this one formatted for input in caret. 
ryan_tune_grid <- expand.grid(nrounds = c(100, 192, 300, 1000)
                         , eta = 0.1       
                         , max_depth = c(4)            
                         , subsample = c(0.5)                            

                         , gamma = 0                   # Default is 0.
                         , colsample_bytree = c(1)     # Default is 1.  
                         , min_child_weight = 1        # Default is 1.                 
                     )
# Tune more thoroughly!
tune_grid <- expand.grid(
    nrounds = c(1189)
  , eta = c(0.03, 0.06, 0.1)       
  , max_depth = c(5, 7)            
  , subsample = c(0.6, 0.7)                            
  , gamma = c(0, 10)                   # Default is 0. Try higher if big difference between train and test metrics. Try 10 or 20.
  , colsample_bytree = c(0.8, 0.9)     # Default is 1.  
  , min_child_weight = c(1, 10)        # Default is 1.                 
)
ya_tune_grid <- expand.grid(
    nrounds = c(400, 800, 1189)
  , eta = c(0.03, 0.06, 0.1)       
  , max_depth = c(5)            
  , subsample = c(0.7)                            
  , gamma = 0                   # Default is 0. Try higher if big difference between train and test metrics. Try 10 or 20.
  , colsample_bytree = c(0.7)     # Default is 1.  
  , min_child_weight = c(1)        # Default is 1.                 
)
# This is the list format. I was using it w/ xgb.cv, but that does not search parameters.
# tune_grid = list(
#   objective = "binary:logistic"        
#   , eval_metric = "auc"       
#   
#   , eta = c(1, 5, 10)*0.01        # default = 0.3
#   , max.depth = c(3, 6, 8)        # default = 6
#   , subsample = c(0.5, 0.6, 0.8)  # default = 1
# 
#   , gamma = c(0, 0.01, 0.02)      # default = 0  
#   , min_child_weight = c(1, 1.5)  # default = 1
#   , colsample_bytree = c(0.8, 1)  # default = 1
# )
###################################################################################################################
## Train using Ryan's original set up: original params and directly into xgboost functions.
ss = Sys.time()
set.seed(1234)
# This is cross validation and picks the best number of rounds, for a single set of model hyperparameters.
cv1 = xgb.cv(params = ryan_orig_parameters,
                  data = sparse_mat_train,
                  label = labs_vitals_random_hour_train[ , any_systemic_mrsa],
                  nrounds = 9000, 
                  nfold = 5,                                                   # number of folds in K-fold
                  prediction = TRUE,                                           # return the prediction using the final model 
                  showsd = TRUE,                                               # standard deviation of loss across folds
                  stratified = TRUE,                                           # sample is unbalanced; use stratified sampling
                  verbose = TRUE,
                  print_every_n = 1, 
                  early_stopping_rounds = 20,
                  missing="NAN")
# Fit the model
# Uses the number of rounds found in the cv step. (See the nrounds arg.)
mod1 = xgboost(data = sparse_mat_train,
                label = labs_vitals_random_hour_train[ , any_systemic_mrsa],
                params = ryan_orig_parameters,
                nrounds = cv1$best_iteration,  
                verbose = TRUE,                                         
                print_every_n = 10,
                early_stopping_rounds = 20                                # stop if no improvement within 20 trees
                )
ss - Sys.time()
save(mod1, cv1, file = 'other_model_objects/original_untuned_labs_only.Rdata')
## Train AUC
mod1$best_score
## Test AUC   
preds_mod1_on_test <- predict(mod1, sparse_mat_test)
test_preds_dt = labs_vitals_all_hrs[testEncounters, .(encounter_id, hour, pred = preds_mod1_on_test, any_systemic_mrsa)]
test_preds_dt[hour == 1, somers2(pred, any_systemic_mrsa)]
# 0.908
test_preds_dt[hour == 6, somers2(pred, any_systemic_mrsa)]
# 0.913

feat_names <- dimnames(sparse_mat_train)[[2]]
importance_matrix <- xgb.importance(feat_names, model = mod1)
xgb.plot.importance(importance_matrix[])
importance_matrix
###################################################################################################################
## Train using more variables
ss = Sys.time()
set.seed(56000)
# This is cross validation and picks the best number of rounds, for a single set of model hyperparameters.
cv_all = xgb.cv(params = a_la_verga_parameters,
                  data = sparse_mat_all_train,
                  label = all_random_hour_train[ , any_systemic_mrsa],
                  nrounds = 7000, 
                  nfold = 4,                                                   # number of folds in K-fold
                  prediction = TRUE,                                           # return the prediction using the final model 
                  showsd = TRUE,                                               # standard deviation of loss across folds
                  stratified = TRUE,                                           # sample is unbalanced; use stratified sampling
                  verbose = TRUE,
                  print_every_n = 10, 
                  early_stopping_rounds = 45,
                  missing="NAN")
save(cv_all, file = 'mod_all/cv_all.Rdata')
Sys.time() - ss
9]  train-auc:0.967640+0.000501     test-auc:0.939874+0.001808

Took 2 hours and stopped at 1189.
ss = Sys.time()
# Fit the model
# Uses the number of rounds found in the cv step. (See the nrounds arg.)
mod_all = xgboost(data = sparse_mat_all_train,
                label = all_random_hour_train[ , any_systemic_mrsa],
                params = a_la_verga_parameters,
                nrounds = 1189, #cv1$best_iteration,  
                verbose = TRUE,                                         
                print_every_n = 100,
                early_stopping_rounds = 100                                # stop if no improvement within 20 rounds
)
Sys.time() - ss 
save(mod_all, cv_all, file = 'mod_all/mod_all_cv_all.Rdata')
plot(cv_all$evaluation_log$train_auc_mean, ylim = c(0, 1))
points(cv_all$evaluation_log$test_auc_mean, col = 'red')
x11()
plot(cv_all$evaluation_log$train_auc_mean, ylim = c(0.75, 1))
points(cv_all$evaluation_log$test_auc_mean, col = 'red')
## Train AUC
mod_all$best_score

## Test AUC   
preds_mod_all_on_test <- predict(mod_all, sparse_mat_all_test)
test_preds_all_df = all_all_hrs[testEncounters, .(encounter_id, hour, pred = preds_mod_all_on_test, any_systemic_mrsa)]
test_preds_all_df[hour == 1, somers2(pred, any_systemic_mrsa)]
           C            Dxy              n        Missing 
     0.9234872      0.8469744 815399.0000000      0.0000000 


test_preds_all_df[hour == 6, somers2(pred, any_systemic_mrsa)]
          C            Dxy              n        Missing 
     0.9257115      0.8514230 815399.0000000      0.0000000 


feat_names <- dimnames(sparse_mat_all_test)[[2]]
importance_matrix_mod_all <- xgb.importance(feat_names, model = mod_all)
xgb.plot.importance(importance_matrix_mod_all[])
importance_matrix

twoClassSummary(data = test_preds_all_df[hour == 1, .(obs = factor(any_systemic_mrsa, levels = c(1, 0), labels = c('yes', 'no')), 
                 pred = factor(as.numeric(pred > 0.01), levels = c(1, 0), labels = c('yes', 'no')),
                 yes =                    pred)],
  lev = c('yes', 'no'))
      ROC      Sens      Spec 
0.9234872 0.6203537 0.9524395 

twoClassSummary(data = test_preds_all_df[hour == 6, .(obs = factor(any_systemic_mrsa, levels = c(1, 0), labels = c('yes', 'no')), 
                 pred = factor(as.numeric(pred > 0.01), levels = c(1, 0), labels = c('yes', 'no')),
                 yes =                    pred)],
  lev = c('yes', 'no'))
 ROC      Sens      Spec 
0.9257115 0.6171057 0.9538042 

twoClassSummary(data = test_preds_all_df[hour == 6, .(obs = factor(any_systemic_mrsa, levels = c(1, 0), labels = c('yes', 'no')), 
                 pred = factor(as.numeric(pred > 0.0491), levels = c(1, 0), labels = c('yes', 'no')),
                 yes =                    pred)],
  lev = c('yes', 'no'))
      ROC      Sens      Spec 
0.9257115 0.2623602 0.9909725 

# For alert rate of 1%, cutoff is 0.0491
test_preds_all_df[hour == 6, .(obs = factor(any_systemic_mrsa, levels = c(1, 0), labels = c('yes', 'no')), 
                 pred = factor(as.numeric(pred > 0.0491), levels = c(1, 0), labels = c('yes', 'no')),
                 yes =                    pred)][
                   , confusionMatrix(pred, obs)
                   ]


# For alert rate of 2%, cutoff is 0.026
test_preds_all_df[hour == 6, .(obs = factor(any_systemic_mrsa, levels = c(1, 0), labels = c('yes', 'no')), 
                 pred = factor(as.numeric(pred > 0.026), levels = c(1, 0), labels = c('yes', 'no')),
                 yes =                    pred)][
                   , confusionMatrix(pred, obs)
                   ]

# For alert rate of 3%, cutoff is 0.017
test_preds_all_df[hour == 6, .(obs = factor(any_systemic_mrsa, levels = c(1, 0), labels = c('yes', 'no')), 
                 pred = factor(as.numeric(pred > 0.017), levels = c(1, 0), labels = c('yes', 'no')),
                 yes =                    pred)][
                   , confusionMatrix(pred, obs)
                   ]

#################################################################################################################
## Train using Ryan's original hyperparams but using caret.
# Object is to confirm that results are consistent between the packages. If this works, it will be useful because of a built-in grid search framework.
ss = Sys.time()
set.seed(1234)
mod_obj_ryan_orig_caret = train(as.formula(paste("mrsa_fac ~ ", ryans_vars_hour_RHS)),
  data = labs_vitals_random_hour_train,
  method = 'xgbTree',
  trControl = trainControl(method = 'cv',
                           summaryFunction = twoClassSummary,
                           number = 5, 
                           savePredictions = 'final',                           
                           verboseIter = TRUE,
                           classProbs = TRUE),
  verbose = TRUE,
  tuneGrid = ryan_tune_grid)
# "tuning parameter grid should have columns nrounds, max_depth, eta, gamma, colsample_bytree, min_child_weight, subsample"
Sys.time() - ss  
# Took 47 mins
 
# Hmm. caret selected a different number of nrounds than xgb.cv found. The possible nrounds I provided were c(100, 192, 300, 1000). xgb.cv picked 192. However, all the nrounds i tested had similar results and small sd. I chalk it up to the randomness inherent in fitting.
Fitting nrounds = 300, max_depth = 4, eta = 0.1, gamma = 0, colsample_bytree = 1, min_child_weight = 1, subsample = 0.5 on full training set
getTrainPerf(mod_obj_ryan_orig_caret)
   TrainROC    TrainSens TrainSpec  method
1 0.9097072 0.0001564945 0.9999942 xgbTree
######################################################################################################
##  Variable importance
var_imp_ryan_orig <- varImp(mod_obj_ryan_orig_caret, scale = TRUE)
ggplot(var_imp_ryan_orig, top = 35)
######################################################################################################
##  Performance with different tuning parameters
mod_obj_ryan_orig_caret[['results']]
  eta max_depth gamma colsample_bytree min_child_weight subsample nrounds       ROC         Sens      Spec       ROCSD       SensSD       SpecSD
1 0.1         4     0                1                1       0.5     100 0.9073690 0.0000000000 0.9999984 0.001536306 0.0000000000 1.444259e-06
2 0.1         4     0                1                1       0.5     192 0.9096273 0.0000000000 0.9999963 0.001626138 0.0000000000 3.006458e-06
3 0.1         4     0                1                1       0.5     300 0.9097072 0.0001564945 0.9999942 0.001776502 0.0003499324 5.072059e-06
4 0.1         4     0                1                1       0.5    1000 0.9045041 0.0015647005 0.9999731 0.001879656 0.0007824728 6.825287e-06
###################################################################################################################
## Tune model with Ryan's variables
ss = Sys.time()
set.seed(1234)
mod_obj_ryans_vars_tuned = train(as.formula(paste("mrsa_fac ~ ", ryans_vars_hour_RHS)),
  data = labs_vitals_random_hour_train,
  method = 'xgbTree',
  trControl = trainControl(method = 'cv',
                           summaryFunction = twoClassSummary,
                           number = 2, 
                           savePredictions = 'final',                           
                           verboseIter = TRUE,
                           classProbs = TRUE),
  verbose = TRUE,
  tuneGrid = tune_grid)
Sys.time() - ss  
# Error in { : task 89 failed - "cannot allocate vector of size 24.2 Mb"
# In addition: There were 50 or more warnings (use warnings() to see the first 50)
# > Sys.time() - ss  
# Time difference of 1.047236 days
 
getTrainPerf(mod_obj_ryan_orig_caret)

? = data.frame(obs = factor(test_preds_dt$any_systemic_mrsa, levels = c(1, 0), labels = c('yes', 'no')), 
                 pred = factor(as.numeric(test_preds_dt$pred > 0.01), levels = c(1, 0), labels = c('yes', 'no')),
                 yes =                    test_preds_dt$pred)   

twoClassSummary(
  data = pred_dat_tuned,
  lev = c('yes', 'no'))


###################################################################################################################
## Tune model with more variables
ss = Sys.time()
set.seed(1234)
mod_obj_more_vars = train(as.formula(paste("mrsa_fac ~ ", ryans_vars_hour_RHS)),
  data = all_random_hour_train,
  method = 'xgbTree',
  trControl = trainControl(method = 'cv',
                           summaryFunction = twoClassSummary,
                           number = 2, 
                           savePredictions = 'final',                           
                           verboseIter = TRUE,
                           classProbs = TRUE),
  verbose = TRUE,
  tuneGrid = tune_grid)
Sys.time() - ss  
# Error in { : task 89 failed - "cannot allocate vector of size 24.2 Mb"
# In addition: There were 50 or more warnings (use warnings() to see the first 50)
# > Sys.time() - ss  
# Time difference of 1.047236 days
 
getTrainPerf(mod_obj_ryan_orig_caret)

? = data.frame(obs = factor(test_preds_dt$any_systemic_mrsa, levels = c(1, 0), labels = c('yes', 'no')), 
                 pred = factor(as.numeric(test_preds_dt$pred > 0.01), levels = c(1, 0), labels = c('yes', 'no')),
                 yes =                    test_preds_dt$pred)   

twoClassSummary(
  data = pred_dat_tuned,
  lev = c('yes', 'no'))
###################################################################################################################
# What did this go to?
# It is output from two class summary, I believe.
      ROC      Sens      Spec 
0.9112544 0.6512089 0.9344103 
   
# Next: try using a recipe object with xgboost function.   
xgb.plot.tree(feature_names = agaricus.train$data@Dimnames[[2]], model = mod1)  
xgb.plot.tree(feature_names = names(all_all_hrs), model = mod_all)
xgb.plot.multi.trees(model = mod1, feature_names = agaricus.train$data@Dimnames[[2]], features_keep = 3
###################################################################################################################
## ROC Curve
# proc package functions
# Make an ROC object
roc_obj_mod1 = roc(
    response = test_preds_dt$any_systemic_mrsa
  , predictor = test_preds_dt$pred
  , direction = "auto"
  , algorithm = 3)
plot(roc_obj_mod1)

plot(roc_obj_mod1
     , print.thres = c(0.001, 0.005, 0.01, 0.5)   # can give this a vector
     , type = "S", 
     print.thres.pattern = "%.3f (Spec = %.2f, Sens = %.2f)", 
     print.thres.cex = .8,
     legacy.axes = TRUE)
#ggroc












########################################################################################################
## Look at different cut off values
# Find the best threshold
# coords() extracts coordinates of the ROC curve
# In packages pROC.

# coords function looks really cool. 
# Ask for the "best" coordinates ant it will return the coordinates plus their stats
#   you can give it some weights so it will take into account the relative cost of fp and fn, as well as the prevalence.
# You can give it a threshold
# You can give it a sensitivity
# you can give it a specificity.

# Not sure if you can give it a vector of say, thresholds

coords(roc_obj_mod1
       , x = c('best')
       , ret = c("threshold", "sens", "ppv", 'spec'))
  threshold sensitivity         ppv specificity 
0.003835119 0.846264886 0.017070250 0.833837254 

#local maximas: upper angles of the roc curve. 
#coords(roc_obj_mod1
#       , x = c('local maximas')
#       , ret = c("threshold", "sens", "ppv", 'spec'))

# If x = 'best', 
# Weights can be supplied if false pos and false neg have different costs. 
# vector of length 2: first element is relative cost of false neg, second element is prevalence. 

# This gave a lot of output. I don't even know what it's doing. I thought it would give me one cut point.
#coords(roc_obj_mod1
#       , x = c('best')
#       , best.weights = c(0.003, 1)   # First element is prevalence. Second element is relative cost of a false negative compared #with false positive classification. 
#       , ret = c("threshold", "sens", "ppv", 'spec'))

# Can experiment with different values of the relative cost:  
prev = 0.003
cost = 1
# This is r from the documentation of coords().
r = function(p, c){(1-p)/(c*p)}
########################################################################################################
## Look at different cut off values, using predictions on the training data.
# preds_mod1_on_test <- predict(mod1, sparse_mat_test)
# test_preds_dt = labs_vitals_all_hrs[testEncounters, .(encounter_id, hour, pred = preds_mod1_on_test, any_systemic_mrsa)]


mod1_prediction_obj = test_preds_dt[ , prediction(pred, any_systemic_mrsa)]
#.... = prediction(pred, outcome)
#.... = performance(...., 'rpp')
mod1_alert_rates_cuttoff_perf_obj = performance(mod1_prediction_obj, 'rpp')

mod_all_prediction_obj = test_preds_all_df[ , prediction(pred, any_systemic_mrsa)]
#.... = prediction(pred, outcome)
#.... = performance(...., 'rpp')
mod_all_alert_rates_cuttoff_perf_obj = performance(mod_all_prediction_obj, 'rpp')
# Make the result into a data table.

# Make the result into a data table.


#pred <- prediction(subset(a, trainset==1)$mrsa_pred, subset(a, trainset==1)$OUTCOME)
#perf <- performance(pred,"rpp")
mod1_alert_rates_cuttoffs <- data.table(rpp = mod1_alert_rates_cuttoff_perf_obj@y.values[[1]]
                                        , cutoff = mod1_alert_rates_cuttoff_perf_obj@x.values[[1]])
mod_all_alert_rates_cuttoffs <- data.table(rpp = mod_all_alert_rates_cuttoff_perf_obj@y.values[[1]]
                                        , cutoff = mod_all_alert_rates_cuttoff_perf_obj@x.values[[1]])


cutoff_ar_01 = mod1_alert_rates_cuttoffs[rpp < 0.01, min(cutoff)]
cutoff_ar_01_mod_all = mod_all_alert_rates_cuttoffs[rpp < 0.01, min(cutoff)]
cutoff_ar_02_mod_all = mod_all_alert_rates_cuttoffs[rpp < 0.02, min(cutoff)]
cutoff_ar_03_mod_all = mod_all_alert_rates_cuttoffs[rpp < 0.03, min(cutoff)]


coords(roc_obj_mod1
       , x = cutoff_ar_01
       , ret = c("threshold", "sens", "ppv", 'spec'))



alert_rates = c(0.01, 0.02, 0.03)
for(i in seq_along(alert_rates)){
  print(
  coords(roc_obj_mod1
       , x = mod1_alert_rates_cuttoffs[rpp < alert_rates[i], min(cutoff)]
       , ret = c("threshold", "sens", "ppv", 'spec'))
)}
 threshold sensitivity         ppv specificity 
 0.05204263  0.25802959  0.08769240  0.99084637 
  threshold sensitivity         ppv specificity 
 0.02980429  0.39227716  0.06665645  0.98127003 
  threshold sensitivity         ppv specificity 
 0.02059692  0.47383616  0.05367619  0.97151403 



cutoffs_1<-subset(cutoffs, rpp>0.01)
c<-max(cutoffs_1$cutoff)
outcomes$c1<-ifelse(outcomes$mrsa_pred>=c, 1, 0)

# Returns a cutoff (scalar).
performance(mod1_prediction_obj, measure = "prbe")






a<-subset(outcomes, select=c('ENCOUNTERID','HOUR','OUTCOME','mrsa_pred','trainset'))


###################################################################################################################
## Thresholds




  
  
## This needs to be for a specified *alert rate*, say 0.01.  
  
coords(rocobj,
       
       )
  
  







###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################

##################
# I thought about filtering out those with missing hour, (because they were not in dataframe wide.) but I need to figure out why they werent in the dataframe wide data 
##$$labs_vitals_subset =  main_w_prev_enc[ , ..labs_vitals_vars]
##################
## Now for a model with both labs and vitals and the previous encounters data.
# Let ks indicate kitchen sink.
ks_vars = c(fewer_vars, setdiff(labs_vitals_vars, 'age_at_admit'))
ks_subset = main_w_prev_enc[ , ..ks_vars]
##################
## Now for a model with *selected* labs and vitals and the previous encounters data.
ks_rd_vars = c(fewer_vars, reduced_lbvs)
ks_rd_subset = main_w_prev_enc[ , ..ks_rd_vars]
##################
## Now for a model with both labs and vitals and the previous encounters data, but no facility_cd.
# Let ks indicate kitchen sink.
ks_vars_sin_fac = c(setdiff(fewer_vars, 'facility_cd'), setdiff(labs_vitals_vars, c('age_at_admit', 'any_systemic_mrsa')))
ks_sin_fac_subset = main_w_prev_enc[ , ..ks_vars_sin_fac]

######################################################################################################
##  Recipe 
######################################################################################################
# I think I may not need to use the subsets. I can just use main_w_prev_enc.
# I can use main_w_prev_enc as the data arg in the recipes
# I need separate recipes to indicate which features to use. (Via the formula arg to the recipefunction.)

# I *think* I can use main_w_prev_enc in the data arg (with the right training index) to caret's train function.
# I will need to supply the specific recipe to train(), though.
# I don't know if I can supply the same data or trained recipes to the prediction functions...definitely need to provide the specific model object.

fewer_vars_RHS = paste(fewer_vars, collapse = " + ")
labs_vitals_RHS = paste(labs_vitals_vars, collapse = " + ")
ks_rd_RHS = paste(ks_rd_vars, collapse = " + ")
ks_RHS = paste(fewer_vars_RHS, labs_vitals_RHS, sep = " + ")
ks_sin_fac_RHS = paste(paste(ks_vars_sin_fac, collapse = " + "), labs_vitals_RHS, sep = " + ")
ryans_vars_hour_RHS = paste(c(ryans_vars, 'hour'), collapse = " + ")
##!! Hope I can use one data set for all three recipes!!
## Then maybe I can get rid of the ks_subset, labs_vitals_subset, and fewer_vars_subset
## I think for these recipes and the training it is ok, but double check running on test data.
# Note that the data supplied to the data arguments are all the same for these three recipes.

# Looks like if there are no variables that need changing, it doesn't run but throws and error saying no vars or terms were selected.
ryans_vars_hour_recipe_obj = recipe(
  as.formula(paste("any_systemic_mrsa ~ ", ryans_vars_hour_RHS)),
             data = labs_vitals_random_hour_train) %>% 
  step_modeimpute(all_predictors(), -all_numeric()) %>%
 # step_bin2factor(all_outcomes()) %>%           # If want to use *classification* with caret, I think you need to use a factor for the outcome. 
  # I could exclude this step if I use the xgboost training function. ## What??
  step_novel(all_nominal(), -all_outcomes()) %>%   # If there are any new levels in the new data, this will make a new level for it.
  step_dummy(all_nominal(), -all_outcomes())

ks_recipe_obj = recipe(
  as.formula(paste("any_systemic_mrsa ~ ", ks_RHS)),
             data = main_w_prev_enc[trainIndex]) %>% 
  step_modeimpute(all_predictors(), -all_numeric()) %>%
  step_bin2factor(all_outcomes()) %>%           # If want to use *classification* with caret, I think you need to use a factor for the outcome. 
  # I could exclude this step if I use the xgboost training function. ## What??
  step_novel(all_nominal(), -all_outcomes()) %>%   # If there are any new levels in the new data, this will make a new level for it.
  step_dummy(all_nominal(), -all_outcomes())
  
ks_sin_fac_recipe_obj = recipe(
  as.formula(paste("any_systemic_mrsa ~ ", ks_sin_fac_RHS)),
             data = main_w_prev_enc[trainIndex]) %>% 
  step_modeimpute(all_predictors(), -all_numeric()) %>%
  step_bin2factor(all_outcomes()) %>%           # If want to use *classification* with caret, I think you need to use a factor for the outcome. 
  # I could exclude this step if I use the xgboost training function. ## What??
  step_novel(all_nominal(), -all_outcomes()) %>%   # If there are any new levels in the new data, this will make a new level for it.
  step_dummy(all_nominal(), -all_outcomes())  
  
labs_vitals_recipe_obj = recipe(
  as.formula(paste("any_systemic_mrsa ~ ", labs_vitals_RHS)),
             data = main_w_prev_enc[trainIndex]) %>% 
  step_modeimpute(all_predictors(), -all_numeric()) %>%
  step_bin2factor(all_outcomes()) %>%           # If want to use *classification* with caret, I think you need to use a factor for the outcome. 
  # I could exclude this step if I use the xgboost training function. ## What??
  step_novel(all_nominal(), -all_outcomes()) %>%   # If there are any new levels in the new data, this will make a new level for it.
  step_dummy(all_nominal(), -all_outcomes())  
  
fewer_vars_recipe_obj = recipe(
  as.formula(paste("any_systemic_mrsa ~ ", fewer_vars_RHS)),
             data = main_w_prev_enc[trainIndex]) %>% 
  step_modeimpute(all_predictors(), -all_numeric()) %>%
  step_bin2factor(all_outcomes()) %>%           # If want to use *classification* with caret, I think you need to use a factor for the outcome. 
  # I could exclude this step if I use the xgboost training function. ## What??
  step_novel(all_nominal(), -all_outcomes()) %>%   # If there are any new levels in the new data, this will make a new level for it.
  step_dummy(all_nominal(), -all_outcomes())   


ks_rd_recipe_obj = recipe(
  as.formula(paste("any_systemic_mrsa ~ ", ks_rd_RHS)),
             data = main_w_prev_enc[trainIndex]) %>% 
  step_modeimpute(all_predictors(), -all_numeric()) %>% # Should be imputing factors (categorical vars?) with the most common values
  step_bin2factor(all_outcomes()) %>%           # If want to use *classification* with caret, I think you need to use a factor for the outcome. 
  # I could exclude this step if I use the xgboost training function, because the label is not supposed to be a factor but rather a 0/1.
  step_novel(all_nominal(), -all_outcomes()) %>%   # If there are any new levels in the new data, this will make a new level for it.
  step_dummy(all_nominal(), -all_outcomes()) %>%    # All nominal vars will be converted from factors to dummies
  check_missing(all_predictors())  %>%              # This should cause prep to break if there are still any missing values.
  check_cols(all_predictors())              # This should cause bake to break if there are still any missing *columns*.

 

######################################################################################################
##  Select hyperparameters to try 
######################################################################################################
# Try setting max_delta_step. Try values of 0, 1
smaller_tune_grid <- expand.grid(nrounds = c(1000), #seq(from = 10, to = 1000, by = 100),  # boosting iterations
                         max_depth = c(3,9),            #  (controls model complexity)
                         eta = c(1, 10)*0.01,           # step size/learning rate (caret's documentation calls this shrinkage?)                         
                         gamma = 0, #c(1, 10)*0.01,     # min loss reduction (controls model complexity)????? Default is 0. “Lagrangian multiplier” 
                         colsample_bytree = c(0.7),     # subsample ratio of columns (adds randomness)    
                         min_child_weight = 1,                         
                         subsample = c(0.5)             # subsample percentage (adds randomness)
                     )          
# Let sg indicate small tuning grid.
tune_grid <- expand.grid(nrounds = c(700, 900)                #seq(from = 10, to = 1000, by = 100),  # boosting iterations
                      , max_depth = c(3, 5, 9)                # (controls model complexity)                  ## Reduced according to Hastie's guidelines
                      , eta = c(0.5, 1)*0.01                    # step size/learning rate (caret's documentation calls this shrinkage?)
                      , gamma = 0               # min loss reduction (controls model complexity)                      
                      , colsample_bytree = 0.7        # subsample ratio of columns (adds randomness)                      
                      , min_child_weight = c(0, 1)         # minimum sum of instance weight (controls model complexity)
                      , subsample = c(0.5, 0.8)               # subsample percentage (adds randomness) . Friedman[3] obtained that {\displaystyle 0.5\leq f\leq 0.8} {\displaystyle 0.5\leq f\leq 0.8} leads to good results for small and moderate sized training sets
) 
######################################################################################################
##  Train models 
######################################################################################################
############################################
## KS teeny grid
set.seed(34)
ss = Sys.time()
mod_obj_ks_sin_fac_teeny_g = train(
  ks_sin_fac_recipe_obj,
  data = main_w_prev_enc[trainIndex],
  method = 'xgbTree',
  trControl = trainControl(method = 'cv',
                           summaryFunction = twoClassSummary,
                           number = 5, 
                           savePredictions = 'final',                           
                           verboseIter = TRUE,
                           classProbs = TRUE),
  verbose = TRUE,
  tuneGrid = teeny_tune_grid)
# "tuning parameter grid should have columns nrounds, max_depth, eta, gamma, colsample_bytree, min_child_weight, subsample"
Sys.time() - ss

Fitting nrounds = 700, max_depth = 3, eta = 0.01, gamma = 0, colsample_bytree = 1, min_child_weight = 1, subsample = 0.7 on full training set
Warning message:
In train.recipe(ks_sin_fac_recipe_obj, data = main_w_prev_enc[trainIndex],  :
  The metric "Accuracy" was not in the result set. ROC will be used instead.
> # "tuning parameter grid should have columns nrounds, max_depth, eta, gamma, colsample_bytree, min_child_weight, subsample"
#Time difference of 35.20257 mins
############################################
## Fewer features
# When using the recipe method of train(), x should be an ***unprepared*** recipe object 
# I guess the train function "applies" the prep step? 
set.seed(34)
ss = Sys.time()
# I want to use a name needs to indicate that it's a model object and that what feaures were used. 
mod_obj_fewer_vars_teeny_w_fac = train(
  fewer_vars_recipe_obj,
  data = main_w_prev_enc[trainIndex],
  method = 'xgbTree',
  #early_stopping_rounds = 100,  # passed to xgb.train through the ... argument
  trControl = trainControl(method = 'cv',
                           summaryFunction = twoClassSummary,
                           number = 3, 
                           savePredictions = 'final',                           
                           # repeats? 
                           verboseIter = TRUE,
                           classProbs = TRUE), #method = 'none', 
  verbose = TRUE,
  #max_delta_step = 1,                # minimum sum of instance weight (controls model complexity) 
  tuneGrid = teeny_tune_grid #,    #tune_grid,           # step size/learning rate (caret's documentation calls this shrinkage?) 
  #, metric = "ROC"
  )
# "tuning parameter grid should have columns nrounds, max_depth, eta, gamma, colsample_bytree, min_child_weight, subsample"
Sys.time() - ss
save(mod_obj_fewer_vars_teeny_w_fac, file = '../model_development/mod_obj_fewer_vars_teeny_w_fac_2018-06-06_13_mo_new_feats.Rdata')
############################################
## Train KS Reduced Model 
# When using the recipe method of train(), x should be an ***unprepared*** recipe object 
# I guess the train function "applies" the prep step? 
set.seed(8000)
ss = Sys.time()
# I want to use a name needs to indicate that it's a model object and that what feaures were used. 
mod_obj_ks_rd_grid_24 = train(
  ks_rd_recipe_obj,
  data = main_w_prev_enc[trainIndex],
  method = 'xgbTree',
  trControl = trainControl(method = 'cv',
                           summaryFunction = twoClassSummary,
                           number = 4, 
                           savePredictions = 'final',                           
                           # repeats? 
                           verboseIter = TRUE,
                           classProbs = TRUE), #method = 'none', 
  verbose = TRUE,
  #max_delta_step = 1,                # minimum sum of instance weight (controls model complexity) 
  tuneGrid = tune_grid_24 #,    #tune_grid,           # step size/learning rate (caret's documentation calls this shrinkage?) 
  )
# "tuning parameter grid should have columns nrounds, max_depth, eta, gamma, colsample_bytree, min_child_weight, subsample"
Sys.time() - ss

save(mod_obj_ks_rd_grid_24, file = '../model_development/mod_obj_ks_rd_grid_24_2018-06-07.Rdata')
############################################
## Labs and vitals model
## with small tuning grid
set.seed(34)
ss = Sys.time()
mod_obj_labs_vitals = train(
  labs_vitals_recipe_obj,
  data = main_w_prev_enc[trainIndex],
  method = 'xgbTree',
  #early_stopping_rounds = 100,  # passed to xgb.train through the ... argument
  
  # Eventually want to try stratifying by outcome
  trControl = trainControl(method = 'cv',
                           summaryFunction = twoClassSummary,
                           number = 3, 
                           # repeats? 
                           savePredictions = 'final',                           
                           verboseIter = TRUE,
                           classProbs = TRUE), #method = 'none', 
  verbose = TRUE,
  tuneGrid = smaller_tune_grid #,    #tune_grid,           # step size/learning rate (caret's documentation calls this shrinkage?) 
  )
Sys.time() - ss
save(less_tuned_labs_vitals_fit_caret_xgb, file = 'less_tuned_labs_vitals_fit_2018-05-17.Rdata')
######################################################################################################
##  Variable importance
######################################################################################################
#caret_xgb_imp_unscaled <- varImp(fit_caret_xgb, scale = FALSE)
var_imp_ks_rd_grid_24 <- varImp(mod_obj_ks_rd_grid_24, scale = TRUE)
ggplot(var_imp_ks_rd_grid_24, top = 35)
var_imp_fewer_vars_teeny_w_fac[[1]]

less_tuned_labs_vitals_fit_caret_xgb_imp <- varImp(less_tuned_labs_vitals_fit_caret_xgb, scale = FALSE)

ggplot(less_tuned_labs_vitals_fit_caret_xgb_imp, top = 10)

######################################################################################################
##  Performance with different tuning parameters 
######################################################################################################

## First look at performance in training data, at different model structures (tuning params.)
results_dt_fewer_vars_teeny = data.table(mod_obj_fewer_vars_teeny[['results']])
  eta max_depth gamma colsample_bytree min_child_weight subsample nrounds       ROC         Sens      Spec       ROCSD      SensSD       SpecSD
1: 0.01         5     0              0.8                1       0.7    1100 0.7134482 0.0004704775 0.9999985 0.008916134 0.000814891 2.652837e-06
2: 0.10         5     0              0.8                1       0.7    1100 0.6458027 0.0004704775 0.9999694 0.027575271 0.000814891 2.769644e-05

getTrainPerf(mod_obj_fewer_vars_teeny)
   TrainROC    TrainSens TrainSpec  method
1 0.7134482 0.0004704775 0.9999985 xgbTree
 



## less_tuned_labs_vitals_fit_caret_xgb
Fitting nrounds = 1200, max_depth = 3, eta = 0.01, gamma = 0, colsample_bytree = 0.7, min_child_weight = 1, subsample = 0.5 on full training set
less_tuned_labs_vitals_fit_caret_xgb[["results"]]
   eta max_depth gamma colsample_bytree min_child_weight subsample nrounds       ROC        Sens      Spec       ROCSD      SensSD       SpecSD
1 0.01         3     0              0.7                1       0.5    1200 0.8953016 0.006472492 0.9999799 0.009565211 0.005605342 1.737803e-05
3 0.10         3     0              0.7                1       0.5    1200 0.7850465 0.000000000 0.9999197 0.006495222 0.000000000 4.597824e-05
2 0.01         9     0              0.7                1       0.5    1200 0.8767344 0.000000000 1.0000000 0.003026087 0.000000000 0.000000e+00
4 0.10         9     0              0.7                1       0.5    1200 0.7947065 0.000000000 0.9999900 0.002552620 0.000000000 1.737855e-05

###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################

######################################################################################################
##  Evaluate Model on ***out of sample*** data.
######################################################################################################

## Prepare the new/test data to be scored by the trained models 
labs_vitals_prod_data = data.table(bake(trained_labs_vitals_rec, newdata = labs_vitals_subset[-trainIndex]))
prod_data_fewer_vars = data.table(bake(fewer_vars_recipe_obj, newdata = main_w_prev_enc[-trainIndex]))
prod_data_ks_rd = data.table(bake(ks_rd_recipe_obj, newdata = main_w_prev_enc[-trainIndex]))

############################################################
## Get predictions from the models

# Find best cutoffs?
resample_stats <- thresholder(mod_obj_fewer_vars_g, 
                              threshold = seq(0, 1, by = 0.05), 
                              final = TRUE)

resample_stats <- thresholder(mod_obj_fewer_vars_g, 
                              threshold = seq(0, 0.05, by = 0.01), 
                              final = TRUE)


preds_ks_mod_g = predict(mod_obj_ks_g, newdata = main_w_prev_enc[-trainIndex], type = 'prob')


ggplot(resample_stats, aes(x = prob_threshold, y = J)) +
geom_point()
ggplot(resample_stats, aes(x = prob_threshold, y = Dist)) +
geom_point()
ggplot(resample_stats, aes(x = prob_threshold, y = Sensitivity)) +
154 train
geom_point() +
geom_point(aes(y = Specificity), col = "red")
confusionMatrix(
  factor(preds_ks_mod_g$yes > 0.00377, levels = c(TRUE, FALSE), labels = c('yes', 'no')), 
  factor(main_w_prev_enc[-trainIndex]$any_systemic_mrsa, levels = c(1, 0), labels = c('yes', 'no'))
)

posPredValue(factor(preds_ks_mod_g$yes > 0.000001, levels = c(TRUE, FALSE), labels = c('yes', 'no')), 
             factor(main_w_prev_enc[-trainIndex]$any_systemic_mrsa, levels = c(1, 0), labels = c('yes', 'no')))

postResample(pred = dumbdataframe$pred, obs = labs_vitals_prod_data$any_systemic_mrsa)




fit_caret_xgb_pred <- predict(fit_caret_xgb, newdata = main_sub[-trainIndex], type = 'prob')
less_tuned_labs_vitals_fit_caret_xgb_pred <- predict(less_tuned_labs_vitals_fit_caret_xgb, newdata = labs_vitals_subset[-trainIndex], type = 'prob')
str(fit_caret_xgb_pred)
#Can use twoClassSummary to look at model stats? This function is in caret.
#after using the predict funtion

twoClassSummary(
  data = data.frame(obs = prod_data$any_systemic_mrsa, 
                    pred = factor(as.numeric(fit_caret_xgb_pred$yes > 0.01), levels = c(1, 0), labels = c('yes', 'no')),
                    yes = fit_caret_xgb_pred$yes),
  
  lev = c('yes', 'no')
  #, class1 = the predicted probabilities? I am going to try and use the predicted probabilities directly in the pred variable.
  )
      ROC      Sens      Spec 
0.8281651 0.3454545 0.9390053 
# What's the distribution of the predictions (probs)?
describe(fit_caret_xgb_pred)
# What's the distribution of the predictions (yes/no)?
> swurf = factor(as.numeric(fit_caret_xgb_pred$yes > 0.01), levels = c(1, 0), labels = c('yes', 'no'))
describe(swurf)
swurf 
       n  missing distinct 
   33326        0        2 
                      
Value        yes    no
Frequency   2064 31262
Proportion 0.062 0.938

The alert rate of 0.06 gives these characteristics.
# What's the distribution of the observed outcomes?  


caret_xgb_imp_unscaled <- varImp(fit_caret_xgb, scale = FALSE)
caret_xgb_imp_scaled <- varImp(fit_caret_xgb, scale = TRUE)

ggplot(caret_xgb_imp, top = 10)
# any prev hiv was ranked number 33.


dumbdataframe = data.frame(obs = labs_vitals_prod_data$any_systemic_mrsa, 
                    pred = factor(as.numeric(less_tuned_labs_vitals_fit_caret_xgb_pred$yes > 0.01), levels = c(1, 0), labels = c('yes', 'no')),
                    yes = less_tuned_labs_vitals_fit_caret_xgb_pred$yes)
twoClassSummary(
  data = ,
  
  lev = c('yes', 'no')
  #, class1 = the predicted probabilities? I am going to try and use the predicted probabilities directly in the pred variable.
  )
      ROC      Sens      Spec 
0.9070445 0.5818182 0.9531551 
twoClassSummary(
  data = data.frame(obs = labs_vitals_prod_data$any_systemic_mrsa, 
                    pred = factor(as.numeric(less_tuned_labs_vitals_fit_caret_xgb_pred$yes > 0.006), levels = c(1, 0), labels = c('yes', 'no')),
                    yes = less_tuned_labs_vitals_fit_caret_xgb_pred$yes),
  
  lev = c('yes', 'no')
  #, class1 = the predicted probabilities? I am going to try and use the predicted probabilities directly in the pred variable.
  )
      ROC      Sens      Spec 
0.9070445 0.7545455 0.9161248 
#Also try defaultSummary. Says it can compute kappa?
# see function getTrainPerf
!! Yay!


thresholder()? Could be helpful!!
also see oneSE?
also see selectionFunction

#####################
## I think this is an approach using a different package. ROCR?
prediction_obj_mod1 = 
  prediction(predictions = test_preds_dt$pred, 
  labels = test_preds_dt$any_systemic_mrsa)
  
x.measure = 'cutoff' by default.


# Should give the rpp for various cutoffs
# Could plot
 <- performance(obj, measure = "rpp")
 
# Returns a cutoff (scalar).
# Could add this as a point or line on another plot.
 <- performance(obj, measure = "prbe")
 
# Plot the calibration error? 
 <- performance(obj, measure = "cal") 
 
# Others to look at: 
auc

rpp
cal: calipration error. 
prbe: precision recal break even point. Cutoff where precision and recall are equal. At this point, positive and negative preds are made at the same rate as their prevalence in the data.
sens
spec
ppv
  
# Get predictions (probabilities)  
preds_fewer_vars_teeny = predict(mod_obj_fewer_vars_teeny, newdata = main_w_prev_enc[-trainIndex], type = 'prob')
 
predRF <- prediction(rf_pred[,1], rm_test$PREGNANT)
perfRF <- performance(predRF, "tpr", "fpr")
'rpp' is the rate of positive predictions.
 
# Make a "prediction object"  
pred_obj_fewer_vars_teeny = 
  prediction(predictions = preds_fewer_vars_teeny$yes, 
  labels = prod_data_fewer_vars$any_systemic_mrsa)  

##!! rpp is rate of positive prediction! I think that's alert rate.
# Try tpr, fpr, cutoff, auc, rpp
x11()
perf_obj = performance(pred_obj_fewer_vars_teeny, measure = 'ppv', 'sens') 
#perf_obj = performance(pred_obj_fewer_vars_teeny, measure = 'prec', 'rec')  # This is the same thing
plot(perf_obj)

x11()
perf_obj = performance(pred_obj_fewer_vars_teeny, measure = 'rpp') 
plot(perf_obj)

x11()
perf_obj = performance(pred_obj_fewer_vars_teeny, measure = 'sens') 
plot(perf_obj)

x11()
perf_obj = performance(pred_obj_fewer_vars_teeny, measure = 'sens', 'spec') 
plot(perf_obj)

confusionMatrix(
  factor(preds_fewer_vars_teeny$yes > 0.00377, levels = c(TRUE, FALSE), labels = c('yes', 'no')), 
  factor(main_w_prev_enc[-trainIndex]$any_systemic_mrsa, levels = c(1, 0), labels = c('yes', 'no'))
)

# Make an ROC object
roc_obj_fewer_vars_teeny = roc(
  response = prod_data_fewer_vars$any_systemic_mrsa
  , predictor = preds_fewer_vars_teeny$yes
  , direction = ">")

# Find the best threshold
# coords() extracts coordinates of the ROC curve
coords(roc_obj_fewer_vars_teeny
       , x = c('best', 'local maximas')
       # , best.weights = c(0.003, 1)   # First element is prevalence. Second element is relative cost of a false negative compared with false positive classification. 
       , ret = c("threshold", "sens", "ppv", 'precision', 'recall'))  # precision and recall should be the same as sens and ppv.
  threshold sensitivity         ppv 
0.002547421 0.703937861 0.999358892   
add F1 score and maybe specificity. Possibly some other F scores, like F2?  
  
  
  
  
# I think I need a cut off for this.
confusionMatrix(
  preds_fewer_vars_teeny
  , factor(prod_data_fewer_vars$any_systemic_mrsa))
  
  
performance(prediction.obj = less_tuned_labs_vitals_predictionObj, measure = 'auc')     
metrixxxs = c('f', 'auc', 'ppv', 'npv', 'sens', 'spec', 'acc', 'err', 'rpp'  #'cutoff')

performanceExtractBinary = function(metric, truth, predictions) {
  #evalCat = classLabelsKeys[penaltyDx]
  predObj = prediction(predictions = predictions, labels = truth)            # as.numeric(binaryPreds > 0.5)
  perfObj = performance(prediction.obj = predObj, measure = metric)
  ctep = unlist(perfObj@y.values)
  names(ctep) = rep(perfObj@y.name, length(ctep))
  ctep = ifelse(length(ctep) == 3, ctep[2], ifelse(length(ctep) == 1, ctep, NA))
  
  return(ctep)
    }

predictionObjRoc = roc(prod_data$any_systemic_mrsa, fit_caret_xgb_pred$yes, direction="auto")
coords(predictionObjRoc, "best", ret=c("threshold", "sens", "ppv"))
  threshold sensitivity         ppv 
0.002547421 0.703937861 0.999358892 

less_tuned_labs_vitals_predictionObjRoc = roc(labs_vitals_prod_data$any_systemic_mrsa, less_tuned_labs_vitals_fit_caret_xgb_pred$yes, direction="auto")
coords(less_tuned_labs_vitals_predictionObjRoc, "best", ret=c("threshold", "sens", "ppv"))
  threshold sensitivity         ppv 
0.003774115 0.873103324 0.999379717 
confusionMatrix(factor(less_tuned_labs_vitals_fit_caret_xgb_pred$yes > 0.00377, levels = c(TRUE, FALSE), labels = c('yes', 'no')), labs_vitals_prod_data$any_systemic_mrsa)
Confusion Matrix and Statistics

          Reference
Prediction   yes    no
       yes    92  4218
       no     18 28998
                                          
               Accuracy : 0.8729          
                 95% CI : (0.8693, 0.8765)
    No Information Rate : 0.9967          
    P-Value [Acc > NIR] : 1               
                                          
                  Kappa : 0.0354          
 Mcnemar's Test P-Value : <2e-16          
                                          
            Sensitivity : 0.836364        
            Specificity : 0.873013        
         Pos Pred Value : 0.021346        
         Neg Pred Value : 0.999380        
             Prevalence : 0.003301        
         Detection Rate : 0.002761        
   Detection Prevalence : 0.129328        
      Balanced Accuracy : 0.854688        
                                          
       'Positive' Class : yes     




confusionMatrix(factor(fit_caret_xgb_pred$yes > 0.0025, levels = c(TRUE, FALSE), labels = c('yes', 'no')), prod_data$any_systemic_mrsa)
Confusion Matrix and Statistics

          Reference
Prediction   yes    no
       yes    95  9962
       no     15 23254
                                          
               Accuracy : 0.7006          
                 95% CI : (0.6957, 0.7055)
    No Information Rate : 0.9967          
    P-Value [Acc > NIR] : 1               
                                          
                  Kappa : 0.0122          
 Mcnemar's Test P-Value : <2e-16          
                                          
            Sensitivity : 0.863636        
            Specificity : 0.700084        
         Pos Pred Value : 0.009446        
         Neg Pred Value : 0.999355        
             Prevalence : 0.003301        
         Detection Rate : 0.002851        
   Detection Prevalence : 0.301776        
      Balanced Accuracy : 0.781860        
                                          
       'Positive' Class : yes             
                                         






#coords(xgboostObjRoc, "local maximas", ret=c("threshold", "sens", "spec", "ppv", "npv"))
    
sapply(metrixxxs, 
  FUN = performanceExtractBinary, 
  truth = testDat$cardiacEvent, 
  predictions = xgbPredProb)    

sapply(metrixxxs, 
  FUN = performanceExtractBinary, 
  truth = testDat$cardiacEvent, 
  predictions =  ifelse (xgbPredProb > 0.012,1,0))






###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
































######################
## This is COPIED from above for reference. This is a bad idea because it should be synced later.
#######################

# Try making a recipe obj without the outcome. Hope I get an error!
##! Yes! The input data must have the outcome. I don't get why it looks like my earlier code ran...
# ks_rd_recipe_obj_sin_y = recipe(
#   as.formula(paste("any_systemic_mrsa ~ ", ks_rd_RHS)),
#              data = main_w_prev_enc[trainIndex, ..ks_rd_vars]) %>% 
#   step_modeimpute(all_predictors(), -all_numeric()) %>%
#   step_bin2factor(all_outcomes()) %>%           # If want to use *classification* with caret, I think you need to use a factor for the outcome. 
#   # I could exclude this step if I use the xgboost training function. ## What??
#   step_novel(all_nominal(), -all_outcomes()) %>%   # If there are any new levels in the new data, this will make a new level for it.
#   step_dummy(all_nominal(), -all_outcomes())  

# Check whether each variable has any missing values.
prepared_training_dat_ks_rd[ , lapply(.SD, anyNA)]
# The only two that had any remaining missing values were previous dialysis and ambulance_or_fire_not_walking. 
# Have since fixed the dialysis one before it goes to recipe.

ks_rd_recipe_obj = recipe(
  as.formula(paste("any_systemic_mrsa ~ ", ks_rd_RHS)),
             data = main_w_prev_enc[trainIndex]) %>% 
  step_modeimpute(all_predictors(), -all_numeric()) %>% # Should be imputing factors (categorical vars?) with the most common values
  step_bin2factor(all_outcomes()) %>%           # If want to use *classification* with caret, I think you need to use a factor for the outcome. 
  # I could exclude this step if I use the xgboost training function. ## What??
  step_novel(all_nominal(), -all_outcomes()) %>%   # If there are any new levels in the new data, this will make a new level for it.
  step_dummy(all_nominal(), -all_outcomes()) %>%    # All nominal vars will be converted from factors to dummies
  check_missing(all_predictors())  %>%              # This should cause prep to break if there are still any missing values.
  check_cols(all_predictors())              # This should cause bake to break if there are still any missing *columns*.
trained_recipe_ks_rd = prep(ks_rd_recipe_obj, 
                    training = main_w_prev_enc[trainIndex], 
                    stringsAsFactors = FALSE,      # You need this...
                    retain = TRUE)
# Apply the recipe to the original data
# Had lots of trouble getting a data.table to work in the xgboost function (withouth using caret's train function)...

# Could also use bake() here, but this (juice) works if you want to apply the recipe to the same data you used to prep it. (Probably training data).
prepared_training_dat_ks_rd = data.table(juice(trained_recipe_ks_rd))
!!##!! Why are there still missing values in this?? There were missing values because I did not handle missing numeric variables (in the recipe function or beforehand). _Still working on figuring out how to do this._
# Must consider reason for missing. Could be that missing is informative, eg, a lab wasn't measured because the doc wasn't suspicious of it --> use extreme value. Could be that missing is an artifact of data collection independent of the true value --> mode impute.

## Prepare the new/test data to be scored by the trained models 
prepared_test_dat_ks_rd = data.table(bake(trained_recipe_ks_rd, newdata = main_w_prev_enc[-trainIndex]))
identical(names(prepared_training_dat_ks_rd), names(prod_data_ks_rd))
FALSE
identical(sort(names(prepared_training_dat_ks_rd)), sort(names(prod_data_ks_rd)))
TRUE
# Check this when testing on out of sample data. Right now the problem is with training performance.








tune_grid_24 = expand.grid(nrounds = c(100, 500, 1000),
                         max_depth = c(4, 7),
                         eta = c(0.01, 0.05),             
                         gamma = 0,              # default
                         colsample_bytree = 0.8,   #?? default   
                         min_child_weight = 1,   # default                         
                         subsample = c(0.7)) 
## Train KS Reduced Model 
# When using the recipe method of train(), x should be an ***unprepared*** recipe object 
# I guess the train function "applies" the prep step? 
set.seed(8000)
ss = Sys.time()

mod_obj_ks_rd_grid_24 = train(
  ks_rd_recipe_obj,
  data = main_w_prev_enc[trainIndex],
  method = 'xgbTree',
  trControl = trainControl(method = 'cv',
                           summaryFunction = twoClassSummary,
                           number = 2, 
                           savePredictions = 'final',                           
                           # repeats? 
                           verboseIter = TRUE,
                           classProbs = TRUE), #method = 'none', 
  verbose = TRUE,
  #max_delta_step = 1,                # minimum sum of instance weight (controls model complexity) 
  tuneGrid = tune_grid_24 #,    #tune_grid,           # step size/learning rate (caret's documentation calls this shrinkage?) 
  )
  
tune_grid_508 = expand.grid(nrounds = c(500, 750, 1000, 1500, 2000),
                         max_depth = c(4),
                         eta = c(0.01),             
                         gamma = 0,              # default
                         colsample_bytree = 0.8,   #?? default   
                         min_child_weight = 1,   # default                         
                         subsample = c(0.7))   
mod_obj_ks_rd_grid_508 = train(
  ks_rd_recipe_obj,
  data = main_w_prev_enc[trainIndex],
  method = 'xgbTree',
  trControl = trainControl(method = 'cv',
                           summaryFunction = twoClassSummary,
                           number = 5, 
                           savePredictions = 'final',                           
                           # repeats? 
                           verboseIter = TRUE,
                           classProbs = TRUE), #method = 'none', 
  verbose = TRUE,
  #max_delta_step = 1,                # minimum sum of instance weight (controls model complexity) 
  tuneGrid = tune_grid_508
  )  
  
Sys.time() - ss
save(mod_obj_ks_rd_grid_508, tune_grid_508, file = '../model_development/mod_obj_ks_rd_grid_508_2018-06-21.Rdata')
# Fitting nrounds = 500, max_depth = 4, eta = 0.01, gamma = 0, colsample_bytree = 0.8, min_child_weight = 1, subsample = 0.7 
# getTrainPerf(mod_obj_ks_rd_grid_24)
#   TrainROC    TrainSens TrainSpec  method
# 1 0.5872194 0.0002072109 0.9999806 xgbTree

############################################################
############################################################
############################################################
## Try same thing in xgboost function without caret.
############################################################
############################################################
############################################################

############################################################
############################################################
############################################################
## Ryan's code
trainset<-data[sample(nrow(data), nrow(data)*0.5), ]
testset <- data[ !(data$ENCOUNTERID %in% trainset$ENCOUNTERID), ]
trainset_ids<-subset(trainset, select=c('ENCOUNTERID'))
testset_ids<-subset(testset, select=c('ENCOUNTERID'))
trainset$ENCOUNTERID<-NULL
testset$ENCOUNTERID<-NULL

sparse_matrix2 = sparse.model.matrix(outcome~.-1, data = testset)
pred <- predict(xgb_1, sparse_matrix2)
somers2(pred, testset$outcome)    #F

hour_one_feat<-head(importance_matrix$Feature, 50)


############################################################
############################################################
############################################################
## The -1 in the formula below means to not include an intercept.
sparse_mat = sparse.model.matrix(any_systemic_mrsa ~ . - 1, data = main_w_prev_enc[trainIndex, c('any_systemic_mrsa', ryans_vars), with = FALSE])

set.seed(1234)
xgb_cv_1 = xgb.cv(params = xgb_params_1,
                  data = sparse_mat,
                  label = main_w_prev_enc[trainIndex, any_systemic_mrsa],
                  nrounds = 9000, 
                  nfold = 5,                                                   # number of folds in K-fold
                  prediction = TRUE,                                           # return the prediction using the final model 
                  showsd = TRUE,                                               # standard deviation of loss across folds
                  stratified = TRUE,                                           # sample is unbalanced; use stratified sampling
                  verbose = TRUE,
                  print.every.n = 1, 
                  early.stop.round = 20,
                  missing="NAN")

# max_auc = max(xgb_cv_1$dt$test.auc.mean)
max_auc_index_ja = xgb_cv_1$best_iteration

# fit the model with the arbitrary parameters specified above
xgb_1_ja = xgboost(data = sparse_mat,
                label = main_w_prev_enc[trainIndex, any_systemic_mrsa],
                params = xgb_params_1,
                nrounds = max_auc_index_ja,                                                 # max number of trees to build                verbose = TRUE,                                         
                print_every_n = 1,
                early_stopping_rounds = 20                                # stop if no improvement within 20 trees
                
)

pred <- predict(xgb_1_ja, sparse_mat)
somers2(pred, trainset$outcome)     #AUC = 0.986 with difftime and COMFORT_MEASURES














############################################################
## This is my older code (circa 2018-06-26)
############################################################
############################################################
# Now need to separate the outcome from the features. !Unless its in that xgb.Dmatrix??
# xgboost package recommends xgb.DMatrix: xgboost's own class
# Also need the label to be numeric.
#feats = setdiff(names(prepared_training_data), 'any_systemic_mrsa')
!!## TEST THIS
# xgb.DMatrix function doesn't accept a data frame. It accepts a dense or sparse matrix.   
prepared_training_dat_ks_rd_xgbdmat = 
  xgb.Dmatrix(
    data = matrix(prepared_training_dat_ks_rd[ , !'any_systemic_mrsa']), # must be a matrix-like object
    label = prepared_training_dat_ks_rd[ , as.numeric(any_systemic_mrsa == 'yes')]
    missing = NA)
!! I have been applying data.table to the baked/juiced data. Probably shouldn't do that for xgboost function, since it wants a matrix?
# what's the difference between objective function
xgb_cv_ks_rd_obj = 
xgb.cv(
  data = prepared_training_dat_ks_rd_xgbdmat
  params = smaller_tune_grid,  # list  
  nrounds = 12000
  early_stopping_rounds = 20,
  
  metrics = list('auc', 'aucpr'),  # list of evaluation metrics to be used in cross validation.
  showsd = TRUE,  # show sd of cross validation 
  
  nfold = 5,
  stratified = TRUE,
  
  verbose = TRUE,   # print statistics during the process.
  print_every_n = 100)
How are metrics = 'auc' and eval_metric = 'auc' different?

set.seed(34)
ss = Sys.time()
fit_xgboost = xgboost(
  data = as.matrix(prepared_training_data[trainIndex, ..feats, with = TRUE]), 
  label = prepared_training_data[trainIndex , as.numeric(any_systemic_mrsa == 'yes')],
  objective = 'binary:logistic',
  nrounds = 4,           # boosting iterations
  eval_metric = 'auc',
 
  # params argument must be a list
  params = list(

  max_depth = 5,         #  (controls model complexity)
  min_child_weight = 0.03,   # minimum sum of instance weight (controls model complexity)                        
  gamma = 0.2,               # min loss reduction (controls model complexity)

  colsample_bytree = 0.3,# subsample ratio of columns (adds randomness)
  subsample = 0.5,       # subsample percentage (adds randomness)  

  eta = 0.01),
  early_stopping_rounds = 20,
  print_every_n = 100
  )
Sys.time() - ss 
# Are these the same model? How do these perform?



############################################################
############################################################
## Predict from the model fit directly by xgboost 
fit_xgboost_pred <- predict(fit_xgboost, as.matrix(prepared_training_data[ , ..feats, with = TRUE]))
str(fit_xgboost_pred)


























predictionObj = 
  prediction(predictions = xgbPredProb, 
  labels = testDat$cardiacEvent)
performance(prediction.obj = predictionObj, measure = 'auc')     
metrixxxs = c('f', 'auc', 'ppv', 'npv', 'sens', 'spec', 'acc', 'err')

performanceExtractBinary = function(metric, truth, predictions) {
  #evalCat = classLabelsKeys[penaltyDx]
  predObj = prediction(predictions = predictions, labels = truth)            # as.numeric(binaryPreds > 0.5)
  perfObj = performance(prediction.obj = predObj, measure = metric)
  ctep = unlist(perfObj@y.values)
  names(ctep) = rep(perfObj@y.name, length(ctep))
  ctep = ifelse(length(ctep) == 3, ctep[2], ifelse(length(ctep) == 1, ctep, NA))
  
  return(ctep)
    }

coords(xgboostObjRoc, "best", ret=c("threshold", "sens", "ppv"))
#  threshold sensitivity         ppv 
# 0.01095007  0.84531340  0.08603535 
#coords(xgboostObjRoc, "local maximas", ret=c("threshold", "sens", "ppv", "spec", "npv"))
    
sapply(metrixxxs, 
  FUN = performanceExtractBinary, 
  truth = testDat$cardiacEvent, 
  predictions = xgbPredProb)    

sapply(metrixxxs, 
  FUN = performanceExtractBinary, 
  truth = testDat$cardiacEvent, 
  predictions =  ifelse (xgbPredProb > 0.012,1,0))
sapply(metrixxxs, 
  FUN = performanceExtractBinary, 
  truth = testDat$cardiacEvent, 
  predictions =  ifelse (xgbPredProb > 0.1,1,0))

sapply(metrixxxs, 
  FUN = performanceExtractBinary, 
  truth = testDat$cardiacEvent, 
  predictions =  ifelse (xgbPredProb > 0.2,1,0))
  sapply(metrixxxs, 
  FUN = performanceExtractBinary, 
  truth = testDat$cardiacEvent, 
  predictions =  ifelse (xgbPredProb > 0.3,1,0))
sapply(metrixxxs, 
  FUN = performanceExtractBinary, 
  truth = testDat$cardiacEvent, 
  predictions =  ifelse (xgbPredProb > 0.4,1,0))
sapply(metrixxxs, 
  FUN = performanceExtractBinary, 
  truth = testDat$cardiacEvent, 
  predictions =  ifelse (xgbPredProb > 0.5,1,0))
sapply(metrixxxs, 
  FUN = performanceExtractBinary, 
  truth = testDat$cardiacEvent, 
  predictions =  ifelse (xgbPredProb > 0.6,1,0))
sapply(metrixxxs, 
  FUN = performanceExtractBinary, 
  truth = testDat$cardiacEvent, 
  predictions =  ifelse (xgbPredProb > 0.978,1,0))  
#         f        auc        ppv        npv       sens       spec        acc        err 
#0.04386952 0.51121334 1.00000000 0.98380582 0.02242668 1.00000000 0.98381184 0.01618816 




































































































  
# recipe_obj = recipe(
#                     any_systemic_mrsa ~ 
#                     facility_cd +
#                     emergency_room_flg +
#                     age +
#                     # los so far. Don't have yet
#                     # time of day of presentation?
#                     do_not_resuscitate_ind +
#                     sex_cd +
#                     race +
#                     marital_status +
#                     financial_class +
#                     ub_admit_source_desc_cat +
#                     #snf +
#                     num_encounters_past_year +
#                     last_ip_los +
#                     time_since_last_ed +
#                     time_since_last_ic +
#                     time_since_last_ip + 
#                     any_prev_hiv_dx +
#                     any_prev_transplant
#                     , 
#                     data = main_sub[trainIndex]) %>% 
#   step_modeimpute(all_predictors(), -all_numeric()) %>%
#   step_bin2factor(all_outcomes()) %>%           # If want to use *classification* with caret, I think you need to use a factor for the outcome. 
#   
#   # I could easily exclude this step if I use the xgboost training function. ## What??
#   step_novel(all_nominal(), -all_outcomes()) %>%   # If there are any new levels in the new data, this will make a new level for it.
#   step_dummy(all_nominal(), -all_outcomes())

  #step_other(all_nominal(), threshold = 0.03)      # I think this makes an "other" category for categories that had small frequencies.
  #step_string2factor(all_nominal()) %>%   # If you do the step_other step, it already changes the strings to factors, so you can't add this.
  # Don't know if i need this...    step_string2factor(all_nominal(), -all_outcomes()) %>%
# ## Make recipe objects for the labs and vitals only model.
# labs_vitals_recipe_obj = recipe(
#                     any_systemic_mrsa ~ .,
#                     data = labs_vitals_subset) %>% 
#   step_modeimpute(all_predictors(), -all_numeric()) %>%
#   step_bin2factor(all_outcomes()) %>%           # If want to use *classification* with caret, I think you need to use a factor for the outcome. 
#   
#   # I could easily exclude this step if I use the xgboost training function. ## What??
#   step_novel(all_nominal(), -all_outcomes()) %>%   # If there are any new levels in the new data, this will make a new level for it.
#   step_dummy(all_nominal(), -all_outcomes())

## Might be able to cut this, but may need it.
# I think you need this for getting the predictions??
# trained_rec = prep(recipe_obj, 
#                    training = main_sub[trainIndex], 
#                    stringsAsFactors = FALSE,      # You need this...
#                    retain = TRUE)
# # Apply the recipe to the original data
# # Had lots of trouble getting a data.table to work in the xgboost function (withouth using caret's train function)...
# prepared_training_data = data.table(juice(trained_rec))
# 
# trained_labs_vitals_rec = prep(labs_vitals_recipe_obj, 
#                    training = labs_vitals_subset[trainIndex], 
#                    stringsAsFactors = FALSE,      # You need this...
#                    retain = TRUE)
# 
# prepared_labs_vitals_training_data = data.table(juice(trained_labs_vitals_rec))
# # Apply the recipe to out of sample or production data.



# The result is still a data.table, and it still has the same keys. The encounter ids are distinct.
# n = 210
# ss = Sys.time()
# ccc = data.table(id = as.character(rep(1:n, each = 3)), hour = rep(c(1, 3, 10), n))
# setkey(ccc, id, hour)
# # Randomly sorts the hours:
# hrs = sample(rep(c(1, 3, 10), n/3))
# eee = ccc[.(unique(id), hrs)]
# Sys.time() - ss

# This ended up taking much longer:
# Use the dplyr package
# This should get one row per encounter
#trainDat = main_w_prev_enc[trainIndex] %>% group_by(encounter_id) %>% sample_n(1)
#lb_train_dat = labs_vitals_dat[trainEncounters] %>% group_by(encounter_id) %>% sample_n(1)


## Gather features to use for test model with only labs and vitals
reduced_lbvs = Cs(
    height_m
  , weight_kg
  , bmi
  , cellulitis
  , abscess
  , sepsis
  , lactic_acid_first
  , lactic_acid_last
  , lactic_acid_change
  , agap_first
  , albumin_lvl_first
  , ast_first  
  , bili_total_first
  , bun_first
  , calcium_lvl_first
  , creatinine_lvl_first
  , creatinine_lvl_last
  , gfr_african_amc_first
  , gfr_non_african_am_first
  , hr_first
  , hr_min 
  , hr_max
  , hr_last
  , monocyte_abs_first
  , neutrophil_rel_lastc_first
  , platelet_countc_first
  , platelet_countc_last
  , rr_first
  , rr_mean
  , rr_max
  , rslt_band_first
  , rslt_glucose_first
  , rslt_glucose_last
  , rslt_map_first
  , spo2_first
  , spo2_mean
  , spo2_change
  , sbp_first
  , sbp_min
  , sbp_mean
  , sbp_max
  , sbp_last
  , sbp_age_start
  , sbp_age_end
  , si_start
  , si_end 
  , sia_start
  , sia_end
  , temp_first
  , temp_mean
  , temp_min
  , temp_max
  , temp_last
  , temp_change
  , wbc_first
  )


