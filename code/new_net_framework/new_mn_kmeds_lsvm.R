# Updated multi net kmeds run

library(tidyverse)
library(tidymodels)
library(here)
library(janitor)
library(gprofiler2)
library(magrittr)
library(igraph)
library(ranger)
library(kernlab)
library(baguette)
library(furrr)
library(R.utils)
library(parallel)
library(doParallel)
library(safejoin)
library(cluster)
library(fpc)
library(themis)

# Source in Functions -----------------------------------------------------

fun_dir <- list.files(here("src"), full.names = TRUE)

for(i in 1:length(fun_dir)) {
  
  source(fun_dir[i])
  
}


# Read in Data ------------------------------------------------------------
new_mn_tib <- read_csv(here("data/new_multi_net_med.csv"))


# Split Data --------------------------------------------------------------
#set.seed(210)
#multi_net_initial_split <- initial_split(multi_net_kmeds_tibble, strata = class)
#multi_net_train <- training(multi_net_initial_split)
#multi_net_test <- testing(multi_net_initial_split)

set.seed(42)

multi_net_cv_folds <- nested_cv(new_mn_tib,
                                outside = vfold_cv(v = 10, strata = class),
                                inside = bootstraps(times = 25, strata = class))


# Set Up Cost Grid and Run Model ------------------------------------------
cost_vals = 10^seq(-5, 3, 1)

dir.create(file.path(here("results"), paste("new_multi_net", "linear", "results", sep = "_")))

res_dir <- here("results/new_multi_net_linear_results")

plan(multisession)
multi_net_trial_results <- future_map2(multi_net_cv_folds$inner_resamples,
                                       multi_net_cv_folds$id,
                                       summarizing_results,
                                       res_dir = res_dir,
                                       learner = "svm_linear",
                                       cost = cost_vals,
                                       .options = furrr_options(seed = 42)) 
write_rds(multi_net_trial_results, here("results/new_multi_net_linear_results/new_multi_net_trial_results.rds"))


# Find Best Hyper-Parameters ----------------------------------------------
multi_net_fold_res <- future_map(multi_net_trial_results, format_folds, learner = "svm_linear")

multi_net_score_res <- future_map(multi_net_fold_res, median_wrapper)

best_hypers <- function(dat) dat[which.max(dat$auc),]

best_multi_net_hyper_df <- map_dfr(multi_net_score_res, best_hypers) %>% 
  mutate(fold_id = multi_net_cv_folds$id) 


# Read in Best Models -----------------------------------------------------
param_list <- list(cost = best_multi_net_hyper_df$cost,
                   fold_id = best_multi_net_hyper_df$fold_id)

best_multi_net_models <- best_multi_net_hyper_df %>% 
  mutate(models =  future_pmap(param_list, 
                               pick_model,
                               learner = "svm_linear",
                               results_dir = here("results/new_multi_net_linear_results"))) %>% 
  mutate(splits = multi_net_cv_folds$splits)


# Set up Final Model Training ---------------------------------------------

dir.create(file.path(here("results"), paste("new_multi_net", "linear", "final", "predictions", sep = "_")))

best_multi_net_mod_list <- list(splits = best_multi_net_models$splits,
                                models = best_multi_net_models$models,
                                fold_id = best_multi_net_models$fold_id)

final_res_dir <- here("results/new_multi_net_linear_final_predictions")

final_multi_net_auc <- best_multi_net_models %>% 
  mutate(assessment_auc = pmap_dbl(best_multi_net_mod_list, 
                                   possibly(bag_wrapper, otherwise = NA_real_),
                                   res_dir = final_res_dir,
                                   tissue = "multi_net")) %>% 
  select(cost, fold_id, auc, assessment_auc)

write_rds(final_multi_net_auc, here("results/new_multi_net_linear_final_predictions/new_final_multi_net_auc.rds"))
