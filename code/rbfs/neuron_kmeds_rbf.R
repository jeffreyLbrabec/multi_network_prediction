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

#Source in Functions
fun_dir <- list.files(here("src"), full.names = TRUE)

for(i in 1:length(fun_dir)) {
  
  source(fun_dir[i])
  
}


# Read in Data ------------------------------------------------------------

neuron_kmeds_tibble <- read_csv(here("data/neuron_kmeds_tibble.csv"))


# Get Folds ---------------------------------------------------------------
set.seed(46)
neuron_cv_folds <- nested_cv(neuron_kmeds_tibble,
                             outside = vfold_cv(v = 10, strata = class),
                             inside = bootstraps(times = 25, strata = class))


# Set Hypers and Run Models -----------------------------------------------
cost_vals = 10^seq(-5, 3, 1)

dir.create(file.path(here("results"), paste("neuron", "rbf", "results", sep = "_")))

res_dir <- here("results/neuron_rbf_results")

plan(multisession)
neuron_trial_results <- future_map2(neuron_cv_folds$inner_resamples,
                                    neuron_cv_folds$id,
                                    summarizing_results,
                                    res_dir = res_dir,
                                    learner = "svm_rbf",
                                    cost = cost_vals,
                                    .options = furrr_options(seed = 42)) 

write_rds(neuron_trial_results, here("results/neuron_rbf_results/neuron_trial_rbf_results.rds"))

# Aggregate Model Performance Across Hyper-Parameters ---------------------

neuron_fold_res <- map(neuron_trial_results, format_folds, learner = "svm_rbf")

neuron_score_res <- map(neuron_fold_res, median_wrapper)

best_hypers <- function(dat) dat[which.max(dat$auc),]

best_neuron_hyper_df <- map_dfr(neuron_score_res, best_hypers) %>% 
  mutate(fold_id = neuron_cv_folds$id) 


# Read in Best Models -----------------------------------------------------

param_list <- list(cost = best_neuron_hyper_df$cost,
                   fold_id = best_neuron_hyper_df$fold_id)

best_neuron_models <- best_neuron_hyper_df %>% 
  mutate(models =  future_pmap(param_list, 
                               pick_model,
                               learner = "svm_rbf",
                               results_dir = res_dir)) %>% 
  mutate(splits = neuron_cv_folds$splits)


# Calculate final AUCs ----------------------------------------------------

dir.create(file.path(here("results"), paste("neuron", "rbf", "final", "predictions", sep = "_")))

best_neuron_mod_list <- list(splits = best_neuron_models$splits,
                             models = best_neuron_models$models,
                             fold_id = best_neuron_models$fold_id)

final_res_dir <- here("results/neuron_rbf_final_predictions")

final_neuron_auc <- best_neuron_models %>% 
  mutate(assessment_auc = pmap_dbl(best_neuron_mod_list, 
                                   bag_wrapper,
                                   res_dir = final_res_dir,
                                   learner = "neuron")) %>% 
  select(cost, fold_id, auc, assessment_auc)

write_rds(final_neuron_auc, here("results/neuron_rbf_final_predictions/final_neuron_rbf_auc.rds"))
