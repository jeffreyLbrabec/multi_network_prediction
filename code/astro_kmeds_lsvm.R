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
astro_kmeds_tibble <- read_csv(here("data/astro_kmeds_tibble.csv"))


# Split Data --------------------------------------------------------------
# set.seed(210)
# multi_net_initial_split <- initial_split(multi_net_kmeds_tibble, strata = class)
# multi_net_train <- training(multi_net_initial_split)
# multi_net_test <- testing(multi_net_initial_split)

set.seed(42)

astro_cv_folds <- nested_cv(astro_kmeds_tibble,
                            outside = vfold_cv(v = 10, strata = class),
                            inside = bootstraps(times = 25, strata = class))


# Set Up Cost Grid and Run Model ------------------------------------------
cost_vals = 10^seq(-5, 3, 1)

dir.create(file.path(here("results"), paste("astro", "linear", "results", sep = "_")))

res_dir <- here("results/astro_linear_results")

plan(multisession)
astro_trial_results <- future_map2(astro_cv_folds$inner_resamples,
                                   astro_cv_folds$id,
                                   summarizing_results,
                                   res_dir = res_dir,
                                   learner = "svm_linear",
                                   cost = cost_vals,
                                   .options = furrr_options(seed = 42)) 
write_rds(astro_trial_results, here("results/astro_linear_results/astro_trial_results.rds"))


# Find Best Hyper-Parameters ----------------------------------------------
astro_fold_res <- future_map(astro_trial_results, format_folds, learner = "svm_linear")

astro_score_res <- future_map(astro_fold_res, median_wrapper)

best_hypers <- function(dat) dat[which.max(dat$auc),]

best_astro_hyper_df <- map_dfr(astro_score_res, best_hypers) %>% 
  mutate(fold_id = astro_cv_folds$id) 

write_rds(best_astro_hyper_df, here("results/astro_linear_results/best_astro_hyper_df.rds"))

# Read in Best Models -----------------------------------------------------
param_list <- list(cost = best_astro_hyper_df$cost,
                   fold_id = best_astro_hyper_df$fold_id)

best_astro_models <- best_astro_hyper_df %>% 
  mutate(models =  future_pmap(param_list, 
                               pick_model,
                               learner = "svm_linear",
                               results_dir = here("results/astro_linear_results"))) %>% 
  mutate(splits = astro_cv_folds$splits)


# Set up Final Model Training ---------------------------------------------

dir.create(file.path(here("results"), paste("astro", "linear", "final", "predictions", sep = "_")))

best_astro_mod_list <- list(splits = best_astro_models$splits,
                            models = best_astro_models$models,
                            fold_id = best_astro_models$fold_id)

final_res_dir <- here("results/astro_linear_final_predictions")

final_astro_auc <- best_astro_models %>% 
  mutate(assessment_auc = pmap_dbl(best_astro_mod_list, 
                                   bag_wrapper,
                                   res_dir = final_res_dir)) %>% 
  select(cost, fold_id, auc, assessment_auc)

write_rds(final_astro_auc, here("results/astro_linear_final_predictions/final_astro_auc.rds"))