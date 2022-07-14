# finish neuron

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
neuron_tibble <- read_csv(here("data/neuron_tibble.csv"))
neuron_trial_results <- read_rds(here("results/neuron_linear_results/neuron_trial_results.rds"))

# Split Data --------------------------------------------------------------
# set.seed(210)
# multi_net_initial_split <- initial_split(multi_net_kmeds_tibble, strata = class)
# multi_net_train <- training(multi_net_initial_split)
# multi_net_test <- testing(multi_net_initial_split)

set.seed(42)

neuron_cv_folds <- nested_cv(neuron_tibble,
                             outside = vfold_cv(v = 10, strata = class),
                             inside = bootstraps(times = 25, strata = class))


# Find Best Hyper-Parameters ----------------------------------------------
neuron_fold_res <- future_map(neuron_trial_results, format_folds, learner = "svm_linear")

neuron_score_res <- future_map(neuron_fold_res, median_wrapper)

best_hypers <- function(dat) dat[which.max(dat$auc),]

best_neuron_hyper_df <- map_dfr(neuron_score_res, best_hypers) %>%
  mutate(fold_id = neuron_cv_folds$id)


# Read in Best Models -----------------------------------------------------
param_list <- list(cost = best_neuron_hyper_df$cost,
                   fold_id = best_neuron_hyper_df$fold_id)

best_neuron_models <- best_neuron_hyper_df %>% 
  mutate(models =  future_pmap(param_list,
                               possibly(pick_model, NA_real_),
                               learner = "svm_linear",
                               results_dir = here("results/neuron_linear_results"))) %>%
  mutate(splits = neuron_cv_folds$splits)


# Set up Final Model Training ---------------------------------------------

dir.create(file.path(here("results"), paste("neuron", "full_linear", "final", "predictions", sep = "_")))

best_neuron_mod_list <- list(splits = best_neuron_models$splits,
                             models = best_neuron_models$models,
                             fold_id = best_neuron_models$fold_id)

final_res_dir <- here("results/neuron_full_linear_final_predictions")

final_neuron_auc <- best_neuron_models %>%
  mutate(assessment_auc = pmap_dbl(best_neuron_mod_list,
                                   possibly(bag_wrapper, NA_real_),
                                   res_dir = final_res_dir,
				   tissue = "neuron")) %>%
  select(cost, fold_id, auc, assessment_auc)

write_rds(final_neuron_auc, here("results/neuron_full_linear_final_predictions/final_neuron_auc.rds"))
