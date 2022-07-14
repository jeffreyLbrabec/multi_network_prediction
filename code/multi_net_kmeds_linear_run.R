#multi_net kmeds run
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
fun_dir <- list.files(here("kmeds_analysis/src"), full.names = TRUE)

for(i in 1:length(fun_dir)) {
  
  source(fun_dir[i])
  
}

#Get kmeds network tibble
multi_net_kmeds_tibble <- read_csv(here("data/multi_net_kmeds_tibble.csv"))

#Splitting Data
set.seed(210)
multi_net_initial_split <- initial_split(multi_net_kmeds_tibble, strata = class)
multi_net_train <- training(multi_net_initial_split)
multi_net_test <- testing(multi_net_initial_split)

write_rds(multi_net_test, here("data/multi_net_linear_test_data.rds"))

set.seed(42)

multi_net_cv_folds <- nested_cv(multi_net_train,
                                outside = vfold_cv(v = 2, strata = class),
                                inside = bootstraps(times = 3, strata = class))

#Prep Data for analysis
multi_net_recipe <- recipe(class ~ ., data = multi_net_train) %>% 
  update_role(genes, new_role = "gene_id") %>%
  step_center(all_predictors()) %>% 
  step_scale(all_predictors()) %>% 
  step_zv(all_predictors()) %>% 
  step_corr(all_predictors(), threshold = 0.4) %>% 
  themis::step_downsample(class) 

multi_train_juiced_obj <- prep(multi_net_recipe) %>% 
  juice() %>% 
  select(class, genes, everything())

multi_net_baked_full <- prep(multi_net_recipe, retain = TRUE) %>% 
  bake(new_data = multi_net_train)



rm(multi_net_kmeds_tibble)
rm(multi_net_recipe)
rm(multi_net_juiced_tib)

cost_vals = 10^seq(-5, 3, 1)

dir.create(file.path(here("results"), paste("multi_net", "linear", "results", sep = "_")))

res_dir <- here("results/multi_net_linear_results")

plan(multisession)
multi_net_trial_results <- future_map2(multi_net_cv_folds$inner_resamples,
                                       multi_net_cv_folds$id,
                                       summarizing_results,
                                       res_dir = res_dir,
                                       learner = "svm_linear",
                                       cost = cost_vals,
                                       .options = furrr_options(seed = 42)) 

write_rds(multi_net_trial_results, here("results/multi_net_linear_results/multi_net_linear_trial_results.rds"))
multi_net_trial_results <- read_rds(here("results/multi_net_linear_results/multi_net_linear_trial_results.rds"))

multi_net_fold_res <- future_map(multi_net_trial_results, format_folds)

multi_net_score_res <- future_map(multi_net_fold_res, median_wrapper)

best_hypers <- function(dat) dat[which.max(dat$auc),]

best_multi_net_hyper_df <- map_dfr(multi_net_score_res, best_hypers) %>% 
  mutate(fold_id = multi_net_cv_folds$id) 

param_list <- list(cost = best_multi_net_hyper_df$cost,
                   fold_id = best_multi_net_hyper_df$fold_id)

best_multi_net_models <- best_multi_net_hyper_df %>% 
  mutate(models =  future_pmap(param_list, 
                               pick_model,
                               learner = "svm_linear",
                               results_dir = here("results/multi_net_linear_results"))) %>% 
  mutate(splits = multi_net_cv_folds$splits)

dir.create(file.path(here("results"), paste("multi_net", "linear", "final", "predictions3", sep = "_")))

best_multi_net_mod_list <- list(splits = best_multi_net_models$splits,
                                models = best_multi_net_models$models,
                                fold_id = best_multi_net_models$fold_id)

final_res_dir <- here("results/multi_net_linear_final_predictions3")

final_multi_net_auc <- best_multi_net_models %>% 
  mutate(assessment_auc = pmap_dbl(best_multi_net_mod_list, 
                                   bag_wrapper,
                                   res_dir = final_res_dir)) %>% 
  select(cost, fold_id, auc, assessment_auc)

write_csv(final_multi_net_auc, 
          here("kmeds_analysis/results/multi_net_linear_results/final_multi_net_linear_aucs.txt"))

final_multi_net_auc <- read_csv(here("kmeds_analysis/results/multi_net_linear_results/final_multi_net_linear_aucs.txt"))

multi_net_test_data <- read_rds(here("kmeds_analysis/results/multi_net_linear_test_data.rds"))


finalModel <- ksvm(class ~ ., 
                   data = multi_train_juiced_obj, 
                   C = 0.001, 
                   scaled = TRUE, 
                   prob.model = TRUE)

multi_things <- multi_net_test_data %>% select(-class, -genes)

large_pred <- predict(finalModel, multi_net_test %>% select(-class), type = "prob") 
large_pred_tib <- as_tibble(large_pred) %>% 
  bind_cols(multi_net_test %>% dplyr::select(class)) %>% 
  mutate(class = as.factor(class))

roc_auc(large_pred_tib, truth = class, P)$.estimate

