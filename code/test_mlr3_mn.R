## mlr3 Benchmark Pipeline

library(mlr3verse)
library(tidyverse)
library(here)


# Set seed and logger messages --------------------------------------------
set.seed(7832)
lgr::get_logger("mlr3")$set_threshold("warn")
lgr::get_logger("bbotk")$set_threshold("warn")


# Read in data and convert to tasks ---------------------------------------

neuron_net <- read_csv(here("data/neuron_kmeds_tibble.csv")) %>% 
  mutate(class = as.factor(class))
glia_net <- read_csv(here("data/glia_kmeds_tibble.csv")) %>% 
  mutate(class = as.factor(class))
astro_net <- read_csv(here("data/astro_kmeds_tibble.csv")) %>% 
  mutate(class = as.factor(class))
multi_net <- read_csv(here("data/multi_net_kmeds_tibble.csv")) %>% 
  mutate(class = as.factor(class))

task_neuron <- as_task_classif(neuron_net, target = "class", id = "genes")
task_neuron$col_roles$stratum <- task_neuron$target_names
task_glia <- as_task_classif(glia_net, target = "class", id = "genes")
task_glia$col_roles$stratum <- task_glia$target_names
task_astro <- as_task_classif(astro_net, target = "class", id = "genes")
task_astro$col_roles$stratum <- task_astro$target_names
task_multi <- as_task_classif(multi_net, target = "class", id = "genes")
task_multi$col_roles$stratum <- task_multi$target_names


# Prep learners -----------------------------------------------------------
print("Prepping Learners")
## SVM
learner_svm <- lrn("classif.svm", 
                   type = "C-classification", 
                   predict_type = "prob", 
                   predict_sets = c("train", "test"))

learner_svm$param_set$values$cost <- to_tune(p_dbl(1e-5, 1e5, logscale = TRUE))
learner_svm$param_set$values$gamma <- to_tune(p_dbl(1e-5, 1e5, logscale = TRUE))
learner_svm$param_set$values$kernel <- to_tune(p_fct(c("polynomial", "radial")))
learner_svm$param_set$values$degree <- to_tune(p_int(1, 4))

## Random Forest
learner_rf <- lrn("classif.ranger",
                  predict_type = "prob",
                  predict_sets = c("train", "test"))

learner_rf$param_set$values$mtry.ratio <- to_tune(p_dbl(0, 1))


# Define Preprocessing ----------------------------------------------------

print("Defining Preprocessors")

preproc <- po("removeconstants", ratio = 0.1) %>>%
  po("classbalancing",
     id = "undersample",
     adjust = "major",
     reference = "major", 
     shuffle = FALSE, 
     ratio = 723/25102)

learner_preproc_svm <- as_learner(preproc %>>% learner_svm)
learner_preproc_rf <- as_learner(preproc %>>% learner_rf)


# Setup model AutoTuners --------------------------------------------------

print("Setting up AutoTuners")

inner_resampling <-  rsmp("bootstrap", repeats = 25)

at_svm <- auto_tuner(
  method = "grid_search",
  learner = learner_preproc_svm,
  resampling = inner_resampling,
  measure = msr("classif.auc"),
  store_models = TRUE)

at_rf <- auto_tuner(
  method = "grid_search",
  learner = learner_preproc_rf,
  resampling = inner_resampling,
  measure = msr("classif.auc"),
  store_models = TRUE)

# Set up benchmark design and Run -----------------------------------------

print("Building design and running benchmark")

outer_resampling <- rsmp("cv", folds = 10)

design <-  benchmark_grid(tasks = list(task_neuron,
                                       task_astro,
                                       task_glia,
                                       task_multi), 
                          learners = list(at_svm, at_rf), 
                          resamplings = outer_resampling)

future::plan(list("multisession", "sequential"))
bmr <- benchmark(design, store_models = TRUE)

write_rds(bmr, here("results/multi_net_bmr.rds"))

print("Done!")
