roc_wrapping <- function(cost = 1,
                         degree = 1,
                         object, 
                         learner = c("svm_poly", "svm_rbf", "svm_linear"),
                         boot_dir) {
  
  if(learner == "svm_poly") {
    modeling_roc(object, cost = cost, degree = degree, learner = learner, boot_dir = boot_dir)
  }else if(learner == "svm_rbf") {
    modeling_roc(object, cost = cost, learner = learner, boot_dir = boot_dir)
  }else if(learner == "svm_linear") {
    modeling_roc(object, cost = cost, learner = learner, boot_dir = boot_dir)
  }
}