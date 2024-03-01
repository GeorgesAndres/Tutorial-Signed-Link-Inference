# A custom function to evaluate the in-sample performance of the predictive model.
# It calculates probabilities, determines the optimal threshold for classification, and computes
# a confusion matrix, sensitivity, specificity, and balanced accuracy to assess model accuracy.
evaluate_insample <- function(datt, model, model_name, positive_level, negative_level) {
  prob <- predict(model, type = c("probs"), newdata = datt)
  p_friend <- prob[, c("1", "2")] %>% apply(1, sum)
  p_nonfriend <- prob[, c("-2", "0")] %>% apply(1, sum)
  
  g <- pROC::roc(friend ~ prob, data = datt %>% mutate(prob = p_friend))
  datt$prediction <- p_friend
  optimal <- pROC::coords(g, "best", ret = "threshold") %>% pull(threshold)
  datt %>% mutate(pred = as.factor(if_else(prediction <= optimal, negative_level, positive_level))) -> datt
  confMat <- caret::confusionMatrix(data = datt$pred, reference = datt$friend, positive = positive_level)
  sens <- caret::sensitivity(data = datt$pred, reference = datt$friend, positive = positive_level, negative = negative_level)
  spec <- caret::specificity(data = datt$pred, reference = datt$friend, positive = positive_level, negative = negative_level)
  balanced_accuracy <- (sens + spec) / 2
  
  list(confusionMatrix = confMat, res = tibble(sensitivity = sens, specificity = spec, balancedAccuracy = balanced_accuracy, model = model_name))
}


# A specialized function to perform Leave-One-Out Cross-Validation (LOOCV) for a three-class prediction task.
# This approach is designed to rigorously test the model's performance by training on all but one observation
# and testing on the left-out observation, iterated across the entire dataset.
make_LOOCV_threeClass <- function(dat_small, dat_big, formula, split, model_name, positive_level = "WITH", negative_level = "NOT", seed = 12) {
  set.seed(seed)
  index <- sample(1:nrow(dat_small), size = round(nrow(dat_small) * split), replace = FALSE)
  train_data <- dat_small[index, ]
  test_data <- dat_big %>% left_join(train_data %>% select(from, to) %>% mutate(xx = from, yy = to)) %>% filter(is.na(xx)) %>% select(-xx, -yy)
  
  ctrl <- trainControl(method = "LOOCV")
  mod <- MASS::polr(as.formula(formula), data = train_data)
  strt <- c(coef(mod), mod$zeta)
  
  model_under <- train(as.formula(formula), data = train_data, trControl = ctrl, method = "polr", metric = "Kappa", start = strt)
  prob <- predict(model_under, type = c("prob"), newdata = test_data)
  p_friend <- prob[, c("1", "2")] %>% apply(1, sum)
  
  g <- pROC::roc(friend ~ prob, data = test_data %>% mutate(prob = p_friend))
  test_data$prediction <- p_friend
  optimal <- pROC::coords(g, "best", ret = "threshold") %>% pull(threshold)
  test_data %>% mutate(pred = as.factor(if_else(prediction <= optimal, negative_level, positive_level))) -> test_data
  confMat <- caret::confusionMatrix(data = test_data$pred, reference = test_data$friend, positive = positive_level)
  sens <- caret::sensitivity(data = test_data$pred, reference = test_data$friend, positive = positive_level, negative = negative_level)
  spec <- caret::specificity(data = test_data$pred, reference = test_data$friend, positive = positive_level, negative = negative_level)
  balanced_accuracy <- (sens + spec) / 2
  
  list(confusionMatrix = confMat, res = tibble(sensitivity = sens, specificity = spec, balancedAccuracy = balanced_accuracy, model = model_name))
}
