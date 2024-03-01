evaluate_insample <- function(datt, model, model_name, positive_level, negative_level) {
  # Predicts the probability of friendships using the specified model and calculates
  # the Receiver Operating Characteristic (ROC) curve to assess model performance.
  #
  # Parameters:
  # - datt: Dataset containing the actual friendship links and features used by the model.
  # - model: The predictive model to be evaluated.
  # - model_name: A descriptive name for the model, used for identification in outputs.
  # - positive_level: The label used in the dataset for positive instances (friendships).
  # - negative_level: The label used for negative instances (absence of friendship).
  #
  # The function identifies the optimal threshold for distinguishing between friendships
  # and non-friendships, creates a confusion matrix, and calculates sensitivity,
  # specificity, and balanced accuracy to quantify model performance.
  
  prob <- predict(model, type = c("response"))
  g <- roc(friend ~ prob, data = datt %>% mutate(prob = prob))
  datt$prediction <- prob
  optimal <- pROC::coords(g, "best", ret = "threshold", quiet = TRUE) %>% pull(threshold)
  datt %>% mutate(pred = as.factor(if_else(prediction <= optimal, negative_level, positive_level))) -> datt
  
  confMat <- caret::confusionMatrix(data = datt$pred, reference = datt$friend, positive = positive_level)
  sens <- caret::sensitivity(data = datt$pred, reference = datt$friend, positive = positive_level, negative = negative_level)
  spec <- caret::specificity(data = datt$pred, reference = datt$friend, positive = positive_level, negative = negative_level)
  
  balanced_accuracy <- (sens + spec) / 2
  list(confusionMatrix = confMat, res = tibble(sensitivity = sens, specificity = spec, balancedAccuracy = balanced_accuracy, model = model_name))
}

# Performs k-fold cross-validation to evaluate the model's performance in predicting
# signed links. This method splits the dataset into training and testing subsets, 
# fits the model on the training data, and evaluates its performance on the test data.
#
# Parameters:
# - dat: The dataset used for model training and evaluation.
# - formula: The formula representing the model to be trained.
# - split: The proportion of data to be used for training.
# - model_name: A descriptive name for the model.
# - seed: A seed for random number generation, ensuring reproducibility.
#
# The function enforces preconditions on the dataset to ensure it's properly formatted
# for model training and evaluation. It then partitions the data, trains the model using
# logistic regression within a repeated cross-validation framework, and calculates
# a confusion matrix, sensitivity, specificity, and balanced accuracy to assess performance.

make_KfoldCV <- function(dat, formula, split, model_name, seed = 44) {
  # Preconditions checking
  stopifnot("There must be a column called 'friend'" = "friend" %in% colnames(dat))
  stopifnot("'friend' must be a factor" = is.factor(dat$friend))
  stopifnot("'friend' must have True, False as levels" = all(levels(dat$friend) %in% c("False", "True")))
  
  set.seed(seed)
  index <- createDataPartition(dat$friend, p = split, list = FALSE)
  train_data <- dat[index, ]
  test_data  <- dat[-index, ]
  
  ctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 20, verboseIter = FALSE, sampling = "down")
  
  set.seed(seed)
  model_under <- train(formula, data = train_data, trControl = ctrl, method = "glm", family = binomial())
  
  final_under <- data.frame(actual = test_data$friend, predict(model_under, newdata = test_data, type = "prob"))
  final_under$predict <- ifelse(final_under$True > 0.5, "True", "False")
  final_under %>% mutate(predict = as.factor(predict)) -> final_under
  
  confMat <- caret::confusionMatrix(final_under$predict, test_data$friend, positive = "True")
  sens <- caret::sensitivity(data = final_under$predict, reference = test_data$friend, positive = "True", negative = "False")
  spec <- caret::specificity(data = final_under$predict, reference = test_data$friend, positive = "True", negative = "False")
  
  balanced_accuracy <- (sens + spec) / 2
  list(confusionMatrix = confMat, res = tibble(sensitivity = sens, specificity = spec, balancedAccuracy = balanced_accuracy, model = model_name))
}