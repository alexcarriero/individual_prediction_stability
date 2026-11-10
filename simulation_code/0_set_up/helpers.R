# GENERATE LINEAR DATA 

generate_data <- function (sample_size, seed){
  
  set.seed(seed)
  
  # variable settings ----------------------------------------------------------
  # none
  
  # fixed settings -------------------------------------------------------------
  n_pred = 12 
  correlation = 0.2
  b = 0.1610
  a = -1.4649
  
  # design matrix --------------------------------------------------------------
  sig <- matrix(correlation, n_pred, n_pred)
  diag(sig) <- rep(1, n_pred)
  mus <- rep(0, n_pred)
  
  x_data <- mvrnorm(n = sample_size, mu = mus, Sigma = sig)
  colnames(x_data) <- paste0("X", 1:n_pred)
  
  # generate linear predictor --------------------------------------------------
  
  # weak: X1..X6
  s1 <- b * (x_data[,1] + x_data[,2] + x_data[,3] + x_data[,4] + x_data[,5] + x_data[,6])
  
  # strong: X7..X9
  s2 <- 3*b * (x_data[,7] + x_data[,8] + x_data[,9])
  
  # zero: X10..X12 omitted
  eta <- a + s1 + s2
  
  # generate outcome -----------------------------------------------------------
  p <- plogis(eta)
  y <- rbinom(sample_size, size = 1, prob = p)
  
  # return data frame ----------------------------------------------------------
  df <- x_data %>% as.data.frame()
  df$true_p = p 
  df$y = y %>% as.factor()
  df$id = c(1:nrow(df))
  
  return (df)
}


# GENERATE NON-LINEAR DATA 

generate_data_nl <- function (sample_size, seed){
  
  set.seed(seed)
  
  # variable settings ----------------------------------------------------------
  # none
  
  # fixed settings -------------------------------------------------------------
  n_pred = 12 
  correlation = 0.2
  b = 0.1610
  a = -1.4649
  
  # design matrix --------------------------------------------------------------
  sig <- matrix(correlation, n_pred, n_pred)
  diag(sig) <- rep(1, n_pred)
  mus <- rep(0, n_pred)
  
  x_data <- mvrnorm(n = sample_size, mu = mus, Sigma = sig)
  colnames(x_data) <- paste0("X", 1:n_pred)
  
  # generate linear predictor --------------------------------------------------
  
  # weak: X1..X6
  s1 <- b * (x_data[,1] * x_data[,2] - x_data[,3] * x_data[,4] + x_data[,5]^4 + x_data[,6]^4)
  
  # strong: X7..X9
  s3 <-  3*b * x_data[,7] * x_data[,8]
  s4 <- -3*b * x_data[,9]^4
  
  # zero: X10..X12 omitted
  
  eta <- a + s1 + s3 + s4
  
  # generate outcome -----------------------------------------------------------
  p <- plogis(eta)
  y <- rbinom(sample_size, size = 1, prob = p)
  
  # return data frame ----------------------------------------------------------
  df <- x_data %>% as.data.frame()
  df$true_p = p 
  df$y = y %>% as.factor()
  df$id = c(1:nrow(df))
  
  return (df)
}

rf_fit <- function(data_set) {
  
  # prepare data
  df <- data_set
  X  <- df[, grep("^X", names(df))]
  y  <- df$y
  
  fold_id <- sample(rep(1:10, length.out = nrow(X)))
  
  # random search grid
  grid_size = 30
  
  grid <- data.frame(
    mtry           = sample(1:ncol(X), grid_size, replace = TRUE),
    min.node.size  = sample(1:20, grid_size, replace = TRUE)
  )
  
  # evaluate each combination via 10-fold CV
  scores <- sapply(1:grid_size, function(g) {
    
    fold_logloss <- sapply(1:10, function(k) {
      
      train <- df[fold_id != k, ]
      test  <- df[fold_id == k, ]
      
      fit <- ranger::ranger(
        y             ~ .,
        data          = train[, c(grep("^X", names(train), value = TRUE), "y")],
        num.trees     = 500,
        mtry          = grid$mtry[g],
        min.node.size = grid$min.node.size[g],
        probability   = TRUE
      )
      
      pred_p <- predict(fit, data = test)$predictions[, "1"]
      
      # make sure pred_p is never exactly 0 or 1
      eps <- 1e-15
      pred_p <- pmax(pmin(pred_p, 1 - eps), eps)
      
      # log loss
      y_test <- as.numeric(test$y) - 1
      log_loss = -mean(y_test * log(pred_p) + (1 - y_test) * log(1 - pred_p))
      return(log_loss)
    })
    
    mean(fold_logloss)
  })
  
  # best hyperparameters
  best_idx <- which.min(scores)
  mtry     <- grid$mtry[best_idx]
  min_node <- grid$min.node.size[best_idx]
  
  # fit final model on full data
  fit <- ranger::ranger(
    y             ~ .,
    data          = df[, c(grep("^X", names(df), value = TRUE), "y")],
    num.trees     = 500,
    mtry          = mtry,
    min.node.size = min_node,
    probability   = TRUE,
  )
  
  
  return(fit)
}

# ------------------------------------------------------------------------------

# FIT MODELS 

# GLM 

do_glm <- function(data_set, seed, iter){
  
  # set seed 
  set.seed(seed)
  
  # prepare data
  df <- data_set
  X  <- as.matrix(df[, grep("^X", names(df))])
  y  <- as.numeric(df$y) - 1
  
  # alpha values 
  alpha_seq <- seq(0, 1, by = 0.1)
  
  # fold id 
  fold_id <- sample(rep(1:10, length.out = nrow(X)))
  
  # cross validation fits
  cv_fits <- lapply(alpha_seq, 
                    function(alpha) {glmnet::cv.glmnet(X, y, alpha = alpha, 
                                                       foldid = fold_id, 
                                                       family = "binomial", 
                                                       type.measure = "deviance")})
  
  # find best hyper-parameters 
  
  # which alpha had the lowest CV error
  min_errors <- sapply(cv_fits, function(fit) min(fit$cvm))
  best_idx   <- which.min(min_errors)
  
  # save alpha and lambda
  alpha  <- alpha_seq[best_idx]
  lambda <- cv_fits[[best_idx]]$lambda.min
  
  # store hyperparameters 
  hyperparameters <- tibble(
    alpha  = alpha, 
    lambda = lambda, 
    iter = iter
  )
  
  # fit model 
  fit <- glmnet(X, y, 
                alpha = alpha, 
                lambda = lambda, 
                family = "binomial", 
                maxit = 100000)
  
  # store predictions
  df$pred_p <- as.numeric(predict(fit, newx = X, s = lambda, type = "response"))
  df$model  <- "glmnet"
  df$iter   <- iter
  
  results <- df %>% dplyr::select("id", "true_p", "y", "pred_p", "model", "iter")
  output <- list("results" = results, "hyperparameters" = hyperparameters, "model" = fit)
  
  return (output)
}

# specify different model for non-linear settings 

do_glm_nl <- function(data_set, seed, iter){
  
  # set seed 
  set.seed(seed)
  
  # prepare data
  df <- data_set
  X  <- model.matrix(~ I(X1*X2) + I(X3*X4) + I(X5^4) + I(X6^4) + I(X7*X8) + I(X9^4) +
                       X10 + X11 + X12, data = df)[, -1]
  y  <- as.numeric(df$y) - 1
  
  # alpha values 
  alpha_seq <- seq(0, 1, by = 0.1)
  
  # fold id 
  fold_id <- sample(rep(1:10, length.out = nrow(X)))
  
  # cross validation fits
  cv_fits <- lapply(alpha_seq, 
                    function(alpha) {glmnet::cv.glmnet(X, y, alpha = alpha, 
                                                       foldid = fold_id, 
                                                       family = "binomial", 
                                                       type.measure = "deviance")})
  
  # find best hyper-parameters 
  
  # which alpha had the lowest CV error
  min_errors <- sapply(cv_fits, function(fit) min(fit$cvm))
  best_idx   <- which.min(min_errors)
  
  # save alpha and lambda
  alpha  <- alpha_seq[best_idx]
  lambda <- cv_fits[[best_idx]]$lambda.min
  
  # store hyperparameters 
  hyperparameters <- tibble(
    alpha  = alpha, 
    lambda = lambda, 
    iter = iter
  )
  
  # fit model 
  fit <- glmnet(X, y, 
                alpha = alpha, 
                lambda = lambda, 
                family = "binomial", 
                maxit = 100000)
  
  # store predictions
  df$pred_p <- as.numeric(predict(fit, newx = X, s = lambda, type = "response"))
  df$model  <- "glmnet"
  df$iter   <- iter
  
  results <- df %>% dplyr::select("id", "true_p", "y", "pred_p", "model", "iter")
  output <- list("results" = results, "hyperparameters" = hyperparameters, "model" = fit)
  
  return (output)
}

# RANDOM FOREST 

do_rf <- function(data_set, seed, iter) {
  
  # prepare data
  df <- data_set
  X  <- df[, grep("^X", names(df))]
  y  <- df$y
  
  # set seed and create folds
  set.seed(seed)
  fold_id <- sample(rep(1:10, length.out = nrow(X)))
  
  # random search grid
  grid_size = 30
  
  grid <- data.frame(
    mtry           = sample(1:ncol(X), grid_size, replace = TRUE),
    min.node.size  = sample(1:20, grid_size, replace = TRUE)
  )
  
  # evaluate each combination via 10-fold CV
  scores <- sapply(1:grid_size, function(g) {
    
    fold_logloss <- sapply(1:10, function(k) {
      
      train <- df[fold_id != k, ]
      test  <- df[fold_id == k, ]
      
      fit <- ranger::ranger(
        y             ~ .,
        data          = train[, c(grep("^X", names(train), value = TRUE), "y")],
        num.trees     = 500,
        mtry          = grid$mtry[g],
        min.node.size = grid$min.node.size[g],
        probability   = TRUE,
        seed          = seed
      )
      
      pred_p <- predict(fit, data = test)$predictions[, "1"]
      
      # make sure pred_p is never exactly 0 or 1
      eps <- 1e-15
      pred_p <- pmax(pmin(pred_p, 1 - eps), eps)
      
      # log loss
      y_test <- as.numeric(test$y) - 1
      log_loss = -mean(y_test * log(pred_p) + (1 - y_test) * log(1 - pred_p))
      return(log_loss)
    })
    
    mean(fold_logloss)
  })
  
  # best hyperparameters
  best_idx <- which.min(scores)
  mtry     <- grid$mtry[best_idx]
  min_node <- grid$min.node.size[best_idx]
  
  # fit final model on full data
  fit <- ranger::ranger(
    y             ~ .,
    data          = df[, c(grep("^X", names(df), value = TRUE), "y")],
    num.trees     = 500,
    mtry          = mtry,
    min.node.size = min_node,
    probability   = TRUE,
    seed          = seed
  )
  
  # store predictions
  df$pred_p <- predict(fit, data = df)$predictions[, "1"]
  df$model  <- "random_forest"
  df$iter   <- iter
  
  results <- df %>% dplyr::select("id", "true_p", "y", "pred_p", "model", "iter")
  
  hyperparameters <- tibble(
    mtry          = mtry,
    min.node.size = min_node,
    iter          = iter
  )
  
  return(list("results" = results, "hyperparameters" = hyperparameters, "model" = fit))
}

# Process Original Results 

calculate_pms <- function(results){
  
  results %>% 
    group_by(model, iter) %>% 
    summarise(
      auc    = as.numeric(pROC::auc(pROC::roc(y, pred_p, quiet = TRUE))), 
      brier  = mean((pred_p - (as.numeric(y) - 1))^2), 
      .groups = "drop"
    )
  
}
