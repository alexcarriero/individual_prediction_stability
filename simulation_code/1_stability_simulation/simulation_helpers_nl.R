# NON-LINEAR DATA 

# save glm formula 
glm_formula <- ~ I(X1*X2) + I(X3*X4) + I(X5^4) + I(X6^4) + I(X7*X8) + I(X9^4) + X10 + X11 + X12

# GENERATE DATA 
generate_data_sim <- function (sample_size, seed){
  
  set.seed(seed)
  
  # variable settings ----------------------------------------------------------
  # none
  
  # fixed settings -------------------------------------------------------------
  n_pred = 12 
  correlation = 0.2
  
  # design matrix --------------------------------------------------------------
  sig <- matrix(correlation, n_pred, n_pred)
  diag(sig) <- rep(1, n_pred)
  mus <- rep(0, n_pred)
  
  x_data <- mvrnorm(n = sample_size, mu = mus, Sigma = sig)
  colnames(x_data) <- paste0("X", 1:n_pred)
  
  # return data frame ----------------------------------------------------------
  df <- x_data %>% as.data.frame() 
  
  # read in data generating model ----------------------------------------------
  rf_dgm <- readRDS("rf_dgm.RDS")
  
  rf_p  = predict(rf_dgm, data = df)$predictions[, "1"]
  rf_y  = rbinom(nrow(df), size = 1, prob = rf_p) %>% as.factor()
  
  nl_data = bind_cols(df, "true_p" = rf_p, "y" = rf_y, "id" = seq_len(sample_size))
  
  return (nl_data)
}

# FIT MODELS ------------------------------------------------------------------

# GLM -------------------------------------------------------------------------

do_glm_sim <- function(training_data, original_data, j){
  
  # prepare data
  X_train <- model.matrix(glm_formula, data = training_data)[, -1]
  y_train <- as.numeric(training_data$y) - 1
  X_pred  <- model.matrix(glm_formula, data = original_data)[, -1]
  
  # alpha values 
  alpha_seq <- seq(0, 1, by = 0.1)
  
  # fold id 
  fold_id <- sample(rep(1:10, length.out = nrow(X_train)))
  
  # cross validation fits
  cv_fits <- lapply(alpha_seq, 
                    function(alpha) {glmnet::cv.glmnet(X_train, y_train, 
                                                       alpha = alpha, 
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
  
  # fit model 
  fit <- glmnet(X_train, y_train, 
                alpha = alpha, 
                lambda = lambda, 
                family = "binomial", 
                maxit = 100000)
  
  # store predictions
results <- original_data %>% 
  dplyr::select(id, true_p, y) %>% 
  mutate(
    pred_p = predict(fit, newx = X_pred, s = lambda, type = "response") %>% as.numeric(), 
    model = "glmnet", 
    re_sample = j
  )
  
  return (results)
}

# RANDOM FOREST ---------------------------------------------------------------

do_rf_sim <- function(training_data, original_data, j) {
  
  # set seed 
  seed = runif(1, min = 0, max = 1000000000) %>% floor()
  set.seed(seed)
  
  # prepare data
  X_train  <- training_data[, grep("^X", names(training_data))]
  
  # create folds
  fold_id <- sample(rep(1:10, length.out = nrow(X_train)))
  
  # random search grid
  grid_size = 30
  
  grid <- data.frame(
    mtry           = sample(1:ncol(X_train), grid_size, replace = TRUE),
    min.node.size  = sample(1:20, grid_size, replace = TRUE)
  )
  
  # evaluate each combination via 10-fold CV
  scores <- sapply(1:grid_size, function(g) {
    
    fold_logloss <- sapply(1:10, function(k) {
      
      train <- training_data[fold_id != k, ]
      test  <- training_data[fold_id == k, ]
      
      fit <- ranger::ranger(
        y             ~ .,
        data          = train[, c(grep("^X", names(train), value = TRUE), "y")],
        num.trees     = 500,
        mtry          = grid$mtry[g],
        min.node.size = grid$min.node.size[g],
        probability   = TRUE,
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
    data          = training_data[, c(grep("^X", names(training_data), value = TRUE), "y")],
    num.trees     = 500,
    mtry          = mtry,
    min.node.size = min_node,
    probability   = TRUE,
  )
  
  # store predictions
results <- original_data %>% 
  dplyr::select(id, true_p, y) %>% 
  mutate(
    pred_p = predict(fit, data = original_data)$predictions[, "1"],
    model = "random forest",
    re_sample = j
    
  )

  return(results)
}


# Extra functions for parametric bootstrap -------------------------------------

glm_mod <- function(data_set){
  
  # prepare data
  df <- data_set
  X  <- model.matrix(glm_formula, data = df)[, -1]
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
  
  # fit model 
  fit <- glmnet(X, y, 
                alpha = alpha, 
                lambda = lambda, 
                family = "binomial", 
                maxit = 100000)
  
  return (fit)
}

rf_mod <- function(data_set) {
  
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

# SIM SEEDS --------------------------------------------------------------------

sim_seeds = list(
  independent  = c(1:200), 
  bootstrap = c(201:400), 
  cross_val = c(401:600), 
  boot_par  = c(601:800),
  boot_635  = c(801:2800)
)

# FUNCTION FOR ALL ------------------------------------------------------------- 

assess_stability <- function(dataset, method, scenario, iter, sim_seeds){
  
  # set up 
  sc = scenario 
  og_data = dataset
  sim_seeds = sim_seeds
  
  results = c()
  
  # independent ----------------------------------------------------------------
  if(method == "independent"){
    
    s_indep = sim_seeds$independent
    
    for (j in c(1:200)){
      
      # generate an independent data set
      data = generate_data_sim(sample_size = nrow(og_data), seed = s_indep[j]+ iter * 100000)
      
      # train models, and predict on original individuals
      res_rf = do_rf_sim(data, og_data, j)
      res_el = do_glm_sim(data, og_data, j)
      
      res_j = bind_rows(res_rf, res_el)
      
      results = bind_rows(results, res_j)
    }
    
    # save results
    saveRDS(bind_cols(results, method = "independent", iter = iter), 
            paste0("results/scenario_", sc, "/scenario_", sc, "_independent_", i,".RDS"))
  }
  
  
  # bootstrapping --------------------------------------------------------------
  if(method == "bootstrap"){
    
    s_boot = sim_seeds$bootstrap
    
    for (j in c(1:200)){
      
      # set seed
      set.seed(s_boot[j])
      
      # re-sample data
      data = sample_n(og_data, size = nrow(og_data), replace = TRUE)
      
      # train models, and predict on original individuals
      res_rf = do_rf_sim(data, og_data, j)
      res_el = do_glm_sim(data, og_data, j)
      
      res_j = bind_rows(res_rf, res_el)
      
      results = bind_rows(results, res_j)
    }
    
    # save results
    saveRDS(bind_cols(results, method = "bootstrap", iter = iter), 
            paste0("results/scenario_", sc, "/scenario_", sc, "_bootstrapping_", i,".RDS"))
  }
  
  # cross validation -----------------------------------------------------------
  
  if(method == "cross_validation"){
    
    s_cv = sim_seeds$cross_val
    
    for (j in c(1:200)){
      
      # set seed
      set.seed(s_cv[j])
      
      # set up for cross validation 
      data = 
        og_data %>% 
        mutate(K = sample(rep(1:10, length.out = nrow(og_data))))
      
      # 10-fold cross validation
      res_cv <- c()
      
      for (k in c(1:10)){
        
        # training set vs predicting set
        training = data %>% filter(K != k)
        predicts = data %>% filter(K == k)
        
        res_rf <- do_rf_sim(training,  predicts, j)
        res_el <- do_glm_sim(training, predicts, j)
        
        res_cv <- bind_rows(res_cv, res_rf, res_el)
        
      }
      
      res_cv = res_cv %>% arrange(id)
      results = bind_rows(results, res_cv)
    }
    
    # save results
    saveRDS(bind_cols(results, method = "cross_validation", iter = iter), 
            paste0("results/scenario_", sc, "/scenario_", sc, "_crossvalidation_", i,".RDS"))
  }
  
  # cross-validation 5 FOLD 
  
  if(method == "cross_validation_5"){
    
    s_cv = sim_seeds$cross_val
    
    for (j in c(1:200)){
      
      # set seed
      set.seed(s_cv[j])
      
      # set up for cross validation 
      data = 
        og_data %>% 
        mutate(K = sample(rep(1:5, length.out = nrow(og_data))))
      
      # 10-fold cross validation
      res_cv <- c()
      
      for (k in c(1:5)){
        
        # training set vs predicting set
        training = data %>% filter(K != k)
        predicts = data %>% filter(K == k)
        
        res_rf <- do_rf_sim(training,  predicts, j)
        res_el <- do_glm_sim(training, predicts, j)
        
        res_cv <- bind_rows(res_cv, res_rf, res_el)
        
      }
      
      res_cv = res_cv %>% arrange(id)
      results = bind_rows(results, res_cv)
    }
    
    # save results
    saveRDS(bind_cols(results, method = "cross_validation5", iter = iter), 
            paste0("results/scenario_", sc, "/scenario_", sc, "_crossvalidation5_", i,".RDS"))
  }
  
  # cross-validation 2 FOLD 
  
  if(method == "cross_validation_2"){
    
    s_cv = sim_seeds$cross_val
    
    for (j in c(1:200)){
      
      # set seed
      set.seed(s_cv[j])
      
      # set up for cross validation 
      data = 
        og_data %>% 
        mutate(K = sample(rep(1:2, length.out = nrow(og_data))))
      
      # 10-fold cross validation
      res_cv <- c()
      
      for (k in c(1:2)){
        
        # training set vs predicting set
        training = data %>% filter(K != k)
        predicts = data %>% filter(K == k)
        
        res_rf <- do_rf_sim(training,  predicts, j)
        res_el <- do_glm_sim(training, predicts, j)
        
        res_cv <- bind_rows(res_cv, res_rf, res_el)
        
      }
      
      res_cv = res_cv %>% arrange(id)
      results = bind_rows(results, res_cv)
    }
    
    # save results
    saveRDS(bind_cols(results, method = "cross_validation2", iter = iter), 
            paste0("results/scenario_", sc, "/scenario_", sc, "_crossvalidation2_", i,".RDS"))
  }
  
  # out-of-bag bootstrap (0.635) -----------------------------------------------
  
  if(method == "bootstrap_635"){
    
    # set up 
    j = 1
    s_b6 <- sim_seeds$boot_635
    
    # storage for predictions 
    oob_count <- integer(nrow(og_data))
    oob_preds <- c()
    
    while (min(oob_count) < 200 & j <= 2000){
      
      set.seed(s_b6[j])
      
      # re-sample data
      training = sample_n(og_data, size = nrow(og_data), replace = TRUE)
      oob = og_data %>% filter(!id %in% training$id)
      
      # train models and predict on oob observations
      res_rf <- do_rf_sim(training,  oob, j)
      res_el <- do_glm_sim(training, oob, j) 
      
      oob_preds <- bind_rows(oob_preds, res_rf, res_el)
      oob_count[oob$id] <- oob_count[oob$id] + 1
      
      j = j + 1
    }
    
    results <- oob_preds %>% 
      group_by(id, model) %>% 
      dplyr::slice(1:200) %>% 
      ungroup()
    
    saveRDS(bind_cols(results, method = "bootstrap_635", iter = iter), 
            paste0("results/scenario_", sc, "/scenario_", sc, "_bootstrap635_", iter, ".RDS"))
  }
  
  # parametric bootstrap --------------------------------------------------------
  
  if (method == "bootstrap_par"){
    
    # set up 
    set.seed(1000)
    s_bpar = sim_seeds$boot_par
    
    rf_fit  = rf_mod(og_data)
    glm_fit = glm_mod(og_data)
    
    X = og_data[, grep("^X", names(og_data))]
    
    rf_p  = predict(rf_fit, data = og_data)$predictions[, "1"]
    glm_p = predict(glm_fit, newx = model.matrix(glm_formula, data = og_data)[, -1], type = "response") %>% as.numeric()
    
    results = c()
    
    for (j in c(1:200)){
      
      set.seed(s_bpar[j])
      
      rf_y  = rbinom(nrow(X), size = 1, prob = rf_p) %>% as.factor()
      glm_y = rbinom(nrow(X), size = 1, prob = glm_p) %>% as.factor()
      
      rf_train  = bind_cols(X, "y" = rf_y)
      glm_train = bind_cols(X, "y" = glm_y)
      
      res_rf  <- do_rf_sim(rf_train,  og_data, j)
      res_glm <- do_glm_sim(glm_train, og_data, j) 
      
      results = bind_rows(results, res_rf, res_glm)
      
    }
    
    saveRDS(bind_cols(results, method = "bootstrap_par", iter = iter), 
            paste0("results/scenario_", sc, "/scenario_", sc, "_bootstrapPAR_", iter, ".RDS"))
  }
}
