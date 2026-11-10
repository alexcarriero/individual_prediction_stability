# Set up 
library(tidyverse)
library(MASS)
library(glmnet)
library(ranger)
library(pROC)

source("simulation_helpers.R")

# Scenario and iteration

sc = 2
i =  as.numeric(commandArgs(trailingOnly = TRUE)[1])

data <- readRDS(paste0("raw_data/scenario_", sc, "/sc_", sc, "_dataset_", i,".RDS"))

# Select method and run (with timing saved)

method = commandArgs(trailingOnly = TRUE)[2]

start <- proc.time()

assess_stability(dataset = data, method = method, scenario = sc, iter = i, sim_seeds)

elapsed <- proc.time() - start
cat("Iteration", i, "took", elapsed["elapsed"], "seconds\n")