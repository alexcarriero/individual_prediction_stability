library(arrow)
library(data.table)
library(dplyr)

# Input: downloaded RDS files; source scenario numbers remain unchanged.
# Output: v2 Parquet datasets and manifests identified by dgm and sample_size.
# Run this script before calculate_evaluation_criteria.R.
# Reruns add new methods. To include new iterations in an existing method,
# select its source scenarios/methods and set refresh_existing = TRUE.

# Processing selection ------------------------------------------------------

source_design <- data.table(
  source_scenario = 1:4,
  dgm = c("linear", "linear", "nonlinear", "nonlinear"),
  sample_size = c(458L, 229L, 458L, 229L)
)

# Restrict these selections when refreshing only part of the simulation.
source_scenarios <- source_design$source_scenario
methods_to_process <- NULL # NULL includes all available methods.
# Otherwise use normalized labels, e.g. c("bootstrap", "cross_validation_5").
refresh_existing <- FALSE # TRUE rebuilds selected existing partitions.

if (!length(source_scenarios) ||
    anyNA(source_scenarios) ||
    anyDuplicated(source_scenarios) ||
    !all(source_scenarios %in% source_design$source_scenario)) {
  stop("Select unique source scenario numbers from 1:4.")
}

# Input roots ---------------------------------------------------------------

simulation_root <- here::here(
  "stability_project-updated_code_AC",
  "scripts"
)

resampling_source_root <- file.path(
  simulation_root,
  "2_process_results",
  "raw_results"
)

original_source_root <- file.path(
  simulation_root,
  "0_set_up",
  "results",
  "predictions"
)

hyperparameter_source_root <- file.path(
  simulation_root,
  "0_set_up",
  "results",
  "hyperparameters"
)

# Output roots --------------------------------------------------------------

resampling_output_root <- here::here("results", "parquet_v2")
original_output_root <- here::here("results", "original_predictions_v2")
hyperparameter_output_root <- here::here("results", "hyperparameters_v2")

resampling_manifest_path <- here::here(
  "results",
  "parquet_v2_manifest.parquet"
)
original_manifest_path <- here::here(
  "results",
  "original_predictions_v2_manifest.parquet"
)
hyperparameter_manifest_path <- here::here(
  "results",
  "hyperparameters_v2_manifest.parquet"
)

# Shared helpers ------------------------------------------------------------

# Replace only a fully written partition; update its manifest immediately.
write_partition <- function(result, destination_file, manifest_row,
                            manifest_path, ...) {
  dir.create(dirname(destination_file), recursive = TRUE, showWarnings = FALSE)
  temporary_file <- tempfile("partition-", tmpdir = dirname(destination_file))
  on.exit(unlink(temporary_file), add = TRUE)
  write_parquet(result, temporary_file, compression = "zstd", ...)

  if (!file.rename(temporary_file, destination_file)) {
    stop("Could not replace partition: ", destination_file)
  }

  manifest <- if (file.exists(manifest_path)) {
    as.data.table(read_parquet(manifest_path))
  } else {
    data.table()
  }

  if (nrow(manifest)) {
    manifest <- manifest[destination != destination_file]
  }
  manifest <- rbindlist(list(manifest, manifest_row), use.names = TRUE)
  setorderv(manifest, intersect(c("dgm", "sample_size", "method", "model"),
                               names(manifest)))

  temporary_manifest <- tempfile("manifest-", tmpdir = dirname(manifest_path))
  on.exit(unlink(temporary_manifest), add = TRUE)
  write_parquet(manifest, temporary_manifest)
  if (!file.rename(temporary_manifest, manifest_path)) {
    stop("Could not replace manifest: ", manifest_path)
  }
}

partition_is_current <- function(destination_file, manifest_path) {
  if (refresh_existing || !file.exists(destination_file)) return(FALSE)
  if (!file.exists(manifest_path) ||
      !destination_file %in% read_parquet(manifest_path)$destination) {
    stop("Partition has no manifest entry; use refresh_existing = TRUE: ",
         destination_file)
  }
  message("Keeping existing partition: ", destination_file)
  TRUE
}

validate_columns <- function(result, required_columns, path) {
  missing_columns <- setdiff(required_columns, names(result))

  if (length(missing_columns) > 0L) {
    stop(
      basename(path),
      " is missing: ",
      paste(missing_columns, collapse = ", ")
    )
  }
}

validate_iteration <- function(result, expected_iter, path) {
  if (!identical(unique(as.integer(result$iter)), expected_iter)) {
    stop("Iteration does not match filename: ", basename(path))
  }
}

normalize_model <- function(x) {
  result <- fifelse(
    x %chin% c("random forest", "random_forest"),
    "random_forest",
    x
  )

  if (!all(result %chin% c("glmnet", "random_forest"))) {
    stop(
      "Unknown model value: ",
      paste(
        setdiff(unique(result), c("glmnet", "random_forest")),
        collapse = ", "
      )
    )
  }

  result
}

normalize_method <- function(x) {
  labels <- c(
    "bootstrap_635" = "bootstrap_635",
    "bootstrap_par" = "bootstrap_par",
    "bootstrap" = "bootstrap",
    "cross_validation" = "cross_validation_10",
    "cross_validation_2" = "cross_validation_2",
    "cross_validation2" = "cross_validation_2",
    "crossvalidation2" = "cross_validation_2",
    "cross_validation_5" = "cross_validation_5",
    "cross_validation5" = "cross_validation_5",
    "crossvalidation5" = "cross_validation_5",
    "independent" = "independent"
  )

  result <- unname(labels[x])

  if (anyNA(result)) {
    stop(
      "Unknown method value: ",
      paste(unique(x[is.na(result)]), collapse = ", ")
    )
  }

  result
}

# Resampling predictions ----------------------------------------------------

parse_resampling_filename <- function(path) {
  filename <- basename(path)
  matched <- regexec(
    "^scenario_([0-9]+)_(.+)_([0-9]+)[.]RDS$",
    filename
  )
  values <- regmatches(filename, matched)[[1]]

  if (length(values) != 4L) {
    stop("Unexpected resampling filename: ", filename)
  }

  data.table(
    path = path,
    source_scenario = as.integer(values[2]),
    filename_method = values[3],
    filename_iter = as.integer(values[4])
  )
}

read_resampling_result <- function(path, expected_iter) {
  result <- as.data.table(readRDS(path))
  required_columns <- c(
    "id",
    "true_p",
    "y",
    "pred_p",
    "model",
    "re_sample",
    "method",
    "iter"
  )

  validate_columns(result, required_columns, path)
  validate_iteration(result, expected_iter, path)

  if (uniqueN(result$method) != 1L) {
    stop("Multiple methods found in: ", basename(path))
  }

  y_character <- as.character(result$y)

  if (!all(y_character %chin% c("0", "1"))) {
    stop("Unexpected outcome values in: ", basename(path))
  }

  result[, `:=`(
    id = as.integer(id),
    true_p = as.numeric(true_p),
    y = as.integer(y_character),
    pred_p = as.numeric(pred_p),
    model = normalize_model(model),
    re_sample = as.integer(re_sample),
    method = normalize_method(method),
    iter = as.integer(iter)
  )]

  result
}

# List only direct children, excluding the duplicate cv2/cv5 folders.
resampling_file_metadata <- rbindlist(lapply(source_scenarios, function(source_scenario_value) {
  scenario_directory <- file.path(
    resampling_source_root,
    paste0("scenario_", source_scenario_value)
  )

  files <- list.files(
    scenario_directory,
    pattern = "[.]RDS$",
    full.names = TRUE,
    recursive = FALSE
  )

  metadata <- rbindlist(lapply(files, parse_resampling_filename))
  if (!nrow(metadata)) {
    message("No resampling files for source scenario ", source_scenario_value)
    return(NULL)
  }
  if (any(metadata$source_scenario != source_scenario_value) ||
      anyDuplicated(metadata, by = c("filename_method", "filename_iter"))) {
    stop("Unexpected source scenario or duplicate iteration in ", scenario_directory)
  }
  metadata
}))

if (!nrow(resampling_file_metadata)) stop("No selected resampling files found.")

setorder(
  resampling_file_metadata,
  source_scenario,
  filename_method,
  filename_iter
)

resampling_file_metadata[,
  batch := paste(source_scenario, filename_method, sep = "::")
]

filename_method_labels <- c(
  bootstrap635 = "bootstrap_635", bootstrapPAR = "bootstrap_par",
  bootstrapping = "bootstrap", crossvalidation = "cross_validation_10",
  crossvalidation2 = "cross_validation_2", crossvalidation5 = "cross_validation_5",
  independent = "independent"
)

for (batch_name in unique(resampling_file_metadata$batch)) {
  batch <- resampling_file_metadata[batch == batch_name]
  setorder(batch, filename_iter)

  message(
    "Converting resampling predictions: source scenario ",
    batch$source_scenario[1],
    ", method ",
    batch$filename_method[1]
  )

  design <- source_design[source_scenario == batch$source_scenario[1]]
  dgm_value <- design$dgm
  sample_size_value <- design$sample_size
  expected_method <- unname(filename_method_labels[batch$filename_method[1]])
  if (is.na(expected_method)) stop("Unknown filename method: ", batch$filename_method[1])
  if (!is.null(methods_to_process) && !expected_method %chin% methods_to_process) next

  destinations <- file.path(
    resampling_output_root, paste0("dgm=", dgm_value),
    paste0("sample_size=", sample_size_value), paste0("method=", expected_method),
    paste0("model=", c("glmnet", "random_forest")), "results.parquet"
  )
  current_partitions <- vapply(destinations, partition_is_current, logical(1),
                               manifest_path = resampling_manifest_path)
  if (all(current_partitions)) next

  pieces <- Map(
    read_resampling_result,
    path = batch$path,
    expected_iter = batch$filename_iter
  )

  combined <- rbindlist(pieces, use.names = TRUE)
  rm(pieces)
  gc()

  if (uniqueN(combined$method) != 1L) {
    stop("Filename batch contains multiple normalized methods.")
  }

  method_value <- unique(combined$method)
  if (!identical(method_value, expected_method)) {
    stop("Filename method does not match result method.")
  }

  for (model_value in c("glmnet", "random_forest")) {
    if (current_partitions[match(model_value, c("glmnet", "random_forest"))]) next
    model_results <- combined[
      model == model_value,
      .(id, true_p, y, pred_p, re_sample, iter)
    ]

    if (any(model_results[, uniqueN(id), by = iter]$V1 != sample_size_value) ||
        anyDuplicated(model_results, by = c("iter", "id", "re_sample"))) {
      stop("Unexpected patient count or duplicate prediction key.")
    }
    setorder(model_results, id, iter, re_sample)

    destination <- file.path(
      resampling_output_root,
      paste0("dgm=", dgm_value),
      paste0("sample_size=", sample_size_value),
      paste0("method=", method_value),
      paste0("model=", model_value)
    )
    destination_file <- file.path(destination, "results.parquet")

    manifest_row <- data.table(
      dgm = dgm_value,
      sample_size = sample_size_value,
      method = method_value,
      model = model_value,
      source_files = nrow(batch),
      rows = nrow(model_results),
      patients = uniqueN(model_results$id),
      iterations = uniqueN(model_results$iter),
      destination = destination_file
    )
    write_partition(
      model_results, destination_file, manifest_row, resampling_manifest_path,
      chunk_size = 200000, use_dictionary = TRUE, write_statistics = TRUE
    )

    rm(model_results)
    gc()
  }

  rm(combined)
  gc()
}

# Original predictions ------------------------------------------------------

parse_original_filename <- function(path) {
  filename <- basename(path)
  matched <- regexec(
    "^scenario_([0-9]+)_iter_([0-9]+)_og_results[.]RDS$",
    filename
  )
  values <- regmatches(filename, matched)[[1]]

  if (length(values) != 3L) {
    stop("Unexpected original-prediction filename: ", filename)
  }

  data.table(
    path = path,
    source_scenario = as.integer(values[2]),
    filename_iter = as.integer(values[3])
  )
}

read_original_prediction <- function(path, expected_iter) {
  result <- as.data.table(readRDS(path))
  required_columns <- c("id", "true_p", "y", "pred_p", "model", "iter")

  validate_columns(result, required_columns, path)
  validate_iteration(result, expected_iter, path)

  y_character <- as.character(result$y)

  if (!all(y_character %chin% c("0", "1"))) {
    stop("Unexpected outcome values in: ", basename(path))
  }

  result[,
    .(
      id = as.integer(id),
      true_p = as.numeric(true_p),
      y = as.integer(y_character),
      original_pred_p = as.numeric(pred_p),
      model = normalize_model(model),
      iter = as.integer(iter)
    )
  ]
}

for (source_scenario_value in source_scenarios) {
  design <- source_design[source_scenario == source_scenario_value]
  dgm_value <- design$dgm
  sample_size_value <- design$sample_size
  destination_file <- file.path(
    original_output_root, paste0("dgm=", dgm_value),
    paste0("sample_size=", sample_size_value), "original_predictions.parquet"
  )
  if (partition_is_current(destination_file, original_manifest_path)) next
  scenario_directory <- file.path(
    original_source_root,
    paste0("scenario_", source_scenario_value)
  )

  files <- list.files(
    scenario_directory,
    pattern = "_og_results[.]RDS$",
    full.names = TRUE,
    recursive = FALSE
  )

  file_metadata <- rbindlist(lapply(files, parse_original_filename))
  if (!nrow(file_metadata)) stop("Missing original predictions: ", scenario_directory)
  if (any(file_metadata$source_scenario != source_scenario_value) ||
      anyDuplicated(file_metadata$filename_iter)) {
    stop("Unexpected source scenario or duplicate original iteration.")
  }
  setorder(file_metadata, filename_iter)

  message("Converting original predictions: ", dgm_value, ", n = ", sample_size_value)

  original_predictions <- rbindlist(Map(
    read_original_prediction,
    path = file_metadata$path,
    expected_iter = file_metadata$filename_iter
  ))

  setorder(original_predictions, model, iter, id)

  if (any(original_predictions[, uniqueN(id), by = .(model, iter)]$V1 != sample_size_value) ||
      anyDuplicated(original_predictions, by = c("model", "iter", "id")) ||
      !setequal(original_predictions$model, c("glmnet", "random_forest"))) {
    stop("Unexpected original patient count, models, or duplicate prediction key.")
  }
  manifest_row <- data.table(
    dgm = dgm_value,
    sample_size = sample_size_value,
    source_files = nrow(file_metadata),
    rows = nrow(original_predictions),
    patients = uniqueN(original_predictions$id),
    iterations = uniqueN(original_predictions$iter),
    destination = destination_file
  )
  write_partition(
    original_predictions, destination_file, manifest_row, original_manifest_path,
    use_dictionary = TRUE, write_statistics = TRUE
  )

  rm(original_predictions)
  gc()
}

# Original-model hyperparameters -------------------------------------------

parse_hyperparameter_filename <- function(path) {
  filename <- basename(path)
  matched <- regexec(
    "^scenario_([0-9]+)_iter_([0-9]+)_hps[.]RDS$",
    filename
  )
  values <- regmatches(filename, matched)[[1]]

  if (length(values) != 3L) {
    stop("Unexpected hyperparameter filename: ", filename)
  }

  data.table(
    path = path,
    source_scenario = as.integer(values[2]),
    filename_iter = as.integer(values[3])
  )
}

hyperparameter_models <- data.table(
  source_model = c("elastic_net", "random_forest"),
  model = c("glmnet", "random_forest")
)

for (source_scenario_value in source_scenarios) {
  design <- source_design[source_scenario == source_scenario_value]
  dgm_value <- design$dgm
  sample_size_value <- design$sample_size

  for (model_number in seq_len(nrow(hyperparameter_models))) {
    source_model <- hyperparameter_models$source_model[model_number]
    model_value <- hyperparameter_models$model[model_number]
    destination_file <- file.path(
      hyperparameter_output_root, paste0("dgm=", dgm_value),
      paste0("sample_size=", sample_size_value), paste0("model=", model_value),
      "hyperparameters.parquet"
    )
    if (partition_is_current(destination_file, hyperparameter_manifest_path)) next

    scenario_directory <- file.path(
      hyperparameter_source_root, source_model,
      paste0("scenario_", source_scenario_value)
    )
    files <- list.files(scenario_directory, pattern = "_hps[.]RDS$",
                        full.names = TRUE, recursive = FALSE)
    file_metadata <- rbindlist(lapply(files, parse_hyperparameter_filename))
    if (!nrow(file_metadata)) stop("Missing hyperparameters: ", scenario_directory)
    if (any(file_metadata$source_scenario != source_scenario_value) ||
        anyDuplicated(file_metadata$filename_iter)) {
      stop("Unexpected source scenario or duplicate hyperparameter iteration.")
    }
    setorder(file_metadata, filename_iter)
    message("Converting hyperparameters: ", dgm_value, ", n = ",
            sample_size_value, ", model ", model_value)

    hyperparameters <- rbindlist(lapply(seq_len(nrow(file_metadata)), function(i) {
      metadata <- file_metadata[i]
      result <- as.data.table(readRDS(metadata$path))
      if (nrow(result) != 1L) {
        stop("Expected one hyperparameter row in: ", basename(metadata$path))
      }
      required_columns <- if (model_value == "glmnet") {
        c("alpha", "lambda", "iter")
      } else {
        c("mtry", "min.node.size", "iter")
      }
      validate_columns(result, required_columns, metadata$path)
      validate_iteration(result, metadata$filename_iter, metadata$path)
      result
    }))

    hyperparameters <- if (model_value == "glmnet") {
      hyperparameters[, .(iter = as.integer(iter), alpha = as.numeric(alpha),
                           lambda = as.numeric(lambda))]
    } else {
      hyperparameters[, .(iter = as.integer(iter), mtry = as.integer(mtry),
                           min_node_size = as.integer(min.node.size))]
    }
    setorder(hyperparameters, iter)
    manifest_row <- data.table(
      dgm = dgm_value, sample_size = sample_size_value, model = model_value,
      source_files = nrow(file_metadata), rows = nrow(hyperparameters),
      iterations = uniqueN(hyperparameters$iter), destination = destination_file
    )
    write_partition(
      hyperparameters, destination_file, manifest_row, hyperparameter_manifest_path,
      use_dictionary = TRUE, write_statistics = TRUE
    )
  }
}

# Conversion summaries -----------------------------------------------------

print(read_parquet(resampling_manifest_path))
print(read_parquet(original_manifest_path))
print(read_parquet(hyperparameter_manifest_path))

# Confirm the expected resampling partition sizes without loading all rows.
resampling_counts <- open_dataset(resampling_output_root) |>
  group_by(dgm, sample_size, method, model) |>
  summarise(
    rows = n(),
    patients = n_distinct(id),
    iterations = n_distinct(iter)
  ) |>
  collect()

print(resampling_counts)
