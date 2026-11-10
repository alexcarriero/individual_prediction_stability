library(arrow)
library(data.table)

# Input:
# - Full resampling distributions in results/parquet_v2.
# - Original in-sample predictions in results/original_predictions_v2.
#
# Output:
# - One row per dgm, sample_size, method, model, iteration, and patient in
#   results/evaluation_v2.
#
# Run parquet_conversion.R first. Reruns add newly converted methods.
# After refreshing an existing conversion partition, select the matching
# dgm/sample_size/method here and set refresh_existing = TRUE.

results_manifest_path <- here::here(
  "results",
  "parquet_v2_manifest.parquet"
)

original_predictions_root <- here::here(
  "results",
  "original_predictions_v2"
)

output_root <- here::here("results", "evaluation_v2")
output_manifest_path <- here::here(
  "results",
  "evaluation_v2_manifest.parquet"
)

# Processing selection ------------------------------------------------------

# NULL selects all partitions in the conversion manifest.
dgms_to_process <- NULL
sample_sizes_to_process <- NULL
methods_to_process <- NULL
refresh_existing <- FALSE # TRUE recalculates selected existing partitions.

results_manifest <- as.data.table(read_parquet(results_manifest_path))
if (!is.null(dgms_to_process)) {
  results_manifest <- results_manifest[dgm %chin% dgms_to_process]
}
if (!is.null(sample_sizes_to_process)) {
  results_manifest <- results_manifest[sample_size %in% sample_sizes_to_process]
}
if (!is.null(methods_to_process)) {
  results_manifest <- results_manifest[method %chin% methods_to_process]
}
if (!nrow(results_manifest)) stop("No selected prediction partitions found.")
if (anyDuplicated(results_manifest, by = c("dgm", "sample_size", "method", "model"))) {
  stop("Duplicate prediction partition in the conversion manifest.")
}
setorder(results_manifest, dgm, sample_size, method, model)

# Original predictions ------------------------------------------------------

# Load only the DGM/sample-size combinations being evaluated.
selected_design <- unique(results_manifest[, .(dgm, sample_size)])
original_predictions <- rbindlist(lapply(seq_len(nrow(selected_design)), function(i) {
  design <- selected_design[i]
  source_file <- file.path(
    original_predictions_root,
    paste0("dgm=", design$dgm),
    paste0("sample_size=", design$sample_size),
    "original_predictions.parquet"
  )
  if (!file.exists(source_file)) stop("Missing original predictions: ", source_file)
  result <- as.data.table(read_parquet(source_file))
  result[, `:=`(dgm = design$dgm, sample_size = as.integer(design$sample_size))]
  result
}))

if (!all(original_predictions$model %chin% c("glmnet", "random_forest"))) {
  stop("Unexpected model in the original predictions.")
}

original_key <- c("dgm", "sample_size", "model", "iter", "id")

if (anyDuplicated(original_predictions, by = original_key)) {
  stop(
    "The original predictions do not have a unique dgm-sample_size-model-iter-id key."
  )
}

setkeyv(original_predictions, original_key)

# Process one DGM/sample-size/method/model partition at a time.
output_manifest <- if (file.exists(output_manifest_path)) {
  as.data.table(read_parquet(output_manifest_path))
} else {
  data.table()
}

for (partition_number in seq_len(nrow(results_manifest))) {
  partition <- results_manifest[partition_number]

  dgm_value <- partition$dgm
  sample_size_value <- partition$sample_size
  method_value <- partition$method
  model_value <- partition$model
  source_file <- partition$destination

  destination <- file.path(
    output_root, paste0("dgm=", dgm_value),
    paste0("sample_size=", sample_size_value), paste0("method=", method_value),
    paste0("model=", model_value)
  )
  destination_file <- file.path(destination, "evaluation.parquet")
  if (!refresh_existing && file.exists(destination_file)) {
    if (!nrow(output_manifest) || !destination_file %in% output_manifest$destination) {
      stop("Partition has no manifest entry; use refresh_existing = TRUE: ",
           destination_file)
    }
    message("Keeping existing evaluation: ", destination_file)
    next
  }

  message(
    "Summarising ",
    dgm_value,
    ", n = ", sample_size_value,
    ", method ",
    method_value,
    ", model ",
    model_value
  )

  # Only one full prediction partition is held in memory.
  predictions <- as.data.table(read_parquet(source_file))

  required_columns <- c(
    "id",
    "true_p",
    "y",
    "pred_p",
    "re_sample",
    "iter"
  )
  missing_columns <- setdiff(required_columns, names(predictions))

  if (length(missing_columns) > 0L) {
    stop(
      source_file,
      " is missing: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  if (
    anyNA(predictions$pred_p) ||
      anyNA(predictions$true_p) ||
      any(predictions$pred_p < 0 | predictions$pred_p > 1) ||
      any(predictions$true_p < 0 | predictions$true_p > 1)
  ) {
    stop("Missing or out-of-range probabilities in: ", source_file)
  }

  if (anyDuplicated(predictions, by = c("iter", "id", "re_sample")) ||
      any(predictions[, uniqueN(id), by = iter]$V1 != sample_size_value)) {
    stop("Duplicate prediction key or unexpected sample size: ", source_file)
  }

  # Add the original in-sample prediction used to calculate MAPE.
  partition_original_predictions <- original_predictions[
    .(dgm_value, sample_size_value, model_value)
  ][,
    .(iter, id, original_pred_p, original_true_p = true_p, original_y = y)
  ]

  predictions[
    partition_original_predictions,
    `:=`(original_pred_p = i.original_pred_p,
         original_true_p = i.original_true_p, original_y = i.original_y),
    on = .(iter, id)
  ]

  if (anyNA(predictions$original_pred_p)) {
    stop("Original predictions could not be matched for: ", source_file)
  }

  if (anyNA(predictions$original_true_p) || anyNA(predictions$original_y) ||
      any(predictions$true_p != predictions$original_true_p) ||
      any(predictions$y != predictions$original_y) ||
      any(!is.finite(predictions$original_pred_p)) ||
      any(predictions$original_pred_p < 0 | predictions$original_pred_p > 1)) {
    stop("Original predictions do not match patient data: ", source_file)
  }

  # Calculate all patient-level summaries in one grouped pass.
  evaluation <- predictions[,
    {
      prediction_quantiles <- quantile(
        pred_p,
        probs = c(0.025, 0.975),
        names = FALSE,
        type = 7
      )

      list(
        true_p = data.table::first(true_p),
        y = data.table::first(y),
        original_pred_p = data.table::first(original_pred_p),
        n_resamples = .N,
        mean_pred = mean(pred_p),
        variance = var(pred_p),
        q025 = prediction_quantiles[1],
        q975 = prediction_quantiles[2],
        pi_width = prediction_quantiles[2] - prediction_quantiles[1],
        tMAPE = mean(abs(pred_p - true_p)),
        MAPE = mean(abs(pred_p - original_pred_p)),
        n_true_p = uniqueN(true_p),
        n_y = uniqueN(y),
        n_original_pred_p = uniqueN(original_pred_p)
      )
    },
    by = .(iter, id)
  ]

  # Every patient should contribute 200 predictions with constant source data.
  invalid_groups <- evaluation[
    n_resamples != 200L |
      n_true_p != 1L |
      n_y != 1L |
      n_original_pred_p != 1L
  ]

  if (nrow(invalid_groups) > 0L) {
    stop(
      "Invalid resample count or non-constant patient data in ",
      source_file,
      ". First affected iteration/id: ",
      invalid_groups$iter[1],
      "/",
      invalid_groups$id[1]
    )
  }

  evaluation[,
    c("n_true_p", "n_y", "n_original_pred_p") := NULL
  ]

  setorder(evaluation, id, iter)

  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  temporary_file <- tempfile("evaluation-", tmpdir = destination)

  # Retain the same Hive partitions as the full prediction dataset.
  write_parquet(
    evaluation,
    temporary_file,
    compression = "zstd",
    use_dictionary = TRUE,
    write_statistics = TRUE
  )

  if (!file.rename(temporary_file, destination_file)) {
    stop("Could not replace evaluation partition: ", destination_file)
  }

  manifest_row <- data.table(
    dgm = dgm_value,
    sample_size = sample_size_value,
    method = method_value,
    model = model_value,
    rows = nrow(evaluation),
    patients = uniqueN(evaluation$id),
    iterations = uniqueN(evaluation$iter),
    min_resamples = min(evaluation$n_resamples),
    max_resamples = max(evaluation$n_resamples),
    destination = destination_file
  )

  if (nrow(output_manifest)) {
    output_manifest <- output_manifest[destination != destination_file]
  }
  output_manifest <- rbindlist(list(output_manifest, manifest_row), use.names = TRUE)
  setorder(output_manifest, dgm, sample_size, method, model)
  temporary_manifest <- tempfile("manifest-", tmpdir = dirname(output_manifest_path))
  write_parquet(output_manifest, temporary_manifest)
  if (!file.rename(temporary_manifest, output_manifest_path)) {
    stop("Could not replace evaluation manifest.")
  }

  rm(
    predictions,
    partition_original_predictions,
    evaluation,
    invalid_groups
  )
  gc()
}

print(output_manifest)
