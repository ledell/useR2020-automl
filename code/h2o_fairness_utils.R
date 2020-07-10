# Code by Tomas Fryda, H2O AutoML Team (H2O.ai): https://github.com/tomasfryda


header <- function(string, level = 1) {
  cat("\n\n", strrep("#", level), " ", string, "\n", sep = "")
  flush.console()
}



# ------------------------------ Confusion Matrix ------------------------------
.calculate_confusion_matrix <-
  function(predicted, actual, favorable_class = "1") {
    t <- predicted == actual
    f <- predicted != actual

    p <- predicted == favorable_class
    n <- predicted != favorable_class

    tp <- sum(t & p)
    fp <- sum(f & p)
    tn <- sum(t & n)
    fn <- sum(f & n)

    total <- tp + fp + tn + fn

    return(
      list(
        TP = tp,
        FP = fp,
        TN = tn,
        FN = fn,
        Accuracy = (tp + tn) / total,
        Precision = tp / (tp + fp),
        Sensitivity = tp / (tp + fn),
        Specificity = tn / (fp + tn),
        F1 = (2 * tp) / (2 * tp + fp + fn),
        FalsePositiveRate = fp / (fp + tn),
        FalseNegativeRate = fn / (fn + tp),
        Selected = tp + fp,
        SelectedRatio = (tp + fp) / total,
        Total = total
      )
    )
  }

.render_confusion_matrix <- function(cm) {
  fmt <-
    function(x)
      sprintf("%d (%.2f %%)", as.integer(x), 100 * x / cm$Total)

  print(
    data.frame(
      "Predicted Positive" = c(fmt(cm$TP), fmt(cm$FP), fmt(cm$TP + cm$FP)),
      "Predicted Negative" = c(fmt(cm$FN), fmt(cm$TN), fmt(cm$FN + cm$TN)),
      "Total" = c(fmt(cm$TP + cm$FN), fmt(cm$FP + cm$TN), fmt(cm$Total)),
      row.names = c("Actual Positive", "Actual Negative", "Total"),
      check.names = FALSE
    )
  )
  cat("\n")

  printf("Accuracy: %.4f\n", cm$Accuracy)
  printf("Precision: %.4f\n", cm$Precision)
  printf("Sensitivity/Recall: %.4f\n", cm$Sensitivity)
  printf("Specificity: %.4f\n", cm$Specificity)
  printf("F1: %.4f\n", cm$F1)
}

make_confusion_matrix <- function(predicted, actual, favorable_class = "1") {
    cm <- .calculate_confusion_matrix(predicted, actual, favorable_class)
    .render_confusion_matrix(cm)
  }

# ------------------------------------------------------------------------------

# For a given model, calculate AIR, ME related cols & binary "adverse_impact" value for each subgroup
calculate_disparate_measures <- function(model,
                                         newdata,
                                         columns,
                                         reference_groups = list(),
                                         criterion = "SelectedRatio",
                                         favorable_class = "1",
                                         thresholds = list(air = c(0.8, 1.25), p_value = 0.05)) {
    if (is(model, "H2OAutoML")) {
      y <- model@leader@allparameters$y
    } else {
      y <- model@allparameters$y
    }

    predictions <- as.data.frame(predict(model, newdata)$predict)
    newdata_df <- as.data.frame(newdata)

    results <- expand.grid(lapply(columns, function(col)
      unlist(c(
        NA, levels(newdata_df[[col]])
      ))),
      stringsAsFactors = FALSE)
    names(results) <- columns

    results <- tail(results, n = -1)
    row.names(results) <- seq_len(nrow(results))

    results$reference <- NA

    for (idx in seq_len(nrow(results))) {
      mask <- rep_len(1, nrow(predictions))
      name <- NULL
      for (col in seq_len(ncol(results[idx, ]))) {
        if (!is.na(results[idx, col])) {
          col <- names(results)[col]
          val <- results[idx, col][[1]]
          name <- c(name, paste0(col, "=", val))
          mask <- mask & (newdata_df[[col]] == val)
        }
      }
      cm <-
        .calculate_confusion_matrix(predictions[mask, "predict"], newdata_df[mask, y],
                                    favorable_class = favorable_class)

      results[idx, "TP"] <- cm$TP
      results[idx, "FP"] <- cm$FP
      results[idx, "TN"] <- cm$TN
      results[idx, "FN"] <- cm$FN

      results[idx, "Accuracy"] <- cm$Accuracy
      results[idx, "Precision"] <- cm$Precision
      results[idx, "Sensitivity"] <- cm$Sensitivity
      results[idx, "Specificity"] <- cm$Specificity
      results[idx, "F1"] <- cm$F1
      results[idx, "Total"] <- cm$Total
      results[idx, "Selected"] <- cm$Selected
      results[idx, "SelectedRatio"] <- cm$SelectedRatio
    }

    groups <-
      combn(c(rep_len(NA, length(columns) - 1), columns), length(columns))

    # Get reference
    for (gid in seq_len(ncol(groups))) {
      cols <- columns[columns %in% groups[, gid]]
      fix_to_na <- (columns[!columns %in% groups[, gid]])

      group_mask <- Reduce("&", Map(function(col) {
        if (col %in% fix_to_na) {
          is.na(results[col])
        } else {
          !is.na(results[col])
        }
      }, columns))

      references_mask <-
        group_mask & Reduce("&", Map(function(col) {
          if (is.null(reference_groups[[col]])) {
            TRUE
          } else {
            results[col] == reference_groups[[col]]
          }
        }, cols))

      refs <- results[references_mask,]

      results[group_mask, "reference"] <-
        row.names(refs[which.max(refs[[criterion]]),])
    }

    # Calculate adverse impact & marginal effect
    for (row_id in row.names(results)) {
      ref <- results[row_id, "reference"]
      results[row_id, "air"] <-
        results[row_id, criterion] / results[ref, criterion]
      results[row_id, "me"] <-
        results[ref, criterion] - results[row_id, criterion]

      m <-
        matrix(
          c(results[row_id, "Total"], results[row_id, "Selected"],
            results[ref, "Total"], results[ref, "Selected"]),
          nrow = 2,
          dimnames = list(c("Total", "Selected"), c("Protected", "Reference"))
        )
      ft <- fisher.test(m)
      results[row_id, "p_value"] <- ft$p.value
    }

    results[["adverse_impact"]] <- FALSE

    if (!is.null(thresholds$air))
      results[["adverse_impact"]] <-
      results[["adverse_impact"]] |
      results[["air"]] < thresholds$air[[1]] |
      results[["air"]] > thresholds$air[[2]]
    if (!is.null(thresholds$me))
      results[["adverse_impact"]] <-
      results[["adverse_impact"]] |
      results[["me"]] < thresholds$me[[1]] |
      results[["me"]] > thresholds$me[[2]]
    if (!is.null(thresholds$p_value))
      results[["adverse_impact"]] <-
      results[["adverse_impact"]] &
      results[["p_value"]] < thresholds$p_value[[1]]

    return(results)
  }

# Add AIR, ME related cols and a binary "adverse_impact" column to an AutoML leaderboard
h2o.get_disparate_analysis <- function(aml,
                                       newdata,
                                       columns,
                                       reference_groups = list(),
                                       criterion = "SelectedRatio",
                                       favorable_class = "1",
                                       thresholds = list(air = c(0.8, 1.25), p_value = 0.05)) {
    da <- cbind(as.data.frame(aml@leaderboard), t(
      sapply(aml@leaderboard$model_id, function(model_id) {
        capture.output({
          dm <-
            calculate_disparate_measures(
              h2o.getModel(model_id),
              newdata,
              columns,
              favorable_class = favorable_class,
              criterion = criterion,
              reference_groups = reference_groups,
              thresholds = thresholds
            )
        })
        c(
          air_min = min(dm$air, na.rm = TRUE),
          air_mean = mean(dm$air, na.rm = TRUE),
          air_median = median(dm$air, na.rm = TRUE),
          air_max = max(dm$air, na.rm = TRUE),
          p_value_min = min(dm[which.min(dm$air), "p_value"], na.rm = TRUE),
          me_min = min(dm$me, na.rm = TRUE),
          me_mean = mean(dm$me, na.rm = TRUE),
          me_median = median(dm$me, na.rm = TRUE),
          me_max = max(dm$me, na.rm = TRUE),
          adverse_impact = any(dm$adverse_impact)
        )
      })
    ))
    da$adverse_impact <- as.logical(da$adverse_impact)
    return(da)
  }
