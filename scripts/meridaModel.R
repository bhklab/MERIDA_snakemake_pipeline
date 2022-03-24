library(data.table)
library(matrixStats)
library(checkmate)
library(caret)
library(mltools)


#' Predict whether a sample is sensitive to a compound based on a MERIDA logical
#'   model
#'
#' @param test_mat `matrix` of `integer`s, where 1 represents the presence of
#'   of that feature in that sample.
#' @param sensitivity_signature `character` Vector of feature names for the
#'   sensitivity signature. Must match the column names of the test matrix.
#' @param resistance_signature `character` Vector fo feature names for
#'   the resistance signature. Must match the column names of the test matrix.
#'
#' @value `matrix` of `integer`s, where 1 represents samples predicted to be
#'   sensitive based on the provide `senstivity_signature` and
#'   `resistance_signature` feature vectors.
#'
#' @importFrom checkmate assertMatrix
#' @importFrom matrixStats rowAnys
#' @export
predictSensitivityMERIDA <- function(test_mat, sensitivity_signature,
        resistance_signature) {

    # input validation
    assertMatrix(test_mat, mode="integer")

    # parse features, warn if signature features are missing
    features <- colnames(test_mat)
    sens_sig <- intersect(sensitivity_signature, features)
    resist_sig <- intersect(resistance_signature, features)

    if (!length(sens_sig) && !length(resist_sig)) {
        warning("No features are present in the data for the sensitivity or ",
            "resistance signature. Please ensure the signature features match",
            "colnames(test_mat)!")
        return(NA_integer_)
    }
    if (length(sens_sig) != length(sensitivity_signature)) {
        warning("Missing sensitivity features: ",
            paste0(setdiff(sensitivity_signature, features), collapse=", "),
            call.=FALSE
        )
    }
    if (length(resist_sig) != length(resistance_signature)) {
        warning("Missing resistance features: ",
            paste0(setdiff(resistance_signature, features), collapse=", "),
            call.=FALSE
        )
    }

    # subset matrices of sensitivty and resistance features
    sens_mat <- test_mat[, sens_sig, drop=FALSE]
    # deal with case when no resistance features are selected
    resist_mat <- if (!length(resist_sig)) matrix(FALSE) else
        test_mat[, resist_sig, drop=FALSE]

    # apply MERIDA algorithm and return results
    predict_mat <- matrix(
        as.integer(rowAnys(sens_mat) & !rowAnys(resist_mat)),
        dimnames=list(rownames(sens_mat), "predicted_sensitive")
    )

    return(predict_mat)
}

#' @param models `data.table` Table with the list columns `sensitivity` and
#'   `resistance` containing `character` vectors of feature names predictive
#'   of sensitivity or resistance to the drug of interest.
#' @param mat `matrix` Of data on which to make MERIDA predictions.
#' @param labels `integer` Vector of true labels for `mat`.
#'
#' @value `data.table` The `models` table with `performance` and `mcc` columns,
#'   containing the results of `caret::confusionMatrix` and `mltools::mcc` for
#'   each row of `models`.
#'
#' @importFrom mltools mcc
#' @importFrom caret confusionMatrix
#' @importFrom data.table data.table copy
#' @export
evaluateMERIDAmodels <- function(models, mat, labels) {

    models <- copy(models)

    # make prediciton for each model
    models[,
        predictions := Map(f=predictSensitivityMERIDA,
            list(mat), sensitivity, resistance)
    ]

    # assess the classification performance with caret::confusionMatrix
    models[, predictions_factor := lapply(predictions, factor, levels=c(1, 0))]
    models[,
        performance := Map(f=caret::confusionMatrix, predictions_factor,
            list(factor(labels, levels=c(1, 0))))
    ]
    models[,
        mcc := vapply(predictions_factor, FUN=mltools::mcc,
            actuals=factor(labels, levels=c(1, 0)), numeric(1))
    ]

    return(models)
}

#' @param models `data.table` Table with columns `M`, `v`, `fold` specifying the
#'   model parameters and `performance`, `mcc` containing model evaluations,
#'   as returned by `evaluateMERIDAmodels`.
#'
#' @value `data.table` With model parameters and unpacked evaluation statistics,
#'   convenient for comparing the performance of various MERIDA model runs for
#'   hyperparameter tuning.
#'
#' @importFrom data.table data.table copy rbindlist setorderv
#' @export
extractMERIDAevaluation <- function(models) {
    class_eval1 <- models[,
        rbindlist(lapply(performance, \(x) as.list(x$byClass)))
    ]
    class_eval2 <- models[,
        rbindlist(lapply(performance, \(x) as.list(x$overall)))
    ]
    keep_cols <- intersect(colnames(models), c("source_file", "M", "v", "fold",
        "mcc", "objective_value", "sensitivity", "resistance"))
    class_eval <- cbind(models[, .SD, .SDcols=keep_cols], class_eval1,
        class_eval2)
    setorderv(class_eval, cols="Balanced Accuracy", order=-1L)

    return(class_eval)
}



# equivalent to if __name__ == "__main__": in Python
if (sys.nframe() == 0) {
    # -- Testing data

    training_set <- "CCLE"
    testset_list <- setdiff(
        unique(gsub("procdata/|_.*$", "", list.files("procdata", "*_Input_Matrix.txt", full.names=TRUE))),
        training_set
    )
    for (testing_set in testset_list) {
        testing_drugs <- unique(
            gsub(
                paste0("procdata/|", testing_set, "_|_Input_Matrix.txt"),
                "",
                list.files("procdata", paste0(testing_set,
                    ".*_Input_Matrix.txt"), full.names=TRUE)
        ))
        for (drug in testing_drugs) {
            print(testing_set)
            print(drug)
            # read in the test data
            test_data <- read.table(paste0("procdata/", testing_set,
                "_", drug, "_Input_Matrix.txt"), sep=" ")
            test_mat <- as.matrix(test_data[, !grepl("w|r", colnames(test_data))])
            test_labels <- matrix(test_data$r, ncol=1,
                dimnames=list(rownames(test_data), "response")
            )

            # read in all models
            models <- fread("results/merida_models.csv")[!(sensitivity %ilike% "resist"), ]  # drop models which predicted no features
            test_models <- copy(models)[source_file %ilike% drug, ]
            if (nrow(test_models) == 0) next()  # skip drugs which don't have any features predicted
            # removing the gene version numbers to make feature names match
            test_models[,
                sensitivity := strsplit(gsub("\\.\\d+", "", sensitivity), split="\\|")
            ]
            test_models[,
                resistance := strsplit(gsub("\\.\\d+", "", resistance), split="\\|")
            ]
            test_models <- evaluateMERIDAmodels(test_models, test_mat, test_labels)
            test_eval <- extractMERIDAevaluation(test_models)[order(fold), ]

            fwrite(test_eval, file=paste0("results/MERIDA_", testing_set, "_", drug, "_model_eval_test.csv"))
        }
    }
    training_drugs <- unique(
        gsub(
            paste0("procdata/|", training_set, "_|_Input_Matrix.txt"),
            "",
            list.files("procdata", paste0(training_set, ".*_Input_Matrix.txt"), full.names=TRUE)
    ))
    # -- Training data
    for (drug in training_drugs) {
        train_data <- read.table(paste0("procdata/", training_set, "_", drug, "_Input_Matrix.txt"))
        train_mat <- as.matrix(train_data[, !grepl("w|r", colnames(train_data))])
        train_labels <- matrix(train_data$r, ncol=1,
            dimnames=list(rownames(train_data), "response")
        )

        train_models <- fread("results/merida_models_by_fold.csv")[!(sensitivity %ilike% "resist"), ] # drop models which predicted no features
        models1 <- copy(models)
        train_models <- rbind(train_models,
            models1[, `:=`(fold="all")]
        )
        train_models <- train_models[source_file %ilike% drug, ]
        if (nrow(train_models) == 0) next()  # skip drug with no predicted features

        train_models[, sensitivity := strsplit(unlist(sensitivity), split="\\|")]
        train_models[, resistance := strsplit(unlist(resistance), split="\\|")]

        train_models <- evaluateMERIDAmodels(train_models, train_mat, train_labels)
        train_eval <- extractMERIDAevaluation(train_models)[order(fold), ]

        fwrite(train_eval, file=paste0("results/MERIDA_", training_set, "_", drug, "_model_eval_train.csv"))
    }


    # train_eval_summary <- train_eval[
    #     fold != "all",
    #     lapply(.SD, mean),
    #     .SDcols=!c("M", "fold"),
    #     by=.(M, v)
    # ]
    # train_eval_summary$fold <- "mean_of_folds"
    # train_eval_summary <- rbind(train_eval[fold == "all", ],
    #     train_eval_summary, fill=TRUE)
    # fwrite(train_eval_summary,
    #     file="results/MERIDA_model_eval_train_folds_vs_all.csv")

}