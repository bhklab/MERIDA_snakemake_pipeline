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
#'   sensitivity signature. Must match the colmn names of the test matrix.
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

    if (!length(sens_sig) || !length(resist_sig)) {
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
    resist_mat <- test_mat[, resist_sig, drop=FALSE]

    # apply MERIDA algorithm and return results
    predict_mat <- matrix(
        as.integer(rowAnys(sens_mat) & !rowAnys(resist_mat)),
        dimnames=list(rownames(sens_mat), "predicted_sensitive")
    )

    return(predict_mat)
}


# equivalent to if __name__ == "__main__": in Python
if (sys.nframe() == 0) {
    # read in the test data
    test_data <- read.table("procdata/test_matrix.txt", sep=" ")
    test_mat <- as.matrix(test_data[, !grepl("w|r", colnames(test_data))])
    test_labels <- matrix(test_data$r, ncol=1,
        dimnames=list(rownames(test_data), "response")
    )

    # read in all models
    models <- fread("results/merida_models.csv")
    setnames(models, "sensitivty", "sensitivity") # fix type in column name
    # removing the gene version numbers to make feature names match
    models[,
        sensitivity := strsplit(gsub("\\.\\d+", "", sensitivity), split="\\|")
    ]
    models[,
        resistance := strsplit(gsub("\\.\\d+", "", resistance), split="\\|")
    ]
    # make prediciton for each model
    models[,
        predictions := Map(f=predictSensitivityMERIDA,
            list(test_mat), sensitivity, resistance)
    ]
    # assess the classification performance with caret::confusionMatrix
    models[, predictions_factor := lapply(predictions, factor, levels=c(1, 0))]
    models[,
        performance := Map(f=confusionMatrix, predictions_factor,
            list(factor(test_labels, levels=c(1, 0))))
    ]
    models[,
        mcc := vapply(predictions_factor, FUN=mltools::mcc, 
            actuals=factor(test_labels, levels=c(1, 0)), numeric(1))
    ]
    class_eval1 <- models[,
        rbindlist(lapply(performance, \(x) as.list(x$byClass)))
    ]
    class_eval2 <- models[,
        rbindlist(lapply(performance, \(x) as.list(x$overall)))
    ]
    class_eval <- cbind(models[, .(M, v, mcc)], class_eval1, class_eval2)
    setorder(class_eval, -`Balanced Accuracy`)

    fwrite(class_eval, file="results/MERIDA_model_eval.csv")
}