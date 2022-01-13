#!/bin/R

library(data.table)

input_dir <- list.files("results", pattern=".*M24.*cubic", full.names=TRUE)
file_path <- list.files(input_dir, pattern="^Result", full.names=TRUE)

#' Read in a MERIDA result file and parse it into a list of results
#'
#' @param file_path `character(1)` Absolute or relative path to the input file
#'
#' @value `list(4)` Results including parameters, best features, feauture result
#'  data.table and sample result data.table
#'
#' @importFrom data.table fread merge.data.table
#' @export
parse_merida_result_file <- function(file_path) {
    # Read in file line by line as character string
    file_txt <- readLines(file_path)

    # Define boundaries of header
    header_end <- 7
    feature_end <- grep(pattern="s_0", file_txt) - 1

    # Extract parameters and best features from the header
    merida_params <- data.table(as.numeric(gsub(".*: ", "", file_txt[3:4])))
    names(merida_params) <- c("M", "objective_value")

    best_features <- gsub("^[^:]*:[^:]*:", "", file_txt[6])
    best_feature_list <- unlist(strsplit(best_features, 
        "\tResistance features: "))
    names(best_feature_list) <- c("sensitivity", "resistance")
    for (i in seq_along(best_feature_list)) {
        best_feature_list[[i]] <- best_feature_list[[i]] |> 
            trimws() |>
            strsplit(split="\t") |>
            unlist() |>
            gsub(pattern='"', replacement='')
    }

    # Parse feature results
    feature_df <- fread(
        paste0(file_txt[seq(header_end, feature_end)], collapse="\n"),
        colClasses="character"
    )
    colnames(feature_df) <- c("feature", "value")

    # Parse sample results
    patterns <- c(sensitive="s_\\d+", resistant="r_\\d+", 
        prediction="sample_\\d+")
    sample_df_list <- vector("list", 3L)
    for (i in seq_along(patterns)) {
        sample_df_list[[i]] <- fread(
            paste0(file_txt[grepl(patterns[i], file_txt)], collapse="\n"),
            colClasses="character"
        )
        sample_df_list[[i]][, V1 := NULL]
        setnames(sample_df_list[[i]], c(V2="sample_name", V3=names(patterns)[i]))
    }
    merge_by <- function(x, y) merge.data.table(x, y, by="sample_name")
    sample_df <- Reduce(merge_by, sample_df_list)

    return(list(
        params=merida_params,
        best_features=best_feature_list,
        feature_df=feature_df,
        sample_df=sample_df
    ))
}