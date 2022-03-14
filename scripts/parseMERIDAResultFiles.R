#!/bin/R

# 0.0 -- Load dependenices
library(data.table)


# 0.1 -- Function definitions

# #' Read in a MERIDA result file and parse it into `data.table`s and write to
# #'   disk
# #'
# #' @param file_path `character(1)` Absolute or relative path to the input file
# #'
# #' @value None, writes to disk.
# #'
# #' @importFrom data.table fread merge.data.table
# #' @export
# parse_merida_result_file <- function(file_path) {
#     # Read in file line by line as character string
#     file_txt <- readLines(file_path)

#     # Define boundaries of header
#     header_end <- 7
#     feature_end <- grep(pattern="s_0", file_txt) - 1

#     # Extract parameters and best features from the header
#     merida_params <- gsub(".*: ", "", file_txt[3:4]) |>
#         as.numeric() |>
#         as.list() |>
#         as.data.table()
#     names(merida_params) <- c("M", "objective_value")
#     merida_params$v <- gsub("^.*_v|_cv.*$", "", file_path)

#     best_features <- gsub("^[^:]*:[^:]*:", "", file_txt[6])
#     best_feature_list <- strsplit(best_features, "\tResistance features: ") |>
#         unlist() |>
#         as.list()
#     names(best_feature_list) <- c("sensitivity", "resistance")
#     for (i in seq_along(best_feature_list)) {
#         best_feature_list[[i]] <- best_feature_list[[i]] |>
#             trimws() |>
#             strsplit(split="\t") |>
#             unlist() |>
#             gsub(pattern='"', replacement='')
#     }
#     merida_params$sensitivty <- paste0(best_feature_list$sensitivity,
#         collapse="|")
#     merida_params$resistance <- paste0(best_feature_list$resistance,
#         collapse="|")

#     # Parse feature results
#     feature_df <- fread(
#         paste0(file_txt[seq(header_end, feature_end)], collapse="\n"),
#         colClasses="character"
#     )
#     colnames(feature_df) <- c("feature", "value")
#     feature_df <- cbind(merida_params, feature_df)

#     # Parse sample results
#     patterns <- c(sensitive="s_\\d+", resistant="r_\\d+",
#         prediction="sample_\\d+")
#     sample_df_list <- vector("list", 3L)
#     for (i in seq_along(patterns)) {
#         sample_df_list[[i]] <- fread(
#             paste0(file_txt[grepl(patterns[i], file_txt)], collapse="\n"),
#             colClasses="character"
#         )
#         sample_df_list[[i]][, V1 := NULL]
#         setnames(sample_df_list[[i]], c(V2="sample_name",
#             V3=names(patterns)[i]))
#     }
#     merge_by <- function(x, y) merge.data.table(x, y, by="sample_name")
#     sample_df <- Reduce(f=merge_by, sample_df_list)
#     sample_df <- cbind(merida_params, sample_df)

#     # Save results to disk
#     results <- list(
#         feature_df=feature_df,
#         sample_df=sample_df
#     )

#     for (i in seq_along(results)) {
#         output_file <- gsub(
#             pattern=basename(file_path),
#             replacement=paste0(names(results)[i], ".csv"),
#             file_path)
#         fwrite(results[[i]], file=output_file)
#     }
# }

#' Read in a MERIDA result file and parse it into `data.table`s and write to
#'   disk
#'
#' @param file_path `character(1)` Absolute or relative path to the input file
#'
#' @value None, writes to disk.
#'
#' @importFrom data.table fread merge.data.table
#' @export
parse_merida_result_file <- function(file_path) {
    print(file_path)
    # Read in file line by line as character string
    file_txt <- readLines(file_path)

    # Define boundaries of header
    header_end <- 7
    feature_end <- grep(pattern="s_0", file_txt) - 1

    # Extract parameters and best features from the header
    merida_params <- gsub(".*: ", "", file_txt[3:4]) |>
        as.numeric() |>
        as.list() |>
        as.data.table()
    names(merida_params) <- c("M", "objective_value")
    merida_params$v <- gsub("^.*_v|_cv.*$", "", file_path)
    hasFold <- grepl("Fold\\d+", file_path)
    merida_params$fold <- if (hasFold) gsub("^.*Fold|Result.*$", "", file_path) else "all"

    best_features <- gsub("^[^:]*:[^:]*:", "", file_txt[6])
    best_feature_list <- strsplit(best_features, "\tResistance features: ") |>
        unlist() |>
        as.list()
    # Add empty string if there are no resistance features
    if (length(best_feature_list) != 2)
        best_feature_list <- append(best_feature_list, list(""))
    names(best_feature_list) <- c("sensitivity", "resistance")
    for (i in seq_along(best_feature_list)) {
        best_feature_list[[i]] <- best_feature_list[[i]] |>
            trimws() |>
            strsplit(split="\t") |>
            unlist() |>
            gsub(pattern='"', replacement='')
    }
    merida_params$sensitivity <- paste0(best_feature_list$sensitivity,
        collapse="|")
    merida_params$resistance <- paste0(best_feature_list$resistance,
        collapse="|")

    # Parse feature results
    feature_df <- fread(
        paste0(file_txt[seq(header_end, feature_end)], collapse="\n"),
        colClasses="character"
    )
    colnames(feature_df) <- c("feature", "value")
    feature_df <- cbind(merida_params, feature_df)

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
        setnames(sample_df_list[[i]], c(V2="sample_name",
            V3=names(patterns)[i]))
    }
    merge_by <- function(x, y) merge.data.table(x, y, by="sample_name")
    sample_df <- Reduce(f=merge_by, sample_df_list)
    sample_df <- cbind(merida_params, sample_df)

    results <- list(
        feature_df=feature_df,
        sample_df=sample_df
    )
    # Add fold identifier if the result is from a cross-validation run
    if (hasFold) {
        names(results) <- paste0("Fold", merida_params$fold, "_",
            names(results))
    }
    for (i in seq_along(results)) {
        output_file <- gsub(
            pattern=basename(file_path),
            replacement=paste0(names(results)[i], ".csv"),
            file_path)
        fwrite(results[[i]], file=output_file)
    }
}

# 0.2 -- Parse snakemake parameters
input_dir <- snakemake@input$results_dir


# 1.0 -- Parse results
result_paths <- list.files(input_dir, pattern="^Result_.*", recursive=TRUE,
    full.names=TRUE)
for (path in result_paths) parse_merida_result_file(path)

# Only read in folds for jobs which have competed
fold_dirs <- gsub("\\/[^\\/]*txt", "", result_paths)
fold_paths <- lapply(fold_dirs, FUN=list.files,
    pattern=".*Fold\\d+Result", full.names=TRUE)
for (paths in fold_paths) {
    for (path in paths) parse_merida_result_file(path)
}

# 2.0 -- Read in and collate results for each set of parameters
# Overall results
file_names <- c(features="^feature_df", samples="^sample_df")
file_paths <- lapply(file_names,
    FUN=\(x) list.files(input_dir, pattern=paste0(x, ".csv"),
        recursive=TRUE, full.names=TRUE)
)
for (i in seq_along(file_paths)) {
    df_ <- lapply(file_paths[[i]], FUN=fread) |>
        setNames(file_paths[[i]]) |>
        rbindlist(idcol="source_file")
    fwrite(df_, file=file.path(input_dir, paste0(names(file_paths)[i], ".csv")))
}

# Cross-validation results
fold_files <- unlist(lapply(fold_dirs, list.files,
    pattern="Fold\\d+_feature_df.csv", full.names=TRUE))
fold_df <- lapply(fold_files, FUN=fread) |>
    setNames(fold_files) |>
    rbindlist(idcol="source_file")
merida_models <- unique(fold_df[, .(source_file, M, v, fold, sensitivity, resistance)])
fwrite(merida_models, file=file.path(input_dir, "merida_models_by_fold.csv"))