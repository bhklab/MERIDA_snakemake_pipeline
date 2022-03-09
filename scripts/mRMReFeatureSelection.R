library(mRMRe)

training_data_all <- read.table("procdata/CCLE_bimodal_genes.txt")
training_data <- training_data_all
training_data <- data.frame(lapply(training_data, factor, ordered=TRUE))
rownames(training_data) <- rownames(training_data_all)

## Removing weight column for feature selection
training_data[, "w"] <- NULL

input_mrmr <- mRMR.data(training_data)

mrmr_res <- mRMR.ensemble(input_mrmr, feature_count=100, solution_count=1,
    target_indices=ncol(training_data))

selected_features <- colnames(training_data)[mrmr_res@filters[[1]]]


filtered_input_data <- training_data_all[,
    c(selected_features, "w", "r")]

write.table(filtered_input_data,
    file=file.path("CCLE/Input_Matrix_filtered.txt"),
    sep=" "
)