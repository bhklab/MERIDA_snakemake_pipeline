library(mRMRe)

# -- read in Snakemake parameters
feature_matrix <- snakemake@input
nFeature <- snakemake@params$nFeature
nSolution <- snakemake@params$nSolution
output <- snakemake@output

# -- mRMRe feature selection
class(feature_matrix)
print(str(feature_matrix))
training_data_all <- read.table(feature_matrix[[1]])
training_data <- training_data_all
training_data <- data.frame(lapply(training_data, factor, ordered=TRUE))
rownames(training_data) <- rownames(training_data_all)

## Removing weight column for feature selection
training_data[, "w"] <- NULL

input_mrmr <- mRMR.data(training_data)

mrmr_res <- mRMR.ensemble(input_mrmr, feature_count=nFeature,
    solution_count=nSolution, target_indices=ncol(training_data))

selected_features <- colnames(training_data)[mrmr_res@filters[[1]]]

filtered_input_data <- training_data_all[,
    c(selected_features, "w", "r")]

print(str(output))

write.table(filtered_input_data,
    file=file.path(output[[1]]),
    sep=" "
)