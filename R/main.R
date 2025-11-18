#' Encode DNA sequences into 5-mer frequency vectors
#'
#' This function converts input DNA sequences into a numerical matrix representation
#' based on 5-mer frequency encoding for m6A site prediction.
#'
#' @param dna_strings A character vector of DNA sequences (e.g., "ATGCGT...")
#'
#' @return A data frame containing 5-mer frequency features.
#' @examples
#' dna_encoding(c("ATGCGTACGA", "CGTACGTAGC"))
#' @export
 dna_encoding <- function(dna_strings){
nn <- nchar(dna_strings[1])
seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol=nn, byrow=TRUE)
colnames(seq_m) <- paste0("nt_pos", 1:nn)
seq_df <- as.data.frame(seq_m)
seq_df[] <- lapply(seq_df, factor, levels=c("A", "T", "C", "G"))
return(seq_df)
}


#' Predict m6A sites for multiple sequences
#'
#' This function performs batch prediction of m6A sites using a pre-trained
#' random forest model included in the package.
#'
 #' @param ml_fit A trained random forest model object.
 #' @param feature_df A data frame containing sequence features for prediction.
 #' @param positive_threshold A numeric cutoff (default 0.5) for classification.
#'
#' @return A data frame containing prediction probabilities for m6A sites.
#'
#' @import randomForest
#' @importFrom stats predict
#' @examples
#' model <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#' input <- read.csv(system.file("extdata", "m6A_input_example.csv", package = "m6APrediction"))
#' prediction_multiple(model, input)
#' @export
prediction_multiple <- function(ml_fit, feature_df, positive_threshold=0.5){
  required_cols <- c("gc_content", "RNA_type", "RNA_region", "exon_length",
                     "distance_to_junction", "evolutionary_conservation", "DNA_5mer")
  if (!all(required_cols %in% colnames(feature_df))) {
    stop("Input data is missing required columns")
  }

  feature_df$RNA_type <- factor(feature_df$RNA_type,
                                levels=c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region,
                                  levels=c("CDS", "intron", "3'UTR", "5'UTR"))

  dna_encoded <- dna_encoding(feature_df$DNA_5mer)

  predict_data <- data.frame(
    gc_content = feature_df$gc_content,
    RNA_type = feature_df$RNA_type,
    RNA_region = feature_df$RNA_region,
    exon_length = feature_df$exon_length,
    distance_to_junction = feature_df$distance_to_junction,
    evolutionary_conservation = feature_df$evolutionary_conservation,
    dna_encoded
  )

  prob_matrix <- predict(ml_fit, newdata = predict_data, type = "prob")
  feature_df$predicted_m6A_prob <- prob_matrix[, "Positive"]

  feature_df$predicted_m6A_status <- ifelse(
    feature_df$predicted_m6A_prob > positive_threshold,
    "Positive",
    "Negative"
  )

  return(feature_df)
}

#' Predict m6A site for a single sequence
#'
#' This function predicts whether a single DNA sequence is likely to contain an m6A modification.
#'
#' @param ml_fit A trained random forest model object.
#' @param gc_content GC content of the sequence (0–1).
#' @param RNA_type RNA type (e.g., "mRNA", "lncRNA").
#' @param RNA_region RNA region (e.g., "CDS", "3'UTR").
#' @param exon_length Numeric length of the exon.
#' @param distance_to_junction Distance from m6A site to splice junction.
#' @param evolutionary_conservation Numeric conservation score.
#' @param DNA_5mer A string representing the 5-mer motif sequence.
#' @param positive_threshold Probability threshold for classification (default 0.5).
#' @return A numeric value between 0 and 1 representing predicted m6A probability.
#' @examples
#' model <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#' prediction_single(model, 0.45, "mRNA", "3'UTR", 10, 8, 0.6, "GGACA")
#' @export
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length,
                              distance_to_junction, evolutionary_conservation,
                              DNA_5mer, positive_threshold = 0.5) {
  feature_df <- data.frame(
    gc_content = gc_content,
    RNA_type = RNA_type,
    RNA_region = RNA_region,
    exon_length = exon_length,
    distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer
  )


  result <- prediction_multiple(ml_fit, feature_df, positive_threshold)

  return(list(
    predicted_m6A_prob = result$predicted_m6A_prob,
    predicted_m6A_status = as.character(result$predicted_m6A_status)
  ))
}
