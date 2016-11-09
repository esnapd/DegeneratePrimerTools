#' Plot the MSA Object
#'
#' @importFrom purrr reduce
#' @export
plot_msa <- function(msa) {
  if (!class(msa) %in% c("DNAStringSet", "DNAMultipleAlignment"))  {
    stop("This function requires a DNAStringSet or DNAMultipleAlignment") 
  }
  
  if (class(msa) == "DNAStringSet") dna <- msa
  if (class(msa) == "DNAMultipleAlignment") dna <- DNAStringSet(msa)
  
  msanames <- names(dna)
  
  dfs <- lapply(seq(dna), function(x){
    dna1 <- dna[[x]]
    name <- names(dna)[[x]]
    df <- data.frame(seq(dna1), strsplit(as.character(dna1), ""))
    names(df) <- c("position", "base")
    df$seqname <- name
    df})
  
  df <- purrr::reduce(dfs, rbind)
  
  # order the sequence to be the same as the MSA
  df$seqname <- factor(df$seqname, levels=msanames)
  
  # filter out sequences
  df <- df %>% filter(base!="-")
  
  ggplot(df, aes(x=position, y=seqname, fill=base)) + 
    geom_tile()  + 
    theme_minimal()
  
}