#' Calculate Rarefaction Curves
#'
#' Generate Rarefaction Curves from Phyloseq.
#' Adopted from @@and3k https://github.com/joey711/phyloseq/issues/143
#' 
#' @param physeq phyloseq object
#' @param measures vector of rarefaction measures
#' @param depths vector of depths to rarefy to
#' @param parallel wether to use ldplyr's parallel features
#' @param ncpus nubmer of cpus
#' 
#' @importFrom plyr ldply
#' @importFrom plyr summarise
#' @importFrom reshape2 melt
#' @import foreach
#' @import doParallel
#' 
#' @export
calculate_rarefaction_curves <- function(physeq, measures, depths, parallel=FALSE, ncpus=1) {
  estimate_rarified_richness <- function(physeq, measures, depth) {
    if(max(sample_sums(physeq)) < depth) return()
    physeq <- prune_samples(sample_sums(physeq) >= depth, physeq)
    
    rarified_physeq <- rarefy_even_depth(physeq, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_physeq, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), 
                                   varnames = c('Sample', 'Measure'), 
                                   value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  if (parallel){
    #if parallel setup the cluster
    library(doParallel)
    print("Running Calculation in Parallel...")
    cl <- makeCluster(ncpus)
    registerDoParallel(cl)
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- plyr::ldply(depths, estimate_rarified_richness, 
                            physeq = physeq, measures = measures, 
                            .id = 'Depth', 
                            .paropts = list(.packages = c('phyloseq', 'vegan','reshape2')),
                            .progress = ifelse(interactive(), 'text', 'none'),
                            .parallel = parallel)
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  #add a standard deviation column
  rarefaction_curve_data_summary <- plyr::ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), 
                                    plyr::summarise, 
                                    #Alpha_diversity_mean = mean(Alpha_diversity), 
                                    #Alpha_diversity_sd = sd(Alpha_diversity)) 
                                    Alpha_diversity_mean = mean(value), 
                                    Alpha_diversity_sd = sd(value))
  
  #add the sample data
  #rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, 
  #                                                data.frame(sample_data(physeq)), 
  #                                                by.x = 'Sample', 
  #                                                by.y = 'row.names')
  #return(rarefaction_curve_data_summary_verbose)
  return(rarefaction_curve_data_summary)
}