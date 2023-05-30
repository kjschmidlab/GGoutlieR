#' calculate genetic similarity from ancestry coefficients.
#' @details
#' Since ancestry coefficients can be interpreted as a propotion of genome derived from a specific ancestral population, this function calculate genetic similarity as the probability of a random genome segment of two individuals derived from the same ancestral population.
#' @param anc_coef a matrix of ancestry coefficients with samples by rows. Each column corresponds to an ancestral population.

anc_coeff_to_GeneticSimilarityMatrix <- function(anc_coef){

  # a function to calculate the probability of drawing alleles from the same ancestral populations
  prob_same_anc_pop <- function(x1 ,x2){
    return(sum(x1*x2))
  }
  n = nrow(anc_coef)
  out <- matrix(0, n,n)

  to_fill <- which(lower.tri(out), arr.ind = T)

  message("calculating probability of shared ancestral alleles...\n")
  sprob <- apply(to_fill,1, function(z){
    prob_same_anc_pop(anc_coef[z[1],], anc_coef[z[2],])
  })

  for (i in 1:length(sprob)) {
    out[to_fill[i,1], to_fill[i,2]] <- sprob[i]
  }
  out <- out + t(out)
  return(out)
} # anc_coeff_to_GeneticSimilarityMatrix end
