phi_element_hyperg_under <- function (edgecount, xi, sum_xi, m) {
  # Calculates the lower-tail probability under the hypergeometric distribution for undirected networks.
  # This function assesses the likelihood that the observed edge count between two nodes is significantly 
  # lower than expected by chance.
  #
  # Parameters:
  # - edgecount: Observed number of edges between two nodes.
  # - xi: Expected number of edges between the nodes based on the model.
  # - sum_xi: Total possible number of edges in the network, adjusted for the model.
  # - m: Total number of observed edges in the network.
  #
  # Returns:
  # A probability value indicating the likelihood of observing a lower edge count than expected.
  
  phyp <- phyper(q = edgecount, m = xi, n = sum_xi - xi, k = m, lower.tail = TRUE) 
  phyp <- phyp - dhyper(x = edgecount, m = xi, n = sum_xi - xi, k = m)
  return(phyp)
}



phi_element_hyperg_over <- function (edgecount, xi, sum_xi, m) {
  # Calculates the upper-tail probability under the hypergeometric distribution for undirected networks.
  # This function assesses the likelihood that the observed edge count between two nodes is significantly 
  # higher than expected by chance.
  #
  # Parameters are the same as `phi_element_hyperg_under`.
  #
  # Returns:
  # A probability value indicating the likelihood of observing a higher edge count than expected.
  
  phyp <- phyper(q = edgecount, m = xi, n = sum_xi - xi, k = m, lower.tail = FALSE) 
  return(phyp)
}


phi_element_CL_under <- function(edgecount, xi, m) {
  # Calculates the lower-tail probability under the Chung-Lu model for undirected networks, adjusted for the model's
  # assumptions about independent edge probabilities proportional to node degrees.
  #
  # Parameters:
  # - edgecount: Observed number of edges between two nodes.
  # - xi: Expected number of edges based on the Chung-Lu model.
  # - m: Total number of observed edges in the network.
  #
  # Returns:
  # A probability value indicating the likelihood of observing a lower edge count than expected under the Chung-Lu model.
  
  CL <- ppois(q = edgecount, lambda = xi / m / 4, lower.tail = TRUE, log.p = FALSE)
  CL - dpois(x = edgecount, lambda = xi / m / 4)
}


phi_element_CL_over <- function(edgecount, xi, m) {
  # Calculates the upper-tail probability under the Chung-Lu model for undirected networks.
  #
  # Parameters are similar to `phi_element_CL_under`.
  #
  # Returns:
  # A probability value indicating the likelihood of observing a higher edge count than expected under the Chung-Lu model.
  
  ppois(q = edgecount, lambda = xi / m / 4, lower.tail = FALSE, log.p = FALSE)
}




phi_element_bccm_under <- function(edgecount, m, prob) {
  # Calculates the lower-tail probability under the Block-Constrained Configuration Model (BCCM) for undirected networks.
  # This function assesses the likelihood that the observed edge count is significantly lower than expected,
  # based on the probability of edge formation in the BCCM.
  #
  # Parameters:
  # - edgecount: Observed number of edges between two nodes.
  # - m: Total number of observed edges in the network.
  # - prob: Probability of edge formation between two nodes in the BCCM.
  #
  # Returns:
  # A probability value indicating the likelihood of observing a lower edge count than expected under the BCCM.
  
  pbin <- pbinom(q = edgecount, size = m, prob = prob, lower.tail = TRUE, log.p = FALSE)
  pbin - dbinom(x = edgecount, size = m, prob = prob)
}

                   

phi_element_bccm_over <- function(edgecount, m, prob) {
  # Calculates the upper-tail probability under the Block-Constrained Configuration Model (BCCM) for undirected networks.
  #
  # Parameters are similar to `phi_element_bccm_under`.
  #
  # Returns:
  # A probability value indicating the likelihood of observing a higher edge count than expected under the BCCM.
  
  pbinom(q = edgecount, size = m, prob = prob, lower.tail = FALSE, log.p = FALSE)
}

