
phi_element_hyperg_under <- function (edgecount, xi, sum_xi, m) {
  # Calculates the probability of observing a number of edges less than or equal to 
  # the observed count under the hypergeometric distribution, adjusted by subtracting 
  # the exact probability of the observed count. This adjustment accounts for the 
  # discrete nature of the distribution and provides a more accurate measure of 
  # significantly lower interactions than expected.
  #
  # Parameters:
  # edgecount: The observed number of edges between two nodes.
  # xi: The expected number of edges between the two nodes based on the model.
  # sum_xi: The total possible number of edges in the network, adjusted for the model.
  # m: The total number of edges observed in the network.
  #
  # Returns:
  # A probability value indicating how likely it is that the observed edge count 
  # is significantly lower than expected under the hypergeometric distribution.
  
  phyp <- phyper(q = edgecount, m = xi, n = sum_xi - xi, k = m, lower.tail = TRUE) 
  phyp <- phyp - dhyper(x = edgecount, m = xi, n = sum_xi - xi, k = m)
  return(phyp)
}



phi_element_hyperg_over <- function (edgecount, xi, sum_xi, m) {
  # Calculates the probability of observing a number of edges greater than the 
  # observed count under the hypergeometric distribution. This probability reflects 
  # how likely it is that the observed edge count is significantly higher than 
  # what would be expected based on the model.
  #
  # Parameters are the same as for phi_element_hyperg_under.
  #
  # Returns:
  # A probability value indicating the likelihood of the observed edge count 
  # being significantly higher than expected under the hypergeometric distribution.
  
  phyp <- phyper(q = edgecount, m = xi, n = sum_xi - xi, k = m, lower.tail = FALSE) 
  return(phyp)
}


phi_element_CL_under <- function(edgecount, xi, m) {
  # Calculates the probability of observing a number of edges less than or equal to 
  # the observed count under the Chung-Lu model, adjusted by subtracting the exact 
  # probability of the observed count. This model assumes edges between nodes are 
  # independent and the probability of an edge is proportional to the product of node 
  # degrees.
  #
  # Parameters:
  # edgecount, xi, and m are defined as in the previous functions but applied in the
  # context of the Chung-Lu model.
  #
  # Returns:
  # A probability value indicating how likely it is that the observed edge count 
  # is significantly lower than expected under the Chung-Lu model.
  
  CL <- ppois(q = edgecount, lambda = xi/m, lower.tail = TRUE, log.p = FALSE)
  CL - dpois(x = edgecount, lambda = xi/m)
}


phi_element_CL_over <- function(edgecount, xi, m) {
  # Calculates the probability of observing a number of edges greater than the 
  # observed count under the Chung-Lu model. This function is analogous to 
  # phi_element_hyperg_over but applies the Chung-Lu model's assumptions.
  #
  # Parameters and returns are similar to phi_element_CL_under but for significantly
  # higher interactions than expected.
  
  ppois(q = edgecount, lambda = xi/m, lower.tail = FALSE, log.p = FALSE)
}

