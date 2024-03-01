# adj2el_2: Convert Adjacency Matrix to Edge List
#
# Description:
# This function converts an adjacency matrix (adj) into a full edge list (reporting disconnected pairs as well), facilitating the analysis
# and visualization of graph data. The edge list format is particularly useful for graph algorithms
# and is widely supported by graph analysis libraries.
#
# Parameters:
# - adj: Adjacency matrix representing the graph. Can be a dense matrix or a sparse matrix (from the Matrix package).
# - directed (logical): Specifies whether the graph is directed. Defaults to TRUE. If set to FALSE, 
#   entries below the diagonal (or above, depending on the include_diag parameter) are ignored.
# - names (character vector): Optional vector of names for the rows and columns of the adjacency matrix. 
#   If provided, it must match the dimensions of the matrix.
# - include_diag (logical): Specifies whether to include diagonal elements in the edge list. Defaults to FALSE.
#   This parameter is relevant for networks where self-loops are meaningful.
#
# Returns:
# A tibble (data frame) representing the edge list of the graph. The tibble has three columns:
# - source: Source node of each edge. For undirected graphs, source is always less than target to avoid duplication.
# - target: Target node of each edge.
# - attr: The value associated with each edge, typically representing weight or capacity.
#
# Details:
# The function first checks if the 'names' parameter is provided and matches the dimensions of the adjacency matrix.
# If 'directed' is FALSE, it sets entries below the diagonal to NA to ignore them in the edge list (depending on 'include_diag').
# It then converts the adjacency matrix to a tibble, filtering out NA values to exclude unused edges.
# For matrices stored using the Matrix package, it specifically renames columns to match the expected output format.
#
# Note:
# If the adjacency matrix format is not recognized or the 'names' parameter does not match the matrix dimensions,
# the function will stop and return an error message.
#
# Example usage:
# adj_matrix <- matrix(c(0, 1, 1, 0), nrow = 2)
# el <- adj2el_2(adj = adj_matrix, directed = FALSE, names = c("A", "B"))

adj2el_2 <- function (adj, directed = TRUE, names = NULL, include_diag = FALSE) 
{
  if ((!is.null(names)) && length(names) != nrow(adj)) 
    stop("length(names) is different from nrow(adj)")
  if (!is.null(names)) 
    rownames(adj) <- colnames(adj) <- names
  if (isFALSE(directed)) 
    adj[lower.tri(adj, !include_diag)] <- NA
  if (is.matrix(adj)) {
    el <- reshape2::melt(adj, value.name = "value") %>% 
      as_tibble()  %>% filter(!is.na(.data$value)) %>% rename(source = "Var1", 
                                                          target = "Var2", attr = "value") %>% mutate_at(c("source", 
                                                                                                           "target"), as.character)
  }
  else {
    if (attributes(class(adj))$package == "Matrix") {
      el <- as_tibble(adj) %>% rename(source = "row", 
                                      target = "column", attr = "value")
    }
    else {
      stop("Wrong adjacency format")
    }
  }
  return(el)
}
