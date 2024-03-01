# Load required libraries for network analysis, data manipulation, visualization, and statistical modeling.
library(tidyverse) # For data manipulation and visualization
library(tidygraph) # For working with graph objects in a tidy manner
library(ghypernet) # For network analysis
library(magrittr) # For advanced piping operations
library(caret)

# Check if the 'adjHelpR' package is installed, as it contains helper functions for adjacency matrix operations.
# If not installed, attempt to install it from GitHub using the 'remotes' package.
if (!require(adjHelpR, quietly = TRUE)) {
  if (!require(remotes, quietly = TRUE)) {
    install.packages("remotes") # Install 'remotes' if not already installed
  }
  remotes::install_github("gi0na/adjHelpR") # Install 'adjHelpR' from GitHub
}
library(adjHelpR) # Load 'adjHelpR'

# Source custom R scripts containing auxiliary functions not available in standard R packages.
source("phi_element_undirected.R") # Contains functions for calculating phi coefficients
source("adj2el_2.R") # Converts adjacency matrices to edge lists



# Load data from ergmcount
# The 'zach' dataset from the 'ergm.count' package is used as an example network for analysis.
# This dataset is loaded and then converted into a tbl_graph object, which is a tidygraph format 
# suitable for network analysis in the tidyverse paradigm.
data('zach', package = "ergm.count")

##################
# IMPORTANT NOTE: 
##################
# The Karate Club dataset captures the affiliations between
# members and two leaders, classified into positive (weak or strong), negative,
# and neutral relations. This dataset comprises 29 positive relations indicating
# varying degrees of association with the leaders, one explicit negative
# relation between the leaders themselves, and six neutral relations where
# individuals express no preference towards either leader.
#
# Analytical Approach: Training Phase: Multiclass Classification For the
# training phase, we employ a multiclass classification approach focused on the
# 36 explicitly declared relationships. This method allows for a nuanced
# understanding of the social dynamics within the dataset by differentiating
# between the varying intensities of positive associations and recognizing the
# presence of neutral or negative stances.
#
# Evaluation Phase: Binary Classification In contrast, the evaluation phase
# adopts a binary classification framework, categorizing relationships simply as
# "Friends" (positive relations) or "Non-Friends" (negative and neutral
# relations). This adjustment ensures consistency with the analysis of other
# datasets, facilitating a standardized assessment of our model's capability to
# discern between affiliative and non-affiliative ties within the network.

# Convert the 'zach' dataset into a tbl_graph object for further analysis.
# The as_tbl_graph function is used to transform the dataset into a graph format, storing it in 'full_graph'.
zach %>%
  as_tbl_graph() ->
  full_graph

# Extract the adjacency matrix from the 'full_graph' object.
# The get.adjacency function from the igraph package is used to obtain the adjacency matrix, 
# specifying 'contexts' as the attribute of interest. This matrix is stored in 'adj' for later use.
# The 'sparse = FALSE' option ensures the matrix is in a dense format, which is easier to work with in subsequent steps.
adj <- igraph::get.adjacency(full_graph, attr = 'contexts', sparse = FALSE)


# Calculate the potential edges (xi) by as the outer product of degrees,
# then adding the transpose to itself. This step doubles the potential edges (undirected network) but avoids 
# double-counting by halving the diagonal values.
xi <- apply(adj, 1, sum) %*% t(apply(adj, 2, sum))
xi <- xi + t(xi)
diag(xi) <- diag(xi) / 2


# Convert the model's omega values and the calculated xi to edge lists and merge them
# with the original adjacency matrix's edge list for comprehensive edge attribute analysis.
el <- adj2el_2(xi, directed = FALSE, names = colnames(adj)) %>% rename(xi = attr) %>% 
  left_join(adj2el_2(adj, directed = FALSE, names = colnames(adj)) %>% rename(count = attr), by = c("source", "target"))

# Calculate the total number of edges and the squared sum of potential edges for normalization.
m <- sum(adj) / 2
sum_xi <- 4 * m^2


# Calculate significance levels for edge under the hypergeometric configuration model (hyperg)
el <- el %>% 
  mutate(
    lower_than = Vectorize(phi_element_hyperg_under, vectorize.args = c('edgecount', 'xi'))(edgecount = count, xi = xi, sum_xi = sum_xi, m = m),
    higher_than = Vectorize(phi_element_hyperg_over, vectorize.args = c('edgecount', 'xi'))(edgecount = count, xi = xi, sum_xi = sum_xi, m = m),
  ) -> el

# Merge the calculated edge attributes back into the full graph and prepare a regression table
# for analyzing the impact of these attributes on network structure, focusing on nodes of interest.
full_graph %>% activate(edges) %>% 
  right_join(el %>% mutate(from = match(source, full_graph %>% as_tibble() %$% name),
                           to = match(target, full_graph %>% as_tibble() %$% name)) %>%
               select(from, to, lower_than, higher_than)) -> dat

# Filter edges for specific analysis, categorize closeness based on node connections, and prepare
# data for evaluation focusing on the declared Friend and Not-Friend relations (see note at the top).
dat %>%
  activate(edges)  %>% 
  filter(from == 1 | to == 34, !(from == 1 & to == 34)) %>% 
  mutate(faction.id = .N()$faction.id[if_else(to == 34, from, to)]) %>% 
  mutate(faction = .N()$faction[if_else(to == 34, from, to)]) %>% 
  mutate(closeness = case_when(from == 1 & faction.id < 0 ~ -faction.id, 
                               from == 1 & faction.id > 0 ~ -faction.id, 
                               to == 34 & faction.id > 0 ~ +faction.id, 
                               to == 34 & faction.id < 0 ~ +faction.id,
                               TRUE ~ 0)) %>% 
  mutate(closeness = if_else(from == 1 & to == 34, -2, closeness)) %>% 
  mutate(contexts = if_else(is.na(contexts), 0, contexts)) %>% 
  select(closeness, contexts, higher_than, lower_than, faction) %>% 
  as_tibble() -> regression_table
regression_table %>% mutate(friend = factor(if_else(closeness > 0, TRUE, FALSE))) -> regression_table


# The dataset 'dat' is filtered to focus on edges connected to nodes 1 or 34,
# which are of particular interest for the analysis. This step narrows down the
# dataset to the declared relations that will be used in the multi-class
# training of the Phi. (see note at the top)
dat %>%
  activate(edges) %>% 
  filter(from == 1 | to == 34) %>%
  # Determine the faction ID for each edge, based on the relationship to nodes 1 or 34.
  # This step categorizes each edge by the faction of the connected node, facilitating 
  # analysis of faction-based relationships within the network.
  mutate(faction.id = .N()$faction.id[if_else(to == 34, from, to)]) %>%
  mutate(faction = .N()$faction[if_else(to == 34, from, to)]) %>%
  # Calculate closeness indicators based on the faction IDs and specific node relations.
  # This categorization helps in understanding the proximity of relationships within the network.
  mutate(closeness = case_when(
    from == 1 & faction.id < 0 ~ -faction.id,
    from == 1 & faction.id > 0 ~ -faction.id,
    to == 34 & faction.id > 0 ~ +faction.id,
    to == 34 & faction.id < 0 ~ +faction.id,
    TRUE ~ 0
  )) %>%
  # Special handling for the direct relation between nodes 1 and 34.
  mutate(closeness = if_else(from == 1 & to == 34, -2, closeness)) %>%
  # Ensure all 'contexts' values are non-NA, setting missing values to 0.
  mutate(contexts = if_else(is.na(contexts), 0, contexts)) %>%
  # Select relevant attributes for regression analysis.
  select(closeness, contexts, higher_than, lower_than, faction) %>%
  as_tibble() %>%
  # Filter to include only edges with positive closeness,
  # or the special case of the direct connection between nodes 1 and 34.
  filter(closeness >= 0 | (from == 1 & to == 34)) ->
  regression_table_2

# The regression table is further processed to categorize 'closeness' into an ordered factor,
# which will serve as the dependent variable in the ordinal regression analysis. This transformation
# facilitates the modeling of ordered relationships based on closeness indicators.
regression_table_2 %>% 
  mutate(ordered_closeness = factor(closeness, levels = c(-2, 0, 1, 2), ordered = TRUE)) ->
  regression_table_2

# Faction identifiers are converted into factors to be used as independent variables or covariates
# in the regression analysis. This conversion allows for the examination of how faction affiliations
# may influence the relationships within the network.
regression_table_2 %>% 
  mutate(faction_factor = as.factor(faction)) ->
  regression_table_2


################################################################################
# multi-class prediction
################################################################################
# Source auxiliary functions for evaluation and cross-validation
source('karate_auxiliary_functions.R')

# Fit an ordinal logistic regression model (using polr from the MASS package) to predict the ordered closeness
# between nodes based on network structure attributes such as higher_than and lower_than metrics.
# This model aims to understand how these attributes influence the closeness of relationships within the network.
ordreg.over_under <- MASS::polr(ordered_closeness ~ higher_than + lower_than, data = regression_table_2)


# Execute the evaluation function on the regression table, categorizing edges based on their closeness
# to predict friendship links within the network.
eva_under_over <- evaluate_insample(datt = regression_table, model = ordreg.over_under, model_name = "over_under", positive_level = "TRUE", negative_level = "FALSE")
eva_under_over$res

spilt <- .7
seed <- 2

# Execute LOOCV for the three-class prediction task, using the specified formula, data split, and seed.
CV_under_over <- make_LOOCV_threeClass(dat_small = regression_table_2, dat_big = regression_table, formula = "ordered_closeness ~ higher_than + lower_than", split = .7, seed = seed, model_name = "over_under", positive_level = "TRUE", negative_level = "FALSE")
CV_under_over$res # given the small data size, this result changes heavily with seeding

