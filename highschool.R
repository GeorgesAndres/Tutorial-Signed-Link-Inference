library(tidyverse)
library(tidygraph)
library(texreg)
library(ghypernet)
library(pROC)
library(caret)

# Check if adjHelpR is installed
if (!require(adjHelpR, quietly = TRUE)) {
  # If not installed, install remotes package if not already installed
  if (!require(remotes, quietly = TRUE)) {
    install.packages("remotes")
  }
  # Use remotes to install adjHelpR from GitHub
  remotes::install_github("gi0na/adjHelpR")
}
library(adjHelpR)

source("functions/adj2el_2.R")
source("functions/phi_element_undirected.R")

###### GET DATA

# Check and create directories if they don't exist
if (!dir.exists('data')) {
  dir.create('data', showWarnings = FALSE)
}
if (!dir.exists('data/highschool')) {
  dir.create('data/highschool', showWarnings = FALSE)
}

# Define file URLs and destination paths
files_info <- list(
  highschool_data = list(
    url = 'http://www.sociopatterns.org/wp-content/uploads/2015/07/High-School_data_2013.csv.gz',
    destfile = 'data/highschool/High-School_data_2013.csv'
  ),
  metadata = list(
    url = 'http://www.sociopatterns.org/wp-content/uploads/2015/09/metadata_2013.txt',
    destfile = 'data/highschool/metadata_2013.txt'
  ),
  facebook_known_pairs = list(
    url = 'http://www.sociopatterns.org/wp-content/uploads/2015/07/Facebook-known-pairs_data_2013.csv.gz',
    destfile = 'data/highschool/Facebook-known-pairs_data_2013.csv'
  ),
  friendship_network = list(
    url = 'http://www.sociopatterns.org/wp-content/uploads/2015/07/Friendship-network_data_2013.csv.gz',
    destfile = 'data/highschool/Friendship-network_data_2013.csv'
  )
)

# Download and decompress files if they do not exist
for (file_info in files_info) {
  if (!file.exists(file_info$destfile)) {
    # Download compressed file
    download.file(url = file_info$url, destfile = paste0(file_info$destfile, ".gz"))
    # Decompress
    R.utils::gunzip(paste0(file_info$destfile, ".gz"), overwrite = TRUE)
  }
}


#######################################################################################
# Load and Preparation of High School Social Interaction Data
#######################################################################################

# Define a threshold to consider subsequent interactions as a single interaction if they occur within a 20-time unit window.
my_threshold <- 20

# Load interaction data from CSV, specifying column types explicitly for efficiency and clarity.
# - 'X1' represents the time of interaction, which is renamed to 't' for clarity.
# - 'X2' and 'X3' represent the IDs of the interacting students, treated as character data for flexibility.
# The data is then processed to:
# 1. Calculate 'from' and 'to' by ordering the IDs to ensure consistency in directionality (min to max).
# 2. Group by these 'from' and 'to' pairs to aggregate interactions.
# 3. Calculate time differences between consecutive interactions to filter based on the defined threshold.
el_contacts <- read_table("data/highschool/High-School_data_2013.csv", 
                          col_names = FALSE, col_types = cols(X2 = col_character(), X3 = col_character())) %>%
  rename(t='X1') %>%
  mutate(from = as.character(apply(cbind(X2,X3),1, function(row) min(as.integer(row)))),
         to = as.character(apply(cbind(X2,X3),1, function(row) max(as.integer(row))))) %>%
  select(-X2,-X3) %>%
  group_by(from,to) %>%
  mutate(t_shifted = c(t[1],head(t,-1))) %>%
  mutate(dif = abs(t - t_shifted)) %>%
  filter(dif == 0 | dif > my_threshold) %>%
  select(t, from, to)

# Load metadata containing additional information about vertices (students).
# - 'X1' is assumed to represent student IDs.
metadata <- read_table("data/highschool/metadata_2013.txt", 
                       col_names = FALSE, col_types = cols(X1 = col_character()))

# Load friendship network data, indicating reported friendships between students.
el_fr <- read_table("data/highschool/Friendship-network_data_2013.csv", 
                    col_names = FALSE, col_types = cols(X1 = col_character(), X2 = col_character()))

# Extract unique nodes from the contact data to ensure consistency across analyses.
nodes <- el_contacts %>% select(from, to) %>% t %>% c %>% unique() %>% sort()

# Generate adjacency matrix from contact data, configuring it for undirected, non-sparse representation.
adj <- get_adjacency(x = el_contacts, select_cols = 2:3, multiedge = TRUE, sparse = FALSE, directed = FALSE, nodes = nodes)

# Reorder metadata to match the adjacency matrix's row order and remove any vertices not present.
# This ensures alignment between the network data and its corresponding metadata.
metadata <- metadata[-which(!(metadata$X1 %in% rownames(adj))),]
maporder <- function(id, v) which(id == v)
metadata <- metadata[Vectorize(FUN = maporder, vectorize.args = 'id')(id = rownames(adj), metadata$X1),]

##################################################################################################################
# Calculate phi (undirected)
# Description: This section calculates the phi coefficient for an undirected network, 
# integrating node attributes and network structure. It involves transforming the adjacency 
# matrix to calculate potential connections (xi), merging network statistics with the edge list, 
# and applying models to infer probabilities of connections.
##################################################################################################################

# Calculate potential connections (xi) by summing the adjacency matrix rows and columns,
# then adding the transpose to itself and halving the diagonal to avoid double-counting.
xi <- apply(adj, 1, sum) %*% t(apply(adj, 2, sum))
xi <- xi + t(xi)
diag(xi) <- diag(xi)/2

# Convert xi to an edge list and merge with the original adjacency matrix edge list,
# renaming attributes for clarity.
el <- adj2el_2(xi, directed = FALSE, names = colnames(adj)) %>% rename(xi = attr) %>% 
  left_join(adj2el_2(adj, directed = FALSE, names = colnames(adj)) %>% rename(count = attr), by = c("source", "target"))

# Calculate the total potential connections (m) and the squared sum of xi (sum_xi) for normalization.
m <- sum(el$count)
sum_xi <- 4 * m^2

# Fit the Block-Constrained Configuration Model (BCCM) using adjacency matrix and metadata,
# specifying model parameters such as homophily and self-loops.
model_bccm <- bccm(adj, labels = metadata$X2, homophily = TRUE, directed = FALSE, selfloops = FALSE, xi = xi)

# Merge BCCM omega values with the edge list and calculate normalized probabilities of connections.
el <- adj2el_2(model_bccm$omega, directed = FALSE, names = colnames(adj)) %>% rename(omega = attr) %>% left_join(el, by = c("source", "target"))
norm <- sum(el$xi * el$omega)
el <- el %>% mutate(prob = omega * xi / norm)

# Calculate indicators for edges being lower or higher than expected under various models
# (hypergeometric, BCCM, and Chung-Lu models) based on count, xi, and probability.
el %>% mutate(
  lower_than = Vectorize(phi_element_hyperg_under, vectorize.args = c('edgecount', 'xi'))(
    edgecount = count, xi = xi, sum_xi = sum_xi, m = m),
  higher_than = Vectorize(phi_element_hyperg_over, vectorize.args = c('edgecount', 'xi'))(
    edgecount = count, xi = xi, sum_xi = sum_xi, m = m),
  lower_than_bccm = Vectorize(phi_element_bccm_under, vectorize.args = c('edgecount', 'prob'))(
    edgecount = count, prob = prob, m = m),
  higher_than_bccm = Vectorize(phi_element_bccm_over, vectorize.args = c('edgecount', 'prob'))(
    edgecount = count, prob = prob, m = m),
  lower_than_CL = Vectorize(phi_element_CL_under, vectorize.args = c('edgecount', 'xi'))(
    edgecount = count, xi = xi, m = m),
  higher_than_CL = Vectorize(phi_element_CL_over, vectorize.args = c('edgecount', 'xi'))(
    edgecount = count, xi = xi, m = m)
) -> el

# Join phi calculations with metadata, adjusting source and target IDs for consistency.
# Also, integrate friendship status from the friendship network data.
left_join(el %>% mutate(source_ok = apply(cbind(source, target), 1, min), 
                        target_ok = apply(cbind(source, target), 1, max)), 
          el_fr %>% mutate(source = apply(cbind(X1, X2), 1, min), 
                           target = apply(cbind(X1, X2), 1, max)) %>% mutate(friend = TRUE) %>%
            group_by(source, target) %>% summarise(friend = sum(friend) > 0, .groups = 'drop'), 
          by = c('source_ok' = 'source', 'target_ok' = 'target')) %>%
  mutate(friend = if_else(is.na(friend), FALSE, friend)) %>%
  select(-xi, -source_ok, -target_ok, -omega, -prob) %>%
  left_join(metadata %>% rename(class_source = 'X2', gender_source = 'X3'), by = c('source' = 'X1')) %>%
  left_join(metadata %>% rename(class_target = 'X2', gender_target = 'X3'), by = c('target' = 'X1')) %>%
  mutate(class = apply(cbind(class_source, class_target), 1, function(vec) paste(c(min(vec), max(vec)), collapse = '_')),
         gender = apply(cbind(gender_source, gender_target), 1, function(vec) paste(c(min(vec), max(vec)), collapse = '_'))) -> dat

# Handle missing values for interaction count
dat <- dat %>% mutate(count = ifelse(is.na(count), 0, count))
dat %>% mutate(same_class = ifelse(class_source == class_target,1,0)) -> dat
###############################################################################
# Fitting to Friendship Data
# Description: This section fits a logistic regression model to predict friendship
# formation based on the calculated lower_than and higher_than phi coefficients.
# It begins by identifying active participants (nodes with outgoing connections)
# and then filters the dataset to include only interactions between these active participants.
# Finally, it fits a logistic regression model to understand the influence of network structure
# on friendship formation and creates a new variable 'phi' based on the model coefficients.
###############################################################################

# Identify active participants based on outgoing degree (deg_out) greater than 0.
# A tbl_graph is created from the nodes and edges, and the centrality_degree function
# calculates both the in-degree and out-degree for each node.
tbl_graph(nodes = tibble(name = nodes), edges = el_fr) %>% 
  mutate(deg_out = centrality_degree(mode = "out"),
         deg_in = centrality_degree(mode = "in")) %>% 
  filter(deg_out > 0) %>% as_tibble() %>% pull(name) -> participants

# Filter the dataset to include only interactions between active participants.
# The 'friend' variable is converted to a factor with levels "No" and "Yes" to 
# represent the absence or presence of a friendship, respectively.
dat %>% 
  filter(source %in% participants & target %in% participants) %>%
  mutate(friend = factor(friend, labels = c("FALSE" = "No", "TRUE" = "Yes"))) -> dat

# Fit a logistic regression model using the glm function with a binomial family,
# predicting the presence of friendship based on the lower_than and higher_than variables.
# This model aims to assess the influence of these network statistics on friendship formation.
logistic.under_over <- glm(friend ~ lower_than + higher_than, data = dat, family = binomial)

# Create a new variable 'phi' in the dataset based on the logistic regression model coefficients.
# This variable combines the effects of being lower_than and higher_than the expected values
# under the network models, providing a measure of the strength or direction of each relationship.
dat %>%
  mutate(phi = logistic.under_over$coefficients["lower_than"]*lower_than +
           logistic.under_over$coefficients["higher_than"]*higher_than) -> dat

###############################################################################
# Evaluation of Signed Links
# This section evaluates the predictive models' performance in identifying 
# friendships within the high school network. It includes both in-sample 
# evaluation and k-fold cross-validation to assess model accuracy and generalizability.
###############################################################################

# Source auxiliary functions for evaluation and cross-validation
source('functions/highschool_auxiliary_functions.R')

# Fit logistic regression models with different predictors
# Model using Chung-Lu derived predictors
logistic.under_over_CL <- glm(friend ~ lower_than_CL + higher_than_CL, data = dat, family = binomial)

# Model using BCCM predictors with interaction terms for class homophily
logistic.under_overbccm <- glm(friend ~ lower_than_bccm*same_class + higher_than_bccm*same_class, data = dat, family = binomial)

# Set a seed for reproducibility in k-fold cross-validation
seed <- 41

# Evaluate in-sample model performance using the evaluate_insample function
# Basic model evaluation
eva_under_over <- evaluate_insample(datt=dat, model=logistic.under_over, model_name="over_under", positive_level="Yes", negative_level="No")

# Chung-Lu model evaluation
eva_under_over_CL <- evaluate_insample(datt=dat, model=logistic.under_over_CL, model_name="over_under_CL", positive_level="Yes", negative_level="No")

# BCCM model with class homophily evaluation
eva_bccm <- evaluate_insample(datt=dat, model=logistic.under_overbccm, model_name="over_under_bccm", positive_level="Yes", negative_level="No")

# Combine results from in-sample evaluations for comparison
bind_rows(eva_under_over$res, eva_under_over_CL$res, eva_bccm$res)

# Perform k-fold cross-validation to evaluate model generalizability
# Basic model cross-validation
CV_under_over <- make_KfoldCV(dat=dat %>% mutate(friend=fct_recode(friend, "True"="Yes", "False"="No")),
                              formula=logistic.under_over$formula, split=.7, seed=seed, model_name="under_over")

# Chung-Lu model cross-validation
CV_under_over_CL <- make_KfoldCV(dat=dat %>% mutate(friend=fct_recode(friend, "True"="Yes", "False"="No")),
                                 formula=logistic.under_over_CL$formula, split=.7, seed=seed, model_name="under_over_CL")

# BCCM model with class homophily cross-validation
CV_bccm <- make_KfoldCV(dat=dat %>% mutate(friend=fct_recode(friend, "True"="Yes", "False"="No")),
                        formula=logistic.under_overbccm$formula, split=.7, seed=seed, model_name="bccm")

# Combine results from k-fold cross-validation for comparison
bind_rows(CV_under_over$res, CV_under_over_CL$res, CV_bccm$res)
