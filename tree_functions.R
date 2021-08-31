
## Some useful functions I use for
## handling trees, specially for plotting


date_edges = function(tree) {
  # Input: phyDat tree
  #        phenotype vector in order: first tips then nodes
  # Output: date of the external node of each edge
  
  t_root = tree$root.time
  n_root = length(tree$tip.label) + 1
  
  
  edges = tree$edge
  time_nodes = c()
  for (i in 1:nrow(edges)) {
    int_node = edges[i, 1]
    ext_node = edges[i, 2]
    
    n_length = tree$edge.length[i]
    if (int_node == n_root) {
      t_node = t_root + n_length
      time_nodes = c(time_nodes, setNames(t_node, ext_node))
    } else {
      t_node = time_nodes[[as.character(int_node)]] + n_length
      time_nodes = c(time_nodes, setNames(t_node, ext_node))
    }
  }
  edges = cbind(edges, unname(time_nodes))
  return(as.data.frame(edges))
}


