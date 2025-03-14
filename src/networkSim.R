# install.packages('igraph'); install.packages('ggraph')
library(ggplot2)
library(dplyr)
library(igraph)
library(ggraph)


# Function to simulate a disease transmission tree
generate_transmission_tree <- function(initial_cases = 1, max_generations = 5, max_infections = 3) {
  edges <- data.frame(from = integer(), to = integer())
  
  case_counter <- 1
  queue <- data.frame(case_num = 1, generation = 0)
  
  while (nrow(queue) > 0) {
    current <- queue[1, ]
    queue <- queue[-1, ]
    
    if (current$generation >= max_generations) next
    
    new_infections <- 0
    
    if(new_infections > 0){
      for (i in 1:new_infections) {
        case_counter <- case_counter + 1
        edges <- rbind(edges, data.frame(from = current$case_num, to = case_counter))
        queue <- rbind(queue, data.frame(case_num = case_counter, generation = current$generation + 1))
      }
    }
  }
  
  graph_from_data_frame(edges, directed = TRUE)
}


transmissiontree_sim <- function(n_reps, parameters){
  
  parvec = c( format(n_reps, scientific = F),
              parameters$model_type, # 0, 
              parameters$R0, # 1, 
              parameters$k, # 2, gamma shape parameter
              parameters$R0_SS, # 3, 
              parameters$prob_SS, # 4, 
              
              parameters$data_type, # 5, 
              parameters$reps, # 6, 
              parameters$N_threshold, # 7, 
              parameters$data_path # 9, 
  )
  
  strvec = format(parvec, digits = 5)
  
  setwd("~/Desktop/Repos/Nipah_localPhi/src") ## call has to be from location of .exe file or parameter read-in fails???
  
  ## Run the model
  ## The path to the bTB cpp binary file must be set correctly in the sys call below:
  nm = paste0("./transmissionTree_model.exe")
  r <- system2( nm, args = strvec, stdout = TRUE)
  
  ## read output from console
  out <- read.table(text = r, header = TRUE, sep = ';', check.names = FALSE) %>% mutate_all(as.numeric)
  
  return( out ) 
}



parameters <- data.frame(model_type = 2, R0 = 1, R0_SS = 15, prob_SS = 0.1, k = 0.5,
                         data_type = 0, reps = 10, N_threshold = 100, 
                         data_path = "na"
)

for(i in 1:1){
  if(i==1){
    tst <- transmissiontree_sim(n_reps = 1, parameters = parameters)
    tst$tree_ID = i 
  }else{
    tmp <- transmissiontree_sim(n_reps = 1, parameters = parameters)
    tmp$tree_ID = i 
    tst <- rbind(tst, tmp)
  }
  
  tst_graph <- graph_from_data_frame(tst, directed = TRUE)
}

# Plot tree with primary case in the center
degreeOut = as.data.frame(degree(tst_graph, mode = 'out'))
names(degreeOut) <- c('degree')
degreeOut$degree[1] = degreeOut$degree[1] - 1 ## reduce degree of first node by 1 to correct for self-connection

tst_graph
ggraph(tst_graph, layout = "dendrogram", circular = TRUE) +
  geom_edge_link(arrow = arrow(length = unit(4, 'mm')), end_cap = circle(3, 'mm')) +
  geom_node_point(size = 5, shape = 21, fill = "white", color = "darkgreen", stroke = 1.5) +
  theme_void() +
  labs(title = "Simulated Disease Transmission tree")



## compute likelihood
tree_LL = 0
for(i in 1:dim(degreeOut)[1]){
  
  d = as.numeric(degreeOut$degree[i])
  if(parameters$model_type == 0){ 
    ## compute likelihood for NB model
    tree_LL = tree_LL + log(dnbinom(x=d, mu = parameters$R0, size = parameters$k))
    
  }else if(parameters$model_type == 1){
    ## compute likelihood for poisson mixture model
    tree_LL = tree_LL + log( as.numeric(1-parameters$prob_SS)*dpois(x=d, lambda = as.numeric(parameters$R0)) + 
                                   as.numeric(parameters$prob_SS)*dpois(x=d, lambda = as.numeric(parameters$R0_SS)) )
  }else{
    ## compute likelihood for poisson-gamma mixture model
    
  }
  
}
tree_LL
exp(tree_LL)




# Generate tree
graph <- generate_transmission_tree()

# Plot tree with primary case in the center
ggraph(ex, layout = "dendrogram", circular = TRUE) +
  geom_edge_link(arrow = arrow(length = unit(4, 'mm')), end_cap = circle(3, 'mm')) +
  geom_node_point(size = 5, shape = 21, fill = "white", color = "darkgreen", stroke = 1.5) +
  theme_void() +
  labs(title = "Simulated Disease Transmission tree")


