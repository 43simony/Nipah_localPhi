# install.packages('igraph'); install.packages('ggraph')
library(ggplot2)
library(dplyr)
library(igraph)
library(ggraph)


# Function to simulate a disease transmission tree
generate_transmission_network <- function(initial_cases = 1, max_generations = 5, max_infections = 3) {
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


transmissionTree_sim <- function(n_reps, parameters){
  
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
  nm = paste0("./transmissionNetwork_model.exe")
  r <- system2( nm, args = strvec, stdout = TRUE)
  
  ## read output from console
  out <- read.table(text = r, header = TRUE, sep = ';', check.names = FALSE) %>% mutate_all(as.numeric)
  
  return( out ) 
}



parameters <- data.frame(model_type = 2, R0 = 1, R0_SS = 15, prob_SS = 0.1, k = 0.5,
                         data_type = 0, reps = 10, N_threshold = 10, 
                         data_path = "na"
)

for(i in 1:1){
  if(i==1){
    tst <- transmissionTree_sim(n_reps = 1, parameters = parameters)
    tst$network_ID = i 
  }else{
    tmp <- transmissionTree_sim(n_reps = 1, parameters = parameters)
    tmp$network_ID = i 
    tst <- rbind(tst, tmp)
  }
  
  tst_graph <- graph_from_data_frame(tst, directed = TRUE)
}

# Plot network with primary case in the center
degreeOut = as.data.frame(degree(tst_graph, mode = 'out'))
names(degreeOut) <- c('degree')
degreeOut$degree[1] = degreeOut$degree[1] - 1 ## reduce degree of first node by 1 to correct for self-connection

tst_graph
ggraph(tst_graph, layout = "dendrogram", circular = TRUE) +
  geom_edge_link(arrow = arrow(length = unit(4, 'mm')), end_cap = circle(3, 'mm')) +
  geom_node_point(size = 5, shape = 21, fill = "white", color = "darkgreen", stroke = 1.5) +
  theme_void() +
  labs(title = "Simulated Disease Transmission Network")



## compute likelihood
network_LL = 0
for(i in 1:dim(degreeOut)[1]){
  
  d = as.numeric(degreeOut$degree[i])
  if(parameters$model_type == 0){ 
    ## compute likelihood for NB model
    network_LL = network_LL + log(dnbinom(x=d, mu = parameters$R0, size = parameters$k))
    
  }else if(parameters$model_type == 1){
    ## compute likelihood for poisson mixture model
    network_LL = network_LL + log( as.numeric(1-parameters$prob_SS)*dpois(x=d, lambda = as.numeric(parameters$R0)) + 
                                   as.numeric(parameters$prob_SS)*dpois(x=d, lambda = as.numeric(parameters$R0_SS)) )
  }else{
    ## compute likelihood for poisson-gamma mixture model
    
  }
  
}
network_LL
exp(network_LL)




# Generate network
graph <- generate_transmission_network()

# Plot network with primary case in the center
ggraph(ex, layout = "dendrogram", circular = TRUE) +
  geom_edge_link(arrow = arrow(length = unit(4, 'mm')), end_cap = circle(3, 'mm')) +
  geom_node_point(size = 5, shape = 21, fill = "white", color = "darkgreen", stroke = 1.5) +
  theme_void() +
  labs(title = "Simulated Disease Transmission Tree")

library(ggplot2)
library(hexbin)
library(dplyr)

# Generate clustered hexagon centers
hex_centers <- data.frame(
  x = c(1, 2, 3, 1.5, 2.5, 3.5, 2, 3, 4),
  y = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
  fill_value = rbeta(9, 0.1, 2) # Fill values between 0 and 1
)

# Create the plot
ggplot(hex_centers, aes(x = x, y = y, fill = fill_value)) +
  geom_hex(stat = "identity") +  # Use stat="identity" to assign values directly
  scale_fill_viridis(option = "plasma") + 
  theme_void() +
  xlim(c(0,5)) +
  ylim(c(0,5)) +
  labs(fill = "Value (0-1)")


install.packages("sf")
install.packages("ggplot2")
install.packages("gridExtra")
library(ggplot2)
library(sf)
library(gridExtra)

bd_data1<- read_sf("~/Desktop/Repos/Nipah_localPhi/model_data/BGD_adm1.shp")
map1 <- ggplot(bd_data1)+
  geom_sf(aes(fill=NAME_1))+labs(title="Divisional Fragments",x="Longitude",y="Lattitude")+
  theme(legend.position = "None")


remotes::install_github("ovirahman/bangladesh")
library(bangladesh)
library(tmap)
population <- bangladesh::pop_district_2011[, c("district", "population")]
district <- get_map("district")

map_data <- dplyr::left_join(district, population, by = c("District" = "district"))

country <- get_map("country")
division <- get_map("division")
district <- get_map("district")
upazila <- get_map("upazila")
union <- get_map("union")
