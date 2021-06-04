#This script generate clustering coefficient table in section 5
#graphs: graph 1 - 6
#column: true SE, EMB NB0.1 NB0.4 NB0.7 NB0.9 NBauto
#settings and libraries
library(parallel)
library(foreach)
library(igraph)
#library(viridis)
source('../graph_utils.r')
source('../graph_boot_funs.R')
#Sys.sleep(3600*5)
size = 200
M = 500
M_true = 1000000
B = 500
node_sampling = TRUE
cl_nodes = 32

#specific function for estimate T
getT <- function(adj.matrix){
  #for mean degree
  #mean(degree(graph_from_adjacency_matrix(adj.matrix,mode = "undirected")))
  
  #for clustering coefficient 
  (transitivity(graph_from_adjacency_matrix(adj.matrix,mode = "undirected")))
  
  #for mean betweenness
  #mean(betweenness(graph_from_adjacency_matrix(adj.matrix,mode = "undirected")))
}

#tool function for parallel
generate_graph <- function(x){
  if(x<=2){
    #graph1-2
    W <- graphon_u_to_p_smoothness2(size=size,smoothness = c(0.01,10)[x],sampling_on_u = node_sampling)
  }else if(x<=4){
    #graph3-4
    W <- graphon_u_to_p(size=size,pattern = "SAS6",smoothness = c(1000,1.5)[x-2],sampling_on_u = node_sampling)
  }else if(x<=6){
    #graph5-6
    W <- graphon_u_to_p(size=size,pattern = "sinsin",smoothness= c(8,10)[x-4],sampling_on_u = node_sampling)
  }
  return(W)
}


#generate true se
pattern_list = 1:6
time_start <- Sys.time()
TrueT_list <- mclapply(pattern_list, function(pattern) {
    true_T <- c() 
    for(m in c(1:M_true)){
      W <- generate_graph(pattern)
      adj.matrix <- gmodel.P(W,symmetric.out = TRUE)
      
      T <- getT(adj.matrix)
      true_T <- c(true_T,T)
    }
    return(var(true_T))
}, mc.cores=cl_nodes) 
time_end <- Sys.time()
print(time_end-time_start)

true_se_column <- unlist(TrueT_list)
true_se_column

#generate bootstrap se
time_start <- Sys.time()
pattern_list = 1:6
neighbor_size_list = seq(0.05,1,0.05)

res_mclapply_boot <- mclapply(pattern_list, function(pattern) {
  se_matrix <- matrix(0,nrow=M,ncol=(length(neighbor_size_list)+1))
  for(m in c(1:M)){
    W <- generate_graph(pattern)
    g.adj <- gmodel.P(W,symmetric.out = TRUE)
    ###boot
    ### 1.naive boot
    method_index=1
    boot.Tlist <- naive_boot(g.adj,B,add_edge = FALSE,return = "getT",getT_function = getT)
    
    se_matrix[m,1] <- var(unlist(boot.Tlist))
    
    for(j in 1:length(neighbor_size_list)){
      ### zhu boot
      boot.glist <- zhu_nb_boot(g.adj,neighbor_size_list[j],B,kowning_u=NULL,induced_sampling = TRUE,weighted=FALSE,getT = getT)
      T.matrix.array.nb <- unlist(boot.glist)
      #se hat
      se_matrix[m,(j+1)] <- var(T.matrix.array.nb)
    }
  }
  res <- apply(se_matrix,2,mean)
  return(res)
}, mc.cores=cl_nodes) 

time_end <- Sys.time()
print(time_end-time_start)
boot_se_column <- ((do.call("rbind",res_mclapply_boot)))
result_table <- cbind(true_se_column,boot_se_column)
result_table


