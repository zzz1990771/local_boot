#This script generate weighted mean degree table in section 5.2
#graphs: graph 7 - 12
#column: true SE, EMB NB0.1 NB0.4 NB0.7 NB0.9 NBauto
#settings and libraries
source("../graph_boot_funs.R")
source("../graph_utils.r")
library(igraph)
library(parallel)
library(foreach)
#library(pushoverr)

#tool function for parallel
#note this is different from binary one
generate_graph <- function(x,u=NULL){
  if(x<=3){
    #D1
    W <- graphon_u_to_p_smoothness(size=size,smoothness = c(8,6,3.5)[x],sampling_on_u = node_sampling, u_input=u)
  }else if(x<=6){
    #graph10-12
    W <- graphon_u_to_p(size=size,pattern = "sinsin",smoothness= c(2,5,8)[x-3],sampling_on_u = node_sampling, u_input=u)
  }
  return(W)
}

getT <- function(adj.matrix){
  #v1 count self-loop once
  sum(adj.matrix)/(NCOL(adj.matrix)-1)
  #v2 count self-loop twice
  #(sum(adj.matrix)+sum(diag(adj.matrix)))/(NCOL(adj.matrix))
  
  #weigthed betweeness
  #mean(betweenness(graph_from_adjacency_matrix(adj.matrix,weighted = TRUE)))
  #mean(betweenness(graph_from_adjacency_matrix(adj.matrix,mode = "undirected")))
}

### parameter setting
cl_nodes <- 8
N <- 100
size <- N
M_true <- 50000
M_boot <- 128
B <- 500
node_sampling <- TRUE
neighbor_size_list = c(0.1,0.4,0.7,0.9)

cl<- makeCluster(cl_nodes)
doParallel::registerDoParallel(cl)

#generate weighted graph by overlay five graphs


for(i in c(23)){
  #true se
  time_start <- Sys.time()
  true.se <- foreach::foreach(m = (1:M_true),.combine = "c",.packages = "igraph") %dopar% {
    u = runif(n = N)
    W <- generate_graph(i,u=u)
    Adj <- Reduce("+",gmodel.P(W,rep=5))
    
    T.vec <- getT(Adj)
    return(T.vec)
  }
  print(var(true.se))
  print(Sys.time()-time_start)
  
  
  
  #boot se
  time_start <- Sys.time()
  se_boot <- foreach::foreach(m = (1:M_boot),.combine = "rbind",.packages = "igraph") %dopar% {
    u = runif(n = N)
    W <- generate_graph(i,u=u)
    Adj <- Reduce("+",gmodel.P(W,rep=5))
    
    se_vector <- rep(NA,length(neighbor_size_list)+1)
    ###boot
    ### 1.naive boot
    method_index=1
    boot.Tlist <- naive_boot(Adj,B,add_edge = FALSE,return = "getT",getT_function = getT)
    
    se_vector[1] <- var(unlist(boot.Tlist))
    
    for(j in 1:length(neighbor_size_list)){
      ### zhu boot
      boot.glist <- zhu_nb_boot(Adj,neighbor_size_list[j],B,kowning_u=NULL,induced_sampling = TRUE,weighted=TRUE)
      T.matrix.array.nb <- array(0,B)
      for(b in 1:B){
        adj.matrix <- boot.glist[[b]]
        T.matrix.array.nb[b] <- getT(adj.matrix)
      }
      #se hat
      se_vector[(j+1)] <- var(T.matrix.array.nb)
    }
    return(se_vector)
  }
  print(apply(se_boot,2,mean))
  print(Sys.time()-time_start)
}


stopCluster(cl)
