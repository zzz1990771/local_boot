#This script generate clustering coefficient table in section 3
#graphs: ER, SMB2, SBM10, MOD
#column: true SE, EMB NB0.1 NB0.4 NB0.7 NB0.9 NBauto
#settings and libraries
library(parallel)
library(foreach)
library(viridis)
source('../graph_utils.r')
source('../graph_boot_funs.R')
#Sys.sleep(3600*5)
size = 200
M = 500
M_true = 100000
B = 500
cl_nodes = 4

#specific function for estimate T
getT <- function(adj.matrix){
  #for mean degree
  #mean(degree(graph_from_adjacency_matrix(adj.matrix,mode = "undirected")))
  
  #for clustering coefficient 
  (transitivity(graph_from_adjacency_matrix(adj.matrix,mode = "undirected")))
}

#tool function for parallel
generate_graph <- function(x){
  if(x<=3){
    #SBM
    W <- graphon_u_to_p(size=size,pattern = "zhu1",smoothness= c(1,2,10)[x],sampling_on_u = TRUE)
  }else if(x<=4){
    #W <- graphon_u_to_p_smoothness2(size=size,smoothness = 100,sampling_on_u = TRUE)
    W <- graphon_u_to_p(size=size,pattern = "mod",smoothness= NULL,sampling_on_u = TRUE)
    
  }
  return(W)
}

#generate true se

cl<- makeCluster(cl_nodes)
doParallel::registerDoParallel(cl)
pattern_list = 1:4
time_start <- Sys.time()
TrueT_list <- foreach::foreach(pattern = pattern_list,.packages=c("igraph")) %dopar% {
    true_T <- c() 
    for(m in c(1:M_true)){
      W <- generate_graph(pattern)
      adj.matrix <- gmodel.P(W,symmetric.out = TRUE)
      
      T <- getT(adj.matrix)
      true_T <- c(true_T,T)
    }
    return(sd(true_T))
}
time_end <- Sys.time()
print(time_end-time_start)

true_se_column <- unlist(TrueT_list)
true_se_column

#generate bootstrap se
time_start <- Sys.time()
pattern_list = 1:4
neighbor_size_list = c(0.1,0.4,0.7,0.9)
BootT_list <- foreach::foreach(pattern = pattern_list,.packages=c("igraph")) %dopar% {
  se_matrix <- matrix(0,nrow=M,ncol=(length(neighbor_size_list)+1))
  for(m in c(1:M)){
    W <- generate_graph(pattern)
    g.adj <- gmodel.P(W,symmetric.out = TRUE)
    
    ### 1.naive boot
    method_index=1
    boot.glist <- naive_boot(g.adj,B)
    T.matrix.array.nb <- array(0,B)
    for(b in 1:B){
      adj.matrix <- boot.glist[[b]]
      T.matrix.array.nb[b] <- getT(adj.matrix)
    }
    se_matrix[m,1] <- sd(T.matrix.array.nb)
    
    for(i in 1:length(neighbor_size_list)){
      ### zhu boot
      method_index=4
      boot.glist <- zhu_nb_boot(g.adj,neighbor_size_list[i],B,method = "own",induced_sampling = TRUE)
      T.matrix.array.nb <- array(0,B)
      for(b in 1:B){
        adj.matrix <- boot.glist[[b]]
        T.matrix.array.nb[b] <- getT(adj.matrix)
      }
      #se hat
      se_matrix[m,(i+1)] <- sd(T.matrix.array.nb)
    }
  }
  res <- apply(se_matrix,2,mean)
  return(res)
}

time_end <- Sys.time()
print(time_end-time_start)
stopCluster(cl)
boot_se_column <- ((do.call("rbind",BootT_list)))
result_table <- cbind(true_se_column,boot_se_column)
result_table                    

