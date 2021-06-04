#######################bootstrap methods##########################

#naive methods
naive_boot <- function(adj_matrix,M,N=NROW(adj_matrix),removeincomplete=FALSE,add_edge=FALSE,resampling=TRUE,return="boot",getT_function=NULL){
  res <- list()
  while(M>0){
    if (resampling==TRUE){
      blist <- sample(1:(dim(adj_matrix)[1]),N, replace = T)
    }else{
      blist <- sample(1:(dim(adj_matrix)[1]),N, replace = FALSE)
    }
    
    #blist <- sort(blist)
    if(add_edge==FALSE){
      bgraph <- adj_matrix[blist,blist]
    }else{
      p<- sum(apply(adj_matrix,1,sum))/(N*(N-1))
      bgraph <- adj_matrix[blist,blist]
      for(i in 1:(N-1)){
        dup_index <- c((rep(FALSE,i)),(blist[(i+1):N]==blist[i]))
        fill.p <- (runif(sum(dup_index))<p)*1
        bgraph[i,dup_index]<- fill.p
        bgraph[dup_index,i]<- fill.p
      }
    }
    
    if(removeincomplete == TRUE){
      if(components(graph_from_adjacency_matrix(bgraph,mode = "undirected"))$no==1){
        M <- M-1
        res <- c(res, list(bgraph))
      }
    }else{
      if(return=="boot"){
        M <- M-1
        res <- c(res, list(bgraph))
      }else{
        M <- M-1
        res <- c(res, getT_function(bgraph))
        bgraph <- NULL
      }
      
    }
  }
  res
}

#embedding bootstrap

embedding_boot <- function(g.adj.matrix,M,d,direct=FALSE,N=NROW(g.adj.matrix),truncatep=FALSE,half=FALSE){
  res <- list()
  
  if(direct==FALSE){
    L=svd(g.adj.matrix)
    ud = cbind(L$u[,1:d],matrix(0,N,N-d))
    dd = rbind(cbind(diag(L$d[1:d]),matrix(0,d,N-d)),matrix(0,N-d,N))
    vd = rbind(t(L$v)[1:d,],matrix(0,N-d,N))
    g.adj.approx <- ud%*%dd%*%vd
    ## trancate to [0,1]
    
    
    #####truncated at 0.0001
    if(truncatep){
      zero <- 0.1
      one <- 0.9
    }else{
      zero <- 0
      one <- 1
    }
    g.adj.approx[g.adj.approx>one] <- 1
    g.adj.approx[g.adj.approx<zero] <- 0
    while (M>0){
      #induced sampling from approx matrix
      if(half){
        first_half <- sample(1:(dim(g.adj.approx)[1]),round(N/2,0), replace = FALSE)
        second_half <- setdiff(1:(dim(g.adj.approx)[1]),first_half)
        blist <- c(first_half,sample(second_half,length(second_half), replace = T) )
      }else{
        blist <- sample(1:(dim(g.adj.approx)[1]),N, replace = T)
      }
      #blist <- sort(blist)
      g.adj.approx_b <- g.adj.approx[blist,blist]
      #random matrix
      random.matrix <- matrix(runif(N*N,0,1),N,N)
      random.matrix[lower.tri(random.matrix)] <- t(random.matrix)[lower.tri(t(random.matrix))]
      ## diagonal set to zero
      diag(random.matrix) <- 100
      #induced sampled realization graph
      g.adj.id <- 1*(random.matrix<g.adj.approx_b)
      if(components(graph_from_adjacency_matrix(g.adj.id,mode = "undirected"))$no==1){
        M <- M-1
        res <- c(res, list(g.adj.id))
      }
    }
  }
  else{
    L=svd(g.adj.matrix)
    ud = cbind(L$u[,1:d],matrix(0,N,N-d))
    dd = rbind(cbind(diag(L$d[1:d]),matrix(0,d,N-d)),matrix(0,N-d,N))
    vd = rbind(t(L$v)[1:d,],matrix(0,N-d,N))
    g.adj.approx <- ud%*%dd%*%vd
    ## trancate to [0,1]
    g.adj.approx[g.adj.approx>1] <- 1
    g.adj.approx[g.adj.approx<0] <- 0
    while (M>0){
      random.matrix <- matrix(runif(N*N,0,1),N,N)
      random.matrix[lower.tri(random.matrix)] <- t(random.matrix)[lower.tri(t(random.matrix))]
      ## diagonal set to zero
      diag(random.matrix) <- 100
      g.adj.dd <- 1*(random.matrix<g.adj.approx)
      if(components(graph_from_adjacency_matrix(g.adj.dd,mode = "undirected"))$no==1){
        M <- M-1
        res <- c(res, list(g.adj.dd))
      }
    }
  }
  res
}

#degree boot
degree_boot <- function(g.adj.matrix){
  blist <- sample(1:(dim(g.adj.matrix)[1]),dim(g.adj.matrix)[1], replace = T)
  g.adj.matrix[blist,]
}

rewire_boot <- function(g.adj.matrix){
  
}

#neibor boot with embedding and l2 distance
emb_l2_nb_boot <- function(g.adj,quantile_n,B,embed_to=NROW(g.adj)){
  #quantile_n <- 0.05
  g.adj.matrix <- g.adj
  
  #embedding
  N=NROW(g.adj.matrix)
  d=embed_to
  L=svd(g.adj.matrix)
  ud = cbind(L$u[,1:d],matrix(0,N,N-d))
  dd = rbind(cbind(diag(L$d[1:d]),matrix(0,d,N-d)),matrix(0,N-d,N))
  vd = rbind(t(L$v)[1:d,],matrix(0,N-d,N))
  uddd <- ud%*%sqrt(dd)
  
  
  #find the neibors for each node
  neibors_matrix <- matrix(0,nrow = NROW(g.adj),ncol = (NROW(g.adj)*quantile_n))
  l2.matrix <- as.matrix(dist(uddd, method = "euclidean",diag = TRUE))
  
  
  for(i in (1:NROW(g.adj))){
    #neibor_index <- (l2.matrix[i,] <= quantile(l2.matrix[i,],0.05))
    neibor_index <- order(l2.matrix[i,], decreasing=FALSE)[1:(quantile_n*NROW(g.adj))]
    neibors_matrix[i,] <- neibor_index
  }
  
  #estimate p for each node pair
  p_hat_matrix <- matrix(0,nrow = NROW(g.adj),ncol = (NROW(g.adj)))
  for(a in 1:(-1+NROW(p_hat_matrix))){
    for(b in (a+1):(NROW(p_hat_matrix))){
      ab_connection <- g.adj[neibors_matrix[a,],neibors_matrix[b,]]
      p_hat_matrix[a,b] <- mean(ab_connection)
    }
  }
  for(a in 1:(NROW(p_hat_matrix))){
    b=a
    ab_connection <- g.adj[neibors_matrix[a,],neibors_matrix[b,]]
    p_hat_matrix[a,b] <- sum(ab_connection)/(NROW(ab_connection)*(NROW(ab_connection)-1))
  }
  p_hat_matrix[lower.tri(p_hat_matrix)] <- t(p_hat_matrix)[lower.tri(t(p_hat_matrix))]
  
  #induced sampling
  nb_boot.list <- list()
  for (i in 1:B){
    blist <- sample(1:(dim(p_hat_matrix)[1]),NROW(p_hat_matrix), replace = T)
    p_hat_matrix_b <- p_hat_matrix[blist,blist]
    N <- NROW(p_hat_matrix)
    random.matrix <- matrix(runif(N*N,0,1),N,N)
    random.matrix[lower.tri(random.matrix)] <- t(random.matrix)[lower.tri(t(random.matrix))]
    diag(random.matrix) <- 100
    g.adj.nb <- 1*(random.matrix<p_hat_matrix_b)
    nb_boot.list[[i]]<- g.adj.nb
  }
  return(nb_boot.list)
}


#neibor boot with zhu's distance
zhu_nb_boot <- function(A,quantile_n=0,B,returns = "boot",method = "own", distance = "zhu",
                        kowning_u=NULL, induced_sampling=TRUE, weighted=FALSE,getT=NULL,...){
  ## construct distance matrix in Zhu 2017
  now <- Sys.time()
  #SIZE
  N <- NROW(A)
  
  #max edge
  max_A <- max(A)

  #quantile_n
  if(quantile_n==0){
    quantile_n <- (log(N)/N)^0.5
  }
  
  #dissimilarity measure
  get_dist_zhu <- function(){
    dist.matrix <- matrix(0,nrow = N,ncol = N)
    A_sq = (A%*%A)/N
    for(i in c(1:(N-1))){
      for(ip in c((i+1):N) ){
        tgtvec = abs(A_sq[i,]-A_sq[ip,])
        tgtvec[i] = 0
        tgtvec[ip] = 0
        max_d = max(tgtvec) # tgtvec2 is Li Chen's
        dist.matrix[i,ip] <- max_d
      }
    }
    dist.matrix[lower.tri(dist.matrix)] <- t(dist.matrix)[lower.tri(t(dist.matrix))]
    dist.matrix
  }
  
  get_dist_test <- function(){
    dist.matrix <- matrix(0,nrow = N,ncol = N)
    dist.matrix <- as.matrix(dist(A, method = "euclidean",diag = TRUE))
    dist.matrix
    #dist.matrix[lower.tri(dist.matrix)] <- t(dist.matrix)[lower.tri(t(dist.matrix))]
  }
  
  if(is.null(kowning_u)){
    dist.matrix <- do.call(switch("zhu","zhu"="get_dist_zhu","test"= "get_dist_test"),args = list())
  }else{
    dist.matrix <- as.matrix(dist(kowning_u))
  }

  
  
  if(method=="own"){
    #get neighbors
    neibors_matrix <- matrix(0,nrow = NROW(A),ncol = (ceiling(quantile_n*N)))
    for(i in (1:N)){
      neibor_index <- order(dist.matrix[i,], decreasing=FALSE)[1:ceiling(quantile_n*N)]
      neibors_matrix[i,] <- neibor_index
    }
    #not weighted graph
    if(weighted==FALSE){
      #estimate p for each node pair
      p_hat_matrix <- matrix(0,nrow = NROW(A),ncol = (NROW(A)))
      for(a in 1:(-1+NROW(p_hat_matrix))){
        for(b in (a+1):(NROW(p_hat_matrix))){
          ab_connection <- A[neibors_matrix[a,],neibors_matrix[b,]]
          p_hat_matrix[a,b] <- mean(ab_connection)
          #p_hat_matrix[a,b] <- sample(as.vector(ab_connection),1)
        }
      }
      for(a in 1:(NROW(p_hat_matrix))){
        b=a
        ab_connection <- A[neibors_matrix[a,],neibors_matrix[b,]]
        #p_hat_matrix[a,b] <- sum(ab_connection)/(N*N-1))
        p_hat_matrix[a,b] <- mean(ab_connection)
        #p_hat_matrix[a,b] <- sample(as.vector(ab_connection),1)
      }
      p_hat_matrix[lower.tri(p_hat_matrix)] <- t(p_hat_matrix)[lower.tri(t(p_hat_matrix))]
    }else{
      #nb_array <- array(0,c(N,N,(ceiling(quantile_n*N)^2)))
      nb_array <- array(0,c(N,N,max_A+1))
      for(a in 1:(-1+N)){
        for(b in (a):(N)){
          ab_connection <- as.vector(A[neibors_matrix[a,],neibors_matrix[b,]])
          freq_table <- (table(ab_connection)/length(ab_connection))
          full_freq_table <- sapply(0:max_A,function(x){ifelse(sum(names(freq_table)==x)>0,freq_table[names(freq_table)==x],0)})
          #print(full_freq_table)
          
          nb_array[a,b,] <- cumsum(full_freq_table)
          nb_array[b,a,] <- nb_array[a,b,]
        }
      }
    }
    
  }else if(method=="zhu"){
    # 3. quantiled as logical
    kernel_mat = matrix(0,N,N)
    for (i in 1:N){
      kernel_mat[i,] = as.double(dist.matrix[i,]<quantile(dist.matrix[i,],quantile_n))
    }
    # 4. L1 normalization of each row
    kernel_mat = kernel_mat/(outer(rowSums(kernel_mat),rep(1,N))+1e-10)
    
    # 5. Compute P
    P = kernel_mat %*% A;
    P = (P+t(P))/2;
    p_hat_matrix <- P
  }
  finish_time <- Sys.time()-now
  if(returns=="p_and_time"){
    return(list(p_hat_matrix=p_hat_matrix,finish_time=finish_time))
  }
  
  #induced sampling
  nb_boot.list <- list()
  for (i in 1:B){
    if(induced_sampling){
      blist <- sample((1:N),N, replace = T)
    }else{
      blist <- 1:N
    }
    if(weighted==FALSE){
      p_hat_matrix_b <- p_hat_matrix[blist,blist]
      random.matrix <- matrix(runif(N*N,0,1),N,N)
      random.matrix[lower.tri(random.matrix)] <- t(random.matrix)[lower.tri(t(random.matrix))]
      #diag(random.matrix) <- 100
      g.adj.nb <- 1*(random.matrix<p_hat_matrix_b)
      if(is.null(getT)){
        nb_boot.list[[i]]<- g.adj.nb 
      }else{
        nb_boot.list[[i]]<- getT(g.adj.nb) 
      }
             
    }else{
      nb_array_b <- nb_array[blist,blist,]
      random.matrix <- matrix(runif(N*N,0,1),N,N)
      random.matrix[lower.tri(random.matrix)] <- t(random.matrix)[lower.tri(t(random.matrix))]
      g.adj.nb <- Reduce("+",lapply((1:(max_A+1)),function(x){
        1*(random.matrix>nb_array_b[,,x])
      }))
      if(is.null(getT)){
        nb_boot.list[[i]]<- g.adj.nb 
      }else{
        nb_boot.list[[i]]<- getT(g.adj.nb)  
      }
       
    }
    
  }
  if(returns=="p_and_boot"){
    return(list(p_hat_matrix=p_hat_matrix,nb_boot.list=nb_boot.list))
  }else{
    return(nb_boot.list)
  }
  
}

# weighted hamming bootstrap
wh_nb_boot <- function(g.adj,quantile_n=0,B,weigted_nerghbors=FALSE,hamming_weight="sd",returns = "boot",method = "own",test=FALSE,induced_sampling=TRUE){
  now <- Sys.time()
  #size
  A <- g.adj
  N <- NROW(A)
  
  #quantile_n
  if(quantile_n==0){
    quantile_n <- (log(N)/N)^0.5
  }
  
  #weighted hamming
  #cap the emp_p
  emp_p <- apply(A,2,sum)/(N-1)
  emp_p[emp_p==0] <- 1/N
  emp_p[emp_p==1] <- (1-(1/N))
  
  if(hamming_weight == "sd"){
    weight_p <- sqrt(1/sapply(emp_p,function(x){x*(1-x)}))
  }else{
    weight_p <- 1/sapply(emp_p,function(x){2*x*(1-x)})
  }
  
  weighted.h2.matrix <- matrix(0,nrow = N,ncol = N)
  #dp <- N*mean(emp_4p)
  #A_weighted <- A%*%diag(weight_p)
  # for(i in c(1:(N-1))){
  #   for(ip in c((i+1):N) ){
  #     #v0
  #     #weighted.h2.matrix[i,ip]<- sum((abs(A_i-A_ip))*weight_p)
  #     #v1
  #     #weighted.h2.matrix[i,ip]<- sum((abs(A_i-A_ip))*weight_p)+(N*(emp_p[i]-emp_p[ip]))
  #     #v2
  #     weighted.h2.matrix[i,ip]<- sum((abs(A_weighted[i,]-A_weighted[ip,])))+((emp_p[i]-emp_p[ip]))*dp 
  #     #weighted.h2.matrix[i,ip]<- sum((abs(A[i,]-A[ip,]))*weight_p)+((emp_p[i]-emp_p[ip]))*dp
  #     
  #   }
  # }
  # weighted.h2.matrix[lower.tri(weighted.h2.matrix)] <- t(weighted.h2.matrix)[lower.tri(t(weighted.h2.matrix))]
  D <- (1 - A) %*% t(A%*%diag(weight_p))
  weighted.h2.matrix <- D + t(D)
  if(test){
    weighted.h2.matrix <- matrix(0,nrow = N,ncol = N)
    for(i in c(1:(N-1))){
      for(ip in c((i+1):N) ){
        weighted.h2.matrix[i,ip]<- sum((abs(A_i-A_ip)))
        
      }
    }
    weighted.h2.matrix[lower.tri(weighted.h2.matrix)] <- t(weighted.h2.matrix)[lower.tri(t(weighted.h2.matrix))]
  }
  
  #estimation:find neighbors and estimate
  if(method=="own"){
    #find the neighbors for each node
    neibors_matrix <- matrix(0,nrow = NROW(g.adj),ncol = (ceiling(quantile_n*N)))
    for(i in (1:N)){
      #neibor_index <- (weighted.h2.matrix[i,] <= quantile(weighted.h2.matrix[i,],0.05))
      neibor_index <- order(weighted.h2.matrix[i,], decreasing=FALSE)[1:ceiling(quantile_n*N)]
      neibors_matrix[i,] <- neibor_index
    }
    
    #estimate p for each node pair
    p_hat_matrix <- matrix(0,nrow = NROW(g.adj),ncol = (NROW(g.adj)))
    for(a in 1:(-1+NROW(p_hat_matrix))){
      for(b in (a+1):(NROW(p_hat_matrix))){
        ab_connection <- g.adj[neibors_matrix[a,],neibors_matrix[b,]]
        if(weigted_nerghbors == FALSE){
          p_hat_matrix[a,b] <- mean(ab_connection)
        }else{
          a_distance <- weighted.h2.matrix[a,neibors_matrix[a,]]
          b_distance <- weighted.h2.matrix[b,neibors_matrix[b,]]
          c=1
          alpha=0.5
          nb_weight <- matrix(0,nrow = length(a_distance),ncol = (length(b_distance)))
          for(i in 1:length(a_distance)){
            for(j in 1:length(b_distance)){
              nb_weight[i,j] <- c*(1+a_distance[i]+b_distance[j])^(-alpha)
            }
          }
          ab_connection <- ab_connection*nb_weight
          p_hat_matrix[a,b] <- sum(ab_connection)/sum(nb_weight)
        }
        
      }
    }
    for(a in 1:(NROW(p_hat_matrix))){
      b=a
      ab_connection <- g.adj[neibors_matrix[a,],neibors_matrix[b,]]
      p_hat_matrix[a,b] <- sum(ab_connection)/(NROW(ab_connection)*(NROW(ab_connection)-1))
    }
    p_hat_matrix[lower.tri(p_hat_matrix)] <- t(p_hat_matrix)[lower.tri(t(p_hat_matrix))]
  }else if(method == "local"){
    #degree  and  wh at the same time by intersect

    #emp_p distance
    emp_p.dist <- as.matrix(dist(emp_p,method = "euclidean",diag = TRUE))
    #find the neighbors for each node
    #method 1, changing nb size for both distance
    # neibors_matrix <- matrix(0,nrow = NROW(g.adj),ncol = (ceiling(quantile_n*N)))
    # for(i in (1:N)){
    #   neibor_index_first <- order(emp_p.dist[i,], decreasing=FALSE)
    #   neibor_index_second <- order(weighted.h2.matrix[i,], decreasing=FALSE)

    #   start_n=target_n <- (ceiling(quantile_n*N))
    #   final_set <- intersect(neibor_index_first[1:start_n],neibor_index_second[1:start_n])
    #   while(length(final_set)< target_n){
    #     start_n = start_n+1
    #     final_set <- intersect(neibor_index_first[1:start_n],neibor_index_second[1:target_n])
    #   }
    #   neibors_matrix[i,] <- final_set

    # }

    #method 2, changing nb size for wh distance, fixing degree
    neibors_matrix <- matrix(0,nrow = NROW(g.adj),ncol = (ceiling(quantile_n*N)))
    target_n <- (ceiling(quantile_n*N))
    degree_n <- ifelse(2*target_n>=N,N,2*target_n)
    for(i in (1:N)){
      neibor_index_first <- order(emp_p.dist[i,], decreasing=FALSE)[1:degree_n]
      neibor_index_second <- order(weighted.h2.matrix[i,neibor_index_first], decreasing=FALSE)[1:target_n]
      neibors_matrix[i,] <- neibor_index_first[neibor_index_second]

    }
    
    #estimate p for each node pair
    p_hat_matrix <- matrix(0,nrow = NROW(g.adj),ncol = (NROW(g.adj)))
    for(a in 1:(-1+NROW(p_hat_matrix))){
      for(b in (a+1):(NROW(p_hat_matrix))){
        ab_connection <- g.adj[neibors_matrix[a,],neibors_matrix[b,]]
        if(weigted_nerghbors == FALSE){
          p_hat_matrix[a,b] <- mean(ab_connection)
        }else{
          a_distance <- weighted.h2.matrix[a,neibors_matrix[a,]]
          b_distance <- weighted.h2.matrix[b,neibors_matrix[b,]]
          c=1
          alpha=0.5
          nb_weight <- matrix(0,nrow = length(a_distance),ncol = (length(b_distance)))
          for(i in 1:length(a_distance)){
            for(j in 1:length(b_distance)){
              nb_weight[i,j] <- c*(1+a_distance[i]+b_distance[j])^(-alpha)
            }
          }
          ab_connection <- ab_connection*nb_weight
          p_hat_matrix[a,b] <- sum(ab_connection)/sum(nb_weight)
        }
        
      }
    }
    for(a in 1:(NROW(p_hat_matrix))){
      b=a
      ab_connection <- g.adj[neibors_matrix[a,],neibors_matrix[b,]]
      p_hat_matrix[a,b] <- sum(ab_connection)/(NROW(ab_connection)*(NROW(ab_connection)-1))
    }
    p_hat_matrix[lower.tri(p_hat_matrix)] <- t(p_hat_matrix)[lower.tri(t(p_hat_matrix))]
  }else if(method == "local1"){
    #degree first and then wh

    #emp_p distance
    emp_p.dist <- as.matrix(dist(emp_p,method = "euclidean",diag = TRUE))
    #find the neighbors for each node
    neibors_matrix <- matrix(0,nrow = NROW(g.adj),ncol = (ceiling(quantile_n*N)))
    for(i in (1:N)){
      # #neibor_index <- (weighted.h2.matrix[i,] <= quantile(weighted.h2.matrix[i,],0.05))
      # neibor_index_first <- order(weighted.h2.matrix[i,], decreasing=FALSE)[1:ceiling(quantile_n*N*2)]
      # #index second is the index of index first
      # neibor_index_second <- order(emp_p.dist[i,neibor_index_first], decreasing=FALSE)[1:ceiling(quantile_n*N)]
      # neibors_matrix[i,] <- neibor_index_first[neibor_index_second]


      neibor_index_first <- order(emp_p.dist[i,], decreasing=FALSE)[1:(2*target_n)]
      neibor_index_second <- order(weighted.h2.matrix[i,neibor_index_first], decreasing=FALSE)[1:target_n]
      neibors_matrix[i,] <- neibor_index_first[neibor_index_second]
    }
    
    #estimate p for each node pair
    p_hat_matrix <- matrix(0,nrow = NROW(g.adj),ncol = (NROW(g.adj)))
    for(a in 1:(-1+NROW(p_hat_matrix))){
      for(b in (a+1):(NROW(p_hat_matrix))){
        ab_connection <- g.adj[neibors_matrix[a,],neibors_matrix[b,]]
        if(weigted_nerghbors == FALSE){
          p_hat_matrix[a,b] <- mean(ab_connection)
        }else{
          a_distance <- weighted.h2.matrix[a,neibors_matrix[a,]]
          b_distance <- weighted.h2.matrix[b,neibors_matrix[b,]]
          c=1
          alpha=0.5
          nb_weight <- matrix(0,nrow = length(a_distance),ncol = (length(b_distance)))
          for(i in 1:length(a_distance)){
            for(j in 1:length(b_distance)){
              nb_weight[i,j] <- c*(1+a_distance[i]+b_distance[j])^(-alpha)
            }
          }
          ab_connection <- ab_connection*nb_weight
          p_hat_matrix[a,b] <- sum(ab_connection)/sum(nb_weight)
        }
        
      }
    }
    for(a in 1:(NROW(p_hat_matrix))){
      b=a
      ab_connection <- g.adj[neibors_matrix[a,],neibors_matrix[b,]]
      p_hat_matrix[a,b] <- sum(ab_connection)/(NROW(ab_connection)*(NROW(ab_connection)-1))
    }
    p_hat_matrix[lower.tri(p_hat_matrix)] <- t(p_hat_matrix)[lower.tri(t(p_hat_matrix))]
  }else if(method == "adaptive_fixN"){
    #emp_p distance
    emp_p.dist <- as.matrix(dist(emp_p,method = "euclidean",diag = TRUE))
    #find the neighbors for each node
    neibors_matrix <- matrix(0,nrow = NROW(g.adj),ncol = (ceiling(quantile_n*N)))
    for(i in (1:N)){
      length_nb = 0
      start_N = ceiling(quantile_n*N) - 10
      while(length_nb < ceiling(quantile_n*N)){
        start_N = start_N + 10
        neibor_index_first <- order(weighted.h2.matrix[i,], decreasing=FALSE)[1:start_N]
        #index second is the index of emp
        neibor_index_second <- order(emp_p.dist[i,], decreasing=FALSE)[1:start_N]
        #intersect
        candidate_neibor <- intersect(neibor_index_first,neibor_index_second)
        length_nb <- length(candidate_neibor)
        #print(length_nb)
      }
      neibors_matrix[i,] <- candidate_neibor[c(1: ceiling(quantile_n*N))]
    }
    
    #estimate p for each node pair
    p_hat_matrix <- matrix(0,nrow = NROW(g.adj),ncol = (NROW(g.adj)))
    for(a in 1:(-1+NROW(p_hat_matrix))){
      for(b in (a+1):(NROW(p_hat_matrix))){
        ab_connection <- g.adj[neibors_matrix[a,],neibors_matrix[b,]]
        if(weigted_nerghbors == FALSE){
          p_hat_matrix[a,b] <- mean(ab_connection)
        }else{
          a_distance <- weighted.h2.matrix[a,neibors_matrix[a,]]
          b_distance <- weighted.h2.matrix[b,neibors_matrix[b,]]
          c=1
          alpha=0.5
          nb_weight <- matrix(0,nrow = length(a_distance),ncol = (length(b_distance)))
          for(i in 1:length(a_distance)){
            for(j in 1:length(b_distance)){
              nb_weight[i,j] <- c*(1+a_distance[i]+b_distance[j])^(-alpha)
            }
          }
          ab_connection <- ab_connection*nb_weight
          p_hat_matrix[a,b] <- sum(ab_connection)/sum(nb_weight)
        }
        
      }
    }
    for(a in 1:(NROW(p_hat_matrix))){
      b=a
      ab_connection <- g.adj[neibors_matrix[a,],neibors_matrix[b,]]
      p_hat_matrix[a,b] <- sum(ab_connection)/(NROW(ab_connection)*(NROW(ab_connection)-1))
    }
    p_hat_matrix[lower.tri(p_hat_matrix)] <- t(p_hat_matrix)[lower.tri(t(p_hat_matrix))]
  }else if(method == "adaptive"){
    #emp_p distance
    emp_p.dist <- as.matrix(dist(emp_p,method = "euclidean",diag = TRUE))
    #find the neighbors for each node
    neibors_matrix <- matrix(0,nrow = NROW(g.adj),ncol = (ceiling(quantile_n*N)))
    for(i in (1:N)){
      length_nb = 0
      start_N = (ceiling(quantile_n*N)) - 10
      while(length_nb < 3){
        start_N = start_N + 10
        neibor_index_first <- order(weighted.h2.matrix[i,], decreasing=FALSE)[1:start_N]
        #index second is the index of emp
        neibor_index_second <- order(emp_p.dist[i,], decreasing=FALSE)[1:start_N]
        #intersect
        candidate_neibor <- intersect(neibor_index_first,neibor_index_second)
        length_nb <- length(candidate_neibor)
        #print(length_nb)
      }
      neibors_matrix[i,c(1:length_nb)] <- candidate_neibor
    }
    

    #estimate p for each node pair
    p_hat_matrix <- matrix(0,nrow = NROW(g.adj),ncol = (NROW(g.adj)))
    for(a in 1:(-1+NROW(p_hat_matrix))){
      for(b in (a+1):(NROW(p_hat_matrix))){
        ab_connection <- g.adj[neibors_matrix[a,],neibors_matrix[b,]]
        if(weigted_nerghbors == FALSE){
          p_hat_matrix[a,b] <- mean(ab_connection)
        }else{
          a_distance <- weighted.h2.matrix[a,neibors_matrix[a,]]
          b_distance <- weighted.h2.matrix[b,neibors_matrix[b,]]
          c=1
          alpha=0.5
          nb_weight <- matrix(0,nrow = length(a_distance),ncol = (length(b_distance)))
          for(i in 1:length(a_distance)){
            for(j in 1:length(b_distance)){
              nb_weight[i,j] <- c*(1+a_distance[i]+b_distance[j])^(-alpha)
            }
          }
          ab_connection <- ab_connection*nb_weight
          p_hat_matrix[a,b] <- sum(ab_connection)/sum(nb_weight)
        }
        
      }
    }
    for(a in 1:(NROW(p_hat_matrix))){
      b=a
      ab_connection <- g.adj[neibors_matrix[a,],neibors_matrix[b,]]
      p_hat_matrix[a,b] <- sum(ab_connection)/(NROW(ab_connection)*(NROW(ab_connection)-1))
    }
    p_hat_matrix[lower.tri(p_hat_matrix)] <- t(p_hat_matrix)[lower.tri(t(p_hat_matrix))]
  }else if(method == "degree_only"){
    #emp_p distance
    emp_p.dist <- as.matrix(dist(emp_p,method = "euclidean",diag = TRUE))
    #find the neighbors for each node
    neibors_matrix <- matrix(0,nrow = NROW(g.adj),ncol = (ceiling(quantile_n*N)))
    for(i in (1:N)){
      neibor_index <- order(emp_p.dist[i,], decreasing=FALSE)[1:ceiling(quantile_n*N)]
      neibors_matrix[i,] <- neibor_index
    }
    
    #estimate p for each node pair
    p_hat_matrix <- matrix(0,nrow = NROW(g.adj),ncol = (NROW(g.adj)))
    for(a in 1:(-1+NROW(p_hat_matrix))){
      for(b in (a+1):(NROW(p_hat_matrix))){
        ab_connection <- g.adj[neibors_matrix[a,],neibors_matrix[b,]]
        if(weigted_nerghbors == FALSE){
          p_hat_matrix[a,b] <- mean(ab_connection)
        }else{
          a_distance <- weighted.h2.matrix[a,neibors_matrix[a,]]
          b_distance <- weighted.h2.matrix[b,neibors_matrix[b,]]
          c=1
          alpha=0.5
          nb_weight <- matrix(0,nrow = length(a_distance),ncol = (length(b_distance)))
          for(i in 1:length(a_distance)){
            for(j in 1:length(b_distance)){
              nb_weight[i,j] <- c*(1+a_distance[i]+b_distance[j])^(-alpha)
            }
          }
          ab_connection <- ab_connection*nb_weight
          p_hat_matrix[a,b] <- sum(ab_connection)/sum(nb_weight)
        }
        
      }
    }
    for(a in 1:(NROW(p_hat_matrix))){
      b=a
      ab_connection <- g.adj[neibors_matrix[a,],neibors_matrix[b,]]
      p_hat_matrix[a,b] <- sum(ab_connection)/(NROW(ab_connection)*(NROW(ab_connection)-1))
    }
    p_hat_matrix[lower.tri(p_hat_matrix)] <- t(p_hat_matrix)[lower.tri(t(p_hat_matrix))]
  }else if(method=="moving_average"){
    #sort the graph by degree
    emp_p_order <- order(emp_p,decreasing = TRUE)
    sorted_adj <- A[emp_p_order,emp_p_order]
    #remove diag
    diag(sorted_adj) <- NA
    sorted_adj_nodiag <- matrix(na.omit(c(t(sorted_adj))),ncol=(N-1),byrow = TRUE)

    #estimate p vector for each node with moving window
    windown_n <- ceiling(quantile_n*N)
    est_mov_p <- matrix(NA,N,(N-windown_n))
    for(i in 1:N){
        est_mov_p[i,] <- zoo::rollmean(sorted_adj_nodiag[i,],windown_n)
    }

    
    #moving_p distance
    mov_p.dist <- as.matrix(dist(est_mov_p,method = "euclidean",diag = TRUE))
    #find the neighbors for each node
    neibors_matrix <- matrix(0,nrow = NROW(g.adj),ncol = (ceiling(quantile_n*N)))
    for(i in (1:N)){
      neibor_index <- order(mov_p.dist[i,], decreasing=FALSE)[1:ceiling(quantile_n*N)]
      neibors_matrix[i,] <- neibor_index
    }
    
    #estimate p for each node pair
    p_hat_matrix <- matrix(0,nrow = NROW(g.adj),ncol = (NROW(g.adj)))
    for(a in 1:(-1+NROW(p_hat_matrix))){
      for(b in (a+1):(NROW(p_hat_matrix))){
        ab_connection <- g.adj[neibors_matrix[a,],neibors_matrix[b,]]
        if(weigted_nerghbors == FALSE){
          p_hat_matrix[a,b] <- mean(ab_connection)
        }else{
          a_distance <- weighted.h2.matrix[a,neibors_matrix[a,]]
          b_distance <- weighted.h2.matrix[b,neibors_matrix[b,]]
          c=1
          alpha=0.5
          nb_weight <- matrix(0,nrow = length(a_distance),ncol = (length(b_distance)))
          for(i in 1:length(a_distance)){
            for(j in 1:length(b_distance)){
              nb_weight[i,j] <- c*(1+a_distance[i]+b_distance[j])^(-alpha)
            }
          }
          ab_connection <- ab_connection*nb_weight
          p_hat_matrix[a,b] <- sum(ab_connection)/sum(nb_weight)
        }
        
      }
    }
    for(a in 1:(NROW(p_hat_matrix))){
      b=a
      ab_connection <- g.adj[neibors_matrix[a,],neibors_matrix[b,]]
      p_hat_matrix[a,b] <- sum(ab_connection)/(NROW(ab_connection)*(NROW(ab_connection)-1))
    }
    p_hat_matrix[lower.tri(p_hat_matrix)] <- t(p_hat_matrix)[lower.tri(t(p_hat_matrix))]
  }else if(method=="zhu"){
    # 3. quantiled as logical
    kernel_mat = matrix(0,N,N)
    for (i in 1:N){
      kernel_mat[i,] = as.double(weighted.h2.matrix[i,]<quantile(weighted.h2.matrix[i,],quantile_n))
    }
    # 4. L1 normalization of each row
    kernel_mat = kernel_mat/(outer(rowSums(kernel_mat),rep(1,N))+1e-10)
    
    # 5. Compute P
    P = kernel_mat %*% A;
    P = (P+t(P))/2;
    p_hat_matrix <- P
  }
  finish_time <- Sys.time()-now
  
  if(returns=="p_and_time"){
    return(list(p_hat_matrix=p_hat_matrix,finish_time=finish_time))
  }
  
  
  


  nb_boot.list <- list()
  for (i in 1:B){
    #induced sampling
    if(induced_sampling){
      blist <- sample(1:(dim(p_hat_matrix)[1]),NROW(p_hat_matrix), replace = T)
    }else{
      blist <- 1:NROW(p_hat_matrix)
    }
    p_hat_matrix_b <- p_hat_matrix[blist,blist]
    N <- NROW(p_hat_matrix)
    random.matrix <- matrix(runif(N*N,0,1),N,N)
    random.matrix[lower.tri(random.matrix)] <- t(random.matrix)[lower.tri(t(random.matrix))]
    diag(random.matrix) <- 100
    g.adj.nb <- 1*(random.matrix<p_hat_matrix_b)
    nb_boot.list[[i]]<- g.adj.nb
  }
  if(returns=="p_and_boot"){
    return(list(p_hat_matrix=p_hat_matrix,nb_boot.list=nb_boot.list))
  }else{
    return(nb_boot.list)
  }
}


#neibor boot with iterative p matrix 

ice_boot<- l2p_nb_boot <- function(g.adj,quantile_n=0,B,returns = "boot",method = "own",kowning_u=FALSE,induced_sampling=TRUE,...){
  ## construct distance matrix in Zhu 2017
  now <- Sys.time()
  #SIZE
  A <- g.adj
  N <- NROW(A)
  
  #quantile_n
  if(quantile_n==0){
    quantile_n <- (log(N)/N)^0.5
  }
  
  #dissimilarity measure on a random p
  #initial_p <- matrix(runif(N^2,0,1),ncol=N,nrow=N)
  #initial_p[lower.tri(initial_p)] <- t(initial_p)[lower.tri(t(initial_p))]
  initial_p <- zhu_nb_boot(g.adj,returns = "p_and_time")$p_hat_matrix
  maxit = 100
  it=1
  old_p_hat_matrix = matrix(0,ncol=N,nrow=N)
  p_hat_matrix = initial_p
  #print(norm(p_hat_matrix-old_p_hat_matrix,type = 'f')/norm(old_p_hat_matrix,type='f'))
  while(it<=maxit & norm(p_hat_matrix-old_p_hat_matrix,type = 'f')/norm(old_p_hat_matrix,type='f')>0.02){
    #print(norm(p_hat_matrix-old_p_hat_matrix,type = 'f')/norm(old_p_hat_matrix,type='f'))
    old_p_hat_matrix = p_hat_matrix
    
    dist.matrix <- as.matrix(dist(p_hat_matrix,diag = TRUE))/N
    
    if(method=="own"){
      #get neighbors
      neibors_matrix <- matrix(0,nrow = NROW(g.adj),ncol = (ceiling(quantile_n*N)))
      for(i in (1:N)){
        neibor_index <- order(dist.matrix[i,], decreasing=FALSE)[1:ceiling(quantile_n*N)]
        neibors_matrix[i,] <- neibor_index
      }
      #estimate p for each node pair
      p_hat_matrix <- matrix(0,nrow = NROW(g.adj),ncol = (NROW(g.adj)))
      for(a in 1:(-1+NROW(p_hat_matrix))){
        for(b in (a+1):(NROW(p_hat_matrix))){
          ab_connection <- g.adj[neibors_matrix[a,],neibors_matrix[b,]]
          p_hat_matrix[a,b] <- mean(ab_connection)
        }
      }
      for(a in 1:(NROW(p_hat_matrix))){
        b=a
        ab_connection <- g.adj[neibors_matrix[a,],neibors_matrix[b,]]
        p_hat_matrix[a,b] <- sum(ab_connection)/(NROW(ab_connection)*(NROW(ab_connection)-1))
      }
      p_hat_matrix[lower.tri(p_hat_matrix)] <- t(p_hat_matrix)[lower.tri(t(p_hat_matrix))]
    }else if(method=="zhu"){
      # 3. quantiled as logical
      kernel_mat = matrix(0,N,N)
      for (i in 1:N){
        kernel_mat[i,] = as.double(dist.matrix[i,]<quantile(dist.matrix[i,],quantile_n))
      }
      # 4. L1 normalization of each row
      kernel_mat = kernel_mat/(outer(rowSums(kernel_mat),rep(1,N))+1e-10)
      
      # 5. Compute P
      P = kernel_mat %*% A;
      P = (P+t(P))/2;
      p_hat_matrix <- P
    }
    
    it = it + 1
  }
  
  finish_time <- Sys.time()-now
  if(returns=="p_and_time"){
    return(list(p_hat_matrix=p_hat_matrix,finish_time=finish_time))
  }
  
  #induced sampling
  nb_boot.list <- list()
  for (i in 1:B){
    if(induced_sampling){
      blist <- sample(1:(dim(p_hat_matrix)[1]),NROW(p_hat_matrix), replace = T)
    }else{
      blist <- 1:NROW(p_hat_matrix)
    }
    p_hat_matrix_b <- p_hat_matrix[blist,blist]
    N <- NROW(p_hat_matrix)
    random.matrix <- matrix(runif(N*N,0,1),N,N)
    random.matrix[lower.tri(random.matrix)] <- t(random.matrix)[lower.tri(t(random.matrix))]
    diag(random.matrix) <- 100
    g.adj.nb <- 1*(random.matrix<p_hat_matrix_b)
    nb_boot.list[[i]]<- g.adj.nb
  }
  if(returns=="p_and_boot"){
    return(list(p_hat_matrix=p_hat_matrix,nb_boot.list=nb_boot.list))
  }else{
    return(nb_boot.list)
  }
  
}