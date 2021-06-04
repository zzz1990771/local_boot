#######################graph utils##########################

#generate p matrix from preset graphon function
graphon_u_to_p <- function(size,pattern,smoothness=NULL,return_u = FALSE,sampling_on_u=TRUE,u_input=NULL){
  
  if(!is.null(u_input)){u = u_input}else{
    if(sampling_on_u){
      u = runif(n = size)
      u = sort(u)
    }else{
      #u = seq(from = 0, to = 1,length.out=size)
      #this is to avoid 0 and 1
      u = seq(from = 0, to = 1,length.out=(size+2))[2:(size+1)]
    }
  }
  
  if(pattern == "zhu1"){
    if(is.null(smoothness)){
      K = 5 
    }else{
      K = smoothness
    }
    p_matrix = matrix(0,nrow=size,ncol=size)
    for(i in 1:size){
      for(j in 1:size){
        if(floor(K*u[i]+1) == floor(K*u[j]+1)){
          p <- floor(K*u[i]+1)/(K+2)
        }else{
          p <- 0.5*(floor(K*u[i]+1)/(K+5) + floor(K*u[j]+1)/(K+5))
        }
        p_matrix[i,j] <- p    
      }
    } 
  }
  
  if(pattern == "sbmk10"){
    K=10
    if(is.null(smoothness)){
      smoothness = 1 
    }
    p_matrix = matrix(0,nrow=size,ncol=size)
    for(i in 1:size){
      for(j in 1:size){
        if(floor(K*u[i]+1) == floor(K*u[j]+1)){
          p <- smoothness
        }else{
          p <- 1-smoothness
        }
        p_matrix[i,j] <- p    
      }
    } 
  }
  
  if(pattern == "sbmk100"){
    K=100
    if(is.null(smoothness)){
      smoothness = 1 
    }
    p_matrix = matrix(0,nrow=size,ncol=size)
    for(i in 1:size){
      for(j in 1:size){
        if(floor(K*u[i]+1) == floor(K*u[j]+1)){
          p <- smoothness
        }else{
          p <- 1-smoothness
        }
        p_matrix[i,j] <- p    
      }
    } 
  }
  
  if(pattern == "mod"){
    if(is.null(smoothness)){
      K = 5 
    }else{
      K = smoothness
    }
    p_matrix = matrix(0,nrow=size,ncol=size)
    for(i in 1:size){
      for(j in 1:size){
        u_class <- floor(u[i]*size*0.5+1)
        v_class <- floor(u[j]*size*0.5+1)
        p <- (sin(u_class)+ sin(v_class))*(100/(u_class+v_class))
        p <- 1/(exp(-p)+1)
        p_matrix[i,j] <- p    
      }
    } 
  }

  if(pattern == "sbm_round"){
    if(is.null(smoothness)){
      K = 5 
    }else{
      K = smoothness
      r = 1/(2*K)
    }
    p_matrix = matrix(0,nrow=size,ncol=size)
    for(i in 1:size){
      for(j in 1:size){
        if(floor(K*u[i]+1) == floor(K*u[j]+1)){
          k <- floor(K*u[i]+1)
          center <- 2*r*(k-1) + r
          p_base <- floor(K*u[i]+1)/(K+2)
          suppressWarnings(p <- (sqrt((r)^2-((u[i]-center)/1.2)^2-((u[j]-center)/1.2)^2)/r)*p_base)
          p <- if (is.na(p)) 0 else p
        }else{
          center_ui <- 2*r*(floor(K*u[i]+1)-1) + r
          center_uj <- 2*r*(floor(K*u[j]+1)-1) + r
          
          p_base <- 0.5*(floor(K*u[i]+1)/(K+5) + floor(K*u[j]+1)/(K+5))
          suppressWarnings(p <- (sqrt((r)^2-((u[i]-center_ui)/1.2)^2-((u[j]-center_uj)/1.2)^2)/r)*p_base)
          p <- if (is.na(p)) 0 else p
        }
        p_matrix[i,j] <- p    
      }
    } 
  }

  if(pattern == "corner_low"){
    if(is.null(smoothness)){
      smoothness = 1 
      er_p = 0
    }else{
      er_p = 0.2
    }
    p_matrix = matrix(0,nrow=size,ncol=size)
    for(i in 1:size){
      for(j in 1:size){
        p1 <- cos(pi*sqrt(u[i]^2+u[j]^2)/sqrt(6))
        p2 <- cos(pi*sqrt((u[i]-1)^2+u[j]^2)/sqrt(6))
        p3 <- cos(pi*sqrt(u[i]^2+(u[j]-1)^2)/sqrt(6))
        p4 <- cos(pi*sqrt((u[i]-1)^2+(u[j]-1)^2)/sqrt(6))
        t=0.4
        #p1 = sign(p1) * pmax(0, abs(p1) - t)
        #p2 = sign(p2) * pmax(0, abs(p2) - t)
        #p3 = sign(p3) * pmax(0, abs(p3) - t)
        #p4 = sign(p4) * pmax(0, abs(p4) - t)
        p = p1+p2+p3+p4
        p_matrix[i,j] <-  1/(exp(-p)+1)   
      }
    }
  }
  
  
  if(pattern == "round_spike"){
    if(is.null(smoothness)){
      smoothness = 20
      er_p = 0
    }else{
      er_p = 0.2
    }
    p_matrix = matrix(0,nrow=size,ncol=size)
    for(i in 1:size){
      for(j in 1:size){
        p = (sin(smoothness*u[i])^2)*sin(smoothness*u[j])^2 + sin(u[j]+u[i])^2
        p_matrix[i,j] <-  p /2  
      }
    }
  }
  
  if(pattern == "peak_valley"){
    if(is.null(smoothness)){
      smoothness = 20
      er_p = 0
    }else{
      er_p = 0.2
    }
    p_matrix = matrix(0,nrow=size,ncol=size)
    for(i in 1:size){
      for(j in 1:size){
        p = (sin(smoothness*u[i])^2)+sin(smoothness*u[j])^2
        p_matrix[i,j] <-  p/2
      }
    }
  }
  
  
  if(pattern == "sinsin"){
    if(is.null(smoothness)){
      smoothness = 20
      er_p = 0
    }else{
      er_p = 0.2
    }
    p_matrix = matrix(0,nrow=size,ncol=size)
    for(i in 1:size){
      for(j in 1:size){
        if(u[i]>=u[j]){
          p = sin((smoothness*u[i]-smoothness/2)*sin(smoothness*u[j]-smoothness/2))
        }else{
          p = sin((smoothness*u[j]-smoothness/2)*sin(smoothness*u[i]-smoothness/2))
        }
        
        p_matrix[i,j] <-  (p+1)/2
      }
    }
  }
  
  
  if(pattern == "zhu2"){
    p_matrix = matrix(0,nrow=size,ncol=size)
    for(i in 1:size){
      for(j in 1:size){
        p <- sin(5*pi*(u[i]+u[j]-1)+1)/2 +0.5
        p_matrix[i,j] <- p    
      }
    }
  }
  
  if(pattern == "zhu3"){
    p_matrix = matrix(0,nrow=size,ncol=size)
    for(i in 1:size){
      for(j in 1:size){
        p <- (1+exp(15*(0.8*abs(u[i]-u[j]))^(4/5)-0.1))^(-1)
        p_matrix[i,j] <- p    
      }
    }
  }
  
  if(pattern == "zhu4"){
    p_matrix = matrix(0,nrow=size,ncol=size)
    for(i in 1:size){
      for(j in 1:size){
        p <- (u[i]^2+u[j]^2)/3*(cos(1/(u[i]^2+u[j]^2))) +0.15
        p_matrix[i,j] <- p    
      }
    }
  }
  
  if(pattern == "SBA1"){
    p_matrix = matrix(0,nrow=size,ncol=size)
    for(i in 1:size){
      for(j in 1:size){
        p <- 1/(1+exp(-50*(u[i]^2+u[j]^2)))
        p_matrix[i,j] <- p    
      }
    }
  }
  
  
#SBA2 is the same with SAS1
  if(pattern == "SBA2" | pattern == "SAS1"){
    p_matrix = matrix(0,nrow=size,ncol=size)
    for(i in 1:size){
      for(j in 1:size){
        p <- u[i]*u[j]
        p_matrix[i,j] <- p    
      }
    }
  }
  
  if(pattern == "SAS2"){
    p_matrix = matrix(0,nrow=size,ncol=size)
    for(i in 1:size){
      for(j in 1:size){
        p <- exp(-(u[i]^0.7+u[j]^0.7))
        p_matrix[i,j] <- p    
      }
    }
  }
  
  if(pattern == "SAS3"){
    p_matrix = matrix(0,nrow=size,ncol=size)
    for(i in 1:size){
      for(j in 1:size){
        p <- 0.25*(u[i]^2 + u[j]^2+ u[i]^0.5+ u[j]^0.5)
        p_matrix[i,j] <- p    
      }
    }
  }
  
  if(pattern == "SAS4"){
    p_matrix = matrix(0,nrow=size,ncol=size)
    for(i in 1:size){
      for(j in 1:size){
        p <- 0.5*(u[i] + u[j])
        p_matrix[i,j] <- p    
      }
    }
  }
  if(pattern == "SAS5"){
    p_matrix = matrix(0,nrow=size,ncol=size)
    for(i in 1:size){
      for(j in 1:size){
        p <- 1/(1 + exp(-10*(u[i]^2 + u[j]^2)))
        p_matrix[i,j] <- p    
      }
    }
  }
  if(pattern == "SAS6"){
    if(is.null(smoothness)){
      smoothness = 1 
      er_p = 0
    }else{
      er_p = 0.2
    }
    p_matrix = matrix(0,nrow=size,ncol=size)
    for(i in 1:size){
      for(j in 1:size){
        p <- abs(u[i] - u[j])/smoothness+er_p
        p_matrix[i,j] <- p    
      }
    }
  }
  if(pattern == "SAS7"){
    p_matrix = matrix(0,nrow=size,ncol=size)
    for(i in 1:size){
      for(j in 1:size){
        p <- 1/(1+exp(-((max(u[i],u[j]))^2+(min(u[i],u[j]))^4)))  
        p_matrix[i,j] <- p    
      }
    }
  }
  if(pattern == "SAS8"){
    if(is.null(smoothness)){
      smoothness = 1 
      er_p = 0
    }else{
      er_p = 0.2
    }
    p_matrix = matrix(0,nrow=size,ncol=size)
    for(i in 1:size){
      for(j in 1:size){
        p <- exp(-(max(u[i],u[j]))^(0.75))/smoothness+er_p  
        p_matrix[i,j] <- p    
      }
    }
  }
  if(pattern == "SAS9"){
    p_matrix = matrix(0,nrow=size,ncol=size)
    for(i in 1:size){
      for(j in 1:size){
        p <- exp(-0.5*(min(u[i],u[j])+u[i]^0.5+u[j]^0.5))  
        p_matrix[i,j] <- p    
      }
    }
  }
  if(pattern == "SAS10"){
    p_matrix = matrix(0,nrow=size,ncol=size)
    for(i in 1:size){
      for(j in 1:size){
        p <- log(1+0.5*(max(u[i],u[j])))
        p_matrix[i,j] <- p    
      }
    }
  }

  if(pattern == "er"){
    p_matrix = matrix(0,nrow=size,ncol=size)
    if(is.null(smoothness)){
      p = 0.5 
    }else{
      p <- smoothness
    }
    for(i in 1:size){
      for(j in 1:size){
        p_matrix[i,j] <- p    
      }
    }
  }
  
  if(pattern == "almost_emp"){
    p_matrix = matrix(0,nrow=size,ncol=size)
    for(i in 1:size){
      for(j in 1:size){
        if(floor(2*u[i]+1) == floor(2*u[j]+1)){
          p <- 0.99
        }else{
          p <- 0.01
        }
        p_matrix[i,j] <- p    
      }
    }
  }
  
  if(pattern == "emp"){
    p_matrix = matrix(0,nrow=size,ncol=size)
    for(i in 1:size){
      for(j in 1:size){
        if(floor(2*u[i]+1) == floor(2*u[j]+1)){
          p <- 1
        }else{
          p <- 0
        }
        p_matrix[i,j] <- p    
      }
    }
  }
  
  if(pattern == "extreme_emp"){
    p_matrix = matrix(0,nrow=size,ncol=size)
    for(i in 1:size){
      for(j in 1:size){
        p_matrix[i,j] <- (abs(u[i]-u[j])*1000)%%2>1    
      }
    }
  }
  
  #smoothness is the value to separate two block in this graph
  if(pattern == "block_emp"){
    smoothness <- ifelse(is.null(smoothness),0.5,smoothness)
    p_matrix = matrix(0,nrow=size,ncol=size)
    for(i in 1:size){
      for(j in 1:size){
        p_matrix[i,j] <-  (u[i]<smoothness) == (u[j]<smoothness)
      }
    }
  }
  
  #smoothness is the value to separate two block in this graph
  if(pattern == "block_emp_reverse"){
    smoothness <- ifelse(is.null(smoothness),0.5,smoothness)
    p_matrix = matrix(0,nrow=size,ncol=size)
    for(i in 1:size){
      for(j in 1:size){
        p_matrix[i,j] <-  !(u[i]<smoothness) == (u[j]<smoothness)
      }
    }
  }
  if(pattern == "50block_emp"){
    u_block = floor((u*50+1)%%50)
    p_matrix = matrix(0,nrow=size,ncol=size)
    for(i in 1:size){
      for(j in 1:size){
        p_matrix[i,j] <-  (u_block[i] == u_block[j]) 
      }
    }
  }
  
  if(pattern == "50block_emp_reverse"){
    u_block = floor((u*50+1)%%50)
    p_matrix = matrix(0,nrow=size,ncol=size)
    for(i in 1:size){
      for(j in 1:size){
        p_matrix[i,j] <-  !(u_block[i] == u_block[j]) 
      }
    }
  }
  
  if(return_u==TRUE){
    return(list(p_matrix=p_matrix,u=u))
  }else{
    return(p_matrix)
  }
  
}

#generate p matrix from preset graphon function
graphon_u_to_p_smoothness <- function(size,smoothness,return_u = FALSE,sampling_on_u=TRUE,u_input=NULL){
  if(!is.null(u_input)){u = u_input}else{
    if(sampling_on_u){
      u = runif(n = size)
    }else{
      u = seq(from = 0, to = 1,length.out=size)
    }
  }
  u = sort(u)
  
  #graphon function
  K = floor(log(size))
  p_matrix = matrix(0,nrow=size,ncol=size)
  for(i in 1:size){
    for(j in 1:size){
      p <- cos(pi*(u[i]-u[j]))/smoothness +0.5
      p_matrix[i,j] <- p    
    }
  }

  if(return_u==TRUE){
    return(list(p_matrix=p_matrix,u=u))
  }else{
    return(p_matrix)
  }

}


graphon_u_to_p_smoothness2 <- function(size,smoothness,return_u = FALSE,sampling_on_u=TRUE,u_input=NULL,p_base=0.5){
  
  if(!is.null(u_input)){u = u_input}else{
    if(sampling_on_u){
      u = runif(n = size)
    }else{
      u = seq(from = 0, to = 1,length.out=size)
    }
  }
  u = sort(u)
  
  #graphon function
  p_matrix = matrix(p_base,nrow=size,ncol=size)
  if(smoothness!=0){
    for(i in 1:size){
      for(j in 1:size){
        #p <- sin((u[i]^2+u[j]^2)*pi)/smoothness +p_base
        p <- sin(((u[i]^2+u[j]^2)^(1/4))*smoothness)/2 +0.5
        p_matrix[i,j] <- p    
      }
    }
  }
  
  
  if(return_u==TRUE){
    return(list(p_matrix=p_matrix,u=u))
  }else{
    return(p_matrix)
  }
  
}

graphon_u_to_p_smoothness3 <- function(size,smoothness,return_u = FALSE,sampling_on_u=TRUE,u_input=NULL,p_base=0.5){
  
  if(!is.null(u_input)){u = u_input}else{
    if(sampling_on_u){
      u = runif(n = size)
    }else{
      u = seq(from = 0, to = 1,length.out=size)
    }
  }
  u = sort(u)
  
  #graphon function
  p_matrix = matrix(p_base,nrow=size,ncol=size)
  if(smoothness!=0){
    for(i in 1:size){
      for(j in 1:size){
        p <- sin((u[i]^2+u[j]^2)*pi)/smoothness +p_base
        #p <- sin(((u[i]^2+u[j]^2)^(1/4))*smoothness)/2 +0.5
        p_matrix[i,j] <- p    
      }
    }
  }
  
  
  if(return_u==TRUE){
    return(list(p_matrix=p_matrix,u=u))
  }else{
    return(p_matrix)
  }
  
}

#this function generate graph adj matrix from p matrix
gmodel.P <- function(P,rep=1,noloop=TRUE,symmetric.out=FALSE){
  ## Check P
  cond1 = ((all(P>=0))&&(all(P<=1)))
  cond2 = (nrow(P)==ncol(P))
  if (!(cond1&&cond2)){
    stop("* gmodel.P : P is not a valid probability matrix.")
  }

  ## Parameter
  n = nrow(P)

  ## Rep 1 case
  if (rep==1){
    tmpmat = matrix(runif(n^2),nrow=n)
    if (symmetric.out){
      tmpmat[lower.tri(tmpmat)] <- t(tmpmat)[lower.tri(t(tmpmat))]
    }
    G = (tmpmat<P)*1
    if (noloop){
      diag(G) = 0
    }
  } else {
    G = list()
    for (i in 1:rep){
      tmpmat = matrix(runif(n^2),nrow=n)
      if (symmetric.out){
        tmpmat[lower.tri(tmpmat)] <- t(tmpmat)[lower.tri(t(tmpmat))]
      }
      tmpG = 1*(tmpmat<P)
      if (noloop){
        diag(tmpG) = 0
      }
      G[[i]] = tmpG
    }
  }
  ## return output
  return(G)
}

#this function generate graph adj matrix from p matrix
#old version
gmodel.P.old <- function(P,rep=1,noloop=TRUE,symmetric.out=FALSE){
  ## Check P
  cond1 = ((all(P>=0))&&(all(P<=1)))
  cond2 = (nrow(P)==ncol(P))
  if (!(cond1&&cond2)){
    stop("* gmodel.P : P is not a valid probability matrix.")
  }

  ## Parameter
  n = nrow(P)

  ## Rep 1 case
  if (rep==1){
    tmpmat = matrix(runif(n^2),nrow=n)
    G = (tmpmat<P)*1
    if (noloop){
      diag(G) = 0
    }
  } else {
    G = list()
    for (i in 1:rep){
      tmpmat = matrix(runif(n^2),nrow=n)
      tmpG = 1*(tmpmat<P)
      if (noloop){
        diag(tmpG) = 0
      }
      G[[i]] = tmpG
    }
  }

  ## Symmetric
  if ((symmetric.out)&&(!isSymmetric(P))){
    stop("* gmodel.P : 'symmetric' option is only valid if where probability P matrix is symmetric.")
  }
  if (symmetric.out){
    if (isSymmetric(P)){
      if (rep==1){
        tmpG = (G+t(G))/2
        G = (tmpG>0)*1
      } else {
        for (i in 1:rep){
          tmptmpG = G[[i]]
          tmpG = tmptmpG+t(tmptmpG)
          G[[i]] = (tmpG>0)*1
        }
      }
    }
  }


  ## return output
  return(G)
}

#replace lower tri with upper tri
to_undirected <- function(matrix){
  matrix[lower.tri(matrix)] <- t(matrix)[lower.tri(t(matrix))]
  matrix
}


#this function estimate smoothness of a probablity matrix
#W is a probability matrix
smoothness <- function(W,version=1){
  
  if(version==1){
    #version 1
    #average l2 dist over all node pairs
    l2_dist_vec <- combn(1:nrow(W),2,function(x){
      node_pair <- W[x,]
      sqrt(sum(((node_pair[1,]-node_pair[2,])^2)))
    })
    res <- mean(l2_dist_vec)
  }else if(version==2){
    #version 2
    #sd of l2 dist over all node pairs
    l2_dist_vec <- combn(1:nrow(W),2,function(x){
      node_pair <- W[x,]
      sqrt(sum(((node_pair[1,]-node_pair[2,])^2)))
    })
    res <- sd(l2_dist_vec)
  }else if(version==3){
    #version 3
    #first remove diag of P matrix and then get sd of P i.i.d
    W_vector <-  W[upper.tri(W, diag = FALSE)]
    res <- sd(W_vector)
  }else if(version==4){
    #version 4
    #mean of hamming distance of P
    l2_dist_vec <- combn(1:nrow(W),2,function(x){
      node_pair <- W[x,]
      sum(abs(node_pair[1,]-node_pair[2,]))
    })
    res <- mean(l2_dist_vec)
  }else if(version==5){
    #version 5
    #se of hamming distance of P
    l2_dist_vec <- combn(1:nrow(W),2,function(x){
      node_pair <- W[x,]
      sum(abs(node_pair[1,]-node_pair[2,]))
    })
    res <- sd(l2_dist_vec)
  }else if(version==6){
    #version 6
    #mean of weighted hamming distance of P

    #weighted hamming
    #cap the emp_p
    N=ncol(W)
    emp_p <- apply(W,2,sum)/N
    p_1p <- sapply(emp_p,function(x){x*(1-x)})
    p_1p[p_1p==0] <- Inf
    weight_p <- (1/p_1p)
    

    #v1
    l2_dist_vec <- combn(1:nrow(W),2,function(x){
      node_pair <- W[x,]
      sum(abs(node_pair[1,]-node_pair[2,])*weight_p)
    })
    res <- mean(l2_dist_vec)

    #v2
    #D <- (1 - W) %*% t(W%*%diag(weight_p))
    #weighted.h2.matrix <- D + t(D)
    #res <- mean(weighted.h2.matrix[upper.tri(weighted.h2.matrix, diag = FALSE)])
  }else if(version==7){
    #version 7
    #se of weighted hamming distance of P

    #weighted hamming
    #cap the emp_p
    N=ncol(W)
    emp_p <- apply(W,2,sum)/N
    p_1p <- sapply(emp_p,function(x){x*(1-x)})
    p_1p[p_1p==0] <- Inf
    weight_p <- (1/p_1p)

    #v1
    l2_dist_vec <- combn(1:nrow(W),2,function(x){
      node_pair <- W[x,]
      sum(abs(node_pair[1,]-node_pair[2,])*weight_p)
    })
    res <- sd(l2_dist_vec)

    #v2
    # D <- (1 - W) %*% t(W%*%diag(weight_p))
    # weighted.h2.matrix <- D + t(D)
    # res <- sd(weighted.h2.matrix[upper.tri(weighted.h2.matrix, diag = FALSE)])
  }
  
  return(res)
}

#this function estimate smoothness of graph adjcency matrix
#W is a graph adjcency matrix
smoothness_adj <- function(W,version=1){
  
  if(version==1){
    #version 1
    #average l2 dist over all node pairs
    l2_dist_vec <- combn(1:nrow(W),2,function(x){
      node_pair <- W[x,]
      sqrt(sum(((node_pair[1,]-node_pair[2,])^2)))
    })
    res <- mean(l2_dist_vec)
  }else if(version==2){
    #version 2
    #sd of l2 dist over all node pairs
    l2_dist_vec <- combn(1:nrow(W),2,function(x){
      node_pair <- W[x,]
      sqrt(sum(((node_pair[1,]-node_pair[2,])^2)))
    })
    res <- sd(l2_dist_vec)
  }else if(version==3){
    #version 3
    #first remove diag of P matrix and then get sd of P i.i.d
    W_vector <-  W[upper.tri(W, diag = FALSE)]
    res <- sd(W_vector)
  }else if(version==4){
    #version 4
    #mean of hamming distance of P
    l2_dist_vec <- combn(1:nrow(W),2,function(x){
      node_pair <- W[x,]
      sum(abs(node_pair[1,]-node_pair[2,]))
    })
    res <- mean(l2_dist_vec)
  }else if(version==5){
    #version 5
    #se of hamming distance of P
    l2_dist_vec <- combn(1:nrow(W),2,function(x){
      node_pair <- W[x,]
      sum(abs(node_pair[1,]-node_pair[2,]))
    })
    res <- sd(l2_dist_vec)
  }else if(version==6){
    #version 6
    #mean of weighted hamming distance of P

    #weighted hamming
    #cap the emp_p
    N=ncol(W)
    emp_p <- apply(W,2,sum)/N
    p_1p <- sapply(emp_p,function(x){x*(1-x)})
    p_1p[p_1p==0] <- Inf
    weight_p <- (1/p_1p)

    #v1
    # l2_dist_vec <- combn(1:nrow(W),2,function(x){
    #   node_pair <- W[x,]
    #   sum(abs(node_pair[1,]-node_pair[2,])*weight_p)
    # })
    # res <- mean(l2_dist_vec)

    #v2
    D <- (1 - W) %*% t(W%*%diag(weight_p))
    weighted.h2.matrix <- D + t(D)
    res <- mean(weighted.h2.matrix[upper.tri(weighted.h2.matrix, diag = FALSE)])
  }else if(version==7){
    #version 7
    #se of weighted hamming distance of P

    #weighted hamming
    #cap the emp_p
    N=ncol(W)
    emp_p <- apply(W,2,sum)/N
    p_1p <- sapply(emp_p,function(x){x*(1-x)})
    p_1p[p_1p==0] <- Inf
    weight_p <- (1/p_1p)
    #v1
    # l2_dist_vec <- combn(1:nrow(W),2,function(x){
    #   node_pair <- W[x,]
    #   sum(abs(node_pair[1,]-node_pair[2,])*weight_p)
    # })
    # res <- sd(l2_dist_vec)

    #v2
    D <- (1 - W) %*% t(W%*%diag(weight_p))
    weighted.h2.matrix <- D + t(D)
    res <- sd(weighted.h2.matrix[upper.tri(weighted.h2.matrix, diag = FALSE)])
  }
  
  return(res)
}



plot_adj <- function(X,...){
  #hmcols<-colorRampPalette(c("red","white","blue"))(256)
  image(t(apply(X, 2, rev)),breaks=seq(0, 1, length.out=13),...)
}





#estimate connecting probability by kmean
cp_estimate <- function(adj.matrix,d=5){
  X_svd = eigen(adj.matrix)
  kmeans_svd <- kmeans(x=(X_svd$vectors)[,1:d], centers=d,iter.max = 100)
  #X_svd = svd(adj.matrix)
  #kmeans_svd <- kmeans(x=(X_svd$u)[,1:5], centers=5,iter.max = 100)
  cp.matrix <- matrix(0,d,d)
  for(i in 1:d){
    for(j in 1:d){
      A<-adj.matrix[kmeans_svd$cluster==i,kmeans_svd$cluster==j]
      if(i==j){
        cp.matrix[i,j]<- sum(apply(A,1,sum))/((dim(A)[1]-1)*dim(A)[2])
      }else{
        cp.matrix[i,j]<- sum(apply(A,1,sum))/(dim(A)[1]*dim(A)[2])
      }
    }
  }
  cp.matrix
}


#generate multiple sbm
generate_sbm_list <- function(M,n, pref.matrix, block.sizes){
  res <- list()
  for (i in 1:M){
    g <- sample_sbm(n=n, pref.matrix=pref.matrix, block.sizes=block.sizes)
    res <- c(res, list(g))
  }
  res
}

#generate multiple pa
generate_pa_list <- function(M,n, power, m){
  res <- list()
  for (i in 1:M){
	g <- sample_pa(n, power = power, m=m, directed = directed)
    res <- c(res, list(g))
  }
  res
}

#generate multiple er
generate_er_list <- function(M,n, p){
  res <- list()
  for (i in 1:M){
	g <- erdos.renyi.game(n, p, type = "gnp", directed = FALSE)
    res <- c(res, list(g))
  }
  res
}

#generate multiple sm
generate_sm_list <- function(M,n,nei, p){
  res <- list()
  for (i in 1:M){
	g <- sample_smallworld(1,n,nei,p)
    res <- c(res, list(g))
  }
  res
}


#adj list to graph list
adj_to_graph <- function(adj_list){
  graph_list <- lapply(adj_list,graph_from_adjacency_matrix,mode="undirected")
  graph_list
}



#######################graph attributes##########################

#1 Order
g.order <- function(g,graph_list = FALSE){
  if(graph_list==TRUE){
    lapply(g,gorder)
  }
  else{
    gorder(g)
  }
}

#2 Size
g.size <- function(g,graph_list = FALSE){
  if(graph_list==TRUE){
    lapply(g,gsize)
  }
  else{
    gsize(g)
  }
}

#3 Number of connected components
g.ncc <- function(g,graph_list = FALSE){
  if(graph_list==TRUE){
    lapply(g,FUN = function(x){
      components(x)$no
    })
  }
  else{
    components(g)$no
  }
}

#4 Circuit rank
#to-do
#not implemented yet

#5 Diameter
g.diameter <- function(g,graph_list = FALSE){
  if(graph_list==TRUE){
    lapply(g,FUN = function(x){
      diameter(x)
    })
  }
  else{
    diameter(g)
  }
}

#6 Vertex connectivity, or group cohesion
g.vconnectivity <- function(g,graph_list = FALSE){
  if(graph_list==TRUE){
    lapply(g,FUN = function(x){
      vertex_connectivity(x)
    })
  }
  else{
    vertex_connectivity(g)
  }
}

#7 Edge connectivity
g.econnectivity <- function(g,graph_list = FALSE){
  if(graph_list==TRUE){
    lapply(g,FUN = function(x){
      edge_connectivity(x)
    })
  }
  else{
    edge_connectivity(g)
  }
}

#11 Independence number
#!!! takes forever for large graph
g.ivs <- function(g,graph_list = FALSE){
  if(graph_list==TRUE){
    lapply(g,FUN = function(x){
      ivs(x)
    })
  }
  else{
    ivs(g)
  }
}

#12 cliques number
g.cliques <- function(g,graph_list = FALSE){
  if(graph_list==TRUE){
    lapply(g,FUN = function(x){
      clique_num(x)
    })
  }
  else{
    clique_num (g)
  }
}

#17 Wiener index
#to do

#20 Clustering coefficient or transitivity
g.transitivity <- function(g,graph_list = FALSE){
  if(graph_list==TRUE){
    lapply(g,FUN = function(x){
      transitivity(x)
    })
  }
  else{
    transitivity (g)
  }
}

#21 betweenness centrality
g.betweenness <- function(g,graph_list = FALSE){
  if(graph_list==TRUE){
    lapply(g,FUN = function(x){
      betweenness(x)
    })
  }
  else{
    betweenness (g)
  }
}

#23 Cheeger constant(bottleneckness)
# to do


#24 Estrada index (protein folding)
#subgraph_centrality(g)

#25 Strength of a graph
#to do

#26 degree ditribution
g.degree <- function(g,which = list(pos="LM", howmany=5),graph_list = FALSE){
  if(graph_list==TRUE){
    lapply(g,FUN = function(x){
      degree(x)
    })
  }
  else{
    degree(g)
  }
}


#27 eigen
g.eigens <- function(g,which = list(pos="LM", howmany=5),graph_list = FALSE){
  if(graph_list==TRUE){
    lapply(g,FUN = function(x){
      spectrum(x,which=which,options=list(maxiter=1000000000))[c("values", "vectors")]
    })
  }
  else{
    spectrum(g,which=which,options=list(maxiter=1000000000))[c("values", "vectors")]
  }
}





## test
#g.degree(g)
#g.degree(boot_graph,graph_list = TRUE)

#plot
#par(mfrow=c(1,2))
#plot(g)
#plot(boot_graph[[3]])



