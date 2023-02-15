
library(tensorA)
library(Matrix)
# Load MASS library
library(MASS)

Gmcmc2<-function(G, Z=NULL, n.iter=1000,alpha=NULL,beta=NULL,cloc=NULL,n.burnin=500)
{
  B<-ncol(G)
  n.edge<-nrow(G)
  p<-(sqrt(1+8*n.edge)+1)/2
  m<-matrix(1:p,ncol=p,nrow = p)
  e1<-t(m)[lower.tri(m)]
  e2<-m[lower.tri(m)]
  if(is.null(cloc))
    cloc<-matrix(rnorm(B*2),ncol=2)
  if(is.null(alpha))
    alpha<-rnorm(B)
  
  dim.cond<-ncol(cloc)
  cloc.save<-array(dim = c(B,ncol(cloc),n.iter-n.burnin))
  alpha.save<-matrix(0,nrow=B,ncol=n.iter-n.burnin)
  # log.graph.prob.save<-rep(0,n.iter-n.burnin)
  
  if(!is.null(Z))
  {
    Z<-as.matrix(Z)
    if(is.null(beta))
      beta<-as.matrix(rep(0,ncol(Z)))
    else
      beta<-as.matrix(beta)
      beta.save<-matrix(0,nrow=ncol(Z),ncol=n.iter-n.burnin)
  }
  
  for (k in 1:n.iter){
    y<-as.vector(G)
    for (b in 1:B){
      # update the latent condition locations
      cloc_tensor <- as.tensor(cloc, c(nrow(cloc), ncol(cloc), n.edge))
      result <- rep(cloc_tensor, times = n.edge, along = 3)
      X = matrix(result,ncol=2) * rep(G[,b],B)
      #X = X  * matrix(G[,b], nrow=n.edge*dim.cond, ncol=1, byrow=TRUE)
      #X<-apply(cloc,2,rep,each=n.edge)*rep(G[,b],B)
      #X <- matrix(cloc[rep(1:n.edge, each=dim.cond), ], nrow=n.edge*dim.cond, byrow=TRUE) * matrix(G[,b], nrow=n.edge*dim.cond, ncol=1, byrow=TRUE)
      #X[(b-1)*n.edge+(1:n.edge),]<-matrix(apply(G,1,function(g,cloc,b){apply(cloc*g,2,sum)-cloc[b,]*g[b]},cloc=cloc,b=b),nrow=n.edge,ncol=dim.cond,byrow=T)
      

      X[(b-1)*n.edge+(1:n.edge),]<-matrix(unlist(mclapply(X = split(G, seq(nrow(G))), 
                                                          FUN = function(g,cloc,b){colSums(cloc * g)-cloc[b,]*g[b]},cloc=cloc,b=b,
                                                          mc.cores = 2),use.names = F),nrow=n.edge,ncol=dim.cond,byrow=T)
      
      hlp2<-NULL
      
      for (bb in 1:B){
        if (bb==b){
          hlp2<-c(hlp2,rep(0,n.edge))
        } else {
          #hlp3<-apply(G,1,function(g,cloc,bb,b){crossprod(apply(cloc*g,2,sum)-cloc[b,]*g[b]-cloc[bb,]*g[bb],cloc[bb,])},cloc=cloc,bb=bb,b=b)
          hlp3<-unlist(mclapply(X = split(G, seq(nrow(G))), 
                   FUN = function(g,cloc,b){crossprod(colSums(cloc * g)-cloc[b,]*g[b]-cloc[bb,]*g[bb],cloc[bb,])},cloc=cloc,b=b,
                   mc.cores = 2),use.names = F)
          hlp2<-c(hlp2,hlp3)
        }
      }
      offset<-hlp2+rep(alpha,each=n.edge)
      if(!is.null(Z))
        offset<-hlp2+rep(alpha,each=n.edge)+ rep(Z%*%beta,B)
        ####cloc[b,]<-blr(y,X,offset,theta = cloc[b,],theta_0 = rep(0,dim.cond))
    }
    dist.cond<-matrix(ncol=B,nrow=n.edge)
    for (b in 1:B){
      #updating condition-specific intercept
     # apply(G,1,function(g,cloc,b){crossprod(apply(cloc*g,2,sum)-cloc[b,]*g[b],cloc[b,])},cloc=cloc,b=b)
      dist.cond[,b]<-unlist(mclapply(X = split(G, seq(nrow(G))), FUN = function(g,cloc,b){crossprod(colSums(cloc * g)-cloc[b,]*g[b],cloc[b,])},cloc=cloc,b=b,
                      mc.cores = 2),use.names = F)
      offset<-dist.cond[,b]
      if(!is.null(Z))
        offset<-dist.cond[,b]+ Z%*%beta
      y<-G[,b]
      X<-as.matrix(rep(1,length(y)))
      #####alpha[b]<-blr(y,X,offset,theta = alpha[b],theta_0 = 0)
    }
    if(!is.null(Z))
    {
      y<-as.vector(G)
      X<-apply(Z,2,rep,B)
      offset<-c(dist.cond)+rep(alpha,each=n.edge)
      ####beta<-blr(y,X,offset,theta = beta,theta_0 = rep(0,length(beta)))
      beta<-t(beta)
    }
    
    if (k>n.burnin){
      cloc.save[,,k-n.burnin]<-cloc
      alpha.save[,k-n.burnin]<-alpha
      if(!is.null(Z))
        beta.save[,k-n.burnin]<-beta
      # calculate graph probability
      #    log.graph.prob<-0
      #    for (b in 1:B){
      #      for (i in 2:p){
      #        for (j in 1:(i-1)){
      #          ind<-e1==j & e2==i
      #          if(!is.null(Z))
      #            edge.prob<-pnorm(alpha[b]+dist.cond[ind,b]+Z[ind,]%*%beta)
      #          else
      #            edge.prob<-pnorm(alpha[b]+dist.cond[ind,b])
      #          if (G[ind,b]==0){
      #            log.graph.prob<-log.graph.prob+log(1-edge.prob)
      #          } else {
      #            log.graph.prob<-log.graph.prob+log(edge.prob)
      #          }
      #        }
      #      }
      #    }
      #    log.graph.prob.save[k-n.burnin]<-log.graph.prob
    }
  }
  if(is.null(Z))
    return(list(alpha=alpha.save,cloc=cloc.save,log.graph.prob=log.graph.prob.save))
  else
    #    return(list(alpha=alpha.save,beta=beta.save,cloc=cloc.save,log.graph.prob=log.graph.prob.save))
    return(list(alpha=alpha.save,beta=beta.save,cloc=cloc.save))
}






Gmcmc3<-function(G, Z=NULL, n.iter=1000,alpha=NULL,beta=NULL,cloc=NULL,n.burnin=500)
{
  B<-ncol(G)
  n.edge<-nrow(G)
  p<-(sqrt(1+8*n.edge)+1)/2
  m<-matrix(1:p,ncol=p,nrow = p)
  e1<-t(m)[lower.tri(m)]
  e2<-m[lower.tri(m)]
  if(is.null(cloc))
    cloc<-matrix(rnorm(B*2),ncol=2)
  if(is.null(alpha))
    alpha<-rnorm(B)
  
  dim.cond<-ncol(cloc)
  cloc.save<-array(dim = c(B,ncol(cloc),n.iter-n.burnin))
  alpha.save<-matrix(0,nrow=B,ncol=n.iter-n.burnin)
  # log.graph.prob.save<-rep(0,n.iter-n.burnin)
  
  if(!is.null(Z))
  {
    Z<-as.matrix(Z)
    if(is.null(beta))
      beta<-as.matrix(rep(0,ncol(Z)))
    else
      beta<-as.matrix(beta)
      beta.save<-matrix(0,nrow=ncol(Z),ncol=n.iter-n.burnin)
  }
  
  for (k in 1:n.iter){
    y<-as.vector(G)
    for (b in 1:B){
      # update the latent condition locations
      X<-apply(cloc,2,rep,each=n.edge)*rep(G[,b],B)
      X[(b-1)*n.edge+(1:n.edge),]<-matrix(apply(G,1,function(g,cloc,b){apply(cloc*g,2,sum)-cloc[b,]*g[b]},cloc=cloc,b=b),nrow=n.edge,ncol=dim.cond,byrow=T)
      hlp2<-NULL
      for (bb in 1:B){
        if (bb==b){
          hlp2<-c(hlp2,rep(0,n.edge))
        } else {
          hlp3<-apply(G,1,function(g,cloc,bb,b){crossprod(apply(cloc*g,2,sum)-cloc[b,]*g[b]-cloc[bb,]*g[bb],cloc[bb,])},cloc=cloc,bb=bb,b=b)
          hlp2<-c(hlp2,hlp3)
        }
      }
      offset<-hlp2+rep(alpha,each=n.edge)
      if(!is.null(Z))
        offset<-hlp2+rep(alpha,each=n.edge)+ rep(Z%*%beta,B)
        cloc[b,]<-blr(y,X,offset,theta = cloc[b,],theta_0 = rep(0,dim.cond))
    }
    dist.cond<-matrix(ncol=B,nrow=n.edge)
    for (b in 1:B){
      #updating condition-specific intercept
      dist.cond[,b]<-apply(G,1,function(g,cloc,b){crossprod(apply(cloc*g,2,sum)-cloc[b,]*g[b],cloc[b,])},cloc=cloc,b=b)
      offset<-dist.cond[,b]
      if(!is.null(Z))
        offset<-dist.cond[,b]+ Z%*%beta
      y<-G[,b]
      X<-as.matrix(rep(1,length(y)))
      ####alpha[b]<-blr(y,X,offset,theta = alpha[b],theta_0 = 0)
    }
    if(!is.null(Z))
    {
      y<-as.vector(G)
      X<-apply(Z,2,rep,B)
      offset<-c(dist.cond)+rep(alpha,each=n.edge)
      ####beta<-blr(y,X,offset,theta = beta,theta_0 = rep(0,length(beta)))
      beta<-t(beta)
    }
    
    if (k>n.burnin){
      cloc.save[,,k-n.burnin]<-cloc
      alpha.save[,k-n.burnin]<-alpha
      if(!is.null(Z))
        beta.save[,k-n.burnin]<-beta
    }
  }
  
  
  if(is.null(Z))
    return(list(alpha=alpha.save,cloc=cloc.save,log.graph.prob=log.graph.prob.save))
  else
    return(list(alpha=alpha.save,beta=beta.save,cloc=cloc.save))
}



############################################
####COMBINED WITH BDgraph ###################
#############################################
lsm.bd<-function(data,Z=NULL,initial.graphs=NULL, D=2, initial.cloc=NULL, initial.alpha=NULL, initial.beta=NULL, bd.iter=100,iter=1000,burnin=NULL, method = c("ggm","gcgm"), gcgm.dwpar=NULL)
{
  p<-ncol(data[[1]]) #number of nodes
  n.edge<-p*(p-1)/2 #number of edges
  B<-length(data) #number of conditions
  pi.edgpost<-matrix(0,n.edge,B)
  
  
  sample.graphs<-array(dim=c(n.edge,B,iter))
  if(is.null(initial.graphs))
  {
    for(i in 1:B)
    {
      g<-huge.select(huge(as.matrix(data[[i]]),method="glasso"),criterion="stars")$refit
      sample.graphs[,i,1]<-g[lower.tri(g)]
    }
  }
  else
    sample.graphs[,,1]<-initial.graphs
  
  sample.cloc<-array(dim = c(B,D,iter))
  sample.alpha<-matrix(0,B,iter)
  
  if(is.null(initial.cloc))
    sample.cloc[,,1]<-matrix(rnorm(B*D),ncol=D)
  else
    sample.cloc[,,1]<-initial.cloc
  
  if(!is.null(Z))
  {
    Z<-as.matrix(Z)
    sample.beta<-matrix(0,ncol(Z),iter)
    if(is.null(initial.beta))
    {
      y<-as.vector(sample.graphs[,,1])
      X<-apply(Z,2,rep,B)
      sample.beta[,1]<-coef(glm(y~X, family=binomial(link = "probit")))[-1]
    }
    if(!is.null(initial.beta))
      sample.beta[,1]<-initial.beta
  }
  log.graph.prob<-NULL
  
  if(is.null(burnin))
    burnin<-floor(0.75*iter)
  
  m<-matrix(1:p,ncol=p,nrow = p)
  e1<-t(m)[lower.tri(m)]
  e2<-m[lower.tri(m)]
  
  
  # Initialize K (precision matrix)
  K=vector(mode="list", length=B)
  for (i in 1:B){
    K[[i]]<-diag(p)
  }
  
  if (method=="gcgm"){
    discrete.data<-data
    #calculate truncated points
    tpoints<-vector("list",B)
    for(i in 1:B)
    {
      tpoints[[i]]<-vector("list",2)
      beta.dw<-gcgm.dwpar[[i]]$beta
      q<-gcgm.dwpar[[i]]$q
      pii<-matrix(rep(gcgm.dwpar[[i]]$pii,each=nrow(q)),nrow(q),ncol(q))
      pdw_lb = BDgraph::pdweibull( data[[i]] - 1, q = q, beta = beta.dw)
      pdw_ub = BDgraph::pdweibull( data [[i]], q = q, beta = beta.dw)
      tpoints[[i]][[1]]<-stats::qnorm( ( 1 - pii)*( data[[i]] != 0 ) + pii*pdw_lb)
      tpoints[[i]][[2]] <- stats::qnorm( (1 - pii) + pii*pdw_ub)
    }
  }
  
  for (k in 1: (iter-1))
  {
    # update data if the Gaussian Copula GM (gcgm) is selected
    if (method=="gcgm"){
      data<-sample.data(data,discrete.data, K, tpoints)
    }
    # update latent node and condition locations
    G<-sample.graphs[,,k]
    if(is.null(Z))
      G.loc<-Gmcmc(G,alpha=sample.alpha[,k],cloc=sample.cloc[,,k],n.iter=1,n.burnin = 0)
    else
      G.loc<-Gmcmc(G,Z=Z,alpha=sample.alpha[,k],beta=sample.beta[,k],cloc=sample.cloc[,,k],n.iter=1,n.burnin = 0)
    
    cloc<-G.loc$cloc[,,1]
    alpha<-G.loc$alpha
    beta<- G.loc$beta
    log.graph.prob<-c(log.graph.prob,G.loc$log.graph.prob)
    
    dist.cond<-matrix(ncol=B,nrow=n.edge)
    for (b in 1:B){
      #updating condition-specific intercept
      dist.cond[,b]<-apply(G,1,function(g,cloc,b){crossprod(apply(cloc*g,2,sum)-cloc[b,]*g[b],cloc[b,])},cloc=cloc,b=b)
    }
    Pi = matrix(ncol=B,nrow=n.edge)
    for (b in 1:B){
      for (i in 2:p){
        for (j in 1:(i-1)){
          ind<-e1==j & e2==i
          if(is.null(Z))
            Pi[ind,b]<-pnorm(alpha[b]+dist.cond[ind,b])
          else
            Pi[ind,b]<-pnorm(alpha[b]+dist.cond[ind,b]+Z[ind,]%*%beta)
        }
      }
    }
    K<-vector("list",B)
    for (j in 1:B)
    {
      pi.post<-Pi[,j]
      g.prior<-matrix(0,nrow=p,ncol=p)
      g.prior[lower.tri(g.prior)] <- pi.post
      g.prior<-g.prior+t(g.prior)
      
      g.start<-matrix(0,nrow=p,ncol=p)
      g.start[lower.tri(g.start)] <- sample.graphs[,j,k]
      g.start<-g.start+t(g.start)
      if(k < (iter-1))
        # update K
        res.bd<-BDgraph::bdgraph(data[[j]], iter = bd.iter, g.start=g.start,  g.prior=g.prior, save=FALSE, burnin=0)
      else
      {
        res.bd<-BDgraph::bdgraph(data[[j]], iter = bd.iter*100, g.start=g.start,  g.prior=g.prior, save=TRUE, burnin=0)
        pp<-plinks(res.bd)
        pi.edgpost[,j]<-t(pp)[lower.tri(pp)]
      }
      g<-res.bd$last_graph
      K[[j]]<-res.bd$last_K
      sample.graphs[,j,k+1]<-g[lower.tri(g)]
    }
    sample.cloc[,,k+1]<-cloc
    sample.alpha[,k+1]<-alpha
    if(!is.null(Z))
      sample.beta[,k+1]<-beta
    
  }
  
  sample.cloc<-sample.cloc[,,-(1:burnin)]
  sample.alpha<-sample.alpha[,-(1:burnin)]
  sample.graphs<-sample.graphs[,,-(1:burnin)]
  if(!is.null(Z))
    sample.beta<-sample.beta[,-(1:burnin),drop=FALSE]
  
  ##probit probabilities from latent space
  n.iter<-dim(sample.cloc)[3]
  pi.probit = array(dim=c(n.edge,B,n.iter))
  for (k in 1: n.iter){
    G<-sample.graphs[,,k]
    alpha<-sample.alpha[,k]
    if(!is.null(Z))
      beta<-sample.beta[,k]
    cloc<-sample.cloc[,,k]
    dist.cond<-matrix(ncol=B,nrow=n.edge)
    for (b in 1:B){
      #updating condition-specific intercept
      dist.cond[,b]<-apply(G,1,function(g,cloc,b){crossprod(apply(cloc*g,2,sum)-cloc[b,]*g[b],cloc[b,])},cloc=cloc,b=b)
    }
    Pi = matrix(ncol=B,nrow=n.edge)
    for (b in 1:B){
      for (i in 2:p){
        for (j in 1:(i-1)){
          ind<-e1==j & e2==i
          if(is.null(Z))
            pi.probit[ind,b,k]<-pnorm(alpha[b]+dist.cond[ind,b])
          if(!is.null(Z))
            pi.probit[ind,b,k]<-pnorm(alpha[b]+dist.cond[ind,b]+Z[ind,]%*%beta)
        }
      }
    }
  }
  if(is.null(Z))
    return(list(sample.alpha=sample.alpha,sample.cloc=sample.cloc,sample.graphs=sample.graphs,pi.edgpost=pi.edgpost,pi.probit=pi.probit))
  else
    return(list(sample.alpha=sample.alpha,sample.beta=sample.beta,sample.cloc=sample.cloc,sample.graphs=sample.graphs,pi.edgpost=pi.edgpost,pi.probit=pi.probit))
}



blr_fast <- function(y, X, offset = 0, theta, theta_0 = c(0, 0, 0), N_sim = 1) {
 
  
  # Dimensions of theta
  D <- ncol(X)
  
  # Number of observations
  n <- length(y)
  N1 <- sum(y)
  N0 <- n - N1
  
  # Conjugate prior on the coefficients theta ~ N(theta_0, Q_0)
  Q_0 <- diag(10, D)
  
  # Initialize parameters
  z <- rep(NA, n)
  
  # Matrix storing samples of the theta parameter
  theta_chain <- matrix(0, nrow = N_sim, ncol = D)
  
  # ---------------------------------
  # Gibbs sampling algorithm
  # ---------------------------------
  
  # Compute posterior variance of theta
  prec_0 <- ginv(Q_0)
  V <- ginv(prec_0 + t(X) %*% X)
  
  for (t in 1:N_sim) {
    # Update Mean of z
    mu_z <- X %*% theta + offset
    # Draw latent variable z from its full conditional: z | theta, y, X
    if (N0 > 0) {
      z[y == 0] <- rtruncnorm(N0, mean = mu_z[y == 0], sd = 1, a = -Inf, b = 0)
    }
    if (N1 > 0) {
      z[y == 1] <- rtruncnorm(N1, mean = mu_z[y == 1], sd = 1, a = 0, b = Inf)
    }
    
    # Compute posterior mean of theta
    M <- V %*% (prec_0 %*% theta_0 + t(X) %*% (z - offset))
    # Draw variable theta from its full conditional: theta | z, X
    theta <- rmvnorm(1, M, V)
    
    # Store the theta draws
    theta_chain[t, ] <- theta
  }
  
  return(list(theta = theta_chain, z = z))
}
