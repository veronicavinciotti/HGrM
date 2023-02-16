da=sim_HGrM()

########
#Analysis
set.seed(123)
iter<-1000 #outer iterations for LSM
bd.iter<-20 #inner iterations for bdgraph
#simulating one covariate
X<-runif(n.edge,-0.5,0.5)
X<-as.matrix(X)


initial.graphs<-matrix(nrow=n.edge,ncol=B)
for(i in 1:B){
    g<-huge.select(huge(as.matrix(data[[i]]),method="glasso"),criterion="stars")$refit
    initial.graphs[,i]<-g[lower.tri(g)]
}

res<-lsm.bd(data = da$data,
            D=2,
            Z=X, 
            bd.iter=bd.iter, 
            iter=iter,
            method="ggm",
            initial.graphs = initial.graphs)

#save.image("LSM-simresults-p100.RData")
