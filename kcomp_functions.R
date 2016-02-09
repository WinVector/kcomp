library(ggplot2)
library(GGally)
library(plot3D)


#----------------------------------------------------
#
#  Functions for within-sum-of-squares calculation.
#  I'm really just using this for the single-cluster
#  situation, since kmeans returns totalWSS already.
#
#-------------------------------------------------

#   Function to calculate squared distance
#   between two vectors.
sqr_edist <- function(x, y) {
  sum((x-y)^2)
}

#
#   Function to calculate the WSS for a single
#   cluster, which is represented as a matrix (one row
#   for every point).
#
wss.cluster <- function(clustermat) {
  #   Calculate the centroid of the cluster (the mean of all the points)
  c0 <- apply(clustermat, 2, FUN=mean)
  #   Calculate the squared difference of every point
  # in the cluster from the centroid, and sum all the distances.
  sum(apply(clustermat, 1, FUN=function(row){sqr_edist(row,c0)}))
}

#----------------------------------------------------
#
#  Functions for generating data sets with random clusterings
#
#---------------------------------------------------

#
# make an arbitrary d-dim cluster
#
make_arbitrary_cluster = function(d, nrows) {
  # first, generate from a gaussian sphere
  sphere = matrix(rnorm(d*nrows), nrow=nrows)

  # make the transformation matrix
  rmat = matrix(rnorm(d*d), nrow=d)

  cloud = sphere %*% rmat

  shift = 5*rnorm(d)
  for(i in seq_len(nrows)) {
    cloud[i, ] = cloud[i, ] + shift
  }
  cloud = as.data.frame(cloud)
  colnames(cloud) = paste("x", seq_len(d), sep="")
  cloud
}

#
# adds an additional character column to mark the cluster
#
make_cluster = function(dim, npts, i) {
  clust = make_arbitrary_cluster(dim, npts)
  clust$gp = as.character(i)
  clust
}

# make dim-dimensional data with k clusters
# each cluster has about 100 pts (size drawn from poisson distribution)
make_data = function(dim, k) {
  npts = rpois(k, 100)
  datax = do.call(rbind, lapply(seq_len(k),
                                function(i) {
                                  make_cluster(dim,npts[i],i)
                                }
  ))
  datax
}

#-----------------------------------------------------------
#
#  Functions for generating the "bootstrap" simulations
#
#---------------------------------------------------------

#
# From a data set (assumed a single gaussian cluster), generate
# a gaussian data set of the same size with the same center and principle components
# (so the same hyperellipsoid)
#
# input and output are matrices
#
generate_similar = function(dmatrix) {
  nr = nrow(dmatrix)
  nc = ncol(dmatrix)

  if(nc==1) {
    mu = mean(dmatrix)
    sigma = sd(dmatrix)
    nd3 = as.data.frame(rnorm(nr, mean=mu, sd=sigma))
    colnames(nd3) = colnames(dmatrix)
    return(as.matrix(nd3))
  }

  # get principle axes
  # princ has a rotation matrix rotation and scaling sdev
  # and a center center
  princ = prcomp(dmatrix)

  # generate data from a axis-aligned gaussian centered at the origin
  # with appropriate standard deviations
  nd1 = as.data.frame(vapply(seq_len(nc),
                             FUN = function(i) {rnorm(nr, sd=princ$sdev[i])},
                             numeric(nr)))
  # rotate it
  nd2 = t(princ$rotation %*% t(nd1))

  # center it
  nd3 = as.data.frame(vapply(seq_len(nc),
                             FUN=function(i) {princ$center[i] + nd2[,i]},
                             numeric(nr)))
  colnames(nd3) = colnames(dmatrix)
  as.matrix(nd3)

}

#
# given a dataset and an existing cluster assignment,
# generate other plausibly similar datasets of the same size
#
generate_bootstrap = function(dmatrix, cluster) {
  k = max(cluster)
  do.call(rbind, lapply(seq_len(k),
                        function(i) {generate_similar(dmatrix[cluster==i, , drop=FALSE])}))
}

#
# "null_clustering' is maybe not a good name. Find a
# k-means clustering for a dataset, assuming there are k clusters
#
null_clustering = function(dmatrix, k, nstart=10, iter.max=100) {
  clustering = kmeans(dmatrix, k, nstart=nstart, iter.max=iter.max)
  # return cluster assignments and total WSS
   list(cluster=clustering$cluster, totWSS = clustering$tot.withinss)
}

#-----------------------------------------------------------
#
#  Functions for finding the best k
#
#---------------------------------------------------------


#
# Compare the null hypothesis that the data has k (gaussian) clusters
# against the alternative that there are k+1 clusters by comparing
# a k+1 clustering of the true data to nboot k+1 clusterings of plausible
# simulated data with k clusters.
#
# Return:
# pval: the p-value of the k+1 clustering of the true data, and
# simdata: the WSS of all the simulations, and the clustering of the true data
#
test_k = function(dmatrix, k, nboot) {
  if(k==1){
    kcluster = list(cluster=rep(1, nrow(dmatrix)))
  } else {
      kcluster = null_clustering(dmatrix, k)
  }
  kplus1cluster = null_clustering(dmatrix, k+1)

  wssdist = vapply(seq_len(nboot),
                   FUN=function(i) {
                     kmeans(generate_bootstrap(dmatrix, kcluster$cluster), k+1)$tot.withinss
                     },
                   numeric(1))


# probability of kplus1cluster$WSS this good under null hypothesis is left tail of the distribution
  pval = sum(wssdist < kplus1cluster$totWSS)/nboot

  simdata = data.frame(wss=wssdist, observed=kplus1cluster$totWSS,
                       label=paste("WSS(",(k+1), ") vs. hypotheses k = ", k, "; pval = ", pval))

#   title = paste("Probability of WSS(",(k+1), ") under null hypothesis k =", k, ":", pval)
#   plot = ggplot(data.frame(wss=wssdist), aes(x=wss)) + geom_density() +
#            geom_vline(xintercept=kplus1cluster$totWSS, color="red") + ggtitle(title)

  list(pval=pval, simdata=simdata)
}

#
# Run test_k until we can't reject the null hypothesis to the p=0.05 level
# Note that data no longer has the cluster column on it.
#
# Return:
#  nclusters: the inferred number of clusters
#  pvals: the p-values of each experiment
#  simplot: a plot of all the experiments (distribution of WSS for simulations, vs WSS of true data)
#
find_clusters = function(data) {
  k = 1
  pval=0
  plist = NULL
  datalist=NULL
  while(pval < 0.05) {
    results = test_k(as.matrix(data), k, 100)
    pval=results$pval
    plist = c(plist, pval)
    datalist[[k]] = results$simdata
    k = k+1
  }

  simdata = do.call(rbind, datalist)
  simplot = ggplot(simdata, aes(x=wss)) + geom_density(trim=TRUE) +
    geom_vline(aes(xintercept=observed), color="red") + facet_wrap(~label, ncol=2, scales="free_x")

  list(nclusters = length(plist), pvals = plist, simplot=simplot)
}


#
# Unlike find_clusters, don't stop at the first failure to reject the null hypothesis,
# force the algo to run all the way to max_clusters and pick 1+ the last time you reject the null hypothesis
#
#
# Return:
#  nclusters: the inferred number of clusters
#  pvals: the p-values of each experiment
#  pvalplot: a plot of the pvals, versus the threshold of p=0.05
#  simplot: a plot of all the experiments (distribution of WSS for simulations, vs WSS of true data)
#
find_clusters_scan = function(data, max_clusters=8) {
  k = 1
  pval=numeric(max_clusters)
  plist = NULL
  datalist=NULL

  for(k in seq_len(max_clusters)) {
    results = test_k(as.matrix(data), k, 100)
    pval=results$pval
    plist = c(plist, pval)
    datalist[[k]] = results$simdata
  }

  pthresh = 0.05
  crossings = which(plist < pthresh)
  if(length(crossings)==0) {
    nclusters=1
    nstop=1
  } else {
    lastreject = max(which(plist < pthresh))
    nclusters = lastreject+1
    nstop = pmin(nclusters, max_clusters)
  }

  tmp = data.frame(k=1:max_clusters, p=plist)
  pvalplot = ggplot(tmp, aes(x=k, y=p)) + geom_point() +
    geom_line() + geom_hline(yintercept=pthresh, color="red", linetype=2)

  simdata = do.call(rbind, datalist[seq_len(nstop)])
  simplot = ggplot(simdata, aes(x=wss)) + geom_density(trim=TRUE) +
    geom_vline(aes(xintercept=observed), color="red") + facet_wrap(~label, ncol=2, scales="free_x")

  list(nclusters = nclusters, pvals = plist, simplot=simplot, pvalplot=pvalplot)
}

#
# prints out the p-values of the experiment, for the UI
#
clustering_report = function(pvals) {
  k = length(pvals)
  labels = paste("k=", 1:k,sep="")
  separators = rep("; ", k)
  separators[k] = ""
  paste(paste(labels," : ", pvals, separators, sep=""), collapse= " ")
}

#
# After inferring a probable k, recluster the data accordingly
# Note that data no longer has the cluster column on it
#
recluster = function(data, k) {
  fcluster = kmeans(as.matrix(data), k, nstart=100, iter.max=100)
  dataplus = cbind(data, gp=as.character(fcluster$cluster))
}


#-----------------------------------------------------------
#
#  Functions for plotting data
#
#---------------------------------------------------------

# assumes last column is the cluster assignment
scatterplot_data = function(datax) {
  D = ncol(datax)-1
  #palette="Dark2"
  palette="Set2"

  # for the case when the data has been read back in from a file
  if(is.numeric(datax$gp)) datax$gp = as.character(datax$gp)

  if(D==1) {
    return(ggplot(datax, aes(x=x1, color=gp, fill=gp)) +
             geom_density(trim=TRUE, alpha=0.5) + geom_rug() +
             scale_color_brewer(palette=palette) +
             scale_fill_brewer(palette=palette))
  }
  if(D==2) {
    return(ggplot(datax, aes(x=x1, y=x2, color=gp)) +
             geom_point() + coord_fixed() +
             scale_color_brewer(palette=palette))
  }

  # yikes. Kind of a hack.
  ps = ggpairs(datax, aes(color=gp), 1:D)
  scaleColor = scale_colour_brewer(palette=palette)
  scaleFill = scale_fill_brewer(palette=palette)
  for (row in seq_len(ps$nrow))
    for (col in seq_len(ps$ncol))
      ps[row, col] <- ps[row, col] + scaleColor + scaleFill
  ps

}

# assumes last column is the cluster assignment
plot3d_data = function(datax) {
  D = ncol(datax)-1
  if(D < 3) return(NULL)
  if(D==3) {
    project = as.matrix(datax[, 1:D])
  } else {
    pmatrix = as.matrix(datax[, 1:D])
    princ = prcomp(pmatrix)
    nComp = 3
    project = predict(princ, newdata=pmatrix)[,1:nComp]
  }
  ngp = length(unique(datax$gp))
  ncol = pmax(3, ngp)
  #palette="RdYlBu"
  palette="Set2"
  onecol = "#66c2a5"  # the first color in Set2
  # there is a bug in the colvar argument
  if(ngp > 1) {
    scatter3D(project[,1], project[,2], project[,3],
              colvar=as.numeric(datax$gp),
              bty="g",
              col=RColorBrewer::brewer.pal(ncol, palette))
  } else {
    scatter3D(project[,1], project[,2], project[,3],
              col=onecol,
              bty="g")
  }
}


