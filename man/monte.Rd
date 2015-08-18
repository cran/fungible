\name{monte}
\alias{monte}
\title{Simulate Clustered Data with User-Defined Properties}
\description{
Function for simulating clustered data with user defined characteristics such as: within cluster 
indicator correlations, within cluster indicator skewness values, within cluster indicator kurtosis 
values, and cluster separations as indexed by each variable (indicator validities).
}
\usage{
monte(seed, nvar, nclus, clus.size, eta2, cor.list, random.cor,
               skew.list, kurt.list, secor, compactness, sortMeans)
}
\arguments{
\item{seed}{Required: An integer to be used as the random number seed.}
\item{nvar}{Required: Number of variables to simulate.}
\item{nclus}{Required: Number of clusters to simulate. \emph{Note} that number of clusters must be equal to or greater than 2.}
\item{clus.size}{Required: Number of objects in each cluster.}
\item{eta2}{Required: A vector of indicator validities that range from 0 to 1. Higher numbers produce clusters with greater separation on that indicator.}
\item{cor.list}{Optional: A list of correlation matrices. There should be one correlation matrix 
 for each cluster.  The first correlation matrix will represent the indicator correlations 
 within cluster 1.  The second correlation matrix will represent the indicator correlations 
 for cluster 2. Etc.}
\item{random.cor}{Optional: Set to TRUE to generate a common within cluster correlation matrix.}
\item{skew.list}{Optional: A list of within cluster indicator skewness values.}
\item{kurt.list}{Optional: A list of within cluster indicator kurtosis values.}
\item{secor}{Optional: If 'random.cor = TRUE' then 'secor' determines the standard error of the simulated within group correlation   matrices.}
\item{compactness}{Optional: A vector of cluster compactness parameters. The meaning of this option is explained Waller et al. (1999). Basically, 'compactness' allows users some control over cluster overlap without changing indicator validities. See the example below for an illustration.}
\item{sortMeans}{Optional: A logical that determines whether the latent means will be sorted by taxon. Default = TRUE}
}
\value{
\item{data}{The simulated data. The 1st column of 'data' denotes cluster membership.}                                     
\item{lmn}{The cluster indicator means.}
\item{fl}{The factor loading matrix as described in Waller, et al. 1999.}
\item{fs}{The unique values of the linearized factor scores.}
\item{call}{The call.}
\item{nclus}{Number of clusters.}
\item{nvar}{Number of variables.}
\item{cor.list}{The input within cluster correlation matrices.}
\item{skew.list}{The input within cluster indicator skewness values.}
\item{kurt.list}{The input within cluster indicator kurtosis values.}
\item{clus.size}{The number of observations in each cluster.}
\item{eta2}{Vector of indicator validities.}
\item{seed}{The random number seed.}
}
\references{
Fleishman, A. I (1978). A method for simulating non-normal distributions. \emph{Psychometrika, 43}, 521-532.

Vale, D. C., & Maurelli, V. A. (1983). Simulating multivariate nonnormal distributions. \emph{Psychometrika, 48}, 465-471.

Waller, N. G., Underhill, J. M., & Kaiser, H. A. (1999).  A method for generating simulated plasmodes and artificial test clusters with user-defined shape, size, and orientation. \emph{Multivariate Behavioral Research, 34}, 123--142.
} 
\author{Niels Waller}
\examples{
## Example 1
## Simulating Fisher's Iris data
# The original data were reported in: 
# Fisher, R. A. (1936) The use of multiple measurements in taxonomic
#     problems. Annals of Eugenics, 7, Part II, 179-188.
#
# This example includes 3 clusters. Each cluster represents
# an Iris species: Setosa, Versicolor, and Virginica.
# On each species, four variables were measured: Sepal Length, 
# Sepal Width, Petal Length, and Petal Width.
#
# The within species (cluster) correlations of the flower
# indicators are as follows:
#
# Iris Type 1: 
#      [,1]  [,2]  [,3]  [,4]
# [1,] 1.000 0.743 0.267 0.178
# [2,] 0.743 1.000 0.278 0.233
# [3,] 0.267 0.278 1.000 0.332
# [4,] 0.178 0.233 0.332 1.000
#
# Iris Type 2
#      [,1]  [,2]  [,3]  [,4]
# [1,] 1.000 0.526 0.754 0.546
# [2,] 0.526 1.000 0.561 0.664
# [3,] 0.754 0.561 1.000 0.787
# [4,] 0.546 0.664 0.787 1.000
#
# Iris Type 3
#      [,1]  [,2]  [,3]  [,4]
# [1,] 1.000 0.457 0.864 0.281
# [2,] 0.457 1.000 0.401 0.538
# [3,] 0.864 0.401 1.000 0.322
# [4,] 0.281 0.538 0.322 1.000
#
# 'monte' expects a list of correlation matrices
# 

#create a list of within species correlations
data(iris)
cormat <- cm <- lapply(split(iris[,1:4], iris[,5]), cor)
 
# create a list of within species indicator 
# skewness and kurtosis
 sk.lst <- list(c(0.120,  0.041,  0.106,  1.254),                     
                c(0.105, -0.363, -0.607, -0.031),
                c(0.118,  0.366,  0.549, -0.129) )
              
              
 kt.lst <- list(c(-0.253, 0.955,  1.022,  1.719),
                c(-0.533,-0.366,  0.048, -0.410),
                c( 0.033, 0.706, -0.154, -0.602) )    


#Generate a new sample of iris data
my.iris <- monte(seed=123, nvar = 4, nclus = 3, cor.list = cormat, 
                clus.size = c(50, 50, 50),
                eta2=c(0.619, 0.401, 0.941, 0.929), 
                random.cor = FALSE,
                skew.list = sk.lst, 
                kurt.list = kt.lst, 
                secor = .3, compactness=c(1, 1, 1),
                sortMeans = TRUE)


summary(my.iris)
plot(my.iris)

# Now generate a new data set with the sample indicator validities 
# as before but with different cluster compactness values.

my.iris2<-monte(seed = 123, nvar = 4, nclus = 3, 
               cor.list = cormat, clus.size = c(50, 50, 50),
               eta2 = c(0.619, 0.401, 0.941, 0.929), random.cor = FALSE,
               skew.list = sk.lst ,kurt.list = kt.lst, 
               secor = .3,
               compactness=c(2, .5, .5), 
               sortMeans = TRUE)


summary(my.iris2)

# Notice that cluster 1 has been blow up whereas clusters 2 and 3 have been shrunk.
plot(my.iris2)


### Now compare your original results with the actual 
## Fisher iris data
library(lattice)
data(iris)
super.sym <- trellis.par.get("superpose.symbol")
splom(~iris[1:4], groups = Species, data = iris,
      #panel = panel.superpose,
      key = list(title = "Three Varieties of Iris",
                 columns = 3, 
                 points = list(pch = super.sym$pch[1:3],
                 col = super.sym$col[1:3]),
                 text = list(c("Setosa", "Versicolor", "Virginica"))))


############### EXAMPLE 2 ##################################

## Example 2
## Simulating data for Taxometric
## Monte Carlo Studies.
##
## In this four part example we will 
## generate two group mixtures 
## (Complement and Taxon groups)
## under four conditions.
##
## In all conditions 
## base rate (BR) = .20
## 3 indicators
## indicator validities = .50 
## (This means that 50 percent of the total
## variance is due to the mixture.)
##
##
## Condition 1:
## All variables have a slight degree
## of skewness (.10) and kurtosis (.10).
## Within group correlations = 0.00.
##
##
##
## Condition 2:
## In this conditon we generate data in which the 
## complement and taxon distributions differ in shape.
## In the complement group all indicators have 
## skewness values of 1.75 and kurtosis values of 3.75.
## In the taxon group all indicators have skewness values
## of .50 and kurtosis values of 0.
## As in the previous condition, all within group
## correlations (nuisance covariance) are 0.00.
##
##
## Conditon 3:
## In this condition we retain all previous 
## characteristics except that the within group
## indicator correlations now equal .80
## (they can differ between groups).
##
##
## Conditon 4:
## In this final condition we retain
## all previous data characteristics except that 
## the variances of the indicators in the complement 
## class are now 5 times the indicator variances
## in the taxon class (while maintaining indicator skewness, 
## kurtosis, correlations, etc.).
 

##----------------------------


library(lattice)


############################
##      Condition 1  
############################
in.nvar <- 3  ##Number of variables
in.nclus <-2  ##Number of taxa
in.seed <- 123                
BR <- .20     ## Base rate of higher taxon

## Within taxon indicator skew and kurtosis
in.skew.list <- list(c(.1, .1, .1),c(.1, .1, .1)) 
in.kurt.list <- list(c(.1, .1, .1),c(.1, .1, .1))          

## Indicator validities
in.eta2 <- c(.50, .50, .50)

## Groups sizes for Population
BigN <- 100000
in.clus.size <- c(BigN*(1-BR), BR * BigN) 
 
## Generate Population of scores with "monte"
sample.data <- monte(seed = in.seed, 
                nvar=in.nvar, 
                nclus = in.nclus, 
                clus.size = in.clus.size, 
                eta2 = in.eta2, 
                skew.list = in.skew.list, 
                kurt.list = in.kurt.list)
               
          
output <- summary(sample.data)

z <- data.frame(sample.data$data[sample(1:BigN, 600, replace=FALSE),])
z[,2:4] <- scale(z[,2:4])
names(z) <- c("id","v1","v2","v3")

#trellis.device()
trellis.par.set( col.whitebg() )
print(
 cloud(v3 ~ v1 * v2,
       groups = as.factor(id),data=z, 
       subpanel = panel.superpose,
       zlim=c(-4, 4),
       xlim=c(-4, 4),
       ylim=c(-4, 4),
       main="",
       screen = list(z = 20, x = -70)),
   position=c(.1, .5, .5, 1), more = TRUE)
               
                 

############################
##      Condition 2  
############################

## Within taxon indicator skew and kurtosis
in.skew.list <- list(c(1.75, 1.75, 1.75),c(.50, .50, .50)) 
in.kurt.list <- list(c(3.75, 3.75, 3.75),c(0, 0, 0))          

## Generate Population of scores with "monte"
sample.data <- monte(seed = in.seed, 
               nvar = in.nvar, 
               nclus = in.nclus, 
               clus.size = in.clus.size, 
               eta2 = in.eta2, 
               skew.list = in.skew.list, 
               kurt.list = in.kurt.list)
               
          
output <- summary(sample.data)

z <- data.frame(sample.data$data[sample(1:BigN, 600, replace=FALSE),])
z[,2:4] <- scale(z[, 2:4])
names(z) <-c("id", "v1","v2", "v3")

print(
 cloud(v3 ~ v1 * v2,
       groups = as.factor(id), data = z, 
       subpanel = panel.superpose,
       zlim = c(-4, 4),
       xlim = c(-4, 4),
       ylim = c(-4, 4),
       main="",
       screen = list(z = 20, x = -70)),
       position = c(.5, .5, 1, 1), more = TRUE)
               
                 
############################
##      Condition 3  
############################ 

## Set within group correlations to .80
cormat <- matrix(.80, 3, 3)
diag(cormat) <- rep(1, 3)
in.cor.list <- list(cormat, cormat)

## Generate Population of scores with "monte"
sample.data <- monte(seed = in.seed, 
               nvar = in.nvar, 
               nclus = in.nclus, 
               clus.size = in.clus.size, 
               eta2 = in.eta2, 
               skew.list = in.skew.list, 
               kurt.list = in.kurt.list,
               cor.list = in.cor.list)
               
output <- summary(sample.data)

z <- data.frame(sample.data$data[sample(1:BigN, 600, 
                replace = FALSE), ])
z[,2:4] <- scale(z[, 2:4])
names(z) <- c("id", "v1", "v2", "v3")

##trellis.device()
##trellis.par.set( col.whitebg() )
print(
  cloud(v3 ~ v1 * v2,
       groups = as.factor(id),data=z, 
       subpanel = panel.superpose,
       zlim = c(-4, 4),
       xlim = c(-4, 4),
       ylim = c(-4, 4),
       main="",
       screen = list(z = 20, x = -70)),
position = c(.1, .0, .5, .5), more = TRUE)
                                

############################
##      Condition 4  
############################

## Change compactness so that variance of
## complement indicators is 5 times
## greater than variance of taxon indicators
                     
 v <-  ( 2 * sqrt(5))/(1 + sqrt(5)) 
 in.compactness <- c(v, 2-v)   
 
## Generate Population of scores with "monte"
sample.data <- monte(seed = in.seed, 
               nvar = in.nvar, 
               nclus = in.nclus, 
               clus.size = in.clus.size, 
               eta2 = in.eta2, 
               skew.list = in.skew.list, 
               kurt.list = in.kurt.list,
               cor.list = in.cor.list,
               compactness = in.compactness)
               
output <- summary(sample.data)

z <- data.frame(sample.data$data[sample(1:BigN, 600, replace = FALSE), ])
z[, 2:4] <- scale(z[, 2:4])
names(z) <- c("id", "v1", "v2", "v3")
print(
  cloud(v3 ~ v1 * v2,
       groups = as.factor(id),data=z, 
       subpanel = panel.superpose,
       zlim = c(-4, 4),
       xlim = c(-4, 4),
       ylim = c(-4, 4),
       main="",
       screen = list(z = 20, x = -70)),
 position = c(.5, .0, 1, .5), more = TRUE)
}
\keyword{datagen}% at least one, from doc/KEYWORDS
