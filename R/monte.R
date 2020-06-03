#' Simulate Clustered Data with User-Defined Properties
#' 
#' Function for simulating clustered data with user defined characteristics
#' such as: within cluster indicator correlations, within cluster indicator
#' skewness values, within cluster indicator kurtosis values, and cluster
#' separations as indexed by each variable (indicator validities).
#' 
#' 
#' @param seed Required: An integer to be used as the random number seed.
#' @param nvar Required: Number of variables to simulate.
#' @param nclus Required: Number of clusters to simulate. \emph{Note} that
#' number of clusters must be equal to or greater than 2.
#' @param clus.size Required: Number of objects in each cluster.
#' @param eta2 Required: A vector of indicator validities that range from 0 to
#' 1. Higher numbers produce clusters with greater separation on that
#' indicator.
#' @param cor.list Optional: A list of correlation matrices. There should be
#' one correlation matrix for each cluster.  The first correlation matrix will
#' represent the indicator correlations within cluster 1.  The second
#' correlation matrix will represent the indicator correlations for cluster 2.
#' Etc.
#' @param random.cor Optional: Set to TRUE to generate a common within cluster
#' correlation matrix.
#' @param skew.list Optional: A list of within cluster indicator skewness
#' values.
#' @param kurt.list Optional: A list of within cluster indicator kurtosis
#' values.
#' @param secor Optional: If 'random.cor = TRUE' then 'secor' determines the
#' standard error of the simulated within group correlation matrices.
#' @param compactness Optional: A vector of cluster compactness parameters. The
#' meaning of this option is explained Waller et al. (1999). Basically,
#' 'compactness' allows users some control over cluster overlap without
#' changing indicator validities. See the example below for an illustration.
#' @param sortMeans Optional: A logical that determines whether the latent
#' means will be sorted by taxon. Default = TRUE
#' @return \item{data}{The simulated data. The 1st column of 'data' denotes
#' cluster membership.} \item{lmn}{The cluster indicator means.} \item{fl}{The
#' factor loading matrix as described in Waller, et al. 1999.} \item{fs}{The
#' unique values of the linearized factor scores.} \item{call}{The call.}
#' \item{nclus}{Number of clusters.} \item{nvar}{Number of variables.}
#' \item{cor.list}{The input within cluster correlation matrices.}
#' \item{skew.list}{The input within cluster indicator skewness values.}
#' \item{kurt.list}{The input within cluster indicator kurtosis values.}
#' \item{clus.size}{The number of observations in each cluster.}
#' \item{eta2}{Vector of indicator validities.} \item{seed}{The random number
#' seed.}
#' @author Niels Waller
#' @references Fleishman, A. I (1978). A method for simulating non-normal
#' distributions. \emph{Psychometrika, 43}, 521-532.
#' 
#' Olvera Astivia, O. L. & Zumbo, B. D. (2018). On the solution 
#' multiplicity of the Fleishman method and its impact in 
#' simulation studies. \emph{British Journal of Mathematical and Statistical Psychology, 71}
#' (3), 437-458. 
#' 
#' Vale, D. C., & Maurelli, V. A. (1983). Simulating multivariate nonnormal
#' distributions. \emph{Psychometrika, 48}, 465-471.
#' 
#' Waller, N. G., Underhill, J. M., & Kaiser, H. A. (1999).  A method for
#' generating simulated plasmodes and artificial test clusters with
#' user-defined shape, size, and orientation. \emph{Multivariate Behavioral
#' Research, 34}, 123-142.
#' @keywords datagen
#' @import lattice
#' @export
#' @examples
#' 
#' ## Example 1
#' ## Simulating Fisher's Iris data
#' # The original data were reported in: 
#' # Fisher, R. A. (1936) The use of multiple measurements in taxonomic
#' #     problems. Annals of Eugenics, 7, Part II, 179-188.
#' #
#' # This example includes 3 clusters. Each cluster represents
#' # an Iris species: Setosa, Versicolor, and Virginica.
#' # On each species, four variables were measured: Sepal Length, 
#' # Sepal Width, Petal Length, and Petal Width.
#' #
#' # The within species (cluster) correlations of the flower
#' # indicators are as follows:
#' #
#' # Iris Type 1: 
#' #      [,1]  [,2]  [,3]  [,4]
#' # [1,] 1.000 0.743 0.267 0.178
#' # [2,] 0.743 1.000 0.278 0.233
#' # [3,] 0.267 0.278 1.000 0.332
#' # [4,] 0.178 0.233 0.332 1.000
#' #
#' # Iris Type 2
#' #      [,1]  [,2]  [,3]  [,4]
#' # [1,] 1.000 0.526 0.754 0.546
#' # [2,] 0.526 1.000 0.561 0.664
#' # [3,] 0.754 0.561 1.000 0.787
#' # [4,] 0.546 0.664 0.787 1.000
#' #
#' # Iris Type 3
#' #      [,1]  [,2]  [,3]  [,4]
#' # [1,] 1.000 0.457 0.864 0.281
#' # [2,] 0.457 1.000 0.401 0.538
#' # [3,] 0.864 0.401 1.000 0.322
#' # [4,] 0.281 0.538 0.322 1.000
#' #
#' # 'monte' expects a list of correlation matrices
#' # 
#' 
#' #create a list of within species correlations
#' data(iris)
#' cormat <- cm <- lapply(split(iris[,1:4], iris[,5]), cor)
#'  
#' # create a list of within species indicator 
#' # skewness and kurtosis
#'  sk.lst <- list(c(0.120,  0.041,  0.106,  1.254),                     
#'                 c(0.105, -0.363, -0.607, -0.031),
#'                 c(0.118,  0.366,  0.549, -0.129) )
#'               
#'               
#'  kt.lst <- list(c(-0.253, 0.955,  1.022,  1.719),
#'                 c(-0.533,-0.366,  0.048, -0.410),
#'                 c( 0.033, 0.706, -0.154, -0.602) )    
#' 
#' 
#' #Generate a new sample of iris data
#' my.iris <- monte(seed=123, nvar = 4, nclus = 3, cor.list = cormat, 
#'                 clus.size = c(50, 50, 50),
#'                 eta2=c(0.619, 0.401, 0.941, 0.929), 
#'                 random.cor = FALSE,
#'                 skew.list = sk.lst, 
#'                 kurt.list = kt.lst, 
#'                 secor = .3, compactness=c(1, 1, 1),
#'                 sortMeans = TRUE)
#' 
#' 
#' summary(my.iris)
#' plot(my.iris)
#' 
#' # Now generate a new data set with the sample indicator validities 
#' # as before but with different cluster compactness values.
#' 
#' my.iris2<-monte(seed = 123, nvar = 4, nclus = 3, 
#'                cor.list = cormat, clus.size = c(50, 50, 50),
#'                eta2 = c(0.619, 0.401, 0.941, 0.929), random.cor = FALSE,
#'                skew.list = sk.lst ,kurt.list = kt.lst, 
#'                secor = .3,
#'                compactness=c(2, .5, .5), 
#'                sortMeans = TRUE)
#' 
#' 
#' summary(my.iris2)
#' 
#' # Notice that cluster 1 has been blow up whereas clusters 2 and 3 have been shrunk.
#' plot(my.iris2)
#' 
#' 
#' ### Now compare your original results with the actual 
#' ## Fisher iris data
#' library(lattice)
#' data(iris)
#' super.sym <- trellis.par.get("superpose.symbol")
#' splom(~iris[1:4], groups = Species, data = iris,
#'       #panel = panel.superpose,
#'       key = list(title = "Three Varieties of Iris",
#'                  columns = 3, 
#'                  points = list(pch = super.sym$pch[1:3],
#'                  col = super.sym$col[1:3]),
#'                  text = list(c("Setosa", "Versicolor", "Virginica"))))
#' 
#' 
#' ############### EXAMPLE 2 ##################################
#' 
#' ## Example 2
#' ## Simulating data for Taxometric
#' ## Monte Carlo Studies.
#' ##
#' ## In this four part example we will 
#' ## generate two group mixtures 
#' ## (Complement and Taxon groups)
#' ## under four conditions.
#' ##
#' ## In all conditions 
#' ## base rate (BR) = .20
#' ## 3 indicators
#' ## indicator validities = .50 
#' ## (This means that 50 percent of the total
#' ## variance is due to the mixture.)
#' ##
#' ##
#' ## Condition 1:
#' ## All variables have a slight degree
#' ## of skewness (.10) and kurtosis (.10).
#' ## Within group correlations = 0.00.
#' ##
#' ##
#' ##
#' ## Condition 2:
#' ## In this conditon we generate data in which the 
#' ## complement and taxon distributions differ in shape.
#' ## In the complement group all indicators have 
#' ## skewness values of 1.75 and kurtosis values of 3.75.
#' ## In the taxon group all indicators have skewness values
#' ## of .50 and kurtosis values of 0.
#' ## As in the previous condition, all within group
#' ## correlations (nuisance covariance) are 0.00.
#' ##
#' ##
#' ## Conditon 3:
#' ## In this condition we retain all previous 
#' ## characteristics except that the within group
#' ## indicator correlations now equal .80
#' ## (they can differ between groups).
#' ##
#' ##
#' ## Conditon 4:
#' ## In this final condition we retain
#' ## all previous data characteristics except that 
#' ## the variances of the indicators in the complement 
#' ## class are now 5 times the indicator variances
#' ## in the taxon class (while maintaining indicator skewness, 
#' ## kurtosis, correlations, etc.).
#'  
#' 
#' ##----------------------------
#' 
#' 
#' library(lattice)
#' 
#' 
#' ############################
#' ##      Condition 1  
#' ############################
#' in.nvar <- 3  ##Number of variables
#' in.nclus <-2  ##Number of taxa
#' in.seed <- 123                
#' BR <- .20     ## Base rate of higher taxon
#' 
#' ## Within taxon indicator skew and kurtosis
#' in.skew.list <- list(c(.1, .1, .1),c(.1, .1, .1)) 
#' in.kurt.list <- list(c(.1, .1, .1),c(.1, .1, .1))          
#' 
#' ## Indicator validities
#' in.eta2 <- c(.50, .50, .50)
#' 
#' ## Groups sizes for Population
#' BigN <- 100000
#' in.clus.size <- c(BigN*(1-BR), BR * BigN) 
#'  
#' ## Generate Population of scores with "monte"
#' sample.data <- monte(seed = in.seed, 
#'                 nvar=in.nvar, 
#'                 nclus = in.nclus, 
#'                 clus.size = in.clus.size, 
#'                 eta2 = in.eta2, 
#'                 skew.list = in.skew.list, 
#'                 kurt.list = in.kurt.list)
#'                
#'           
#' output <- summary(sample.data)
#' 
#' z <- data.frame(sample.data$data[sample(1:BigN, 600, replace=FALSE),])
#' z[,2:4] <- scale(z[,2:4])
#' names(z) <- c("id","v1","v2","v3")
#' 
#' #trellis.device()
#' trellis.par.set( col.whitebg() )
#' print(
#'  cloud(v3 ~ v1 * v2,
#'        groups = as.factor(id),data=z, 
#'        subpanel = panel.superpose,
#'        zlim=c(-4, 4),
#'        xlim=c(-4, 4),
#'        ylim=c(-4, 4),
#'        main="",
#'        screen = list(z = 20, x = -70)),
#'    position=c(.1, .5, .5, 1), more = TRUE)
#'                
#'                  
#' 
#' ############################
#' ##      Condition 2  
#' ############################
#' 
#' ## Within taxon indicator skew and kurtosis
#' in.skew.list <- list(c(1.75, 1.75, 1.75),c(.50, .50, .50)) 
#' in.kurt.list <- list(c(3.75, 3.75, 3.75),c(0, 0, 0))          
#' 
#' ## Generate Population of scores with "monte"
#' sample.data <- monte(seed = in.seed, 
#'                nvar = in.nvar, 
#'                nclus = in.nclus, 
#'                clus.size = in.clus.size, 
#'                eta2 = in.eta2, 
#'                skew.list = in.skew.list, 
#'                kurt.list = in.kurt.list)
#'                
#'           
#' output <- summary(sample.data)
#' 
#' z <- data.frame(sample.data$data[sample(1:BigN, 600, replace=FALSE),])
#' z[,2:4] <- scale(z[, 2:4])
#' names(z) <-c("id", "v1","v2", "v3")
#' 
#' print(
#'  cloud(v3 ~ v1 * v2,
#'        groups = as.factor(id), data = z, 
#'        subpanel = panel.superpose,
#'        zlim = c(-4, 4),
#'        xlim = c(-4, 4),
#'        ylim = c(-4, 4),
#'        main="",
#'        screen = list(z = 20, x = -70)),
#'        position = c(.5, .5, 1, 1), more = TRUE)
#'                
#'                  
#' ############################
#' ##      Condition 3  
#' ############################ 
#' 
#' ## Set within group correlations to .80
#' cormat <- matrix(.80, 3, 3)
#' diag(cormat) <- rep(1, 3)
#' in.cor.list <- list(cormat, cormat)
#' 
#' ## Generate Population of scores with "monte"
#' sample.data <- monte(seed = in.seed, 
#'                nvar = in.nvar, 
#'                nclus = in.nclus, 
#'                clus.size = in.clus.size, 
#'                eta2 = in.eta2, 
#'                skew.list = in.skew.list, 
#'                kurt.list = in.kurt.list,
#'                cor.list = in.cor.list)
#'                
#' output <- summary(sample.data)
#' 
#' z <- data.frame(sample.data$data[sample(1:BigN, 600, 
#'                 replace = FALSE), ])
#' z[,2:4] <- scale(z[, 2:4])
#' names(z) <- c("id", "v1", "v2", "v3")
#' 
#' ##trellis.device()
#' ##trellis.par.set( col.whitebg() )
#' print(
#'   cloud(v3 ~ v1 * v2,
#'        groups = as.factor(id),data=z, 
#'        subpanel = panel.superpose,
#'        zlim = c(-4, 4),
#'        xlim = c(-4, 4),
#'        ylim = c(-4, 4),
#'        main="",
#'        screen = list(z = 20, x = -70)),
#' position = c(.1, .0, .5, .5), more = TRUE)
#'                                 
#' 
#' ############################
#' ##      Condition 4  
#' ############################
#' 
#' ## Change compactness so that variance of
#' ## complement indicators is 5 times
#' ## greater than variance of taxon indicators
#'                      
#'  v <-  ( 2 * sqrt(5))/(1 + sqrt(5)) 
#'  in.compactness <- c(v, 2-v)   
#'  
#' ## Generate Population of scores with "monte"
#' sample.data <- monte(seed = in.seed, 
#'                nvar = in.nvar, 
#'                nclus = in.nclus, 
#'                clus.size = in.clus.size, 
#'                eta2 = in.eta2, 
#'                skew.list = in.skew.list, 
#'                kurt.list = in.kurt.list,
#'                cor.list = in.cor.list,
#'                compactness = in.compactness)
#'                
#' output <- summary(sample.data)
#' 
#' z <- data.frame(sample.data$data[sample(1:BigN, 600, replace = FALSE), ])
#' z[, 2:4] <- scale(z[, 2:4])
#' names(z) <- c("id", "v1", "v2", "v3")
#' print(
#'   cloud(v3 ~ v1 * v2,
#'        groups = as.factor(id),data=z, 
#'        subpanel = panel.superpose,
#'        zlim = c(-4, 4),
#'        xlim = c(-4, 4),
#'        ylim = c(-4, 4),
#'        main="",
#'        screen = list(z = 20, x = -70)),
#'  position = c(.5, .0, 1, .5), more = TRUE)
#' 
monte<-function(seed=123,nvar=4, nclus=3, clus.size=c(50,50,50), 
               eta2=c(0.619, 0.401, 0.941, 0.929), cor.list=NULL, random.cor=FALSE,
               skew.list=NULL ,kurt.list=NULL,secor=NULL,compactness=NULL,sortMeans=TRUE)
{

    set.seed(seed)  #
    if(is.null(compactness)) compactness<-rep(1,nclus)
    
    if(!is.null(cor.list) || random.cor==TRUE)
        orientation<-TRUE
    else
        orientation<-FALSE    
        
    if(!is.null(skew.list) && !is.null(kurt.list))
       transform<-TRUE
    else
       transform<-FALSE 
 
       
    seed.vec <- c(25,43,47,21,38,42,18,49,41,22,30,3,33,27,45,16,35,37,2,5,20,26,31,29,12,44,
                  11,36,46, 8,10,50,15,34,39, 4,48,19,13,23,9,28,6,17,14,24,7,1,40,32)

      
        
#=================================================================#
# *****************Arguement definitions**************************#
 
#****************************************************************#
#****************************************************************#
#
#   Beware: Traveling beyond this point is dangerous!
#****************************************************************#
#
#*****************************************************************#
#
# DEFINE external functions: mkclusb, mkclusc, mkclusd
#
#
#*****************************************************************#
#
########################################################################
    call<-match.call()
    #library(stats)
    mkclusb.prg <- function(nsub, nvar, skewvec, kurtvec, seed, desired.cor, orientation, callnum)
    {
#=====================================================================#
# mkclusb.prg                                                         #
#                                                                     #
# This program calls external function mkclusc.prg and mkclusd.prg    # 
#                                                                     #
# a program for generating nonnormal multivariate data with           # 
# specified correlation structure using equations described in        #
# Vale, D. C. & Maurelli, V. A. (1983). Simulating multivariate       # 
#    nonnormal distributions.  Psychometrika, 48, 465-471.            #
#                                                                     #
# nvar=number of variables                                            #
# skewvec=vector of skewness values for nvar variables                #
# kurtvec=vector of kurtosis values for nvar variables                #
# seed = seed number for data generation                              #
# desired.cor is a matrix containing the target correlation matrix    #
# for nonnormal data                                                  #
#                                                                     #
# upon completion the function returns X (nsub x nvar) with non-normal#
# data sampled from a population with correlation: desired.cor        # 
#=====================================================================# 

        mkclusc.prg <- function(x, skew = 0, kurt = 0)
        {

#======================================================================#
# mkclusc.prg                                                          #
# This function is minimized to determine the weights needed to        # 
# transform the data to desired skewness and kurtosis                  #
# See Vale and Maurelli (1983)                                         #
#                                                                      #
# b, c, and d are weights used to transform normal data                #
# f, g, and h are the three nonlinear equations that must be solved    #
#    to find b, c, and d                                               #
# skew and kurt are the desired skewness and kurtosis values           #
#======================================================================# 
            b <- x[1]
            c <- x[2]
            d <- x[3]
   
            f <- (b^2 + 6 * b * d + 2 * c^2 + 15 * d^2 - 1)
            g <- 2 * c * (b^2 + 24 * b * d + 105 * d^2 + 2) - skew
            h <- 24 * (b * d + c^2 * (1 + b^2 + 28 * b * d) + d^2 * (12 + 48 * b * d + 141 * c^2 + 225 * d^2)) - kurt
            return((f^2 + g^2 + h^2))
        }
#----------------------------------------------------------------------#
        mkclusd.prg <- function(p, r, matr)
        {
            f <- (p * (matr[1, 2] * matr[2, 2] + 3 * matr[1, 2] * matr[2, 4] + 3 * matr[1, 4] * matr[2, 2] + 9 * matr[1, 4] * matr[
                2, 4]) + p^2 * (2 * matr[1, 3] * matr[2, 3]) + p^3 * (6 * matr[1, 4] * matr[2, 4])) - r
            
                    return(f^2)
        }
#
        #print(cat("\n", "Generating transformed data for cluster", callnum, "\n")) #   
#=====================================================================#
#  Generating weights for nonnormal transformation by minimizing      #
#  a set of nonlinear equations                                       #
#                                                                     #
#  mkclusc.prg is called here                                         #
#=====================================================================#
# bcdvec is a vector of the b,c, and d weights from Vale and Maurelli (1983) #

        bcdvec <- matrix(0, nrow = nvar, ncol = 3)  #
        for(i in 1:nvar) {
            bcdvec[i,] <-  optim(par = c(1.0, .0, .0),fn = mkclusc.prg,method="L-BFGS-B", lower=-2,upper=2, skew = skewvec[i], 
                kurt = kurtvec[i],,control=list(ndeps=rep(1e-7,3)))$par
                        bcdvec[i,1]<-bcdvec[i,1]*-1
                        bcdvec[i,3]<-bcdvec[i,3]*-1          
                          
        }
# (avec=a) a = -c in the Vale and Maurelli equation          
        avec <-  - bcdvec[, 2]  #
#matrix of weights (the regression constants) 
#used for transformation; (a,b,c,d) is a nvar X 4 (abcd) matrix#
        constant.mat <- as.matrix(cbind(avec, bcdvec))  #
#
#================================================================#
# intermediate correlation matrix (rxx)for normal data generation#
# needed to determine the required correlation structure that    #
# will result in the desired correlation structure after the data#
# been transformed to desired skewness and kurtosis              #
# This is done by minimizing the function found in mkclusd.prg   #
#                                                                #
# matr contains the a b c d weights for the two variables of     #
# desired cor                                                    #
#================================================================#
        rxx <- desired.cor
        if(orientation == TRUE) {
            rxx <- matrix(c(rep(0, (nvar * nvar))), nrow = nvar, ncol = nvar)
            for(r in 2:nvar) {
                for(col in 1:(r - 1)) {

                  rxx[r, col] <- optim(method="L-BFGS-B", par = 0.1, fn = mkclusd.prg,control=list(ndeps=1e-7),lower=-.99,upper=.99,r = desired.cor[r, col], 
                    matr = matrix(c(constant.mat[r,  ], constant.mat[col,  ]), nrow = 2, ncol = 4, byrow = TRUE))$par
                            rxx[col, r] <- rxx[r, col]
                }
            }
            I <- diag(nvar) #
# I is an identity matrix (nvar X nvar) used to place ones on the diagonal
# of rxx
            rxx <- rxx + I  #
        }
#======================================================================#
#  The intermediate correlation matrix is decomposed using a           #
#  Cholesky decomposition.  The resulting triangular matrix and stand- #
#  ardized random variates are used to create a matrix of              #
#  normally distributed data with the desired pattern of               #
#  correlations                                                        #
#======================================================================#
        set.seed(seed)  # within.clus.dev is an nsub X nvar data matrix used to develop clustered data#
        within.clus.dev <- matrix(rnorm(nvar * nsub), ncol = nvar, nrow = nsub) #
        within.clus.dev <- apply(within.clus.dev, 2, scale)
        floadt <- (chol(rxx))   #
# floadt is the triangular matrix of loadings resulting from a Cholesky
# decomposition of the intermediate correlation matrix
        X <- within.clus.dev %*% floadt #
# X is data which has desired correlation structure (based on rxx) 
#===============================================================#
# The weight matrix is now used to transform the data to have   #
# desired skewness, kurtosis and correlation structure          #
# using Fleischman and Vale & Maurelli's power method           #
# Fleishman, A.I (1978). A method for simulating non-normal     #
# distributions, Psychometrika, 43, 521-532.                    #
#===============================================================#
        one <- rep(1, nsub) # a vector of ones #
        X2 <- X^2
        X3 <- X^3
        for(i in 1:nvar) {
            ytemp <- cbind(one, X[, i], X2[, i], X3[, i])
            w <- matrix(c(constant.mat[i,  ]), ncol = 1)
            X[, i] <- ytemp %*% w
        }
        X <- apply(X, 2, scale) #
# X contains the nonnormal, standardized data #
        return(X)
    }
#
#
#==============================================================#
# mkclusd.prg                                                  #
# function based on equation found in Vale and Maurelli(1983)  #
# this function is minimized to obtain an intermediate         #
# correlation matrix (used for normally distributed data)      #
# which will work with certain weights to achieve the desired  #
# correlation structure for the transformed (nonnormal) data   #
#                                                              #
# p is the correlation that is being solved for (the intermed- #
#   iate correlation).  Start values are provided.             #
# r is provided. It is the correlation desired for the non-    #
#   normal data.                                               #
# matr is a matrix of transformation weights (a,b,c,d; from    #
#   constant.mat in mkclusb.prg                                #
#==============================================================#
#*****************************************************************
#
#                        MAIN PROGRAM STARTS HERE
#=================================================================#
# use this transformation if random cluster volumes are desired
#       compactness <- runif(length(clus.size))
#       compactness <- sqrt(length(clus.size) * (compactness/sum(compactness)))   #
#=================================================================# 
   
    if(transform==TRUE){
        skewvec <- matrix(0, nclus, nvar)
        kurtvec <- matrix(0, nclus, nvar)
         for(i in 1:nclus) {
            skewvec[i,  ] <- skew.list[[i]]
            kurtvec[i,  ] <- kurt.list[[i]]
         }
   }      
       
    
    callnum <- 1    #used by mkclusb.prg to inform user of current cluster#
    cornum <- round(1/secor^2 + 3, 0)
    nsub <- sum(clus.size)    #number of subjects #   
    sum.size <- c(0, cumsum(clus.size))   #
#=================================================================#
# sum.size is a nclus+1 vector with element j equal to  
# the cummulative cluster sample sizes 
#
#================See Equation 15 in Waller, Underhill, & Kaiser=================#
    eta2 <- 1 + (nsub * eta2 - nsub)/(eta2 * sum(clus.size * compactness^2) - (nsub * eta2) + nsub)   #
#==================================================================
# 
#=======make design.mat===and dummy factor scores===============#   
    design.mat <- matrix(0, nrow = nclus, ncol = nclus - 1)
    for(i in 1:(nclus - 1)) {
        design.mat[i, i] <- 1
    }
    dm <- design.mat    #
#
# dm or design.mat is a nclus X (nclus-1) matrix with all elements#
# 0 except for those where column and row number are equal (these are 1) #
# 
# initiate factor score matrix #
    fs <- matrix(0, nrow = nsub, ncol = nclus - 1)
    k <- 1  # put in appropriate dummy scores #
    for(i in 1:nclus) {
        for(j in 1:clus.size[i]) {
            fs[k,  ] <- dm[i,  ]
            k <- k + 1
        }
    }
#
#
# fs is the "dummy" factor score matrix (nsub X nclus-1) with 0's and 1's #
# used to demarcate cluster membership#
#
    rawfs <- fs # save backup of dummy factor scores
#===================standardize fs and orthogonalize=========#
    if(ncol(fs) == 1) {
        fs <- scale(fs)
    }
    if(ncol(fs) >= 2) {
        fs <- apply(fs, 2, scale)   #
        fs <- apply(princomp(fs, cor = TRUE)$scores, 2, scale)
    }
#
#===generate dummy factor loading matrix=====================#
#
    fl <- matrix(runif(nvar * (nclus - 1), min = -1, max = 1), nrow = nvar, ncol = nclus - 1)   #
# fl is a factor loading matrix(nvar X nclus-1) from random uniform variates# 
#
    sigmahat <- fl %*% t(fl)    #
#
#=================================================================#
# sigmahat is the sum of squares and crossproducts of the loadings
# 
# Note, at this point the communalities can be greater than 1, thus
# the factor loading matrix must be normalized by rows
#=================================================================#
    h2 <- diag(sigmahat)    #
#
# h2 are the sum of squared factor loadings #
#===========normalize rows of dummy factor loading matrix=====#
    signfl <- sign(fl)  #
#
# signfl is a matrix containing +1 and -1 for sign of loadings#
#
# this normalizes rows #
    for(i in 1:nvar) {
        fl[i,  ] <- sqrt(eta2[i]) * sqrt(fl[i,  ]^2/h2[i])
    }
# putting the correct sign back on the altered  fl#
    fl <- signfl * fl   #
#   
#=================================================================#
    sigmahat <- fl %*% t(fl)    #
    h2 <- diag(sigmahat)    #
#
#==============lmn=latent means==================================#
# the latent means (cluster centroids) are computed by post multiplying
# the normalized factor loading matrix (fl) by the standardized (and 
# orthogonalized) dummy factor scores
#================================================================#
    lmn <- matrix(0, nrow = nvar, ncol = nclus)
    for(k in 1:nvar) {
        j <- 1
        for(i in 1:nclus) {
            lmn[k, i] <- sum(fl[k,  ] * fs[j,  ])
            j <- j + clus.size[i]
        }
    }
    lmn <- t(lmn)
    if(sortMeans) lmn <- apply(lmn, 2, sort)    #
#================================================================#
#===================For zero nuisance covariance=================#
    if(orientation == FALSE) {
        within.clus.dev <- matrix(rnorm(nsub * nvar), nrow = nsub, ncol = nvar) #
#
# within.clus.dev is an nsub X nvar matrix containing the simulated data #
# without the latent means, i.e., the within cluster deviation scores #
#
        for(i in 1:nclus) {
            if(transform == TRUE) {                   
                       
#==============nonlinear transformation of data=======================#
#              external prog. mkclusb.prg called here                 #
#              the desired correlation matrix is an identity matrix   #
#              when orientation = FALSE and thus no transformation of   # 
#              cors is necessary                                      #
#=====================================================================#
                desired.cor <- diag(nvar)
                within.clus.dev[(1 + sum.size[i]):sum.size[i + 1],  ] <- mkclusb.prg(nsub = clus.size[i], nvar = nvar, skewvec
                   = skewvec[i,  ], kurtvec = kurtvec[i,  ], seed = seed + i, desired.cor = desired.cor, orientation = 
                  orientation, callnum = callnum)
                callnum <- callnum + 1
            }
            if(nvar > 1) {
                within.clus.dev[(1 + sum.size[i]):sum.size[i + 1],  ] <- compactness[i] * apply(within.clus.dev[(1 + sum.size[
                  i]):sum.size[i + 1],  ], 2, scale) %*% diag(sqrt(1 - h2)) #
            }
            if(nvar == 1) {
                within.clus.dev[(1 + sum.size[i]):sum.size[i + 1]] <- compactness[i] * scale(within.clus.dev[(1 + sum.size[i]):
                  sum.size[i + 1]]) * sqrt(1 - h2)
            }
        }
    }
#========================for nuisance covariance =================#
# Using the Cholesky decomposition of the desired correlation     #
# matrix we generate a loading matrix (possibly for each cluster) #
# for the within cluster scores.                                  #
# This is done cluster by cluster below  with either randomly     #
# chosen within cluster correlation structure (random.cor=TRUE)      #
# or with desired correlations set by the user (cor.list,         #
# which is renamed desired.cor).                                  #
#=================================================================#
    if(orientation == TRUE) {
       
#
# uscor is a matrix of random normal deviates that will equal 
# within.clus.dev after appropriate weighting by uf (uf=loading matrix 
# from random cor structure
        uscor <- matrix(rnorm(nsub * nvar), nrow = nsub, ncol = nvar)
        for(i in 1:nclus) {               
            uscor[(1 + sum.size[i]):sum.size[i + 1],  ] <- apply(uscor[(1 + sum.size[i]):sum.size[i + 1],  ], 2, scale)
        }
#
#============For within group correlations that are chosen randomly with
#============a standard error = secor
        if(random.cor == TRUE) {
            detx <- -99 #used for determinant of x
            while(detx < 0.05) {
                x <- matrix(rnorm(cornum * nvar), nrow = cornum, ncol = nvar)
                detx <- prod(eigen(abs(cor(x)))$values)
            }
#
# a matrix of x values, standard error of the absolute value of 
# correlations of x=secor
            uf <- chol(abs(cor(x)))
            uf <- t(uf) #
# uf is a triangular matrix from the Cholesky decomp. of abs(correlation of x)#
# (the target correlation matrix) #
            print(cat("Expected within group correlation matrix", "\n"))
            expect.cor <- round(uf %*% t(uf), 3)
            print(expect.cor)
            within.clus.dev <- matrix(0, nrow = nvar, ncol = nsub)
            uscor <- t(uscor)
            for(i in 1:nclus) {
#
#=====================================================================#
                if(transform == TRUE) {                             
                  if(transform == TRUE & i == 1) {
                    within.clus.dev <- t(within.clus.dev)
                  }
#==============================nonlinear transformation of data=====#
#==============================external prog. mkclusb.prg called here=#
                  within.clus.dev[(1 + sum.size[i]):sum.size[i + 1],  ] <- mkclusb.prg(nsub = clus.size[i], nvar = nvar, skewvec
                     = skewvec[i,  ], kurtvec = kurtvec[i,  ], seed = seed + i, desired.cor = expect.cor, orientation = 
                    orientation, callnum = callnum)
                  callnum <- callnum + 1
                }
#
                if(transform != TRUE) {
                  within.clus.dev[, (1 + sum.size[i]):sum.size[i + 1]] <- uf %*% uscor[, (1 + sum.size[i]):sum.size[i + 1]]
                }
            }
        }
#
#============If within cluster correlations are specified by the user==========#
#============via  (cor.list)===============================================#
        if(random.cor == FALSE) {
            within.clus.dev <- matrix(0, nrow = nvar, ncol = nsub)
            uscor <- t(uscor)
            for(i in 1:nclus) {
                if(transform == TRUE) {
                  print(i)
                  if(transform == TRUE & i == 1) {
                    within.clus.dev <- t(within.clus.dev)
                  }
                  desired.cor <- matrix(unlist(cor.list[i]), nrow = nvar)   # unlisting cor.list
#==============================nonlinear transformation of data================#
#==============================external prog. mkclusb.prg called here==========#
                  within.clus.dev[(1 + sum.size[i]):sum.size[i + 1],  ] <- mkclusb.prg(nsub = clus.size[i], nvar = nvar, skewvec
                     = skewvec[i,  ], kurtvec = kurtvec[i,  ], seed = seed + seed.vec[i], desired.cor = desired.cor, orientation = 
                    orientation, callnum = callnum)
                  callnum <- callnum + 1
                }
#
                if(transform != TRUE) {
                  uf <- chol(matrix(unlist(cor.list[i]), nrow = nvar, byrow = T))
                  uf <- t(uf)
                  within.clus.dev[, (1 + sum.size[i]):sum.size[i + 1]] <- uf %*% uscor[, (1 + sum.size[i]):sum.size[i + 1]]
                }
            }
            if(transform == FALSE) {
                within.clus.dev <- t(within.clus.dev)
            }
        }
#
        for(i in 1:nclus) {
            if(orientation == TRUE & random.cor == TRUE & transform == FALSE & i == 1) {
                within.clus.dev <- t(within.clus.dev)
            }
            within.clus.dev[(1 + sum.size[i]):sum.size[i + 1],  ] <- apply(within.clus.dev[(1 + sum.size[i]):sum.size[i + 1],  ], 
                2, scale)
        }
        within.clus.dev <- within.clus.dev %*% diag(sqrt(1 - h2))
        for(i in 1:nclus) {
            within.clus.dev[(1 + sum.size[i]):sum.size[i + 1],  ] <- compactness[i] * within.clus.dev[(1 + sum.size[i]):sum.size[i +
                1],  ]
        }
    }
#=====End "if orientation =TRUE loop=====================================#
# now add in the mean structure ======================================= #
    data.means <- fs %*% t(fl)  
    id <- rep(seq(1, nclus, 1), clus.size)
 
 
 ## SORT MEANS == TRUE
    if(sortMeans){
        data.means <- apply(data.means, 2, sort)   
        dm<-apply(apply(data.means,2,unique),2,sort)
        for(i in 1:nclus){
           for(j in 1:nvar){ 
              data.means[id==i,j]<-dm[i,j]
           }   
        }
    }    
 
              
    data <- data.means + within.clus.dev    #fs <- apply(fs, 2, unique)

    data <- cbind(id, data)
    result<-list(data = data, lmn = lmn, fl = round(fl, 4), fs = fs,call=call,nclus=nclus,nvar=nvar,cor.list=cor.list,
                 skew.list=skew.list,kurt.list=kurt.list,clus.size=clus.size, eta2=eta2,seed=seed)
    class(result)<-"monte"
    result
}
