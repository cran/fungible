# Tips:
# to debug code -- run in R console:    debugonce(simFA)
# then step through calculations until error occurs

## simFA
##
## Purpose:
##    Generate Factor Analysis Models and Data Sets for Simulation Studies 
## Date:
##    October 16, 2019 return W
##    Jan 21 2019 added specific factor correlations
##    Dec 31 IRT parameters
##    Dec 20 fixed user Phi for bifactor models
##    November 11 2018 
##    October 16, 2018
##    July 9, 2018
##    Time 11:00

# Search '2DO' for sections in need of updates



#' Generate Factor Analysis Models and Data Sets for Simulation Studies
#' 
#' A function to simulate factor loadings matrices and Monte Carlo data sets
#' for common factor models and bifactor models.
#' 
#' \code{simFA} was specifically designed to simplify the process of running
#' Monte Carlo studies of factor analysis models. Thus, \code{simFA} can save
#' all relevant output for a user-specified model. Saved output can be accessed
#' by calling one or more of the following object names.
#' 
#' @param Model (list) 
#'      \itemize{ 
#'              \item \code{NFac} (scalar) Number of common or
#'                            group factors; defaults to \code{NFac = 3}.  
#'              \item \code{NItemPerFac}
#'                   \itemize{ 
#'                      \item (scalar) All factors have the same number
#'                         of primary loadings.  
#'                       \item (vector) A vector of length \code{NFac} specifying the
#'                             number of primary loadings for each factor; defaults 
#'                             to \code{NItemPerFac = 3}.  } 
#'                      \item \code{Model} (character) \code{"orthogonal"} or
#'                            \code{"oblique"}; defaults to 
#'                            \code{Model = "orthogonal"}.  }
#'                               
#' @param Loadings (list) 
#'      \itemize{ 
#'              \item \code{FacPattern} (\code{NULL} or matrix).        
#'                  \itemize{ 
#'                           \item \code{FacPattern = M} where \code{M} is 
#'                                 a user-defined factor pattern matrix.  
#'                           \item \code{FacPattern = NULL};
#'                           \code{simFA} will generate a factor pattern based on 
#'                            the arguments specified  under other keywords 
#'                           (e.g., \code{Model}, \code{CrossLoadings}, etc.);
#'                           defaults to \code{FacPattern = NULL}.  } 
#'               \item \code{FacLoadDist}  (character) Specifies the 
#'                   sampling distribution for the common factor
#'                   loadings. Possible values are \code{"runif"}, 
#'                   \code{"rnorm"}, \code{"sequential"}, and \code{"fixed"}; 
#'                   defaults to \code{FacLoadDist = "runif"}.  
#'               \item \code{FacLoadRange} (vector of length \code{NFac}, 
#'                       2, or 1); defaults to \code{FacLoadRange = c(.3, .7)}.  
#'                       \itemize{ 
#'                         \item If \code{FacLoadDist = "runif"} the vector 
#'                         defines the bounds of the uniform distribution; 
#'                         \item If \code{FacLoadDist = "rnorm"} the vector defines the
#'                          mean and standard deviation of the normal distribution from 
#'                          which loadings are sampled.  
#'                        \item If \code{FacLoadDist = "sequential"} the vector
#'                         specifies the lower and upper bound of the loadings sequence.  
#'                        \item If \code{FacLoadDist = "fixed"} and 
#'                                 \code{FacLoadRange} is a vector of length 1
#'                                 then all common loadings will equal the constant 
#'                                 specified in \code{FacLoadRange}. 
#'                                 If \code{FacLoadDist = "fixed"} and 
#'                                 \code{FacLoadRange} is a vector of length 
#'                                 \code{NFac} then each factor will have fixed loadings
#'                                 as specified by the associated element in 
#'                                 \code{FacLoadRange}.  }
#' 
#'                \item \code{h2} (vector) An optional vector of communalities 
#'                     used to constrain the population communalities to user-defined 
#'                     values; defaults to \code{h2 = NULL}.  }
#'                        
#' @param CrossLoadings (list) 
#'      \itemize{ 
#'          \item \code{ProbCrossLoad} (scalar) A value in the (0,1) interval 
#'              that determines the probability that a cross loading will be 
#'              present in elements of the loadings matrix that do not have
#'              salient (primary) factor loadings.  
#'              If set to \code{ProbCrossLoad = 1}, a
#'              single cross loading will be added to each factor;  defaults to
#'              \code{ProbCrossLoad = 0}.
#'          \item \code{CrossLoadRange} (vector of length 2) Controls 
#'              size of the crossloadings; defaults to 
#'              \code{CrossLoadRange= c(.20, .25)}.  
#'          \item \code{CrossLoadPositions} (matrix) Specifies the row and column 
#'              positions of (optional) cross-loadings; defaults to 
#'              \code{CrossLoadPositions = NULL}.
#'          \item \code{CrossLoadValues} (vector) If 
#'                \code{CrossLoadPositions} is specified then 
#'                \code{CrossLoadValues} is a vector of user-supplied
#'                 cross-loadings; defaults to 
#'                 \code{CrossLoadValues = NULL}.  
#'          \item \code{CrudFactor} (scalar) Controls the size of tertiary 
#'             factor loadings. If \code{CrudFactor != 0} then elements of 
#'             the loadings matrix with neither primary nor secondary 
#'             (i.e., cross) loadings will be sampled from a
#'             [-(CrudFactor), (CrudFactor)] uniform distribution; 
#'             defaults to \code{CrudFactor = 0}.  }
#'             
#'             
#' @param Phi (list) 
#'      \itemize{
#'          \item \code{MaxAbsPhi} (scalar) Upper (absolute) bound on factor
#'              correlations; defaults to \code{MaxAbsPhi = .5}.  
#'          \item \code{EigenValPower} (scalar) Controls the skewness of the 
#'             eigenvalues of Phi. Larger values of \code{EigenValPower} result 
#'             in a Phi spectrum that is more right-skewed (and
#'             thus closer to a unidimensional model); 
#'             defaults to \code{EigenValPower = 2}.  
#'          \item \code{PhiType} (character); defaults to 
#'             \code{PhiType = "free"}.
#'             \itemize{ 
#'                \item If \code{PhiType = "free"} factor correlations 
#'                   will be randomly generated under the constraints of 
#'                   \code{MaxAbsPhi} and \code{EigenValPower}.  
#'                 \item If \code{PhiType = "fixed"} all factor
#'                   correlations will equal the value specified in 
#'                   \code{MaxAbsPhi}. A fatal error will be produced if 
#'                   \code{Phi} is not positive semidefinite.  
#'                 \item If \code{PhiType = "user"} the factor correlations are 
#'                    defined by the matrix specified in \code{UserPhi} (see below).} 
#'            \item \code{UserPhi} (matrix) A positive semidefinite (PSD) matrix 
#'              of user-defined factor correlations;defaults to \code{UserPhi = NULL}.}
#'              
#'              
#' @param ModelError (list) 
#'      \itemize{ 
#'         \item \code{ModelError} (logical) If
#'               \code{ModelError = TRUE} model error will be introduced into 
#'               the factor pattern via the method described by Tucker, Koopman, 
#'               and Linn (TKL, 1969); defaults to \code{ModelError = FALSE}.  
#'         \item \code{NMinorFac} (scalar) Number of minor factors in the TKL model; 
#'             defaults to \code{NMinorFac = 150}.  
#'         \item \code{ModelErrorType} (character) If \code{ModelErrorType =
#'             "U"} then \code{ModelErrorVar} is the proportion of uniqueness 
#'             variance that is due to model error.  
#'             If \code{ModelErrorType = "V"} then
#'             \code{ModelErrorVar} is the proportion of total variance that 
#'             is due to model error; defaults to \code{ModelErrorType = "U"}.  
#'         \item \code{ModelErrorVar} (scalar [0,1]) The proportion of 
#'             uniqueness (U) or total (V) variance that is due to model error; 
#'             defaults to \code{ModelErrorVar = .10}.  
#'         \item \code{epsTKL} (scalar [0,1]) Controls the size of the factor 
#'            loadings in successive minor factors; defaults to \code{epsTKL = .20}.  
#'         \item \code{RSpecific} (matrix) Optional correlation
#'             matrix for specific factors; defaults to \code{RSpecific = NULL}. }
#'             
#'             
#' @param Bifactor (list) 
#'      \itemize{ 
#'         \item Bifactor (logical) If \code{Bifactor = TRUE} parameters for 
#'             the bifactor model will be generated; defaults to
#'             \code{Bifactor = FALSE}.  
#'         \item Hierarchical (logical) If \code{Hierarchical
#'              = TRUE} then a hierarchical Schmid Leiman (1957) bifactor model 
#'              will be generated; defaults to \code{Hierarchical = FALSE}.  
#'         \item 
#'            \code{F1FactorDist} (character) Specifies the sampling distribution 
#'               for the general factor loadings.  Possible values are 
#'            \code{"runif"}, 
#'            \code{"rnorm"}, 
#'            \code{"sequential"}, and 
#'            \code{"fixed"}; defaults to
#'            \code{F1FactorDist = "sequential"}.  
#'         \item \code{F1FactorRange} (vector of length 1 or 2) Controls the 
#'            sizes of the general factor loadings in
#'            nonhierarchical bifactor models; defaults to 
#'            \code{F1FactorRange = c(.4, .7)}.  
#'              \itemize{ 
#'                  \item If \code{F1FactorDist = "runif"}, the vector of
#'                      length 2 defines the bounds of the uniform distribution, 
#'                      c(lower, upper);
#'                   \item If \code{F1FactorDist = "rnorm"}, the vector defines 
#'                      the mean and standard deviation of the normal distribution 
#'                      from which loadings are sampled, c(MN, SD).  
#'                   \item If \code{F1FactorDist = "sequential"}, the vector
#'                      specifies the lower and upper bound of the loadings 
#'                      sequence, c(lower, upper).  }
#'         }
#'         
#'         
#' @param MonteCarlo (list) 
#'      \itemize{ 
#'          \item \code{NSamples} (integer) Defines
#'             number of Monte Carlo Samples; defaults to \code{NSamples = 0}.  
#'          \item \code{SampleSize} (integer) Sample size for each Monte 
#'             Carlo sample; defaults to \code{SampleSize = 250}.  
#'          \item \code{Raw} (logical) If \code{Raw = TRUE}, simulated data sets 
#'             will contain raw data.  If \code{Raw = FALSE}, simulated data sets 
#'             will contain correlation matrices; defaults to
#'             \code{Raw = FALSE}.  
#'          \item \code{Thresholds} (list) List elements contain
#'              thresholds for each item. Thresholds are required when generating 
#'              Likert variables.  }
#'              
#'              
#' @param FactorScores (list) 
#'      \itemize{ 
#'          \item \code{FS} (logical) If \code{FS = TRUE} (true) factor scores 
#'             will be simulated; defaults to \code{FS = FALSE}.
#'           \item \code{CFSeed} (integer) Optional starting seed for the common 
#'              factor scores; defaults to \code{CFSeed = NULL} in which case a 
#'              random seed is used.  
#'           \item \code{SFSeed} (integer) Optional starting seed for the 
#'              specific factor scores; defaults to \code{SFSeed = NULL} in 
#'              which case a random seed is used.  
#'           \item \code{EFSeed} (integer) Optional starting seed for the error
#'               factor scores; defaults to 
#'                 \code{EFSeed = NULL} in which case a 
#'                   random seed is used. Note that 
#'                 \code{CFSeed}, 
#'                 \code{SFSeed}, and 
#'                 \code{EFSeed} must be
#'                 different numbers (a fatal error is produced when two or more 
#'                 seeds are specified as equal).  
#'            \item \code{VarRel} (vector) A vector of manifest
#'                variable reliabilities.  The specific factor variance for 
#'                variable \emph{i} will equal \eqn{VarRel[i] - h^2[i]} 
#'                (the manifest variable reliability minus
#'                 its commonality). By default, \eqn{VarRel = h^2} (resulting 
#'                 in uniformly zero specific factor variances).  
#'            \item \code{Population} (logical) If \code{Population = TRUE}, 
#'               factor scores will fit the correlational
#'               constraints of the factor model exactly (e.g., the common factors 
#'               will be orthogonal to the unique factors); defaults to 
#'               \code{Population = FALSE}.
#'             \item \code{NFacScores} (scalar) Sample size for the factor scores; 
#'                 defaults to \code{NFacScores = 250}.  
#'             \item \code{Thresholds} (list) A list of quantiles used to polychotomize 
#'                the observed data that will be generated from the factor scores.  }
#'                
#'                
#' @param Missing (list) 
#'      \itemize{ 
#'         \item Missing (logical) If \code{Missing = TRUE} all data sets will 
#'            contain missing values; defaults to \code{Missing = FALSE}.  
#'         \item \code{Mechanism} (character) Specifies the missing data
#'             mechanism.  Currently, the program only supports missing completely at
#'             random (MCAR): \code{Missing = "MCAR"}.  
#'         \item \code{MSProb} (scalar or vector of length \code{NVar}) Specifies 
#'            the probability of missingness for each variable; defaults to 
#'            \code{MSprob = 0}.  }
#'            
#'            
#' @param Control (list) 
#'      \itemize{ 
#'         \item \code{Maxh2} (scalar) Rows of the
#'             loadings matrix will be rescaled to have a maximum communality of
#'             \code{Maxh2}; defaults to \code{Maxh2 = .98}.  
#'             \code{itemReflect} (logical) If \code{Reflect = TRUE} loadings on 
#'             the common factors will be randomly reflected; defaults to 
#'             \code{Reflect = FALSE}.  }
#'             
#'             
#' @param Seed (integer) Starting seed for the random number generator;
#'          defaults to \code{Seed = NULL}. When no seed is specified by the 
#'          user, the program will generate a random seed.
#'          
#' @return 
#'      \itemize{ 
#'          \item \code{loadings} A common factor or bifactor loadings
#'              matrix.  
#'          \item \code{Phi} A factor correlation matrix.  
#'          \item \code{urloadings} The unrotated loadings matrix.  
#'          \item \code{h2} A vector of item commonalities.  
#'          \item\code{h2PopME} A vector item commonalities that may include
#'                 model approximation error.
#'          \item \code{Rpop} The model-implied population correlation matrix.  
#'          \item \code{RpopME} The model-implied population
#'                 correlation matrix with model error.  
#'          \item \code{CovMatrices} A list
#'                 containing: 
#'                 \itemize{ 
#'                    \item \code{CovMajor} The model implied covariances
#'                       from the major factors.  
#'                    \item \code{CovMinor} The model implied covariances
#'                       from the minor factors.  
#'                    \item \code{CovUnique} The model implied variances
#'                       from the uniqueness factors.  } 
#'                \code{Bifactor} A list containing: 
#'                  \itemize{
#'                      \item \code{loadingsHier} Factor loadings of the 1st 
#'                         order solution of a hierarchical bifactor model.  
#'                      \item \code{PhiHier} Factor correlations of
#'                         the 1st order solution of a hierarchical bifactor model.  } 
#'                  \item \code{Scores} A list containing: 
#'                   \itemize{ 
#'                       \item\code{FactorScores} Factor scores for the common 
#'                            and uniqueness factors.  
#'                       \item \code{FacInd} Factor
#'                           indeterminacy indices for the error free population model.  
#'                       \item \code{FacIndME} Factor score indeterminacy indices for 
#'                          the population model with model error.  
#'                       \item \code{ObservedScores} A matrix of model implied
#'                             \code{ObservedScores}. If 
#'                             \code{Thresholds} were supplied under Keyword
#'                             \code{FactorScores}, 
#'                             \code{ObservedScores} will be transformed into Likert
#'                                  scores.  } 
#'                 \item\code{Monte} A list containing output from the Monte Carlo
#'                      simulations if generated.  
#'                 \item\code{IRT}Factor loadings expressed in the
#'                     normal ogive IRT metric. If \code{Thresholds} were given 
#'                     then IRT difficulty values will also be returned.  
#'                 \item\code{Seed}The initial seed for the
#'                      random number generator.  
#'                 \item\code{call}A copy of the function call.
#'                 \item\code{cn}A list of all active and nonactive function arguments.}
#'                 
#'                 
#' @author Niels G. Waller
#' 
#' @references
#' 
#' Schmid, J. and Leiman, J. M.  (1957).  The development of hierarchical
#' factor solutions.  Psychometrika, 22(1), 53--61.
#' 
#' Tucker, L. R., Koopman, R. F., and Linn, R. L.  (1969).  Evaluation of
#' factor analytic research procedures by means of simulated correlation
#' matrices.  Psychometrika, 34(4), 421--459.
#' 
#' @keywords stats
#' 
#' @examples
#' 
#' #  Ex 1. Three Factor Simple Structure Model with Crossloadings and 
#' #  Ideal Nonsalient Loadings
#'    out <-  simFA(Seed = 1)
#'    print( round( out$loadings, 2 ) )
#' 
#' # Ex 2. Non Hierarchical bifactor model 3 group factors
#' # with constant loadings on the general factor
#'    out <- simFA(Bifactor = list(Bifactor = TRUE,
#'                                 Hierarchical = FALSE,
#'                                 F1FactorRange = c(.4, .4),
#'                                 F1FactorDist = "runif"),
#'                 Seed = 1)
#'    print( round( out$loadings, 2 ) ) 
#'    
#' @importFrom stats cov2cor lm rchisq resid rnorm runif
#' @export simFA
#' 
simFA <- function(Model = list(),
                  Loadings = list(),
                  CrossLoadings = list(),
                  Phi = list(),
                  ModelError = list(),
                  Bifactor = list(), 
                  MonteCarlo = list(),
                  FactorScores = list(),
                  Missing = list(),
                  Control = list(),
                  Seed = NULL){  
  
  # If a fatal script-error is detected I print a custom 
  # error msg and suppress R's built in msg. 
  # This turns error messages back on
  options( show.error.messages = TRUE ) 
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #                ---- SEC 1: Read in Argument Lists----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
    # ---- _____Model (MD) Default values ----
    # NFac          scalar default = 3
    # NItemPerFac   scalar if all factors have the 
    #                      same number of primary loadings
    #               vector if each factor has a 
    #                      possibly different number of 
    #                      primary loadings
    # Model         character  "orthogonal" (Default)  
    #                      or "oblique"
    #
    cnMD <- list(NFac = 3,            
                 NItemPerFac = 3, 
                 Model = "orthogonal") 
    cnMD_Len = 3
    # legitimate names for Model
    MD_Arg_Names = sort(names(cnMD))
  
    cnMD[(namMD <- names(Model))] <- Model


   # ----_____Loadings (LD) Default values ----
   # FacPattern   NULL or matrix. defaults to FacPattern = NULL.
   #                 A user-defined factor pattern matrix can be 
   #                 specified with this argument.  
   # FacLoadRange  vector of length 2. Default(FacLoadRange = c(.3, .7)).
   #              If FacLoadDist = "runif"
   #                 the vector defines the bounds of the
   #                 uniform distribution;
   #              If FacLoadDist = "rnorm" 
   #                 the vector defines the mean and 
   #                 standard deviation of the 
   #                 the normal distribution from which loadings
   #                 are sampled.
   #              If FacLoadDist = "sequential" 
   #                 the vector specifies the lower and upper 
   #                 bound of the loadings sequence.
   #              If FacLoadDist = "fixed" all common loadings
   #                 will equal the constant specified in FacLoadRange
   # FacLoadDist  character: Specifies the sampling distribution for 
   #                 the common factor loadings. Possible values are
   #                 "runif", "rnorm", "sequential". 
   #                 Default( FacLoadDist = "runif").
   # h2           vector of user-defined commonalities
   #              Default(h2 = NULL)
   #
   cnLD <- list(FacPattern = NULL,
                FacLoadDist = "runif",
                FacLoadRange = c(.3, .7), 
                h2 = NULL)
   cnLD_Len = 4
   LD_Arg_Names = sort(names(cnLD))
   cnLD[(namLD <- names(Loadings))] <- Loadings
  
 
   # ----_____CrossLoadings (CL) Default values  ----
   # ProbCrossLoad   scalar: A value in the (0,1) interval
   #                         that determines the probability that a 
   #                         cross loading will be present in
   #                         elements of the loadings matrix that do not 
   #                         have salient (primary) factor loadings.
   #                         If set to 1, a single cross loading is added
   #                         to each factor. 
   # CrossLoadRange vector of length 2 Controls size sizes of the
   #                          crossloadings.  
   #                          Default( CrossLoadRange= c(.20, .25) ).
   # CrossLoadPositions matrix specifying the row and column positions of 
   #                          secondary loadings (Default = NULL)
   # CrossLoadValues   vector  if !is.null(CrossLoadPositions) then 
   #                          user-specified cross-loadings taken from 
   #                          CrossLoadValues
   # CrudFactor    scalar     Controls the size of tertiary factor loadings.
   #                          Default(CrudFactor = 0).
   #
   cnCL = list(ProbCrossLoad = .0,
               CrossLoadRange = c(.20, .25), 
               CrossLoadPositions = NULL,
               CrossLoadValues = NULL,
               CrudFactor = 0)
   cnCL_Len = 5
   CL_Arg_Names = sort(names(cnCL))
   cnCL[(namCL <- names(CrossLoadings))] <- CrossLoadings
   
  

   #  ----_____Phi (PH) Default values ----
   # PhiType       character  free, fixed, user
   #                     If PhiType = free factor correlations
   #                     will be randomaly generated 
   #                     under the constraints of MaxAbsPhi 
   #                     and EigenValPower
   # UserPhi      optional input Phi Matrix
   #
   # MaxAbsPhi  scalar:  upper (absolute) bound on factor
   #                     correlations. Default( MaxAbsPhi = .5 ).
   # EigenValPower scalar: Controls the skewness of
   #                     the eigenvalues of Phi. Larger values of
   #                     EigenValPower result in a Phi spectrum
   #                     that is more right-skewed (and thus 
   #                     closer to a unidimensional model).
   #                     Default( EigenValPower = 2 )
   #
   cnPH = list(PhiType = "free",
               UserPhi = NULL,
               MaxAbsPhi = .5,
               EigenValPower = 2
              )
   cnPH_Len = 4
   PH_Arg_Names = sort(names(cnPH))
   cnPH[(namPhi <- names(Phi))] <- Phi
  
  
   # ----_____Model Error (ME) Default values ----
   # ModelError   logical: If ModelError = TRUE 
   #                       model error will be introduced into 
   #                       the factor pattern via the method 
   #                       described by Tucker, Koopman, and 
   #                       Linn (TKL). Default( ModelError = FALSE). 
   # NMinorFac    scalar:  Number of minor factors in the TKL model. 
   #                       Default( NMinorFac = 150 ).
   # ModelErrorType        character If ModelErrorType = "U" then apply 
   #                       ModelErrorVar to
   #                       uniqueness variances, otherwhise apply to 
   #                       total variances
   # ModelErrorVar scalar  [0,1]: The proportion of uniqueness 
   #                       variance that is due to model error.
   #                       Default( ModelErrorVar = .10 ).
   # epsTKL        scalar[0,1]:  Controls the size of the 
   #                       factor loadings in sucessive
   #                       minor factors.
   #                       Default( epsTKL = .20 ).
   # RSpecific    Optional correlation matrix for specific factors
   #                       Default = NULL.
   
  
   cnME = list(     ModelError = FALSE,
                    NMinorFac = 150,
                    ModelErrorType = "U",
                    ModelErrorVar = .10,
                    epsTKL = .20,
                    RSpecific = NULL) 
   cnME_Len = 5
   ME_Arg_Names = sort(names(cnME))
   cnME[(namME <- names(ModelError))] <- ModelError
 
   if(!is.null(cnME$RSpecific)){
      if(!is.matrix(cnME$RSpecific)){
        # 2DO January 21, 2019  
        # ADD CODE TO TEST FOR PD R MATRIX  
        stop("\nRSpecific must be NULL or a positive definite correation matrix\n")
      } 
    }
  
  # 2DO CORRELATED RESIDUALS AND MODEL ERROR NOT CURRENTLY ALLOWED
  if( !is.null(cnME$RSpecific) & isTRUE(cnME$ModelError) ){
    stop("\nCurrently you can not specify correlated residuals and model error.")
  }
  
  
   # ----_____Bifactor (BF) Models Default values ----
   # Bifactor     logical:  If Bifactor TRUE parameters 
   #                      for the bifactor model will be generated.
   #                      Default( Bifactor = FALSE ).
   # Hierarchical logical: If TRUE then a model
   #                      consistent with a hierarchical Schmid Leiman
   #                      bifactor model will be generated.
   #                      Default( Hierarchical = FALSE ). 
   # F1FactorRange vector of length 1 or 2: Controls the sizes of the 
   #                      general factor loadings in nonhierarchical 
   #                      bifactor models.
   #                      Default( F1FactorRange = c(.4, .7) ).
   #              If F1FactorDist = "runif"
   #                 Vector of length 2 defines the bounds of the
   #                 uniform distribution,  c(lower, upper);
   #              If F1FactorDist = "rnorm" 
   #                 the vector defines the mean and 
   #                 standard deviation of the 
   #                 the normal distribution from which loadings
   #                 are sampled, c(MN, SD).
   #              If F1FactorDist = "sequential" 
   #                 the vector specifies the lower and upper 
   #                 bound of the loadings sequence,  c(lower, upper).
   # F1FactorDist  character: Specifies the sampling distribution for 
   #                 the general factor loadings. Possible values are
   #                 "runif", "rnorm", "sequential".
   #                Default( F1FactorDist = "sequential" ).
   #
   cnBF = list(  Bifactor = FALSE,
                 Hierarchical = TRUE,
                 F1FactorRange = c(.4, .7),
                F1FactorDist = "sequential")
   cnBF_Len = 4
   BF_Arg_Names = sort(names(cnBF))
   cnBF[(namc <- names(Bifactor))] <- Bifactor
  
   
   # ----_____MonteCarlo (MC) Default values ----
   # NSamples    integer: Defines number of Monte Carlo Samples. 
   #                      Default( NSamples = 0) .
   # SampleSize  integer: Sample size for each Monte Carlo sample.
   #                      Default( SampleSize = 250 )
   # Raw         logical: If Raw = TRUE, simulated data sets will be 
   #                      contain raw data. If Raw = FALSE simulated 
   #                      data sets will contain correlation matrices. 
   #                      Default( Raw = FALSE ).
   # Thresholds list:    List elements contain thresholds for each item.
   #                      
   cnMC = list(NSamples = 0, 
               SampleSize = 250,
               Raw = FALSE,
               Thresholds = NULL)
   cnMC_Len = 4
   MC_Arg_Names = sort(names(cnMC))
   cnMC[(namMC <- names(MonteCarlo))] <- MonteCarlo
  
   # Initialize MC return values
   MCData = NULL            
   MCDataME = NULL
  
  
 
   # ----_____FactorScores (FS) Default values ----
   # FS         logical:  If FS = TRUE factor scores will be simulated.
   #                      Default( FS = FALSE ).
   # CFSeed   Optional seed for the common factors
   #
   # SFSeed   Optional seed for the specific factors
   #
   # EFSeed   Optional seed for the error factors
   #
   # VarRel  Vector of manifest variable reliabilities
   #         By Default VarRel = H^2 values (communalities)
   #
   # Population logical:  If Population = TRUE, the factor scores will
   #                      fit the correlational constraints of the factor
   #                      model exactly (e.g., the common factors will be
   #                      orthogonal to the unique factors). 
   #                      Default( Population = FALSE ).
   # NFacScores scalar:   Sample size for the factor scores.
   #                      Default( NFacScores = 250 ).
   # Thresholds vector:   A list of quantiles used to polychotomize the
   #                      observed data to produce Likert scores.
   #
   cnFS = list(FS = FALSE, 
               CFSeed = NULL,
               SFSeed = NULL,
               EFSeed = NULL,
               VarRel = NULL,
               Population = FALSE,
               NFacScores = 250,
               Thresholds = NULL)
   cnFS_Len = 8
   FS_Arg_Names = sort(names(cnFS))
   cnFS[(namFS <- names(FactorScores))] <- FactorScores
   
   # If user supplies non integer value for NFacScores then round to
   # nearest integer
   if(  round(cnFS$NFacScores) !=  cnFS$NFacScores){
     cat("\nWarning: Non integer value for NFacScores. \nNFacScores was rounded to the nearest integer.\n")
     cnFS$NFacScores <- round(cnFS$NFacScores)
   }
     

  
   #CHECK: that the correct number of thresholds 
   # has been specified
   is.scalar <- function(x) is.atomic(x) && length(x) == 1L
    if( !is.null(cnFS$Thresholds) ){
    if( is.scalar(cnMD$NItemPerFac) ){
       if( length(cnFS$Thresholds) != (cnMD$NFac * cnMD$NItemPerFac)){
         stop ("\n\nFATAL ERROR: You seem to have specified the wrong number of Thresholds")
       } #END if( length(cnFS$Thresholds)
    } #END if ( is.scalar(cnMD$NItemPerFac)
    
    if( !is.scalar(cnMD$NItemPerFac) ){
       if( length(cnFS$Thresholds) != sum(cnMD$NItemPerFac) ){
         stop ("\n\nFATAL ERROR: You seem to have specified the wrong number of Thresholds")
       }
    } #END if (!is.scalar(cnMD$NItemPerFac )
  } #END if( !is.null(cnFS$Thresholds) ) 
  
  
   # CHECK if FS seeds are unique 
   if(  !is.null(cnFS$CFSeed) &&  !is.null(cnFS$UFSeed) ){
      if( cnFS$CFSeed == cnFS$UFSeed ){
        stop("\n\n *** FATAL ERROR: CFSeed must differ from UFSeed ***\n" )
      } 
   }

   
   # ----_____Missing (MS) Default values ----
   # Missing  logical
   # Mechanism character
   # MSProb   scalar or vector of length NVar
   cnMS = list(Missing = FALSE, 
               Mechanism = "MCAR",
               MSProb = .10)
   cnMS_Len = 3
   MS_Arg_Names = sort(names(cnMS))
   cnMS[(namMS <- names(Missing))] <- Missing
  
   
   # ----_____Control (CT) Default values ----
   # Maxh2       scalar:  maximum allowable communality (Default = .98)
   # Reflect    logical: If TRUE factors will be randomally 
   #                     reflected
   cnCT = list(Maxh2 = .98,
               Reflect = FALSE)
   cnCT_Len = 2
   CT_Arg_Names = sort(names(cnCT))
   cnCT[(namCT <- names(Control))] <- Control
  
   # ---- Initialize Variables ---- #
   RpopME  = CovMajor = CovMinor = CovUnique = discrim  =
     FacInd = FacIndME =  NULL
  
   ## Number of Variables
   # NItemPerFac is either a scalar (same number items 
   # for all factors) or a vector of length NFac.  
   # If scalar then create vector 
   if (length(cnMD$NItemPerFac) == 1) { 
     cnMD$NItemPerFac <- rep(cnMD$NItemPerFac, cnMD$NFac)
   } # END if (length(cnMD$NItemPerFac) == 1){ 
   NVar <- sum(cnMD$NItemPerFac)
  
   if( length(cnMS$MSProb) == 1){
     cnMS$MSProb <- rep(cnMS$MSProb, NVar) 
   }
  
  
   # Initialize Phi for orthogonal models
   if(cnMD$Model == "orthogonal"){
      Phi <- diag(cnMD$NFac)
   }   
   
   
   
   # ----_____Set Seed ---- 
   ## Generate random seed if not supplied

   if (is.null(Seed)) { 
      Seed <- as.integer((as.double(Sys.time()) * 1000 + Sys.getpid()) %% 2^31) 
   } # END if (is.null(Seed)) { 
   set.seed(Seed)
   
   
#----------------END SEC 1 ---------------------#
 
  
  #  ----_____Record arguments call ----
  cl <- match.call()
  
  

 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                ---- SEC 2: Begin Error Checking ----
# Error checking for each keyword
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

   ## ----_____Check Model Arguments----
   ##     Check that cnMD$Model[1] is specified as either 
   ##     orthogonal or oblique
   if (cnMD$Model[1] %in% c("oblique", "orthogonal") == FALSE) {
       stop("The first element of 'cnMD$Model' must be specified as either 'oblique' or 'orthogonal'")
   } # END if (cnMD$Model[1] %in% c("oblique", "orthogonal") == FALSE) {
    
   # Check for valid arguments if supplied
   if(! length(names(Model)) == 0 ){
   # is number of args less than max number of args
   errorLength <- length(names(Model)) <= cnMD_Len
   # are all arg names valid
   errorName <- min((sort(names(Model)) %in%  MD_Arg_Names))
     if(!errorLength | !errorName){
         options( show.error.messages=FALSE) 
        stop(cat("\n***FATAL ERROR***\nUnrecognized argument in 'Model'.",
        "\nPlease see the help page for a list of valid arguments.\n")) 
      }
   } #END if(! length(names(Model)) == 0 )  
   # END Check Model
  
  
   ## ----_____Check Loadings Arguments----
   ## If cnLD$FacLoadDist argument is mis-specified, stop
     if (cnLD$FacLoadDist %in% c("sequential", "runif", "rnorm", "fixed") == FALSE) {
        stop("Argument 'cnLD$FacLoadDist' is not properly specified")
     }# END if (cnLD$FacLoadDist %in% c("sequential", "runif", "rnorm", "fixed") == FALSE) {
  
     # if only one fixed loading is specified repeat for each factor
     if(cnLD$FacLoadDist == "fixed"){
        if( (cnMD$NFac > 1) && (length(cnLD$FacLoadRange) == 1) ){
        cnLD$FacLoadRange <- rep(cnLD$FacLoadRange, cnMD$NFac )
        }
     }#END if(cnLD$FacLoadDist == "fixed")
    
    # Check for valid arguments if supplied
    if(! length(names(Loadings)) == 0 ){
      # is number of args less than max number of args
      errorLength <- length(names(Loadings)) <= cnLD_Len
      # are all arg names valid
      errorName <- min((sort(names(Loadings)) %in%  LD_Arg_Names))
      if(!errorLength | !errorName){
        options( show.error.messages=FALSE) 
        stop(cat("\n***FATAL ERROR***\nUnrecognized argument in 'Loadings'.",
                 "\nPlease see the help page for a list of valid arguments.\n")) 
      }
    }  
   # END Check Loadings
    
  
   ## ----_____Check CrossLoadings Arguments----
     if (length(cnCL$CrossLoadRange) != 2) {
       stop("Argument 'CrossLoadRange' is not properly specified")
     }# END if (length(CrossLoadRange) != 2) 
  
     if (cnCL$CrossLoadRange[1] < 0 | cnCL$CrossLoadRange[2] > 1) {
       stop("Argument 'CrossLoadRange' is out of bounds")
     } # END if (CrossLoadRange[1] < 0 | CrossLoadRange[2] > 1) 
  
     if (cnCL$CrudFactor < 0 | cnCL$CrudFactor > .25) {
       stop("Argument 'CrudFactor' is out of bounds (maximum CrudFactor = .25")
     } # END if (cnCL$CrudFactor < 0 | cnCL$CrudFactor > .25) 
 
     # Check for valid arguments if supplied
     if(! length(names(CrossLoadings)) == 0 ){
        # is number of args less than max number of args
        errorLength <- length(names(CrossLoadings)) <= cnCL_Len
        # are all arg names valid
        errorName <- min((sort(names(CrossLoadings)) %in%  CL_Arg_Names))
       if(!errorLength | !errorName){
        options( show.error.messages=FALSE) 
       stop(cat("\n***FATAL ERROR***\nUnrecognized argument in 'CrossLoadings'.",
                "\nPlease see the help page for a list of valid arguments.\n")) 
       }
     }  
    #END Check CrossLoadings
  

   ## ----_____Check Phi Arguments----
   ##  PhiType  
    if (cnPH$PhiType %in% c("free", "fixed", "user") == FALSE) {
       stop("Argument 'PhiType' is not properly specified")
    }# END if (cnPH$PhiType %in% c("Free", "Fixed", "User") == FALSE)
     
   # Check for valid arguments if supplied
    if(! length(names(Phi)) == 0 ){
      # is number of args less than max number of args
      errorLength <- length(names(Phi)) <= cnPH_Len
      # are all arg names valid
      errorName <- min((sort(names(Phi)) %in%  PH_Arg_Names))
      if(!errorLength | !errorName){
        options( show.error.messages=FALSE) 
        stop(cat("\n***FATAL ERROR***\nUnrecognized argument in 'Phi'.",
                 "\nPlease see the help page for a list of valid arguments.\n")) 
    }
   }  
   # END Check Phi
    
    
  ## ----_____Check ModelError Arguments----
    # Check for valid arguments if supplied
    if(! length(names(ModelError)) == 0 ){
      # is number of args less than max number of args
      errorLength <- length(names(ModelError)) <= cnME_Len
      # are all arg names valid
      errorName <- min((sort(names(ModelError)) %in%  ME_Arg_Names))
      if(!errorLength | !errorName){
        options( show.error.messages=FALSE) 
        stop(cat("\n***FATAL ERROR***\nUnrecognized argument in 'ModelError'.",
                 "\nPlease see the help page for a list of valid arguments.\n")) 
      }
    }  
   # END Check ModelError 
    
    
  
   ## ----_____Check Bifactor Arguments----
   # Check for valid arguments if supplied
    if(! length(names(Bifactor)) == 0 ){
      # is number of args less than max number of args
      errorLength <- length(names(Bifactor)) <= cnBF_Len
      # are all arg names valid
      errorName <- min((sort(names(Bifactor)) %in%  BF_Arg_Names))
      if(!errorLength | !errorName){
        options( show.error.messages=FALSE) 
        stop(cat("\n***FATAL ERROR***\nUnrecognized argument in 'Bifactor'.",
                 "\nPlease see the help page for a list of valid arguments.\n")) 
      }
    } 
    # END Check Bifactor
    
    
   ## ----_____Check MonteCarlo Arguments----
   # Check for valid arguments if supplied
    if(! length(names(MonteCarlo)) == 0 ){
      # is number of args less than max number of args
      errorLength <- length(names(MonteCarlo)) <= cnMC_Len
      # are all arg names valid
      errorName <- min((sort(names(MonteCarlo)) %in%  MC_Arg_Names))
      if(!errorLength | !errorName){
        options( show.error.messages=FALSE) 
        stop(cat("\n***FATAL ERROR***\nUnrecognized argument in 'MonteCarlo'.",
                 "\nPlease see the help page for a list of valid arguments.\n")) 
      }
    }  
    # END Check MonteCarlo
  
   
   ## ----_____Check FactorScores Arguments----
   # Check for valid arguments if supplied
   if(! length(names(FactorScores)) == 0 ){
     # is number of args less than max number of args
     errorLength <- length(names(FactorScores)) <= cnFS_Len
     # are all arg names valid
     errorName <- min((sort(names(FactorScores)) %in%  FS_Arg_Names))
     if(!errorLength | !errorName){
      options( show.error.messages=FALSE) 
      stop(cat("\n***FATAL ERROR***\nUnrecognized argument in 'FactorScores'.",
               "\nPlease see the help page for a list of valid arguments.\n")) 
     }
     
     if(isTRUE(cnFS$FS) && cnFS$NFacScores <= NVar && isTRUE(cnFS$Population) ){
       cat("\n\n*** Fatal Error: NFacScores < number of observed variables. ***
           \n*** The (full) population factor score matrix cannot have full column rank. ***\n
           \n*** Set Population = FALSE ***\n")
       stop()
     }
   } 
   #END Check FactorScores 
  
  
   ## ----_____Check Missing Arguments----
   # Check for valid arguments if supplied
   if(! length(names(Missing)) == 0 ){
    # is number of args less than max number of args
    errorLength <- length(names(Missing)) <= cnMS_Len
    # are all arg names valid
    errorName <- min((sort(names(Missing)) %in%  MS_Arg_Names))
    if(!errorLength | !errorName){
      options( show.error.messages=FALSE) 
      stop(cat("\n***FATAL ERROR***\nUnrecognized argument in 'Missing'.",
               "\nPlease see the help page for a list of valid arguments.\n")) 
    }
   }
   # END Check Missing
  
  
  ## ----_____Check Control Arguments----
  ## Constrain max communality to be lower than 1, give warning, 
  ## then change to default
  if (cnCT$Maxh2 > 1) {
    warning("Maximum communality value ('Maxh2' argument) is greater than 1, reverted to the default of .98")
    cnCT$Maxh2 <- .98
  } # END if (Maxh2 > 1) {
  # END Check Control
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
  ##           End SEC 2: Error Checking          ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
  
 
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#          ---- SEC 3: Define Internal Functions ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
    ##---------------------------------------------#
    #  corSample
    #  Generate sample correlation matrices
    #  see Browne 1968
    #  Browne, M. (1968). A comparison of factor analytic techniques. 
    #  Psychometrika, 33(3):267-334.
    #
    # Output
    #    a sampled correlation matrix
    #----------------------------------------------#
  
    corSample<- function(R,n, seed){
      set.seed(seed)
      Nvar<-ncol(R)
      Tmat<-matrix(0,Nvar,Nvar)
      Tmat[lower.tri(Tmat)]<-rnorm(n=(Nvar*(Nvar-1)/2))
      for(i in 1:Nvar){
        # Note that Browne 1968 contains a typo for the df -- 
        # see the following for (n-i+1)
        # Kshirsagar, A. (1959). Bartlett decomposition and 
        # Wishart distribution.The Annals of Mathematical Statistics, 
        # 30(1) 239-241.
        Tmat[i,i]<-sqrt(rchisq(n=1,df=(n-i+1)))
      }
    
      H<- Tmat %*% t(Tmat)
      Omega <-t(chol(R))
      A <-  Omega %*% H %*% t(Omega)
      S <- (1/n) * A
      Dmat<-diag(1/sqrt(diag(A)))
      R.samp<-Dmat%*%A%*%Dmat
      R.samp <- .5*(R.samp + t(R.samp))
      R.samp #return value
    } #END  corSample
  
  
  
    #-----------------------------------------------------------------#
    # rGivens: Generate R matrices with user-specified eigenvalues via
    # Givens rotations based on theory described in:
    #  References:
    #   [1] R. B. Bendel and M. R. Mickey, Population correlation matrices
    #      for sampling experiments, Commun. Statist. Simulation Comput.,
    #      B7 (1978), pp. 163-182.
    #   [2] P. I. Davies and N. J. Higham, Numerically stable generation of
    #       correlation matrices and their factors, BIT, 40 (2000), pp. 640-651.
    #-----------------------------------------------------------------#
    
    rGivens <- function(eigs, Seed = NULL){
    
      eps<-1e-12
      n <- length(eigs)
    
      if(abs(sum(eigs)-n)/n > 100*eps){
       stop("Elements of eigs must sum to n,\n")
      }
    
      A = diag(eigs);
      
      ## generate a random orthogonal matrix
      if (is.null(Seed)) 
        Seed <- as.integer((as.double(Sys.time()) * 1000 + Sys.getpid())%%2^31)
      set.seed(Seed)
    
      M <- matrix(rnorm(n * n), nrow = n, ncol = n)
      Q <- qr.Q(qr(M))
      S0 <- A <- Q %*% A %*% t(Q);  # Not exploiting symmetry here.
    
     a = diag(A);
     y = which(a < 1);
     z = which(a > 1);
    
     while( length(y) > 0 && length(z) > 0){
       
      i = y[ceiling(runif(1)*length(y))];
      j = z[ceiling(runif(1)*length(z))];
      
      if( i > j){
        temp = i; 
        i = j; 
        j = temp
      }
      
      alpha = sqrt(A[i,j]^2 - (a[i]-1)*(a[j]-1));
      
      t<-c(0,0)
      t[1] = (A[i,j] + sign(A[i,j])*alpha)/(a[j]-1);
      
      t[2] = (a[i]-1)/((a[j]-1)*t[1]);
      
      t = t[ceiling(runif(1)*2)];  # Choose randomly from the two roots.
      
      c = 1/sqrt(1 + t^2);
      s = t*c;
      
      A[, c(i,j)] =  A[, c(i,j)] %*% matrix(c(c, s, -s, c),2,2, byrow=TRUE);
      A[c(i,j),] = matrix(c(c, -s, s,c), 2, 2, byrow=TRUE) %*% A[c(i,j), ];
      
      # Ensure (i,i) element is exactly 1.
      A[i,i] = 1;
      
      a = diag(A);
      y = which(a < 1);
      z = which(a > 1);
      
     } ## End while loop
    
     # Setlast diagonal to 1:
     diag(A) = 1
    
     A = (A + t(A))/2; # Force symmetry.
     convergence <- TRUE
     if(max(abs(A))>1) convergence <- FALSE
    
     list(R = A, convergence = convergence)
    } # END rGivens
   
  
    #-----------------------------------------------------------------#
    # genPhi  fnc to generate random factor correlation matrices 
    #-----------------------------------------------------------------#
    
    genPhi <- function(NFac, 
                      MaxAbsPhi = .5, 
                      EigenValPower = 6,
                      PhiSeed){
    # generate random positive numbers that sum to cnMD$NFac
    
    ## CG Create descending vector of values
    eigs <- sort(abs(rnorm(NFac))^EigenValPower, decreasing=TRUE)
    
    ## CG Convert above values to be eigenvalues (i.e., sum to cnMD$NFac)
    eigs <- eigs * NFac/sum(eigs)
    
    # generate R matrix with eigenvalues in eigs via Givens rotations
    Phi <- rGivens(eigs, Seed = PhiSeed)$R
    
    # Rescale to bound Phi_ij 
    # make hollow matrix
    Phi <- Phi - diag(NFac)
    maxPhi <- max(abs(Phi))
    S <- diag(NFac) * sqrt(cnPH$MaxAbsPhi)/sqrt(maxPhi)
    Phi <- S %*% Phi %*% S
    diag(Phi) <- 1
    Phi <- round(Phi, 2)
    Phi
   } ## End genPhi
  
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
   # END: Define Internal Functions #
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#  
#         ---- SEC 4: Fncs to Generate Common Factor Loadings ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#   
  
   cnMD$NItemPerFacIndex <- 0
  
   # sequential loadings
   if (cnLD$FacLoadDist == "sequential") {
       FncGenLoadings <- function(){
         cnMD$NItemPerFacIndex <<- cnMD$NItemPerFacIndex + 1
                 seq(from   = cnLD$FacLoadRange[1],
                     to     = cnLD$FacLoadRange[2],
                     length = cnMD$NItemPerFac[cnMD$NItemPerFacIndex])
       }
   } #End if (cnLD$FacLoadDist == "sequential") 
  
   # uniform loadings
   if (cnLD$FacLoadDist == "runif") {
      FncGenLoadings <- function(){
        cnMD$NItemPerFacIndex <<- cnMD$NItemPerFacIndex + 1
                runif(n   = cnMD$NItemPerFac[cnMD$NItemPerFacIndex],
                      min = min(cnLD$FacLoadRange),
                      max = max(cnLD$FacLoadRange))
      }
   } #END if (cnLD$FacLoadDist == "runif") 
  
   # normally distributed loadings
   if (cnLD$FacLoadDist == "rnorm") {
     FncGenLoadings <- function(){
       cnMD$NItemPerFacIndex <<- cnMD$NItemPerFacIndex + 1
                rnorm(n    = cnMD$NItemPerFac[cnMD$NItemPerFacIndex],
                      mean = cnLD$FacLoadRange[1],
                      sd   = cnLD$FacLoadRange[2])
     }
   } # END if (cnLD$FacLoadDist == "rnorm") 
    
  # fixed loadings
  fixIndex <- 0  # initialize
  if (cnLD$FacLoadDist == "fixed") {
       FncGenLoadings <- function(){
        cnMD$NItemPerFacIndex <<- cnMD$NItemPerFacIndex + 1   
        fixIndex <<- fixIndex + 1 
        rep(cnLD$FacLoadRange[fixIndex], cnMD$NItemPerFac[cnMD$NItemPerFacIndex])
       }   
  } #END if (cnLD$FacLoadDist == "fixed") 
  
  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ##             ---- Load User-Supplied Factor Loadings ----
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # ---- Loadings Matrix Supplied by User
  if( !is.null(cnLD$FacPattern) ){
    NVar <- nrow(cnLD$FacPattern)
    cnMD$NFac <- NFac <- ncol(cnLD$FacPattern)
    Fl <- cnLD$FacPattern
    
    # initialize Phi for orthogonal model
    # if not supplied by user
    if(cnPH$PhiType != "user"){
        Phi <- diag(NFac)
    }
    if(cnPH$PhiType == "user"){
      Phi <- cnPH$UserPhi
    }
    
    # Add column labels
    colnames(Fl) <- paste0("F", 1:ncol(Fl))
    rownames(Fl) <- paste0("V", 1:NVar)
  }#if( !is.null(cnLD$FacPattern) )
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#  
  #         ---- SEC 5: Generate Common Factor Loadings ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 
 
  # Add crossloadings to user supplied factor pattern
   if( !is.null(cnLD$FacPattern) ){
       if ( (cnCL$ProbCrossLoad == 0) & !is.null(cnCL$CrossLoadPositions) ){
          if( nrow(cnCL$CrossLoadPositions) != length(cnCL$CrossLoadValues)){
          stop("\nInsufficient number of cross loadings supplied in CrossLoadValues\n")
          }
       }
     
     Fl[cnCL$CrossLoadPositions] <- cnCL$CrossLoadValues
   }# END if( !is.null(cnLD$FacPattern) ){   
 
  if( is.null(cnLD$FacPattern) ){
     # This code allows different number of variables
     # per factor
     FLoadVec <- 999
      for(j in 1:cnMD$NFac){
        FLoadVec <- c(FLoadVec, FncGenLoadings(), rep(0, NVar))
      }
     FLoadVec <- FLoadVec[-1]
     FLoadVec <- FLoadVec[1: (NVar * cnMD$NFac)]
     Fl <- matrix(FLoadVec, nrow = NVar, ncol = cnMD$NFac)
  
     # FsZeros is a logical matrix that flags elements with non primary
     # loadings
     FsZeros <- Fl == 0
  
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
    ##            ----_____Add Cross Loadings ---- 
     
     if ( (cnCL$ProbCrossLoad == 0) & !is.null(cnCL$CrossLoadPositions) ){
       if( nrow(cnCL$CrossLoadPositions) != length(cnCL$CrossLoadValues)){
         stop("\nInsufficient number of cross loadings supplied in CrossLoadValues\n")
       }
        
       Fl[cnCL$CrossLoadPositions] <- cnCL$CrossLoadValues
     }     
     
  
    if ( (cnCL$ProbCrossLoad > 0) & (cnCL$ProbCrossLoad < 1) ) {
      for (i in 1:NVar) {
          for (j in 1:cnMD$NFac) {
            
             if (FsZeros[i, j]) {
                  CoinFlip <- sample(c(1, 0), size = 1, 
                  prob = c(cnCL$ProbCrossLoad, 1 - cnCL$ProbCrossLoad))
                  
                if (CoinFlip) { #If heads, add cross loading
                    Fl[i, j] <- runif(1, min(cnCL$CrossLoadRange), max(cnCL$CrossLoadRange))
                 }
             }# End if (FsZeros[i, j])
         } #End (j in 1:cnMD$NFac) 
      }#End for (i in 1:NVar)
    }#End if (cnCL$ProbCrossLoad < 1)  
  
    ## CG If ProbCrossLoad == 1, add 1 cross loading per factor
    if (cnCL$ProbCrossLoad == 1) {
       Fl[1, cnMD$NFac] <- cnCL$CrossLoadRange[1]
       for (j in 1:(cnMD$NFac - 1)){
          Row <- sum(cnMD$NItemPerFac[1:j])+1
          Col <- j
          Fl[Row, Col] <- cnCL$CrossLoadRange[1]
       }#END for j in 1:(cnMD$NFac - 1)
    }#END if (cnCL$ProbCrossLoad == 1)
     
     
     
     
     ## --- _____Add Crud Factor ----
     if (cnCL$CrudFactor > 0) { 
        numZeroLoadings <- sum(Fl == 0)
        Fl[Fl == 0] <- runif(numZeroLoadings, -cnCL$CrudFactor, cnCL$CrudFactor)
     } #End Add Crud Factor
     
  
    # Randomly reflect factor orientations
    if (cnCT$Reflect) {
      Fl <- Fl %*% diag(sample(c(-1, 1), size = cnMD$NFac, replace = TRUE))
    }
  
    # Add column labels
    colnames(Fl) <- paste0("F", 1:ncol(Fl))
    rownames(Fl) <- paste0("V", 1:NVar)
    
  }#END if( is.null(cnLD$FacPattern) )
  
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#  
#              ---- SEC 6: Create Oblique Model ----  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#  
  
  # cnMD$Model specifies either an 
  #  "orthogonal" or "oblique" factor.  
  # If cnMD$Model == "oblique" then Phi will not be
  # an identy matrix.  If cnPH$Free = FALSE  then 
  # all off diagonal elements of Phi will equal
  # the value specified in MaxAbsPhi.  If 
  # cnPH$Free = TRUE then a random Phi matrix will be generated.
  
 
  
  if (cnMD$Model[1] == "oblique") {
  
  #---- _____Generate fixed Phi  ----  
  # ----_____PhiType == "Fixed" ----  
    # Elements of Phi are fixed to user defined constant
    if (cnPH$PhiType == "fixed") {
      
      Phi <- matrix(cnPH$MaxAbsPhi, cnMD$NFac, cnMD$NFac)
      
      diag(Phi) <- 1
   
        # Check that Phi is PSD
        if( eigen(Phi)$values[cnMD$NFac] < 0){
           stop("\n *** Fatal Error: User Defined Phi is not PSD *** \n")
        } 
    }

    
  #----_____Generate random Phi with desired constraints ----
  #  ----_____PhiType == "Free" ----
    if (cnPH$PhiType == "free") {
      NPDFlag <- TRUE
      
      # loop until PSD Phi is generated
      iPhiSeed <- Seed
        
      while (NPDFlag) {
        Phi <- genPhi(NFac = cnMD$NFac, 
                     MaxAbsPhi = cnPH$MaxAbsPhi, 
                    EigenValPower = cnPH$EigenValPower,
                    PhiSeed = iPhiSeed)
        if (min(eigen(Phi, symmetric = TRUE)$values) > 0) NPDFlag <- FALSE
        iPhiSeed <- iPhiSeed + 1
     } #End While loop 
    } # End generate random Phi
    
    #----_____PhiType == "user"  ----
    if (cnPH$PhiType == "user") {
      Phi <- cnPH$UserPhi
    }  
    
  }# End if cnMD$Model == "oblique"  
  
  #---- When Model defaults to "orthogonal but 
  #     PhiType == "user"  reset model to "oblique"
  if (cnMD$Model == "orthogonal" && cnPH$PhiType == "user") {
    Phi <- cnPH$UserPhi
    cnMD$Model == "oblique"
  } 
  
  
  # ----_____Compute Communalities----
  h2 <- diag(Fl %*% Phi %*% t(Fl))
  
  #Check for communalities ge Maxh2 and rescale loadings if necessary
  if (max(h2) >= cnCT$Maxh2) {
    s <- sqrt(cnCT$Maxh2)/sqrt(h2)
    # if communality i less than Maxh2 do not rescale
    s[h2 <= cnCT$Maxh2] <- 1
    Fl <- diag(s) %*% Fl  #rescale rows of Fl
    h2 <- diag(Fl %*% Phi %*% t(Fl))
  } # END if (max(h) >= cnCT$Maxh2) {
  
  
  #~~~~~~~~~~~~ Generate Rpop~~~~~~~~~~~~~~~~
  # ----Generate Rpop ----
    # No specific factor correlations
  if( is.null(cnME$RSpecific) ){
      Rpop <- Fl %*% Phi %*% t(Fl)
      diag(Rpop) <- 1
  } 
  

  # Specific factors allowed to correlate
  if( !is.null(cnME$RSpecific) ){
 
    # _____**2TEST add R resid #1** ---- 
    SFSstndevs <- diag(sqrt(cnFS$VarRel - h2))
    CovSFS <-   SFSstndevs %*% cnME$RSpecific %*% SFSstndevs
    diag(CovSFS) <- 0
    Rpop <- Fl %*% Phi %*% t(Fl) + CovSFS
    if( max( abs(Rpop)) > 1) stop("\nNon valid entry for RSpecific (|r_ij| > 1")
    diag(Rpop) <- 1
  }  
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#  
#       ---- SEC 7: Create Bifactor Model ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#  
    
  ## Return null unless creating a hierarchical solution
  FlHierarchical <- PhiHierarchical <- NULL
  
  
   if(cnBF$Bifactor == TRUE) { 
      # ----_____Non Hierarchical Bifactor Model ----
      if (cnBF$Hierarchical == FALSE) {
         # If cnMD$Model = Bifactor then add loadings on a general factor  
      
          if (cnBF$F1FactorDist == "sequential") { 
             f1 <- seq(cnBF$F1FactorRange[1], 
                   cnBF$F1FactorRange[2],
                   length = NVar)
             Fl <- cbind(f1, Fl)
          } # END if (F1FactorDist == "sequential") { 
      
          if (cnBF$F1FactorDist == "runif") { 
             f1 <- runif(NVar, min = min(cnBF$F1FactorRange), 
                          max = max(cnBF$F1FactorRange))
             Fl <- cbind(f1, Fl)
          } # END if (F1FactorDist == "runif") { 
      
          if (cnBF$F1FactorDist == "rnorm") { 
             f1 <- rnorm(NVar, mean = cnBF$F1FactorRange[1], 
                          sd = cnBF$F1FactorRange[2])
             Fl <- cbind(f1, Fl)
          } # END if (F1FactorDist == "rnorm") { 
      
          if (cnBF$F1FactorDist == "fixed") {
             f1 <- rep(cnBF$F1FactorRange[1], NVar)
             Fl <- cbind(f1, Fl)
          } # END if (F1FactorDist == "fixed") {
      
          colnames(Fl)[1] <- "G"
          
          
          #If user-supplied communalities are provided, then scale loadings
          if(!is.null(cnLD$h2)){
             h2 <- apply(Fl^2, 1, sum)
             h2ScalingConstants <- sqrt(cnLD$h2/h2)
             Fl <- diag( h2ScalingConstants) %*% Fl
             h2 <- apply(Fl^2, 1, sum)
          }    
          
        
         # ----_____Gen Model implied R for Nonhier BF Model----
 
        
         # No specific factor correlations
         if( is.null(cnME$RSpecific) ){
           Rpop <- Fl %*%  t(Fl)
           diag(Rpop) <- 1
         } #END if( is.null(cnME$RSpecific) ){
         
         # ----_____Specific factors allowed to correlate----
         if( !is.null(cnME$RSpecific) ){
            # _____**2TEST Add R Resid #2 **  ---- 
           SFSstndevs <- diag(sqrt(cnFS$VarRel - h2))
           CovSFS <-   SFSstndevs %*% cnME$RSpecific  %*% SFSstndevs
           diag(CovSFS) <- 0
           Rpop <- Fl %*% t(Fl) + CovSFS
           if( max( abs(Rpop)) > 1) stop("\nNon valid entry for RSpecific (|r_ij| > 1")
           diag(Rpop) <- 1
         }#END if( !is.null(cnME$RSpecific) ){  
         
        
     } # END IF HIERARCHICAL == FALSE
    
    
    # ----_____Hierarchical Bifactor Model ----
    if (cnBF$Hierarchical == TRUE) {
      # Step 1: Let Lambda1 = Fl
      Lambda1 <- Fl
      
        if(cnPH$PhiType != "user"){
          
          if(cnBF$F1FactorDist == "runif"){
             # Step 2: Generate Phi from random loadings on 
             # general factor
             Lambda2 <- runif(ncol(Lambda1), 
                       min = min(cnBF$F1FactorRange), 
                       max = max(cnBF$F1FactorRange))
      
             # at this point Phi is not a cor matrix
             Phi <- Lambda2 %*% t(Lambda2)
      
             # Step 3 Compute sqrt of uniqueness for hier order cnMD$Model
             Psi2 <- diag(sqrt(1 - diag(Phi)))
          } #END  if(cnBFF1FactorDist == "runif")
          
          
          
          if(cnBF$F1FactorDist == "sequential"){
             # gload: loadings of the lower factors on the gen
             gload <- Lambda2 <- seq(cnBF$F1FactorRange[1],
                                     cnBF$F1FactorRange[2],
                                     length = ncol(Fl))
             # at this point Phi is not a cor matrix
             Phi <- gload %*% t(gload)
             
             # Step 3 Compute sqrt of uniqueness for hier order cnMD$Model
             Psi2 <- diag(sqrt(1 - diag(Phi)))
          }   
      
          
          
          if(cnBF$F1FactorDist == "fixed"){
            # gload: loadings of the lower factors on the gen
            gload <- Lambda2 <- seq(cnBF$F1FactorRange[1],
                                    cnBF$F1FactorRange[1],
                                    length = ncol(Fl))
            # at this point Phi is not a cor matrix
            Phi <- gload %*% t(gload)
            
            # Step 3 Compute sqrt of uniqueness for hier order cnMD$Model
            Psi2 <- diag(sqrt(1 - diag(Phi)))
          }  
          
          # make Phi a corr matrix
          diag(Phi) <- 1
          
          
        } #END if(cnPH$PhiType !="user")
      
      
      ## **2PROCESS ADD CODE FOR !="USER AND FIXED F1FACDIST**-----
      
      
      
        if(cnPH$PhiType == "user"){
           # calculate loadings on the general factor
           Lambda2 <- fungible::fals(R = Phi, nfactors = 1)$loadings
         
           PhiFromLambda2 <- Lambda2 %*% t(Lambda2)
         
           # Compute sqrt of uniqueness for hier order cnMD$Model
           Psi2 <- diag(sqrt(1 - diag( PhiFromLambda2)))
         }# END if(cnPH$PhiType =="user")
      
   
      PhiHierarchical <- Phi
      # Compute reproduced R maxtrix and record communalities
      Rhat <- Lambda1 %*% PhiHierarchical %*% t(Lambda1)
      h2 <- diag(Rhat)
      
      # _____Step 4 Scale Lambda1 to insure max communality = Maxh2----
      s <- sqrt(cnCT$Maxh2) / sqrt(h2)
      s[h2 < cnCT$Maxh2] <- 1
      Lambda1 <- diag(s) %*% Lambda1
      
      # Model implied R for Hierarchical Bifactor cnMD$Model
      # No specific factor correlations
      if( is.null(cnME$RSpecific) ){
        Rpop <- Lambda1 %*% PhiHierarchical %*% t(Lambda1)
      } #END if( is.null(cnME$RSpecific) ){
      
      # Specific factors allowed to correlate
      if( !is.null(cnME$RSpecific) ){
        # _____**TEST add R resid #3** ----  
        SFSstndevs <- diag(sqrt(cnFS$VarRel - h2))
        CovSFS <-   SFSstndevs %*% cnME$RSpecific %*% SFSstndevs
        diag(CovSFS) <- 0
        Rpop <- Lambda1 %*% PhiHierarchical %*% t(Lambda1) + CovSFS
        if( max( abs(Rpop)) > 1) stop("\nNon valid entry for RSpecific (|r_ij| > 1")
      }#END if( !is.null(cnME$RSpecific) ){
      
 
      
      # Communalities for hierarchical Bifactor cnMD$Model 
      h2 <- diag(Rpop)
      diag(Rpop) <- 1
      
      
      # Step 5 Create S-L loadings matrix
      b1 <- Lambda1 %*% Lambda2
      B2 <- Lambda1 %*% Psi2
      B <- cbind(b1,B2)
      
      # Hierarchica Bifactor Solution
      colnames(B) <- c("G", paste0("F", 1:cnMD$NFac))
      
      ## Save the 1st order factor loadings
      FlHierarchical <- Fl
      Fl <- B
    } # END if (Hierarchical == TRUE) 
    

      Phi <- diag(ncol(Fl))
      colnames(Phi) <- rownames(Phi) <- c("G", paste0("F", 1:(cnMD$NFac)))
} #END if(cnBF$Bifactor == TRUE) 
  
  if(cnBF$Bifactor == FALSE){
     colnames(Phi) <- rownames(Phi) <- paste0("F", 1:cnMD$NFac)
  }
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#  
  #       ---- SEC 7: IRT PARAMETERS ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#  

  ## IRT scale
  ## Parameters are reported in the logistic metric (D is multipled in the
  ## params)
  
  ## Lookup page numbers in Recakse mult IRT book for 
  ## all relevant equations
     IRTdiscrim = NULL
     IRTdifficulty = NULL
     IRTintercept = NULL
     Dparam <- 1.702
  
     # _____Report loadings in logistic IRT metric-----
     IRTdiscrim <- Dparam * diag(1 / sqrt(1 - diag( Fl %*% Phi %*% t(Fl) ) )) %*% Fl

     
     
     ## _____IRT: Binary Item Models ----
     if( !is.null(cnFS$Thresholds) && length(unlist(cnFS$Thresholds)) == nrow(Fl) ){
     
         # 1 dimensional models 
         if( !is.null(cnFS$Thresholds) && cnMD$NFac == 1 ){
           IRTintercept <- -Dparam * unlist(cnFS$Thresholds)/sqrt(1 - diag( Fl %*% Phi %*% t(Fl) ) )
           IRTdifficulty <- IRTintercept/-IRTdiscrim
         }
    
     
         # > 1 dimensional models only
         if( !is.null(cnFS$Thresholds) && cnMD$NFac > 1 ){
           # multivariate discrimination
           # see Reckase (MIRT book) EQ 4.9
           # Check if this is true with correlated factors
           MVDISC <- sqrt( apply(IRTdiscrim^2, 1, sum) )
       
           IRTintercept <- -Dparam * unlist(cnFS$Thresholds)/sqrt(1 - diag( Fl %*% Phi %*% t(Fl) ) )
           IRTdifficulty <- -IRTintercept/MVDISC
         }## END if( !is.null(cnFS$Thresholds) && cnMD$NFac > 1 )
     } ## END if binary items
  
     
     ## _____IRT: Graded Response Model----
    
     ## _____IRT: GRM 1-Factor----
     if( !is.null(cnFS$Thresholds) && length(unlist(cnFS$Thresholds)) > nrow(Fl) ){
       
       PsiVec <- sqrt(1 - diag( Fl %*% Phi %*% t(Fl) ) ) #uniqueness std devs
       # tau = thresholds 
       tau <- cnFS$Thresholds
       IRTintercept <- IRTdifficulty <- vector("list", length = nrow(Fl))
       for( i in 1:nrow(Fl) ){
         IRTintercept[[i]] <- -Dparam * tau[[i]]/PsiVec[i]
         IRTdifficulty[[i]] <- IRTintercept[[i]]/-IRTdiscrim[i]
       }
    
     } #END IRT GRADED RESPONSE MODEL (1FAC)-
     
     
     # 
     # if(cnBF$Bifactor == FALSE){
     #    min_unique_var <-  min(1 - diag( Fl %*% Phi %*% t(Fl) ) ) 
     #    # in BG models min_unique-var will equal zero
     #    if(min_unique_var > 0){
     #       discrim <- diag(1 / sqrt(1 - diag( Fl %*% Phi %*% t(Fl) ) )) %*% Fl
     #    } 
     # }# END if(cnBF$Bifactor == FALSE)
 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#  
#       ---- SEC 8: Scale communalities to user values ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#  

   if(!is.null(cnLD$h2)){
     if(length(cnLD$h2) != NVar) stop("Commonalities of incorrect length")
     h2ScalingConstants <- sqrt(cnLD$h2/h2)
     Fl <- diag( h2ScalingConstants) %*% Fl 
     # recalculate commonalities
     h2 <- diag (Fl %*% Phi %*% t(Fl))
     
    
     # Generate the Model implied correlation matrix
    
     # No specific factor correlations
     if( is.null(cnME$RSpecific) ){
       Rpop <- Fl %*% Phi %*% t(Fl)
       diag(Rpop) <- 1
     } #END if( is.null(cnME$RSpecific) ){
    
     # Specific factors allowed to correlate
     if( !is.null(cnME$RSpecific) ){
       # _____**TEST add R resid #4** ---- 
       SFSstndevs <- diag(sqrt(cnFS$VarRel - h2))
       CovSFS <-   SFSstndevs %*% cnME$RSpecific %*% SFSstndevs
       diag(CovSFS) <- 0
       Rpop <- Fl %*% Phi %*% t(Fl) + CovSFS
       if( max( abs(Rpop)) > 1) stop("\nNon valid entry for RSpecific (|r_ij| > 1")
       diag(Rpop) <- 1
     }#END if( !is.null(cnME$RSpecific) ){
    
    
  }#END  if(!is.null(cnLD$h2)){
  
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#  
#       ---- SEC 9: Add Model Error ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#      
   # if no model error requested then return regular h2 values
    h2PopME <- h2
   ## Tucker Koopman Linn approach
   if(cnME$ModelError == TRUE){
      # Create matrix of minor factors
      # MacCallum and Tucker 1991 chose
      # 50 minor factors.
      # Briggs and MacCallum chose 150 minor factors
        NMinorFactors <- cnME$NMinorFac
        W <- matrix(0, NVar, NMinorFactors)
     
        # s is a scaling constant
        s <- 1 - cnME$epsTKL
     
       for (i in 0:(NMinorFactors - 1)) {
           W[, i+1] <- rnorm(NVar, 
                         mean = 0, 
                         sd   = s^i)  #TKL EQ 13, p. 431
       }# END for (i in 0:(NMinorFactors - 1)) {
     
       # constrain variance of minor factors to 
       # ModelError
       CovMajor <- Fl %*% Phi %*% t(Fl)
     
       u <- 1 - diag( CovMajor )
     
       wsq <- diag(W %*% t(W))
       
  
       # _____**Check these next two lines**---- 
       ModelErrorVar <- rep(cnME$ModelErrorVar, NVar)
       
       
       # Should ModelErrorVar apply to total or unique variances 
       # "U" (the default) applies to a proprtion of the uniqueness 
       #                   variances (See Briggs and MacCallum 2003 page 32)
       # "V" applies to total variances
       if(cnME$ModelErrorType == "U" ){
         ModelErrorVar <- ModelErrorVar * u
       }   
                                        
       W <- diag(sqrt(ModelErrorVar/wsq)) %*% W
       
     
       # cnMD$Model Error Covariance Matrix for minor factors
       CovMinor <- (W %*% t(W))
       
       # print(W)
       # print(CovMinor)
      
       # Create Pop R with cnMD$Model Error
       RpopME <- Fl %*% Phi %*% t(Fl) +  CovMinor
       
       # These are the true population h2 values
       h2PopME <- diag(RpopME)  
       
       # Cov matrix for Uniqueness factors
       CovUnique <- diag( 1 - diag(RpopME) )
 
       diag(RpopME) <- 1
       ## END: ADD Model ERROR
   }#END if(cnME$ModelError == TRUE)
  
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#  
#       ---- SEC 10: Compute Factor Indeterminancy ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#      
     
     # FacStrc = Factor structure matrix
     # Check if Rpop can be inverted
     EigsRpop <- eigen(Rpop)$values
     if(min(EigsRpop) < 1E-08)
        {
        warning("\n     Rpop is numerically singular\n\n")
        FacInd <- "Rpop not invertible"
        }
     else{
     
        FacStrc <- Fl %*% Phi
        FacInd <- sqrt( diag(  t(FacStrc) %*% solve(Rpop) %*% FacStrc  )  )
     
        FacIndME <- NULL
           if(cnME$ModelError == TRUE){
              FacIndME   <- sqrt( diag(  t(FacStrc) %*% solve(RpopME) %*% FacStrc ) )
           }
     } ## END if else     
 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#  
#       ---- SEC 11: Generate Monte Carlo Samples ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#      
     
    # Preliminary code for Monte Carlo data generation
     
     if(cnMC$NSamples > 0 && cnMC$Raw == FALSE){
        MCData <- MCDataME <- vector("list", length = cnMC$NSamples)
       
        for(i in 1:cnMC$NSamples){
          
         
           MCData[[i]]   <- corSample(Rpop,   
                                      n = cnMC$SampleSize,
                                      seed = i)
         # model error
         if(cnME$ModelError == TRUE){
            MCDataME[[i]] <- corSample(RpopME, 
                                       n = cnMC$SampleSize,
                                       seed = i)
         }#End if(cnME$ModelError == TRUE)
       }#END for(i in 1:cnMC$NSamples)
     }#END if(cnMC$NSamples > 0 && cnMC$Raw == FALSE)
       
     if(cnMC$NSamples > 0 && cnMC$Raw == TRUE){
       
       
       MCData <- MCDataME <- vector("list", length = cnMC$NSamples)
       
       set.seed(Seed)
       
       for(i in 1:cnMC$NSamples){
          
            MCData[[i]] <- MASS::mvrnorm(n = cnMC$SampleSize,
                               mu = rep(0, NVar),
                               Sigma = round(Rpop,12),  
                               empirical = FALSE)
           
            # model error
            if(cnME$ModelError == TRUE){
               MCDataME[[i]] <-   MASS::mvrnorm(n = cnMC$SampleSize, 
                                      mu = rep(0, NVar),
                                      Sigma = round(RpopME,12),
                                      empirical = FALSE)
            }
            
          } # END  for(i in 1:cnMC$NSamples)
       
       
       
       #  ----_____Monte Carlo Likert Variables----
       if( !is.null(cnMC$Thresholds) ) { # thresholds supplied
     
          # Initialize MCLikertData (to all zeros)
          MCLikertData <- vector(mode="list", length = cnMC$NSamples )
          MCLikertDataME <- MCLikertData
          ZeroMatrix <- matrix(0, nrow = cnMC$SampleSize, ncol = NVar)
          
          for(i in 1:cnMC$NSamples){ 
            MCLikertData[[i]] <- MCLikertDataME[[i]] <- ZeroMatrix 
          }#END for(i in 1:cnMC$NSamples)
         
         
           ## cut data at thresholds 
           for( iReps in 1:cnMC$NSamples){   #loop of MC samples
              for(i in 1:NVar){
                 for(j in 1:length(cnMC$Thresholds[[i]])){
              
                    MCLikertData[[iReps]][( MCData[[iReps]][, i] >= 
                                           cnMC$Thresholds[[i]][j]), i] <- j
                    
                    # model error
                    if(cnME$ModelError == TRUE){
                    MCLikertDataME[[iReps]][( MCDataME[[iReps]][, i] >= 
                                              cnMC$Thresholds[[i]][j]), i] <- j
                    }
                    
                 }#END for(j in 1:length(cnMC$Thresholds[[i]]))
              }#END for(i in 1:NVar)
           }#END for( iReps in 1:cnMC$NSamples) 
          
       MCData <- MCLikertData 
       if(cnME$ModelError == TRUE) MCDataME <- MCLikertDataME 
       
       } # END if( !is.null(cnFS$Thresholds) )
     } # END of (cnMC$NSamples > 0 && Raw == TRUE)
     #  END GENERATE MONTE CARLO DATA 
     
    
   
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#  
#       ---- SEC 12: Generate Factor Scores ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#      
 
     
     FSscores <- NULL
     ObservedScores <- NULL

     if( cnFS$FS == TRUE){
      
       NFac <- cnMD$NFac
       
       # Note that in sampled factor scores, common, specific, and unique
       # factors are not exactly orthogonal due to 
       # sampling errors
       
       # Researchers can fix common factor seed 
       # to fix factor scores across samples.
       if( !is.null(cnFS$CFSeed) ) set.seed(cnFS$CFSeed)
       
       
       # ----_____Generate Common Fac Scores ----
       
       # For non bifactor models
       if(cnBF$Bifactor == FALSE){
          CFS <-    MASS::mvrnorm(n = cnFS$NFacScores,
                            mu = rep(0, NFac),
                            Sigma =  Phi,
                            empirical = cnFS$Population) 
       } # END if(cnBF$Bifactor == FALSE)
       
       # For bifactor models
       if(cnBF$Bifactor == TRUE){
         CFS <-    MASS::mvrnorm(n = cnFS$NFacScores,
                           mu = rep(0, NFac + 1),
                           Sigma =  diag( (NFac + 1) ),
                           empirical = cnFS$Population) 
       } # END  if(cnBF$Bifactor == TRUE)
       
       
       # ----_____Generate Specific Fac Scores ----
       # Reset random number seed
       if( !is.null(cnFS$SFSeed) ) set.seed(cnFS$SFSeed)
       
       # No correlated specific factors
       if( is.null(cnME$RSpecific) ){
            SFS <-     MASS::mvrnorm(n = cnFS$NFacScores,
                               mu = rep(0, NVar),
                               Sigma =  diag(NVar),
                               empirical = cnFS$Population)
       }#END if( is.null(cnME$RSpecific) )
         
      # Correlated specific factors present
      if( !is.null(cnME$RSpecific) ){
           SFS <-     MASS::mvrnorm(n = cnFS$NFacScores,
                              mu = rep(0, NVar),
                              Sigma =  cnME$RSpecific,
                              empirical = cnFS$Population)
       } #END if( !is.null(cnME$RSpecific) )
           
         

       # ----_____Generate Error Factor Scores----
       # Reset random number seed
       if( !is.null(cnFS$EFSeed) ) set.seed(cnFS$EFSeed)
       EFS <-       MASS::mvrnorm(n = cnFS$NFacScores,
                            mu = rep(0, NVar),
                            Sigma =  diag(NVar),
                            empirical = cnFS$Population) 
       
     
       #----_____Population Factor Scores ----
        if(cnFS$Population == TRUE){
          # make specific fac scores orthgonal to common fs
          SFS <- resid(lm(SFS ~ CFS))
          
          # No correlated specific factor
          if( is.null(cnME$RSpecific) ){
             # specfic fac scores orthogonal to each other
             SFS <- scale( svd(SFS)$u )
          }#END if( is.null(cnME$RSpecific) )
          
          # Correlated specific factors present
          if( !is.null(cnME$RSpecific) ){
            K <- chol(cnME$RSpecific)
            # specfic fac scores allowed to correlate
            SFS <-  scale( svd(SFS)$u %*% K )
          }#END if( !is.null(cnME$RSpecific) )
          
         
          # make error fac scores orthgonal to common and specific fs
          EFS <- resid(lm(EFS ~ cbind(CFS, SFS) ))
          # unique scores orthogonal to each other
          EFS <- scale( svd(EFS)$u )
       }#END  if(cnFS$Population == TRUE)
      
       
       # If user DID NOT supply indicator reliabilities (VarRel)  
       # the specific factor scores are uniformly zero
       # If Alphas == NULL 
       if(is.null(cnFS$VarRel)){ 
          SFS <- matrix(0, cnFS$NFacScores, NVar)
          ObservedScores <- CFS %*% t(Fl) + 
                            EFS %*% diag(sqrt (1 - h2))
       }    
          
       # If user DID supply indicator reliabilities (Alphas)                       
       if(!is.null(cnFS$VarRel)){
           if(  min(cnFS$VarRel - h2) < 0 ){
              badItems <- (1:NVar)[(cnFS$VarRel - h2) < 0 ]
              stop("\n\n *** FATAL ERROR (Reliabilities too small) ***",
              "\nItems ", paste0(badItems," "), "have reliabilities smaller than their communalities\n\n")
           } 
         
          # No specific factor correlations
           if( is.null(cnME$RSpecific) ){
             ObservedScores <- CFS %*% t(Fl) + 
                            SFS %*% diag(sqrt(cnFS$VarRel - h2))  +
                            EFS %*% diag(sqrt( 1 - cnFS$VarRel) ) 
           }# END if( is.null(cnME$RSpecific) ) 
         
         # Specific factors allowed to correlate
         if( !is.null(cnME$RSpecific) ){
           
           SFSstndevs <- diag(sqrt(cnFS$VarRel - h2))
           ObservedScores <- CFS %*% t(Fl) + 
             SFS %*% SFSstndevs  +
             EFS %*% diag(sqrt( 1 - cnFS$VarRel) ) 
         }# END if( !is.null(cnME$RSpecific) ) 
         
       }# END if(!is.null(cnFS$Alphas))
       
       FactorScores = cbind(CFS, SFS, EFS)
      
       
       # Add names to factor scores
       if(cnBF$Bifactor == FALSE){
       colnames(FactorScores) <- c(  paste0("F", 1:NFac), 
                                 paste0("S", 1:NVar), # specific factor names
                                 paste0("E", 1:NVar)) # error factor names
       }
       if(cnBF$Bifactor == TRUE){
         colnames(FactorScores) <- c(  "G",   # general factor
                                       paste0("F", 1:NFac), 
                                       paste0("S", 1:NVar), # specific factor names
                                       paste0("E", 1:NVar)) # error factor names
       }# END if(cnBF$Bifactor == TRUE)
       
   
       # Create Likert Items
       if( !is.null(cnFS$Thresholds) ) { # thresholds supplied
         # Declare matrix for Likert data
          LikertData <- matrix(0, nrow = nrow(ObservedScores),
                                  ncol = NVar) 
      
         ## cut data at thresholds   
         for(i in 1:NVar){
           for(j in 1:length(cnFS$Thresholds[[i]])){
            LikertData[(ObservedScores[, i] >= cnFS$Thresholds[[i]][j]), i] <- j
           } 
         }   
      ObservedScores <- LikertData    
      } # END if( !is.null(cnFS$Thresholds) )
       
    }# END if( cnFS$FS == TRUE)
     
     
   
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#  
#       ---- SEC 13: Generate Missing Data ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#      
     
     ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~## 
     ##               Generate Missing Data 
     ## for Monte Carlo Samples
     ## To generate MAR see
     ## https://stats.stackexchange.com/questions/109737/
     ##   how-to-generate-mar-data-with-a-fixed-proportion-of-missing-values
     ##
     ## https://stats.stackexchange.com/questions/184741/
     ##   how-to-simulate-the-different-types-of-missing-data
     ##
     if(cnMS$Missing == TRUE){
        ## Missing Completely at Random
       
        if(cnMS$Mechanism == "MCAR"){
          
           for(iSamp in 1:cnMC$NSamples){
             
              for(iVar in 1:NVar){
                
                  MCData[[iSamp]][
                       runif(cnMC$SampleSize, 0, 1) <= cnMS$MSProb[iVar],
                       iVar] <- NA
                  
                  # model error
                  if(cnME$ModelError == TRUE){
                    MCDataME[[iSamp]][
                      runif(cnMC$SampleSize, 0, 1) <= cnMS$MSProb[iVar],
                      iVar] <- NA
                  }#END if(cnME$ModelError == TRUE)
                  
              }#END for(iVar in 1:NVar)
             
           }#END for(iSamp in 1:cnMC$NSamples)  
          
        }#END if(cnMS$Mechanism == "MCAR")
       
     }#END if(cnMS$Missing == TRUE) 
     
     ## END Generate Missing Data 
     
     
     
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#  
#       ---- SEC 14: Organize Return Values ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#      
     
   cn <- list("Model" = cnMD, 
              "Loadings" = cnLD, 
              "CrossLoadings" = cnCL, 
              "Phi" = cnPH, 
              "ModelError" = cnME, 
              "MonteCarlo" = cnMC, 
              "FactorScores" = cnFS,
              "Missing" = cnMS, 
              "Control" = cnCT)
     
   CovMats <- list("CovMajor" = CovMajor,     # major factor covariance matrix
                  "CovMinor" = CovMinor,     # model error covariance matrix
                  "CovUnique" = CovUnique)   # uniqueness factors covariance matrix  
  
  # _____Output related to Factor Scores----
  FacInd <- FacInd            # Factor Indeterminancy
  FacIndME <-  FacIndME       # Factor Indeterminancy Model Error
  Scores <- list(  "FactorScores" = FactorScores,
                    "FacInd" = FacInd,
                    "FacIndME" = FacIndME,
                    "ObservedScores" = ObservedScores)
  
  # _____Output related to bifactor models----
  Bifactor <- list("loadingsHier" = FlHierarchical,  # Factor loadings of the 
                                                     # 1st order solution
                   "PhiHier"  = PhiHierarchical)     # Factor correlations of the 
                                                     #1st order solution)
  
  # _____Output related to Monte Carlo Simulations----
   Monte = list("MCData" = MCData,             # MonteCarlo raw data
                "MCDataME" = MCDataME)         # MonteCarlo raw data
     
   # _____Output related to IRT parameters----                
   IRT = list("discrimination" = IRTdiscrim,
              "difficulty" = IRTdifficulty,
              "intercept"  = IRTintercept,
              "metric" = "logistic (D is included in a)" )
   

 
   # _____Compute Unrotated Loadings ----
   # Given any factor model, find the unrotated loadings
   # Note that this will not work in BG models in which
   # Fl has more columns than rows
   urloadings = NULL
   if( ncol(Fl) > 1 & (nrow(Fl) > ncol(Fl)) ) {
      Rhat <- Fl %*% Phi %*% t(Fl)
      VLV <- eigen(Rhat)
   #A simple PCA of Rhat gives us what we want
   # because Rhat has the true communalities
      V <- VLV$vectors[,1:cnMD$NFac]
      L <- diag(sqrt(zapsmall(VLV$values[1:cnMD$NFac])))
      urloadings <- V %*% L
   }# END  if(nrow(Fl) > ncol(Fl)) 
   
  ## ---- _____Return Values ----
  list(loadings    = Fl,              # Factor loadings
       Phi         = Phi,             # Factor correlations
       urloadings  = urloadings,      # Unrotated loadings
       h2          = h2,              # communalities for ideal Model
       h2PopME     = h2PopME,
       Rpop        = Rpop,            # no Model error
       RpopME      = RpopME,          # with Model error
       CovMatrices = CovMats,    
       Bifactor    = Bifactor,
       Scores      = Scores,          # Output related to Factor Scores (ObservedScores)
       Monte       = Monte,
       IRT         = IRT,             # Multidimentional IRT
       Seed        = Seed,
       call        = cl,              # program call
       cn          = cn)     
 } #END simFA
