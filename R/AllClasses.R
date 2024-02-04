setClass("Outlier", representation(call = "language",
                              counts = "numeric",
                              grp = "factor",
                              wt = "Uvector",
                              flag = "Uvector",
                              method = "character",
                              singularity = "Ulist",
                              "VIRTUAL")
)

setClass("OutlierMahdist", representation(covobj="Ulist"),
                                    contains="Outlier"
)

setClass("OutlierPCOut", representation(covobj="Ulist"),
                                    contains="Outlier"
)

setClass("OutlierPCDist", representation(covobj="Ulist",
                                    k="numeric"),
                                    contains="Outlier"
)

setClass("OutlierSign1", representation(covobj="Ulist"),
                                    contains="Outlier"
)

setClass("OutlierSign2", representation(covobj="Ulist"),
                                    contains="Outlier"
)

setClass("SPcaGrid", representation(),
                    contains="PcaGrid")

###################### SIMCA ####################################
setClass("Simca", representation(call = "language",
                               prior = "vector",
                               counts = "vector",
                               pcaobj="Ulist",
                               k = "Uvector",
                               flag = "Uvector",
                               X = "Umatrix",
                               grp = "factor",
                               "VIRTUAL"))

setClass("CSimca", contains="Simca")
setClass("RSimca",  contains="Simca")
setClass("PredictSimca", representation(classification = "factor",
                                      odsc = "matrix",
                                      sdsc = "matrix",
                                      ct="Utable"))
setClass("SummarySimca", representation(simcaobj = "Simca"))

###################### SosDisc ####################################
setClass("SosDisc", representation(
                               call = "language",
                               prior = "vector",
                               counts = "vector",

                               beta = "matrix",               # Q coefficient vectors of the predictor matrix from optimal scoring
                               theta = "matrix",              # Q coefficient vectors of the dummy matrix for class coding from optimal scoring
                               lambda = "numeric",            # L1 norm penaly parameter
                               varnames = "character",        # vector of names of selected predictor variables
##                               varIndex = "integer",          # We will remove this (use varnames instead)
##                               origP = "numeric",             # We will remove this (get it as ncol(X))
##                               rss = "numeric",               # We will remove this: vector of length Q with weighted residual sum of squares plus L1 penaltyterm (weighted lasso cost) from optimal scoring

                               ## centering and scaling: later we can add parameters center and scale which can be (TRUE or FALSE, any function or vectors with lengthp)
                               ##   see Pca
                               center = "vector",             # centering vector of the predictors (coordinate wise median)
                               scale = "vector",              # scaling vector of the predictors (mad)

                               fit = "Linda",                 # Linda model (robust LDA model) from the low dimensional subspace
                               mahadist2="vector",            # fixme: These will go later to Linda object: squared robust Mahalanobis distance (calculated with estimates from Linda) to the group center in the low dimensional subspace
                               wlinda = "vector",             # fixme: These will go later to Linda object: weights derived from mahadist2

##                               wlts = "matrix",             # We remove these (not needed): Q weighting vectors from initial sparseLTS estimate
##                               wrew = "matrix",             # We remove these (not needed):Q weighting vectors from reweighting step

                               X = "Umatrix",                 # can be missing
                               grp = "factor",
                               "VIRTUAL"))

setClass("SummarySosDisc", representation(obj = "SosDisc"))
setClass("SosDiscClassic", contains="SosDisc")
setClass("SosDiscRobust", contains="SosDisc")

##setClass("PredictSosDisc", representation(classification = "factor",
##                                      posterior = "matrix",
##                                      x = "matrix",
##                                      ct="Utable"))
##
setClass("PredictSosDisc", representation(classification = "factor",
                                        w = "vector",
                                        mahadist2="matrix"))
