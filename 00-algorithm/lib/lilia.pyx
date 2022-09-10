#cython: boundscheck=False, wraparound=False, initializedcheck=False

#COMPILE WITH: cython3 lilia.pyx && gcc -shared -ffast-math -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -lm -I/usr/include/python3.5 -o lilia.so lilia.c

import numpy as np
from itertools import combinations
cimport numpy as np
cimport cython
from libc.math cimport exp, log

cdef extern from "math.h":
    long double logl(long double x)
    long double expl(long double x)


############################################################################################################################################


cdef inline double lp_alphas(double Aui, double alpha_u, double beta_i):
    return ((Aui * exp(-alpha_u * beta_i) * beta_i) / (1 - exp(-alpha_u * beta_i))) - ((1 - Aui) * beta_i)


cdef inline double lp_betas(double Aui, double alpha_u, double beta_i):
    return ((Aui * exp(-alpha_u * beta_i) * alpha_u) / (1 - exp(-alpha_u * beta_i))) - ((1 - Aui) * alpha_u)


@cython.profile(True)
def ClikelihoodPartial(np.ndarray pars, np.ndarray[np.long_t, ndim=2] paMatrix):
    """
    Calculate the partial derivatives of the likelihood function depending a set of taxon-specific parameters (alpha_u)
    and sample-specific parameters (beta_i) according to Pascual-Garcia et al. (2014).
    Parameters are passed in a single list, with the taxon-specific parameters first.        
    """

    cdef double [:] alphas = pars[:paMatrix.shape[0]]
    cdef double [:] betas = pars[paMatrix.shape[0]:]
    cdef double result
    cdef Py_ssize_t u, i
    cdef double alpha_u, beta_i
    cdef long Aui
    
    cdef double [:] results = np.zeros(paMatrix.shape[0] + paMatrix.shape[1])

    for u, alpha_u in enumerate(alphas):
        result = 0
        for i, beta_i in enumerate(betas):
            Aui = paMatrix[u, i]
            result += lp_alphas(Aui, alpha_u, beta_i)
        results[u] = result

    for i, beta_i in enumerate(betas):
        result = 0
        for u, alpha_u in enumerate(alphas):
            Aui = paMatrix[u, i]
            result += lp_betas(Aui, alpha_u, beta_i)
        results[paMatrix.shape[0] + i] = result
    
    return results


############################################################################################################################################


cdef inline double jacLP_alpha_alpha(double Aui, double alpha_u, double beta_i):
    return (-Aui * beta_i * beta_i * exp(-alpha_u * beta_i)) / ((1 - exp(-alpha_u * beta_i)) * (1 - exp(-alpha_u * beta_i)))


cdef inline double jacLP_beta_beta(double Aui, double alpha_u, double beta_i):
    return (-Aui * alpha_u * alpha_u * exp(-alpha_u * beta_i)) / ((1 - exp(-alpha_u * beta_i)) * (1 - exp(-alpha_u * beta_i)))


cdef inline double jacLP_alpha_beta(double Aui, double alpha_u, double beta_i):
    num = (Aui * exp(-alpha_u * beta_i) - Aui * alpha_u * beta_i * exp(-alpha_u * beta_i) - Aui * exp(-alpha_u * beta_i) * exp(-alpha_u * beta_i)) 
    den = (1 - exp(-alpha_u * beta_i)) * (1 - exp(-alpha_u * beta_i))
    return (num / den) - (1 - Aui)


@cython.profile(True)
def CjacLP(np.ndarray pars, np.ndarray[np.long_t, ndim=2] paMatrix):
    """
    Calculate the jacobian of the equation system described in self.likelihoodPartial.
    Our matrix is in the form:

                         -----------------------------------------------------------------------------------------
                         |      d(alpha1)      |      d(alpha2)      |      d(beta1)       |      d(beta2)       |
    --------------------------------------------------------------------------------------------------------------
    |                    | d(lnL) / d(alpha1)  |                     | d(lnL) / d(alpha1)  | d(lnL) / d(alpha1)  |
    | d(lnL) / d(alpha1) | __________________  |          0          | __________________  | __________________  |
    |                    |     d(alpha1)       |                     |     d(beta1)        |     d(beta2)        |
    --------------------------------------------------------------------------------------------------------------
    |                    |                     | d(lnL) / d(alpha2)  | d(lnL) / d(alpha2)  | d(lnL) / d(alpha2)  | 
    | d(lnL) / d(alpha2) |          0          | __________________  | __________________  | __________________  |
    |                    |                     |     d(alpha2)       |     d(beta1)        |     d(beta2)        |
    --------------------------------------------------------------------------------------------------------------
    |                    | d(lnL) / d(beta1)   | d(lnL) / d(beta1)   | d(lnL) / d(beta1)   |                     |
    | d(lnL) / d(beta1)  | __________________  | __________________  | __________________  |          0          |
    |                    |     d(alpha1)       |     d(alpha2)       |     d(beta1)        |                     |
    --------------------------------------------------------------------------------------------------------------
    |                    | d(lnL) / d(beta2)   | d(lnL) / d(beta2)   |                     | d(lnL) / d(beta2)   |
    | d(lnL) / d(beta2)  | __________________  | __________________  |          0          | __________________  |
    |                    |     d(alpha1)       |     d(alpha2)       |                     |     d(beta2)        |
    --------------------------------------------------------------------------------------------------------------
    """

    cdef double [:] alphas = pars[:paMatrix.shape[0]]
    cdef double [:] betas = pars[paMatrix.shape[0]:]
    cdef Py_ssize_t lenAlphas = len(alphas)
    cdef double result
    cdef Py_ssize_t u, i
    cdef double alpha_u, beta_i, Aui

    cdef double [:,:] jacobian = np.zeros([pars.shape[0], pars.shape[0]])

    for u, alpha_u in enumerate(alphas):
        #Second partial derivative of the likelihood function with respect to alpha_u
        result = 0
        for i, beta_i in enumerate(betas):
            Aui = paMatrix[u, i]
            result += jacLP_alpha_alpha(Aui, alpha_u, beta_i)
        jacobian[u, u] = result
        #Derivatives of d(lnL)/d(alpha_u) with respect to other alphas are zero.
        for i, beta_i in enumerate(betas):
            #Derivatives of d(lnL)/d(alpha_u) with respect to betas.
            Aui = paMatrix[u, i]
            jacobian[u, lenAlphas + i] = jacLP_alpha_beta(Aui, alpha_u, beta_i)
    for i, beta_i in enumerate(betas):
        #Second partial derivative of the likelihood function with respect to beta_i
        result = 0
        for u, alpha_u in enumerate(alphas):
            Aui = paMatrix[u, i]
            result += jacLP_beta_beta(Aui, alpha_u, beta_i)
        jacobian[lenAlphas + i, lenAlphas + i] = result
        #Derivatives of d(lnL)/d(beta_i) with respect to other betas are zero.
        for u, alpha_u in enumerate(alphas):
            #Derivatives of d(lnL)/d(beta_u) with respect to alphas. This is actually the same expression as d(lnL)/d(alpha_i) with respect to betas.
            Aui = paMatrix[u, i]
            jacobian[lenAlphas + i, u] = jacLP_alpha_beta(Aui, alpha_u, beta_i)

    return jacobian


############################################################################################################################################


@cython.profile(True)
def CbuildMinCosmoMatrix(np.ndarray [np.long_t, ndim = 2] paMatrix, np.ndarray[np.long_t, ndim=1] cosmo):
    cdef Py_ssize_t i, j
    cdef np.ndarray [np.long_t, ndim = 2] minCosmoMatrix = np.zeros([paMatrix.shape[0], paMatrix.shape[0]], dtype=np.int)
    for i in range(paMatrix.shape[0]):
        for j in range(paMatrix.shape[0]):
            minCosmoMatrix[i, j] = min(cosmo[i], cosmo[j])

    return minCosmoMatrix


@cython.profile(True)
def CbuildCoocMatrix(np.ndarray [np.long_t, ndim = 2] paMatrix):
    cdef Py_ssize_t i, j, cooc, z
    cdef np.ndarray[np.long_t, ndim = 2] coocMatrix = np.zeros([paMatrix.shape[0], paMatrix.shape[0]], dtype=np.int)
    for i in range(paMatrix.shape[0]):
        for j in range(paMatrix.shape[0]):
            cooc = 0
            for z in range(paMatrix.shape[1]):
                cooc += paMatrix[i,z] * paMatrix[j,z]
            coocMatrix[i, j] = cooc
    return coocMatrix


############################################################################################################################################


cdef inline double aggScore(Py_ssize_t i, Py_ssize_t j, long[:, :] paMatrix, long[:, :] coocMatrix, long double[:, :] probMatrix, long double[:, :] likes):
    """
    Calculate aggregation score, defined as -log of the conditional probability of taxa i y j co-occurring in N over M total samples
    the observed times or more.

                                                         ⌈                 nij_obs - 1                   ⌉   
    - Score_ij = log[P(cond_ij)] - log| P(cond_ij) - Σ(Pij(nij = n, M, cond_ij) |
                                                         ⌊                     0                          ⌋
                                
    - P(cond_ij) is the likelihood of the probability matrix realizing into the observed presence absence value for the conditioning
      taxon in each sample, with the conditioning taxon in a given sample being the taxon with the smallest likelihood in that sample.

      like(i,m) = paMatrix[i,m] *  probMatrix[i,m] + (1 - paMatrix[i,m]) * (1 - probMatrix[i,m])
      like(j,m) = paMatrix[j,m] *  probMatrix[j,m] + (1 - paMatrix[j,m]) * (1 - probMatrix[j,m])
      cond(m) = c = i if like(i,m) < like(j,m) else j
      nonCond(m) = z = j if like(i,m) < like(j,m) else i
      log[P(cond_ij)] = logPCond = sum([log(like(cond(m), m)) for m in range(samples)])
            
    - Pij(n,m,cond_ij) is recursively calculated as follows:
            
        Pij(n,m,cond_ij) {if Acm == 1} = (Pij(n-1, m-1,cond_ij) * probMatrix[c][m] * probMatrix[z][m]
                                        + Pij(n, m-1, cond_ij) *  probMatrix[c][m] * (1 - probMatrix[z][m]))
                         {if Acm == 0} = (Pij(n, m-1, cond_ij) * (1 - probMatrix[c][m]))

    In practice we will get Pij(n,m,cond_ij) iteratively, starting for n = 0.
    """

    #Note that we used memoryviews instead of numpy arrays for coocMatrix, probMatrix and likes.
    #This allows us to inline the function and increase the speed.
    #https://jakevdp.github.io/blog/2012/08/16/memoryview-benchmarks-2/
    #Memoryviews could also be used for the scores matrix, but it did not affect execution speed.
    
    cdef Py_ssize_t n, m, c, z 
    #Calculate the conditioning taxon for each sample and P(cond:ij).
    cdef long double logPCond = 0
    cdef np.ndarray[Py_ssize_t, ndim = 1] condTaxon = np.zeros(likes.shape[1], dtype = np.int)
    cdef np.ndarray[Py_ssize_t, ndim = 1] nonCondTaxon = np.zeros(likes.shape[1], dtype = np.int)
    for m in range(likes.shape[1]):
        condTaxon[m] = i if likes[i,m] < likes[j,m] else j
        nonCondTaxon[m] = j if likes[i,m] < likes[j,m] else i
        logPCond += logl(likes[condTaxon[m],m])
    #Initialize result variables.
    cdef long double pCond = expl(logPCond)
    cdef long double prob
    cdef long double pCooc = 0
    cdef long double pCoocN
    cdef double score
    cdef np.ndarray[long double, ndim = 2] coocProbs = np.zeros([coocMatrix[i, j], probMatrix.shape[1]], dtype = np.float128)
    #Calculate aggScore.
    if not coocMatrix[i, j]:
        score = 0 #Since the probability of having 0 or more co-occurrences over M samples is 1.
    else:
        for n in range(coocProbs.shape[0]): #This actually goes from 0 to the number of co-occurrences - 1 (N-1), which is fine.
            for m in range(coocProbs.shape[1]):
                c, z = condTaxon[m], nonCondTaxon[m]
                if n == 0 and m == 0:
                    prob = paMatrix[c, m] * probMatrix[c, m] * (1 - probMatrix[z, m])
                    prob += (1 - paMatrix[c, m]) * (1 - probMatrix[c, m])
                    coocProbs[n, m] = prob
                elif n == 0 and m > 0: #The first part of the expression for Acm == 1 depends on P(-1|m), which is 0.
                    prob = paMatrix[c, m] * coocProbs[n, m - 1] * probMatrix[c, m] * (1 - probMatrix[z, m])
                    prob += (1 - paMatrix[c, m]) * coocProbs[n, m - 1] * (1 - probMatrix[c, m])
                    coocProbs[n, m] = prob
                elif n == 1 and m == 0: #Probability of getting a co-occurrence in the first sample.
                    prob = paMatrix[c, m] * probMatrix[c, m] * probMatrix[z, m]
                    #If paMatrix[c, m] == 0 we can't get a co-occurrence.
                    coocProbs[n, m] = prob
                elif n > m + 1: #We can't get more co-occurrences than samples. Note that m is sample index, thus m + 1 is the number of samples.
                    coocProbs[n,m] = 0
                else:
                    prob = paMatrix[c, m] * (coocProbs[n - 1, m - 1] * probMatrix[c, m] * probMatrix[z, m] +
                                              coocProbs[n, m - 1] * probMatrix[c, m] * (1 - probMatrix[z, m]))
                    prob += (1 - paMatrix[c, m]) * coocProbs[n, m - 1] * (1 - probMatrix[c, m])
                    coocProbs[n, m] = prob
        
        for n in range(coocProbs.shape[0]):
            pCooc += coocProbs[n, coocProbs.shape[1] - 1]

        if pCond > pCooc:
            score = logPCond - logl(pCond - pCooc)
        else:
            #This could actually happen due to rounding errors in corner cases.
            #In this case we define the score as P(n_ij_obs,M,cond_ij)
            coocProbs = np.vstack([coocProbs, np.zeros(probMatrix.shape[1], dtype = np.float128)])
            n = coocMatrix[i, j]
            for m in range(coocProbs.shape[1]):
                c, z = condTaxon[m], nonCondTaxon[m]
                prob = paMatrix[c, m] * (coocProbs[n - 1, m - 1] * probMatrix[c, m] * probMatrix[z, m] +
                                         coocProbs[n, m - 1] * probMatrix[c, m] * (1 - probMatrix[z, m]))
                prob += (1 - paMatrix[c, m]) * coocProbs[n, m - 1] * (1 - probMatrix[c, m])
                coocProbs[n, m] = prob
            pCoocN = coocProbs[n, coocProbs.shape[1]-1]
            score = logPCond - logl(pCoocN)
            assert score > 0
            print('\tThere were rounding errors in score calculation for pair: {}-{}'.format(i, j))
            print('\tA score of {:.2f} was asigned by an approximate method.'.format(score))


    return score


@cython.profile(True)
def CbuildAggScoreMatrix(np.ndarray [np.long_t, ndim=2] paMatrix, np.ndarray [np.long_t, ndim=2] coocMatrix, np.ndarray [long double, ndim=2] probMatrix, taxon = False, Py_ssize_t idx = 0):
    """
    Calculate aggregation scores for a given probability matrix. Aggregation score of a taxon to itself is left at zero.
    If the taxon and idx parameters are provided, calculate scores only for this taxon with the others (should be a numerical index).
    Else, calculate all possible pairs.
    """
    cdef Py_ssize_t i, m, t1, t2
    cdef double score
    #The aggregation score of a taxon with itself will be zero.
    cdef np.ndarray [double, ndim=2] aggScoreMatrix = np.zeros([coocMatrix.shape[0], coocMatrix.shape[1]])
    #Calculate the likelihood of the probability matrix realizing into the observed presence/absence for each independent cell.
    cdef np.ndarray [long double, ndim=2] likes = np.zeros([probMatrix.shape[0], probMatrix.shape[1]], dtype=np.float128)
    for i in range(likes.shape[0]):
        for m in range(likes.shape[1]):
            likes[i,m] = paMatrix[i,m] *  probMatrix[i,m] + (1 - paMatrix[i,m]) * (1 - probMatrix[i,m])
    #Create memoryviews.
    cdef long[:, :] paMatrixV = paMatrix
    cdef long[:, :] coocMatrixV = coocMatrix
    cdef long double [:, :] probMatrixV = probMatrix
    cdef long double [:, :] likesV = likes
    #Calculate aggregation scores for the requested pairs of taxa.
    if not taxon:
        for t1, t2 in combinations(range(coocMatrix.shape[0]), 2):
            score = aggScore(t1, t2, paMatrixV, coocMatrixV, probMatrixV, likesV)
            aggScoreMatrix[t1, t2] = score
            aggScoreMatrix[t2, t1] = score
    else:
        for t2 in range(coocMatrix.shape[0]):
            if idx != t2:
                score = aggScore(idx, t2, paMatrixV, coocMatrixV, probMatrixV, likesV)
                aggScoreMatrix[idx, t2] = score
                aggScoreMatrix[t2, idx] = score
    return aggScoreMatrix

