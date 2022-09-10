from multiprocessing import Pool
from itertools import combinations

import numpy as np
from pandas import DataFrame, read_csv, concat

from scipy.optimize import fsolve

import warnings
import sys

from rpy2 import robjects

from copy import copy

from lib.lilia import ClikelihoodPartial, CjacLP, CbuildMinCosmoMatrix, CbuildCoocMatrix, CbuildAggScoreMatrix




class Model:
    """
    Class containing the presence-absence, min-cosmopolitanism, co-occurrence, probability and aggregation scores matrices for a given
    presence-absence matrix, and the methods to calculate them.
    The presence-absence matrix will be partitioned into sub-matrices, one for each environment, a different probability matrix will be calculated
    for each of those sub-matrices. The resulting matrices will then be concatenated and reordered again to generate a complete probability matrix.
    By computing the model independently for each environment, we relax the assumption that there is no preferential association between taxa and environments.
    If the replicates parameter is set to a non-zero value, also calculate a corrected aggregation scores matrix at the required pVal.
    """
    def __init__(self, paMatrix, environments = None, pVal = 0.001, replicates = 0, processors=1, null = False):
        self.paMatrix = paMatrix
        #Make sure paMatrix and probMatrix will be sorted the same way.
        self.paMatrix.sort_index(inplace = True)
        self.paMatrix.sort_index(axis = 1, inplace = True)
        #Get cosmopolitanism for all taxa.
        self.cosmo = np.sum(self.paMatrix.values, axis = 1)
        self.minCosmoMatrix = DataFrame(CbuildMinCosmoMatrix(self.paMatrix.values, self.cosmo),
                                        index = self.paMatrix.axes[0], columns = self.paMatrix.axes[0], dtype = np.int)
        self.iCosmoMatrix, self.jCosmoMatrix = self.buildCosmoMatrices(self.paMatrix.values, self.cosmo)
        self.null = null
        #Calculate probMatrix
        if not null: #If this is the original matrix, raise warnings as exceptions. Else, raise them normally and let them be catched by self.generateNullMatrix.
            warnings.simplefilter('error')
        else:
            warnings.simplefilter('always')
        if not environments:
           self.probMatrix = self.buildProbabilityMatrix(copy(self.paMatrix))
        else:
            #environments = {envName: set([samples])} etc -> Can be supertypes, types or subtypes, depending on what was passed when creating the object.
            partialMatrices = []
            #One submatrix per environment.
            for environment, samples in environments.items():
                partialMatrix = self.buildProbabilityMatrix(self.paMatrix[list(samples)])
                partialMatrices.append(partialMatrix)
            self.probMatrix = concat(partialMatrices, axis = 1)
        #Make sure paMatrix and probMatrix are sorted the same way.
        self.probMatrix.sort_index(inplace = True)
        self.probMatrix.sort_index(axis = 1, inplace = True)
        #Move on.
        self.coocMatrix = DataFrame(CbuildCoocMatrix(self.paMatrix.as_matrix()),
                                    index = self.paMatrix.axes[0], columns = self.paMatrix.axes[0], dtype = np.int)
        self.aggScoreMatrix = DataFrame(CbuildAggScoreMatrix(self.paMatrix.as_matrix(), self.coocMatrix.as_matrix(), self.probMatrix.as_matrix()),
                                        index = self.coocMatrix.index, columns = self.coocMatrix.columns)
        #Calculate corrected aggregation scores, if required.
        if replicates:
            pool = Pool(processors)
            results = pool.starmap_async(self.generateNullMatrix, [(i,environments) for i in range(replicates)])
            pool.close()
            pool.join()
            self.nullMatrices = results.get()
            self.cutoff, self.meanSpline, self.stdSpline = self.fitZscoreFunction(pVal)
            self.correctScores() #This creates self.correctedAggScoreMatrix, which is a pandas dataframe similar to self.aggScoreMatrix


    def write(self, matrix, name = None):
        if matrix not in ('pa', 'minCosmo', 'cooc', 'prob', 'agg', 'corAgg'):
            raise ValueError('matrix must be "pa", "minCosmo", "cooc", "prob", "agg" or "corAgg"')
        if not name:
            name = '{}.tsv'.format(matrix)
        if matrix == 'pa':
            self.paMatrix.astype(np.int).to_csv(name, sep='\t')
        elif matrix == 'minCosmo':
            self.minCosmoMatrix.astype(np.int).to_csv(name, sep='\t')
        elif matrix == 'cooc':
            self.coocMatrix.astype(np.int).to_csv(name, sep='\t')
        elif matrix == 'prob':
            self.probMatrix.to_csv(name, sep='\t')
        elif matrix == 'agg':
            self.aggScoreMatrix.to_csv(name, sep='\t')
        else: #matrix == 'corAgg'
            self.correctedAggScoreMatrix.to_csv(name, sep='\t')


    def buildCoocMatrix(self, paMatrix):
        coocMatrix = np.zeros((paMatrix.shape[0], paMatrix.shape[0]))
        for i, profile_i in enumerate(paMatrix):
            for j, profile_j in enumerate(paMatrix):
                coocMatrix[i, j] = sum(profile_i * profile_j)

        return coocMatrix


    def buildCosmoMatrices(self, paMatrix, cosmo):
        iCosmoMatrix = np.zeros((paMatrix.shape[0], paMatrix.shape[0]))
        jCosmoMatrix = np.zeros((paMatrix.shape[0], paMatrix.shape[0]))
        for i in range(paMatrix.shape[0]):
            for j in range(paMatrix.shape[0]):
                iCosmoMatrix[i, j] = cosmo[i]
                jCosmoMatrix[i, j] = cosmo[j]
        iCosmoMatrix = DataFrame(iCosmoMatrix, index = self.paMatrix.axes[0], columns = self.paMatrix.axes[0])
        jCosmoMatrix = DataFrame(jCosmoMatrix, index = self.paMatrix.axes[0], columns = self.paMatrix.axes[0])

        return iCosmoMatrix, jCosmoMatrix



    def buildMinCosmoMatrix(self, paMatrix, cosmo):
        minCosmoMatrix = np.zeros((paMatrix.shape[0], paMatrix.shape[0]))
        for i in range(paMatrix.shape[0]):
            for j in range(paMatrix.shape[0]):
                minCosmoMatrix[i, j] = min(cosmo[i], cosmo[j])

        return minCosmoMatrix

    
    def buildProbabilityMatrix(self, paMatrix):
        """
        Build a matrix with the probability of finding taxon u in sample i, according to the null model presented in Pascual-Garcia et al. (2014)
        and Navarro-Alberto & Manly (2009).
        """
        #First remove all-zero and all-one columns and rows, since they will result in divisions by zero. They will be inserted back afterwards.
        allOneRows = [name for (name, rowTotal) in zip(paMatrix.sum(axis=1).index, paMatrix.sum(axis=1)) if rowTotal == paMatrix.shape[1]]
        allOneCols = [name for (name, colTotal) in zip(paMatrix.sum(axis=0).index, paMatrix.sum(axis=0)) if colTotal == paMatrix.shape[0]]
        paMatrix = paMatrix.drop(allOneRows)
        paMatrix = paMatrix.drop(allOneCols, axis = 1)
        allZeroRows = [name for (name, rowTotal) in zip(paMatrix.sum(axis=1).index, paMatrix.sum(axis=1)) if not rowTotal]
        allZeroCols = [name for (name, colTotal) in zip(paMatrix.sum(axis=0).index, paMatrix.sum(axis=0)) if not colTotal]
        paMatrix = paMatrix.drop(allZeroRows)
        paMatrix = paMatrix.drop(allZeroCols, axis = 1)
        #Calculate initial guesses for the alpha_u and beta_i parameters. Parameters are passed in a single list, with the taxon-specific parameters first.
        guesses = []
        for u in range(paMatrix.shape[0]):
            guesses.append(paMatrix.sum(1).iat[u] / np.sqrt(sum(paMatrix.sum(1)))) #Cosmopolitanism of taxon u / sqrt(total observations)
        for i in range(paMatrix.shape[1]):
            guesses.append(paMatrix.sum(0).iat[i] / np.sqrt(sum(paMatrix.sum(1)))) #Occupancy of sample i / sqrt(total observations)
        #Find the set of alpha_u and beta_i parameters that would result in a probability matrix with a maximum likelihood of generating our
        #presence absence matrix. Note that we are calling the cython implementation, not the python implementation shown below.
        pars = fsolve(ClikelihoodPartial, np.array(guesses), args = paMatrix.as_matrix(), fprime = CjacLP, diag=1/np.array(guesses))
        alphas = pars[:paMatrix.shape[0]]    #The 'diag' parameter in fsolve is meant to include scaling factors so all the variables are 
        betas = pars[paMatrix.shape[0]:]     #roughly in the same scale. It helps to reach convergence in some corner cases.
        #Generate the probability matrix. Calculating products on the probability matrix requires quadruple precision floating point numbers. 
        probMatrix = DataFrame(np.zeros(paMatrix.shape), index = paMatrix.axes[0], columns = paMatrix.axes[1], dtype = np.float128)
        for u, alpha_u in enumerate(alphas):
            for i, beta_i in enumerate(betas):
                probMatrix.iat[u, i] = 1 - np.exp(-alpha_u * beta_i)
        #Insert all-zero and all-one rows and columns.
        for row in allZeroRows:
            probMatrix.ix[row] = np.zeros(probMatrix.shape[1])
        for column in allZeroCols:
            probMatrix.insert(0, column, np.zeros(probMatrix.shape[0]))
        for row in allOneRows:
            probMatrix.ix[row] = np.ones(probMatrix.shape[1])
        for column in allOneCols:
            probMatrix.insert(0, column, np.ones(probMatrix.shape[0]))
        probMatrix.sort_index(inplace = True)
        probMatrix.sort_index(axis = 1, inplace = True)
        return probMatrix


    def likelihoodPartial(self, pars, paMatrix):
        """
        Calculate the partial derivatives of the likelihood function depending a set of taxon-specific parameters (alpha_u)
        and sample-specific parameters (beta_i) according to Pascual-Garcia et al. (2014).
        Parameters are passed in a single list, with the taxon-specific parameters first.        
        """
        alphas = pars[:paMatrix.shape[0]]
        betas = pars[paMatrix.shape[0]:]
        results = []
        for u, alpha_u in enumerate(alphas):
            result = 0
            for i, beta_i in enumerate(betas):
                Aui = paMatrix[u, i]
                result += ((Aui * np.exp(-alpha_u * beta_i) * beta_i) / (1 - np.exp(-alpha_u * beta_i))) - ((1 - Aui) * beta_i)
            results.append(result)

        for i, beta_i in enumerate(betas):
            result = 0
            for u, alpha_u in enumerate(alphas):
                Aui = paMatrix[u, i]
                result += ((Aui * np.exp(-alpha_u * beta_i) * alpha_u) / (1 - np.exp(-alpha_u * beta_i))) - ((1 - Aui) * alpha_u)
            results.append(result)
        index = results.index(max(results))
        return results


    def jacLP(self, pars, paMatrix):
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

        alphas = pars[:paMatrix.shape[0]]
        betas = pars[paMatrix.shape[0]:]
        jacobian = np.zeros([len(pars), len(pars)])

        for u, alpha_u in enumerate(alphas):
            #Second partial derivative of the likelihood function with respect to alpha_u
            result = 0
            for i, beta_i in enumerate(betas):
                Aui = paMatrix[u, i]
                result += ((-Aui * (beta_i ** 2)) * np.exp(-alpha_u * beta_i)) / ((1 - np.exp(-alpha_u * beta_i)) ** 2)
            jacobian[u, u] = result
            #Derivatives of d(lnL)/d(alpha_u) with respect to other alphas are zero.
            for i, beta_i in enumerate(betas):
                #Derivatives of d(lnL)/d(alpha_u) with respect to betas.
                Aui = paMatrix[u, i]
                num = (Aui * np.exp(-alpha_u * beta_i) - Aui * alpha_u * beta_i * np.exp(-alpha_u * beta_i) - Aui * (np.exp(-alpha_u * beta_i)) ** 2)
                den = ((1 - np.exp(-alpha_u * beta_i)) ** 2)
                jacobian[u, len(alphas) + i] = (num / den) - (1 - Aui)

        for i, beta_i in enumerate(betas):
            #Second partial derivative of the likelihood function with respect to beta_i
            result = 0
            for u, alpha_u in enumerate(alphas):
                Aui = paMatrix[u, i]
                result += ((-Aui * (alpha_u ** 2)) * np.exp(-alpha_u * beta_i)) / ((1 - np.exp(-alpha_u * beta_i)) ** 2)
            jacobian[len(alphas) + i, len(alphas) + i] = result
            #Derivatives of d(lnL)/d(beta_i) with respect to other betas are zero.
            for u, alpha_u in enumerate(alphas):
                #Derivatives of d(lnL)/d(beta_u) with respect to alphas. This is actually the same expression as d(lnL)/d(alpha_i) with respect to betas.
                Aui = paMatrix[u, i]
                num = (Aui * np.exp(-alpha_u * beta_i) - Aui * alpha_u * beta_i * np.exp(-alpha_u * beta_i) - Aui * (np.exp(-alpha_u * beta_i)) ** 2)
                den = ((1 - np.exp(-alpha_u * beta_i)) ** 2)
                jacobian[len(alphas) + i, u] = (num / den) - (1 - Aui)

        return jacobian


    def buildAggScoreMatrix(self, paMatrix, coocMatrix, probMatrix, taxon = False, idx = 0):
        """
        Calculate aggregation scores for a given probability matrix. Aggregation score of a taxon to itself is left at zero.
        If the taxon and idx parameters are provided, calculate scores only for this taxon with the others (should be a numerical index).
        Else, calculate all possible pairs.
        """
        #The aggregation score of a taxon with itself will be zero.
        aggScoreMatrix = np.zeros(coocMatrix.shape)
        #Calculate the likelihood of the probability matrix realizing into the observed presence/absence for each independent cell.
        likes = np.array([[paMatrix[i,m] *  probMatrix[i,m] + (1 - paMatrix[i][m]) * (1 - probMatrix[i,m]) for m in range(paMatrix.shape[1])]
                 for i in range(paMatrix.shape[0])], dtype = np.float128)
        #Calculate aggregation scores for the requested pairs of taxa.
        if not taxon:
            for t1, t2 in combinations(range(coocMatrix.shape[0]), 2):
                score = self.aggScore(t1, t2, paMatrix, coocMatrix, probMatrix, likes)
                aggScoreMatrix[t1, t2] = score
                aggScoreMatrix[t2, t1] = score
        else:
            for t2 in range(coocMatrix.shape[0]):
                if idx != t2:
                    score = self.aggScore(idx, t2, paMatrix, coocMatrix, probMatrix, likes)
                    aggScoreMatrix[idx, t2] = score
                    aggScoreMatrix[t2, idx] = score

        return aggScoreMatrix


    def aggScore(self, i, j, paMatrix, coocMatrix, probMatrix, likes): #Note that here we use i,j for taxa, and n for samples. In the rest of the program, u is used for taxa and i for samples.
        """
        Calculate aggregation score, defined as -log of the conditional probability of taxa i y j co-occurring in N over M total samples
        the observed times or more.

                                             ⌈           nij_obs - 1                   ⌉   
        - aggScore_ij = log[P(cond_ij)] - log| P(cond_ij) - Σ(Pij(nij = n, M, cond_ij) |
                                             ⌊              0                          ⌋

        
                                             ⌈                 nij_obs - 1                   ⌉   
        - segScore_ij = log[P(cond_ij)] - log| P(cond_ij) - [1 - Σ(Pij(nij = n, M, cond_ij)] |
                                             ⌊                    0                          ⌋
  
        - zero_ij = log[P(cond_ij)] - log[P(cond_ij) - Pij(nij = 0, M, cond_ij)]
        
        - P(cond_ij) is the likelihood of the probability matrix realizing into the observed presence absence value for the conditioning
          taxon in each sample, with the conditioning taxon in a given sample being the taxon with the smallest likelihood in that sample.

                  like(i,m) = paMatrix[i][m] *  probMatrix[i][m] + (1 - paMatrix[i][m]) * (1 - probMatrix[i][m])
                  like(j,m) = paMatrix[j][m] *  probMatrix[j][m] + (1 - paMatrix[j][m]) * (1 - probMatrix[j][m])
                  cond(m) = c = i if like(i,m) < like(j,m) else j
                  nonCond(m) = z = j if like(i,m) < like(j,m) else i
                  log[P(cond_ij)] = logPCond = sum([log(like(cond(m), m)) for m in range(samples)])
        
        - Pij(n,m,cond_ij) is recursively calculated as follows:
        
            Pij(n,m,cond_ij) {if Acm == 1} = (Pij(n-1, m-1,cond_ij) * probMatrix[c][m] * probMatrix[z][m]
                                            + Pij(n, m-1, cond_ij) *  probMatrix[c][m] * (1 - probMatrix[z][m]))
                             {if Acm == 0} = (Pij(n, m-1, cond_ij) * (1 - probMatrix[c][m]))

        In practice we will get Pij(n,m,cond_ij) iteratively, starting for n = 0.
        """
        if not coocMatrix[i, j]:
            score = 0 #Since the probability of having 0 or more co-occurrences over M samples is 1.  
        else:
            condTaxon = [i if likes[i,m] < likes[j,m] else j for m in range(paMatrix.shape[1])]
            nonCondTaxon = [j if likes[i,m] < likes[j,m] else i for m in range(paMatrix.shape[1])]
            logPCond = sum(np.log([likes[i,m] if condTaxon[m] == i else likes[j,m] for m in range(paMatrix.shape[1])]))
            coocProbs = np.zeros([coocMatrix[i, j], probMatrix.shape[1]], dtype = np.float128)
            for n in range(coocProbs.shape[0]): #This actually goes from 0 to the number of co-occurrences - 1; (N-1), which is fine.
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

            if np.exp(logPCond) > np.sum(np.sum(coocProbs, axis = 0)[-1]):
                score = logPCond - np.log(np.exp(logPCond) - (np.sum(coocProbs, axis = 0)[-1]))
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
                score = logPCond - np.log(pCoocN)
                assert score > 0
                print('\tThere were rounding errors in the score calculation for pair: {}-{}'.format(self.paMatrix.index[i], self.paMatrix.index[j]))
                print('\tA score of {:.2f} was asigned by an approximate method.'.format(score))

        return score


    def generateNullMatrix(self, i, environments):
        """
        Generate a realization of the probability matrix.
        Return a Model object.
        """
        print('Generating random matrix {:d}'.format(i))
        np.random.seed() #Since worker process inherits the random number generator from the parent process, every worker would have the same output.
                         #I manually checked that this was actually happening here.
        w = True
        ntries = 1
        while w:
            with warnings.catch_warnings(record = True) as w:
                random = np.random.rand(*self.probMatrix.shape)
                mod =  Model(DataFrame(data = (random < self.probMatrix).astype(float),
                                       index = self.probMatrix.index, columns = self.probMatrix.columns, dtype = np.int),
                                       environments, null = True)
                if w:
                    ntries += 1
                    if ntries < 5:
                        print('\tNo solution, trying again... (#{})'.format(ntries))
                    else:
                        print('\t' + warnings.formatwarning(w[0].message, w[0].category, w[0].filename, w[0].lineno).strip())
        return mod
                

    def fitZscoreFunction(self, pVal):
        """
        Agregation scores between a given pair of taxa depend on their minimum cosmopolitanism, that is,
        the number of samples in which the less cosmopolitan taxon appears. Therefore, in order to make them
        comparable between pairs with different minCosmo, we will normalize them to Zscores,
        with the average and standard deviations of the aggregation being fit with respect to min cosmo.
        - For all the null model matrices and all the possible pairs of taxa obtain aggScore and minCosmo.
        - Sort the data in increasing order of min(cosmo).
        - Divide the data in chunks of 10/p_val pairs.
        - For each chunk:
            - Calculate mean(minCosmo).
            - Calculate mean(aggScore).
            - Calculate std(aggScore).
        - Fit log(mean(aggScore)) ~ mean(minCosmo)
        - Fit log(std(aggScore)) ~ mean(minCosmo)
        - For all of our minCosmo-aggScore pairs:
            - Calculate Zscore:
                  Zscore_ij = (aggScore_ij - scoreMean(minCosmo_ij)) / scoreStd(minCosmo_ij)
        - Sort the resulting Zscores, and obtain the Zscore cutoff that would result in the requested pVal.
        - Return Zscore_cutoff, slope and intercept of mean(aggScore) ~ mean(minCosmo), slope and intercept of std(aggScore) ~ mean(minCosmo).
        """
        minCosmo = np.concatenate([null.minCosmoMatrix for null in self.nullMatrices]).flatten()
        scores = np.concatenate([null.aggScoreMatrix for null in self.nullMatrices]).flatten()
        cooc = np.concatenate([null.coocMatrix for null in self.nullMatrices]).flatten()
        cosmo_i = np.concatenate([null.iCosmoMatrix for null in self.nullMatrices]).flatten()
        cosmo_j = np.concatenate([null.jCosmoMatrix for null in self.nullMatrices]).flatten()
        chunkSize = 1000
        #chunkSize = 10/pVal
        #Sort by minCosmo and split into equal chunks.
        sortIndices = np.array([index for index in minCosmo.argsort() if scores[index]]) #Also avoid non-interacting pairs.
        minCosmo = minCosmo[sortIndices]
        scores = scores[sortIndices]
        cooc = cooc[sortIndices]
        cosmo_i = cosmo_i[sortIndices]
        cosmo_j = cosmo_j[sortIndices]
        numChunks = len(minCosmo) // chunkSize
        assert numChunks > 10
        minCosmoChunks = np.array_split(minCosmo, numChunks)
        scoresChunks = np.array_split(scores, numChunks)
        coocChunks = np.array_split(cooc, numChunks)
        cosmo_iChunks = np.array_split(cosmo_i, numChunks)
        cosmo_jChunks = np.array_split(cosmo_j, numChunks)
        #For each chunk, get mean(minCosmo), mean(aggScore) and std(aggScore).
        meanMinCosmo = []
        meanAggScore = []
        meanCosmo_i = []
        meanCosmo_j = []
        meanCooc = []
        stdAggScore = []
        for chunkCosmo, chunkScores, chunkCooc, chunk_i, chunk_j in zip(minCosmoChunks, scoresChunks, coocChunks, cosmo_iChunks, cosmo_jChunks):
            meanMinCosmo.append(np.mean(chunkCosmo))
            meanAggScore.append(np.mean(chunkScores))
            meanCooc.append(np.mean(chunkCooc))
            meanCosmo_i.append(np.mean(chunk_i))
            meanCosmo_j.append(np.mean(chunk_j))
            stdAggScore.append(np.std(chunkScores, ddof = 1))
        #Fit mean(aggScore) ~ mean(minCosmo) and std(aggScore) ~ mean(minCosmo).
        print('Fitting mean null aggregation scores...') #Use rpy2 for R smooth.spline function, which works... smoothlier than scipy's equivalent.
        meanSpline = robjects.r['smooth.spline'](robjects.FloatVector(meanMinCosmo),robjects.FloatVector(np.log(meanAggScore)))
        fittedMeans = np.exp(robjects.r['predict'](meanSpline, robjects.FloatVector(minCosmo))[1])
        print('Fitting std null aggregation scores...')
        stdSpline = robjects.r['smooth.spline'](robjects.FloatVector(meanMinCosmo),robjects.FloatVector(np.log(stdAggScore)))
        fittedStd = np.exp(robjects.r['predict'](stdSpline, robjects.FloatVector(minCosmo))[1])
        #Calculate Zscores for all our null model aggScores.
        Zscores = (scores - fittedMeans)/fittedStd
        meanScoreFit = robjects.r['predict'](meanSpline, robjects.FloatVector(range(1, int(max(meanMinCosmo)))))[1]
        stdScoreFit = robjects.r['predict'](stdSpline, robjects.FloatVector(range(1, int(max(meanMinCosmo)))))[1]
        #Get Zscore cutoff for the desired pVal.
        cutoffIndex = round(len(Zscores) * pVal)
        cutoff = sorted(Zscores, reverse = True)[cutoffIndex]
##        with open('minCosmo', 'w') as o:
##            for c in minCosmo:
##                o.write('{}\n'.format(c))
##        with open('cosmo_i', 'w') as o:
##            for c in cosmo_i:
##                o.write('{}\n'.format(c))
##        with open('cosmo_j', 'w') as o:
##            for c in cosmo_j:
##                o.write('{}\n'.format(c))
##        with open('scores', 'w') as o:
##            for s in scores:
##                o.write('{}\n'.format(s))
##        with open('cooc', 'w') as o:
##            for c in cooc:
##                o.write('{}\n'.format(c))
##        with open('meanMinCosmo', 'w') as o:
##            for c in meanMinCosmo:
##                o.write('{}\n'.format(c))
##        with open('meanAggScore', 'w') as o:
##            for c in meanAggScore:
##                o.write('{}\n'.format(c))
##        with open('meanCooc', 'w') as o:
##            for c in meanCooc:
##                o.write('{}\n'.format(c))
##        with open('meanCosmoi', 'w') as o:
##            for c in meanCosmo_i:
##                o.write('{}\n'.format(c))
##        with open('meanCosmoj', 'w') as o:
##            for c in meanCosmo_j:
##                o.write('{}\n'.format(c))
##        with open('stdAggScore', 'w') as o:
##            for c in stdAggScore:
##                o.write('{}\n'.format(c))
##        with open('Zscores','w') as o:
##            for c in Zscores:
##                o.write('{}\n'.format(c))
##        with open('meanScoreFit','w') as o:
##            for c in meanScoreFit:
##                o.write('{}\n'.format(c))
##        with open('stdScoreFit','w') as o:
##            for c in stdScoreFit:
##                o.write('{}\n'.format(c))

        return cutoff, meanSpline, stdSpline
    


    def correctScores(self, taxon = False, idx = 0):
        """
        - Correct aggregation scores to eliminate their dependency from minCosmo.
        - aggScores are biased towards pairs with a low mininum cosmopolitanism (see doc of self.fitZscoreFunction). Null model aggScores (which should be
          unbiased) decrease exponentially for increasing values of minCosmo.
        - We will normalize the aggScore to a Zscore, calculating the null-model mean and the null-model standard deviation from the minimum cosmopolitanism
          of the pair, using the fits calculated when running self.fitZscoreFunction.
        - We will then substract the cutoff Zscore for the requested pVal, as calculated when running self.fitZscoreFunction.
        - If called with taxon == False (initial call by the self.__init__ method after fitting the threshold function) run for the whole matrix and store the
          resulting dataframe in the self.correctedAggScoreMatrix attribute.
        - If called with taxon == True (by the self.mergeTaxa method) run only for the requested and return a numpy array with the corrected aggScores
          for that taxon (and zeros for the rest).
        """
        
        cutoff, meanSpline, stdSpline = self.cutoff, self.meanSpline, self.stdSpline
        aggScoreMatrix = self.aggScoreMatrix.values
        minCosmoMatrix = self.minCosmoMatrix.values
        correctedScores = np.zeros(aggScoreMatrix.shape)

        if not taxon:
            print('Correcting aggregation scores...')
            minCosmoV = minCosmoMatrix.flatten() #Flatten the minCosmo array so it can be passed to R in one go.
            means = np.exp(robjects.r['predict'](meanSpline, robjects.FloatVector(minCosmoV))[1]) #Calculate the means with R.
            means.shape = minCosmoMatrix.shape #And reshape back to a 2d array.
            stds = np.exp(robjects.r['predict'](stdSpline, robjects.FloatVector(minCosmoV))[1])
            stds.shape = minCosmoMatrix.shape
            for i in range(correctedScores.shape[0]):
                for j in range(correctedScores.shape[1]):
                    if not correctedScores[i,j] and aggScoreMatrix[i,j]:
                        aggScore = aggScoreMatrix[i,j]
                        minCosmo = minCosmoMatrix[i,j]
                        mean = means[i,j] 
                        std = stds[i,j]
                        Z = (aggScore - mean) / std
                        correctedScore = Z - cutoff
                        correctedScores[i,j] = correctedScore
                        correctedScores[j,i] = correctedScore
            self.correctedAggScoreMatrix = DataFrame(correctedScores, index = self.aggScoreMatrix.index, columns = self.aggScoreMatrix.columns)

        else:
            minCosmoV = minCosmoMatrix[idx]
            #Precalculate means and stds with R.
            means = np.exp(robjects.r['predict'](meanSpline, robjects.FloatVector(minCosmoV))[1]) 
            stds = np.exp(robjects.r['predict'](stdSpline, robjects.FloatVector(minCosmoV))[1])
            for j in range(correctedScores.shape[0]):
                if aggScoreMatrix[idx,j]:
                    aggScore = aggScoreMatrix[idx,j]
                    minCosmo = minCosmoMatrix[idx,j]
                    mean = means[j]
                    std = stds[j]
                    Z = (aggScore - mean) / std
                    correctedScore = Z - cutoff
                    correctedScores[idx,j] = correctedScore
                    correctedScores[j,idx] = correctedScore
            return correctedScores


    def clearNullMatrices(self):
        """Clear the null matrices to free memory."""
        del self.nullMatrices


    def mergeTaxa(self, taxon1, taxon2):
        """
        Merge two taxa (or groups of taxa) u, v according to the following rules:
        - A new taxon uv is created. For each sample i:
            - Auvi = Aui * Avi
            - πuvi = πui * πvi
        - Taxa u, v are substituted by u', v'. u' contains the samples in which u did appear but v didn't, and vice-versa.
            - Au'i = Aui * (1 - Avi)
            - πu'i = πui * (1 - πvi)
            - Av'i = Avi * (1 - Aui)
            - πv'i = πvi * (1 - πui)

        Recalculate the co-occurrence and aggregation score matrices according to the new taxa.    
        """
        mergedName = '-'.join(sorted({taxon for taxon in (taxon1 + '-' + taxon2).split('-')}))
        allTaxa = list(self.paMatrix.index) + [mergedName]
        samples = self.paMatrix.columns
        u = self.coocMatrix.index.get_loc(taxon1)
        v = self.coocMatrix.index.get_loc(taxon2)
        uv = self.paMatrix.shape[0] #We will place the new combined taxon at the end of all matrices.
        #Get dereferenced copies of the Pandas Dataframes values.
        paMatrix = np.array(self.paMatrix.values, dtype = np.int)
        probMatrix = np.array(self.probMatrix.values, dtype = np.float128)
        coocMatrix = np.array(self.coocMatrix.values, dtype = np.int)
        minCosmoMatrix = np.array(self.minCosmoMatrix.values, dtype = np.int)
        aggScoreMatrix = np.array(self.aggScoreMatrix.values, dtype = np.double)
        correctedAggScoreMatrix = np.array(self.correctedAggScoreMatrix.values, dtype = np.double)
        #Resize matrices.
        paMatrix = np.vstack([paMatrix, np.zeros(paMatrix.shape[1], dtype = np.int)])
        probMatrix = np.vstack([probMatrix, np.zeros(probMatrix.shape[1], dtype = np.float128)])
        coocMatrix = np.vstack([coocMatrix, np.zeros(coocMatrix.shape[1], dtype = np.int)])
        coocMatrix = np.hstack([coocMatrix, np.zeros(shape = [coocMatrix.shape[0], 1], dtype = np.int)])
        minCosmoMatrix = np.vstack([minCosmoMatrix, np.zeros(minCosmoMatrix.shape[1], dtype = np.int)])
        minCosmoMatrix = np.hstack([minCosmoMatrix, np.zeros(shape = [minCosmoMatrix.shape[0], 1], dtype = np.int)])
        aggScoreMatrix = np.vstack([aggScoreMatrix, np.zeros(aggScoreMatrix.shape[1], dtype = np.double)])
        aggScoreMatrix = np.hstack([aggScoreMatrix, np.zeros(shape = [aggScoreMatrix.shape[0], 1], dtype = np.double)])
        correctedAggScoreMatrix = np.vstack([correctedAggScoreMatrix, np.zeros(correctedAggScoreMatrix.shape[1], dtype = np.double)])
        correctedAggScoreMatrix = np.hstack([correctedAggScoreMatrix, np.zeros(shape = [correctedAggScoreMatrix.shape[0], 1], dtype = np.double)])
        #Update presence/absence and probability matrices.
        A_uvi = paMatrix[u] * paMatrix[v]
        pi_uvi = probMatrix[u] * probMatrix[v]
        A_uprimei = paMatrix[u] * (1 - paMatrix[v])
        pi_uprimei = probMatrix[u] * (1 - probMatrix[v])
        A_vprimei = paMatrix[v] * (1 - paMatrix[u])
        pi_vprimei = probMatrix[v] * (1 - probMatrix[u])
        paMatrix[uv] = A_uvi
        probMatrix[uv] = pi_uvi
        paMatrix[u] = A_uprimei
        probMatrix[u] = pi_uprimei
        paMatrix[v] = A_vprimei
        probMatrix[v] = pi_vprimei
        #Add new rows and columns to co-occurrence matrix.
        cooc_uvi = np.zeros(paMatrix.shape[0])
        for k in range(paMatrix.shape[0]):
            cooc_uvi[k] = sum(paMatrix[uv] * paMatrix[k])
        cooc_uprimei = np.zeros(paMatrix.shape[0])
        for k in range(paMatrix.shape[0]):
            cooc_uprimei[k] = sum(paMatrix[u] * paMatrix[k])
        cooc_vprimei = np.zeros(paMatrix.shape[0])
        for k in range(paMatrix.shape[0]):
            cooc_vprimei[k] = sum(paMatrix[v] * paMatrix[k])
        coocMatrix[uv] = cooc_uvi
        coocMatrix[:,uv] = cooc_uvi
        coocMatrix[u] = cooc_uprimei
        coocMatrix[:,u] = cooc_uprimei
        coocMatrix[v] = cooc_vprimei
        coocMatrix[:,v] = cooc_vprimei
        #Update cosmopolitanism.
        cosmo = np.sum(paMatrix, axis = 1)
        #Add new rows and columns to minCosmoMatrix.
        minCosmo_uvi = np.zeros(paMatrix.shape[0])
        for k in range(paMatrix.shape[0]):
            minCosmo_uvi[k] = min(cosmo[uv], cosmo[k])
        minCosmo_uprimei = np.zeros(paMatrix.shape[0])
        for k in range(paMatrix.shape[0]):
            minCosmo_uprimei[k] = min(cosmo[u], cosmo[k])
        minCosmo_vprimei = np.zeros(paMatrix.shape[0])
        for k in range(paMatrix.shape[0]):
            minCosmo_vprimei[k] = min(cosmo[v], cosmo[k])
        minCosmoMatrix[uv] = cooc_uvi
        minCosmoMatrix[:,uv] = cooc_uvi
        minCosmoMatrix[u] = cooc_uprimei
        minCosmoMatrix[:,u] = cooc_uprimei
        minCosmoMatrix[v] = cooc_vprimei
        minCosmoMatrix[:,v] = cooc_vprimei
        #Add new rows and columns to aggScoreMatrix.
        aggScore_uvi = CbuildAggScoreMatrix(paMatrix, coocMatrix, probMatrix, True, uv)
        aggScore_uprimei = CbuildAggScoreMatrix(paMatrix, coocMatrix, probMatrix, True, u)
        aggScore_vprimei = CbuildAggScoreMatrix(paMatrix, coocMatrix, probMatrix, True, v)
        #Note that buildAggScoreMatrix returns complete matrices, even if only the requested row has been filled.
        aggScoreMatrix[uv] = aggScore_uvi[uv]
        aggScoreMatrix[:,uv] = aggScore_uvi[uv]
        aggScoreMatrix[u] = aggScore_uprimei[u]
        aggScoreMatrix[:,u] = aggScore_uprimei[u]
        aggScoreMatrix[v] = aggScore_vprimei[v]
        aggScoreMatrix[:,v] = aggScore_vprimei[v]
        #Recast pandas DataFrames (except for correctedAggScoreMatrix).
        self.paMatrix = DataFrame(paMatrix, index = allTaxa, columns = samples)
        self.probMatrix = DataFrame(probMatrix, index = allTaxa, columns = samples)
        self.coocMatrix = DataFrame(coocMatrix, index = allTaxa, columns = allTaxa)
        self.minCosmoMatrix = DataFrame(minCosmoMatrix, index = allTaxa, columns = allTaxa)
        self.aggScoreMatrix = DataFrame(aggScoreMatrix, index = allTaxa, columns = allTaxa)
        #Correct aggregation scores.
        correctedAggScore_uvi = self.correctScores(True, uv)
        correctedAggScore_uprimei = self.correctScores(True, u)
        correctedAggScore_vprimei = self.correctScores(True, v)
        #Note that correctScores returns complete matrices, even if only the requested row has been filled.
        correctedAggScoreMatrix[uv] = correctedAggScore_uvi[uv]
        correctedAggScoreMatrix[:,uv] = correctedAggScore_uvi[uv]
        correctedAggScoreMatrix[u] = correctedAggScore_uprimei[u]
        correctedAggScoreMatrix[:,u] = correctedAggScore_uprimei[u]
        correctedAggScoreMatrix[v] = correctedAggScore_vprimei[v]
        correctedAggScoreMatrix[:,v] = correctedAggScore_vprimei[v]
        #Recast as panda DataFrame.
        self.correctedAggScoreMatrix = DataFrame(correctedAggScoreMatrix, index = allTaxa, columns = allTaxa)
