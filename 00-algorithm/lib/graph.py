import sys

import datetime
from random import seed, shuffle, uniform

import MySQLdb
import MySQLdb.cursors

import networkx as nx

import numpy as np

from lib.node import Node


class Graph:
    """
    """
    def __init__(self, baseModel, baseNodes, samples, args, i, verbose = False):
        print('Generating graph {}'.format(i))
        #seed() #Since worker process inherits the random number generator from the parent process, every worker would have the same output.
                #Turns out this is not happening here, it happened however when generating null matrices in the Model class.
        #print(uniform(0,1))
        self.envData = baseModel #Model object
        self.nodes = baseNodes #{'nodeName': Node object}
        self.samples = samples
        self.args = args
        self.graph = nx.DiGraph(**vars(self.args))
        self.simplifiedGraph = nx.DiGraph(**vars(self.args))
        self.verbose = verbose
        self.build()


    def build(self):
        while True:
            #Return indices of the flattened coocMatrix, ordered from more to less samples.
            coocMatrix = self.envData.coocMatrix.values
            corrAggScoreMatrix = self.envData.correctedAggScoreMatrix.values
            cutoff = self.envData.cutoff
            sortMatrix = corrAggScoreMatrix if not self.args.cluster_by_size else coocMatrix
            sortedIndices = np.argsort(sortMatrix, axis = None)[::-1]

            #Keep only indices and scores of significant aggregations in the upper diagonal.
            upper_diagonal = [i for i in zip(*np.triu_indices(sortMatrix.shape[0]))]
            def isSig(index):
                coords = np.unravel_index(index, sortMatrix.shape)
                numSamples = int(coocMatrix[coords])
                score = float(corrAggScoreMatrix[coords])
                return (numSamples >= self.args.size_cutoff and score > 0 and coords in upper_diagonal)
            
            sigIndices = [index for index in sortedIndices if isSig(index)] #This was sorted in desdending order previously.
            sigScores = [float(corrAggScoreMatrix[np.unravel_index(index, sortMatrix.shape)]) + cutoff for index in sigIndices]

            #Break if there are no pairs co-occurring in enough samples and with enough score.
            if not sigIndices:
                if self.verbose:
                    print()
                break
            
            #Otherwise:
            if self.args.best_scores:
                #Keep the best N scores and sort them randomly.
                sigIndices = sigIndices[:self.args.best_scores]
                shuffle(sigIndices)
                selectedIndex = sigIndices[0]

            else:
                #Randomly select an index, weighted by the score associated to it.
                #Note that we added back the cutoff when building sigScores, since our corrected scores had it substracted so corrScore > 0
                #was considered as significant, but this changes the relative proportion between scores.
                #Also, we are weighting by e**score, since scores are in a logarithmic scale.
                sigScores = np.exp(sigScores) 
                accumTarget = uniform(0, sum(sigScores))
                accumWeight = 0
                for score, index in zip(sigScores, sigIndices):
                    accumWeight += score
                    if accumWeight >= accumTarget:
                        selectedIndex = index
                        break
                else:
                    raise Exception('Wtf?')

            #Join the selected pair of nodes.
            coords = np.unravel_index(selectedIndex, sortMatrix.shape)
            bestCorr = (self.envData.correctedAggScoreMatrix.index[coords[0]], self.envData.correctedAggScoreMatrix.columns[coords[1]])
            numSamples = int(coocMatrix[coords])
            bestCorrScore = float(corrAggScoreMatrix[coords])
            percentSig = ((np.exp(bestCorrScore+cutoff))/sum(sigScores))*100
            if self.verbose:
                print('\t'.join([str(bestCorr), '{:.2f}'.format(bestCorrScore), str(numSamples), '{:.2f}'.format(percentSig)]))
            #Update the PA, prob, cooc and aggScore matrices.
            self.envData.mergeTaxa(*bestCorr)
            #Create a new node.
            mergedName = '-'.join(sorted({taxon for taxon in (bestCorr[0] + '-' + bestCorr[1]).split('-')}))
            modelCoocSamples = set(self.envData.paMatrix.columns[self.envData.paMatrix.loc[mergedName]>0]) #The samples used to calculate the aggregation score.
            assert len(modelCoocSamples) == numSamples
            parents = (self.nodes[bestCorr[0]], self.nodes[bestCorr[1]])
            self.nodes[mergedName] = Node(mergedName, self, parents, modelCoocSamples)
            #Populate networkx graph.           
            self.graph.add_edge(parents[0].name, mergedName, numSamples = numSamples, aggScore = bestCorrScore)
            self.graph.add_edge(parents[1].name, mergedName, numSamples = numSamples, aggScore = bestCorrScore)
            #First parent.
            self.graph.node[parents[0].name]['sources'] = sorted(['%s | %d'%(source, numSamples) for source, numSamples in parents[0].sources.items()],
                                                                 key = lambda x: int(x.split('| ')[1]), reverse = True)
            self.graph.node[parents[0].name]['titles'] = sorted(['%s | %d'%(title, numSamples) for title, numSamples in parents[0].titles.items()],
                                                                 key = lambda x: int(x.split('| ')[1]), reverse = True)
            self.graph.node[parents[0].name]['pmids'] = sorted(['%s | %d'%(pmid, numSamples) for pmid, numSamples in parents[0].pmids.items()],
                                                                 key = lambda x: int(x.split('| ')[1]), reverse = True)
            self.graph.node[parents[0].name]['envSupertypes'] = sorted(['%s | %d'%(env, numSamples) for env, numSamples in parents[0].envSupertypes.items()],
                                                                       key = lambda x: int(x.split('| ')[1]), reverse = True)
            self.graph.node[parents[0].name]['envTypes'] = sorted(['%s | %d'%(env, numSamples) for env, numSamples in parents[0].envTypes.items()],
                                                                  key = lambda x: int(x.split('| ')[1]), reverse = True)
            self.graph.node[parents[0].name]['envSubtypes'] = sorted(['%s | %d'%(env, numSamples) for env, numSamples in parents[0].envSubtypes.items()],
                                                                     key = lambda x: int(x.split('| ')[1]), reverse = True)
            self.graph.node[parents[0].name]['samples'] = sorted(map(int, parents[0].inSamples))
            self.graph.node[parents[0].name]['numSamples'] = len(parents[0].inSamples)
            self.graph.node[parents[0].name]['numTaxa'] = len(parents[0].components)
            self.graph.node[parents[0].name]['consensusEnvSupertype'] = parents[0].consensusEnvSupertype
            self.graph.node[parents[0].name]['consensusEnvType'] = parents[0].consensusEnvType
            self.graph.node[parents[0].name]['consensusEnvSubtype'] = parents[0].consensusEnvSubtype
            self.graph.node[parents[0].name]['taxa'] = sorted([node.name for node in parents[0].components])
            #Second parent.
            self.graph.node[parents[1].name]['sources'] = sorted(['%s | %d'%(source, numSamples) for source, numSamples in parents[1].sources.items()],
                                                                 key = lambda x: int(x.split('| ')[1]), reverse = True)
            self.graph.node[parents[1].name]['titles'] = sorted(['%s | %d'%(title, numSamples) for title, numSamples in parents[1].titles.items()],
                                                                 key = lambda x: int(x.split('| ')[1]), reverse = True)
            self.graph.node[parents[1].name]['pmids'] = sorted(['%s | %d'%(pmid, numSamples) for pmid, numSamples in parents[1].pmids.items()],
                                                                 key = lambda x: int(x.split('| ')[1]), reverse = True)
            self.graph.node[parents[1].name]['envSupertypes'] = sorted(['%s | %d'%(env, numSamples) for env, numSamples in parents[1].envSupertypes.items()],
                                                                       key = lambda x: int(x.split('| ')[1]), reverse = True)
            self.graph.node[parents[1].name]['envTypes'] = sorted(['%s | %d'%(env, numSamples) for env, numSamples in parents[1].envTypes.items()],
                                                                  key = lambda x: int(x.split('| ')[1]), reverse = True)
            self.graph.node[parents[1].name]['envSubtypes'] = sorted(['%s | %d'%(env, numSamples) for env, numSamples in parents[1].envSubtypes.items()],
                                                                     key = lambda x: int(x.split('| ')[1]), reverse = True)
            self.graph.node[parents[1].name]['samples'] = sorted(map(int, parents[1].inSamples))
            self.graph.node[parents[1].name]['numSamples'] = len(parents[1].inSamples)
            self.graph.node[parents[1].name]['numTaxa'] = len(parents[1].components)
            self.graph.node[parents[1].name]['consensusEnvSupertype'] = parents[1].consensusEnvSupertype
            self.graph.node[parents[1].name]['consensusEnvType'] = parents[1].consensusEnvType
            self.graph.node[parents[1].name]['consensusEnvSubtype'] = parents[1].consensusEnvSubtype
            self.graph.node[parents[1].name]['taxa'] = sorted([node.name for node in parents[1].components])
            #Child node.
            self.graph.node[mergedName]['sources'] = sorted(['%s | %d'%(source, numSamples) for source, numSamples in self.nodes[mergedName].sources.items()],
                                                            key = lambda x: int(x.split('| ')[1]), reverse = True)
            self.graph.node[mergedName]['titles'] = sorted(['%s | %d'%(title, numSamples) for title, numSamples in self.nodes[mergedName].titles.items()],
                                                            key = lambda x: int(x.split('| ')[1]), reverse = True)
            self.graph.node[mergedName]['pmids'] = sorted(['%s | %d'%(pmid, numSamples) for pmid, numSamples in self.nodes[mergedName].pmids.items()],
                                                            key = lambda x: int(x.split('| ')[1]), reverse = True)
            self.graph.node[mergedName]['envSupertypes'] = sorted(['%s | %d'%(env, numSamples) for env, numSamples in self.nodes[mergedName].envSupertypes.items()],
                                                                  key = lambda x: int(x.split('| ')[1]), reverse = True)
            self.graph.node[mergedName]['envTypes'] = sorted(['%s | %d'%(env, numSamples) for env, numSamples in self.nodes[mergedName].envTypes.items()],
                                                             key = lambda x: int(x.split('| ')[1]), reverse = True)
            self.graph.node[mergedName]['envSubtypes'] = sorted(['%s | %d'%(env, numSamples) for env, numSamples in self.nodes[mergedName].envSubtypes.items()],
                                                                key = lambda x: int(x.split('| ')[1]), reverse = True)
            self.graph.node[mergedName]['samples'] = sorted(map(int, self.nodes[mergedName].inSamples))
            self.graph.node[mergedName]['numSamples'] = len(self.nodes[mergedName].inSamples)
            self.graph.node[mergedName]['numTaxa'] = len(self.nodes[mergedName].components)
            self.graph.node[mergedName]['consensusEnvSupertype'] = self.nodes[mergedName].consensusEnvSupertype
            self.graph.node[mergedName]['consensusEnvType'] = self.nodes[mergedName].consensusEnvType
            self.graph.node[mergedName]['consensusEnvSubtype'] = self.nodes[mergedName].consensusEnvSubtype
            self.graph.node[mergedName]['taxa'] = sorted([node.name for node in self.nodes[mergedName].components])

        ###Count the number of outgoing edges.
        for node in self.graph.nodes():
            children = self.graph.successors(node)
            for child in children:
                self.graph[node][child]['parentOutEdges'] = len(children)
            self.graph.node[node]['outEdges'] = len(children)

        ###Free memory.
        del self.envData
        del self.nodes
     
