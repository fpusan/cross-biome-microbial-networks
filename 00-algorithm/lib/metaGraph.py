import sys
from multiprocessing import Pool

import datetime
import cProfile
from copy import copy, deepcopy

import MySQLdb
import MySQLdb.cursors

import networkx as nx

from pandas import DataFrame
import numpy as np

from lib.model import Model
from lib.node import Node
from lib.graph import Graph
from lib.XGMML import write


class MetaGraph:
    """
    """
    def __init__(self, args):
        self.args = args
        self.samples = dict()
        self.environments = {'supertype': dict(), 'type': dict(), 'subtype': dict()} #self.environments['supertype']:{'supertype1': set([samples1])} etc...
        self.baseNodes = dict() #{'nodeName': Node object} . Will get filled with self.loadData().
        self.metaGraph = nx.DiGraph(**vars(self.args))
        self.loadData()
        self.subTypeNames = sorted(self.environments['subtype'].keys())
        ###Run the null model.
        self.baseModel = self.runModel()
        ###Write matrices.
        self.baseModel.write('pa', '{}.pa.tsv'.format(self.args.output))
##        self.baseModel.write('minCosmo')
##        self.baseModel.write('cooc')
##        self.baseModel.write('prob')
##        self.baseModel.write('agg')
##        self.baseModel.write('corAgg')
##        for i, nullMat in enumerate(self.baseModel.nullMatrices):
##            nullMat.write('pa', 'nullMat{}.tsv'.format(i))

        self.baseModel.clearNullMatrices()
        ###Generate the random-path graphs.
        pool = Pool(self.args.processors)
        #Deep copy the model object, as it will get modified independently during each clustering process.
        #Shallow copy the node baseNode list: new nodes will be added during each clustering process, but the original ones will remain unchanged.
        #self.graphs = [Graph(deepcopy(self.baseModel), copy(self.baseNodes), self.samples, args, i) for i in range(self.args.graphs)]
        #Probably this is no longer after switching to multiprocessing, everything should be getting deepcopied now. But meh.
        results = pool.starmap_async(Graph, [(deepcopy(self.baseModel), copy(self.baseNodes), self.samples, args, i) for i in range(self.args.graphs)])
        pool.close()
        pool.join()
        self.graphs = results.get()
        ###Merge them into a metagraph.
        for graphObj in self.graphs:
            graph = graphObj.graph
            for node in graph.nodes():
                if node in self.metaGraph.node:
                    self.metaGraph.node[node]['support'] += 1
                else:
                    self.metaGraph.add_node(node, support = 1)
            for source, target in graph.edges():
                if (source, target) not in self.metaGraph.edges():
                    self.metaGraph.add_edge(source, target, support = 1)
                else:
                    self.metaGraph[source][target]['support'] += 1
        ###Generate the metadata for the nodes in the metagraph.
        for node in self.metaGraph.node.keys():
            components = node.split('-')
            nodeObj = self.baseNodes[components[0]]
            for comp in components[1:]:
                prevName = nodeObj.name
                mergedName = '-'.join(sorted({taxon for taxon in (prevName + '-' + comp).split('-')}))
                newNode = self.baseNodes[comp]
                nodeObj = Node(mergedName, self, (nodeObj, newNode))
            #Add metadata to the metagraph. It's oh so meta.
            self.metaGraph.node[node]['sources'] = sorted(['%s | %d'%(source, numSamples) for source, numSamples in nodeObj.sources.items()],
                                                                 key = lambda x: int(x.split('| ')[1]), reverse = True)
            self.metaGraph.node[node]['titles'] = sorted(['%s | %d'%(title, numSamples) for title, numSamples in nodeObj.titles.items()],
                                                                 key = lambda x: int(x.split('| ')[1]), reverse = True)
            self.metaGraph.node[node]['pmids'] = sorted(['%s | %d'%(pmid, numSamples) for pmid, numSamples in nodeObj.pmids.items()],
                                                                 key = lambda x: int(x.split('| ')[1]), reverse = True)
            self.metaGraph.node[node]['envSupertypes'] = sorted(['%s | %d'%(env, numSamples) for env, numSamples in nodeObj.envSupertypes.items()],
                                                                       key = lambda x: int(x.split('| ')[1]), reverse = True)
            self.metaGraph.node[node]['envTypes'] = sorted(['%s | %d'%(env, numSamples) for env, numSamples in nodeObj.envTypes.items()],
                                                                  key = lambda x: int(x.split('| ')[1]), reverse = True)
            self.metaGraph.node[node]['envSubtypes'] = sorted(['%s | %d'%(env, numSamples) for env, numSamples in nodeObj.envSubtypes.items()],
                                                                     key = lambda x: int(x.split('| ')[1]), reverse = True)
            self.metaGraph.node[node]['samples'] = sorted(map(int, nodeObj.inSamples))
            self.metaGraph.node[node]['numSamples'] = len(nodeObj.inSamples)
            self.metaGraph.node[node]['numTaxa'] = len(nodeObj.components)
            self.metaGraph.node[node]['consensusEnvSupertype'] = nodeObj.consensusEnvSupertype
            self.metaGraph.node[node]['consensusEnvType'] = nodeObj.consensusEnvType
            self.metaGraph.node[node]['consensusEnvSubtype'] = nodeObj.consensusEnvSubtype
            self.metaGraph.node[node]['taxa'] = sorted([node.name for node in nodeObj.components])
            self.metaGraph.node[node]['pieSubtypeNames'] = self.subTypeNames #Pie chart metadata.
            self.metaGraph.node[node]['pieSubtype'] = [nodeObj.envSubtypes[s] if s in nodeObj.envSubtypes else 0 for s in self.subTypeNames]

        ###Count the number of outgoing edges.
        for node in self.metaGraph.nodes():
            children = self.metaGraph.successors(node)
            for child in children:
                self.metaGraph[node][child]['parentOutEdges'] = len(children)
            self.metaGraph.node[node]['outEdges'] = len(children)

        ###Add metadata on sample sources and environments.
        sampleEnvs = [' | '.join((str(sample), self.samples[sample]['envs'][2])) for sample in self.samples]
        self.metaGraph.graph['sampleEnvs'] = sampleEnvs
        sampleSources = [' | '.join((str(sample), str(self.samples[sample]['source']))) for sample in self.samples]
        self.metaGraph.graph['sampleSources'] = sampleSources

        #Write output.
        write(self.metaGraph, self.args.output + '.xml')
        nx.write_gpickle(self.metaGraph, self.args.output + '.gpickle')



    def loadData(self):
        """Load input data"""
        if self.args.input:
            self.loadTaxaFromInputTable(self.args.input)
        else:
            self.loadTaxaFromMICRODB()      
        print('Loaded {:d} taxa and {:d} samples'.format(len(self.baseNodes), len(self.samples)))


    def buildOutputName(self):
        if self.args.output:
            return self.args.output
        else:
            if self.args.input:
                return '{}.xml'.format('.'.join(self.args.input.split('.')[:-1]))
            else:
                return 'comm-builder_MICRODB_{:%Y-%m-%d_%H:%M:%S}.xml'.format(datetime.datetime.now())
            return '.'.join(outputName)
             
        
    def loadTaxaFromMICRODB(self):
        """Parse the data contained in the microdb database, using the provided constraints."""
        microdb = MySQLdb.connect(host = 'localhost', user = 'user', db = 'microdb_frozen', cursorclass = MySQLdb.cursors.SSCursor)
        c = microdb.cursor()
        args = self.args
        args.tax_level = 'ordertx' if args.tax_level == 'order' else args.tax_level #Since order is a mysql reserved word.
        ###Build whitelists for envsupertypes, envtypes and envsubtypes, to avoid SQL injections.
        c.execute("SELECT envsupertype FROM AMBIENTE")
        envSupertypeWhiteList = {i[0] for i in c}
        c.execute("SELECT envtype FROM AMBIENTE")
        envTypeWhiteList = {i[0] for i in c}
        c.execute("SELECT envsubtype FROM AMBIENTE")
        envSubtypeWhiteList = {i[0] for i in c}
        ###Get sequence clusters belonging to the selected samples.
        sqlstring = """SELECT TAXON.{}, SECUENCIA.ID_muestra, AMBIENTE.isource, AMBIENTE.envsupertype, AMBIENTE.envtype, AMBIENTE.envsubtype, MUESTRA.titulo, MUESTRA.pmid \
                       FROM SECUENCIA, MUESTRA, AMBIENTE, TAXON WHERE \
                       SECUENCIA.ID_muestra = AMBIENTE.ID_muestra AND MUESTRA.ID_muestra = AMBIENTE.ID_muestra AND \
                       TAXON.ID_cluster = SECUENCIA.ID_cluster""".format(args.tax_level)   
        if args.env_supertype:
            for esupertype in args.env_supertype:
                if esupertype not in envSupertypeWhiteList:
                    raise ValueError('{} is not a valid environmental supertype'.format(esupertype))
            else:
                sqlstring += ' AND AMBIENTE.envsupertype IN ("{}")'.format('", "'.join(args.env_supertype))
        if args.env_type:
            for etype in args.env_type:
                if etype not in envTypeWhiteList:
                    raise ValueError('{} is not a valid environmental type'.format(etype))
            else:
                sqlstring += ' AND AMBIENTE.envtype IN ("{}")'.format('", "'.join(args.env_type))
        if args.env_subtype:
            for esubtype in args.env_subtype:
                if esubtype not in envSubtypeWhiteList:
                    raise ValueError('{} is not a valid environmental subtype'.format(esubtype))
            else:
                sqlstring += ' AND AMBIENTE.envsubtype IN ("{}")'.format('", "'.join(args.env_subtype))
        #sql injection can not happen in the above lines, since tax_level can only take fixed values, and env_supertype, env_type and env_subtype are whitelisted.
        c.execute(sqlstring)
        taxonSampleSourcesEnv = [taxonSampleSourceEnv for taxonSampleSourceEnv in c if taxonSampleSourceEnv[0] not in ('Unresolved', 'Inconsistent', None)]
        c.close()
        ###Limit our search to the selected taxa, if required.
        if args.use_taxon:
            sqlstring = """SELECT {} FROM TAXON WHERE {} = """.format(args.tax_level, args.use_taxon_tax_level)
            c.execute(sqlstring + '%s', (args.use_taxon,))
            goodTaxa = {taxon[0] for taxon in c}
            taxonSampleSourcesEnv = [taxonSampleSourceEnv for taxonSampleSourceEnv in taxonSampleSources if taxonSampleSourceEnv[0] in goodTaxa]
        ###Remove samples with less than --min_richness taxa.
        sampleRichness = dict()
        for taxon, sample, source, *envs, title, pmid in taxonSampleSourcesEnv:
            if sample not in sampleRichness:
                sampleRichness[sample] = {taxon}
            else:
                sampleRichness[sample].add(taxon)
        taxonSampleSourcesEnv = [(taxon, sample, source, *envs, title, pmid) for taxon, sample, source, *envs, title, pmid in taxonSampleSourcesEnv
                                 if len(sampleRichness[sample]) >= args.min_richness]
        ###Filter out taxa appearing in less than --min_cosmopolitanism samples.
        taxaCosmo = dict()
        for taxon, sample, source, *envs, title, pmid in taxonSampleSourcesEnv:
            if taxon not in taxaCosmo:
                taxaCosmo[taxon] = {sample}
            else:
                taxaCosmo[taxon].add(sample)
        taxonSampleSourcesEnv = [(taxon, sample, source, *envs, title, pmid) for taxon, sample, source, *envs, title, pmid in taxonSampleSourcesEnv
                                 if len(taxaCosmo[taxon]) >= args.min_cosmopolitanism]  
        ###To which distinct environments do those samples belong?
        tempBuffer = list()
        for taxon, sample, source, envSupertype, envType, envSubtype, title, pmid in taxonSampleSourcesEnv:
            #Translate NULL mysql values into strings.
            envSupertype = envSupertype if envSupertype else 'NA'
            envType = envType if envType else 'NA'
            envSubtype = envSubtype if envSubtype else 'NA'
            #Add full hierarchy info to children environmental levels.
            envSubtype = '-'.join([envSupertype, envType, envSubtype])
            envType = '-'.join([envSupertype, envType])
            #Classify the sample.
            if envSupertype not in self.environments['supertype']:
                self.environments['supertype'][envSupertype] = {sample}
            else:
                self.environments['supertype'][envSupertype].add(sample)
            if envType not in self.environments['type']:
                self.environments['type'][envType] = {sample}
            else:
                self.environments['type'][envType].add(sample)
            if envSubtype not in self.environments['subtype']:
                self.environments['subtype'][envSubtype] = {sample}
            else:
                self.environments['subtype'][envSubtype].add(sample)
            #Store corrected info into an ugly temporary buffer.
            tempBuffer.append((taxon, sample, source, envSupertype, envType, envSubtype, title, pmid))
        taxonSampleSourcesEnv = [i for i in tempBuffer]
        ###Filter out environments with less than --min_samples samples. We do this only for the level that we will use for sample partitioning. 
        filteredEnvs = [env for env, samples in self.environments[args.env_partition].items() if len(samples) < args.min_samples]
        for env in filteredEnvs:
            del self.environments[args.env_partition][env]
        if args.env_partition == 'supertype':
            taxonSampleSourcesEnv = [(taxon, sample, source, *envs, title, pmid) for taxon, sample, source, *envs, title, pmid in taxonSampleSourcesEnv
                                     if envs[0] in self.environments['supertype']]
        elif args.env_partition == 'type':
            taxonSampleSourcesEnv = [(taxon, sample, source, *envs, title, pmid) for taxon, sample, source, *envs, title, pmid in taxonSampleSourcesEnv
                                     if envs[1] in self.environments['type']]
        else:
            taxonSampleSourcesEnv = [(taxon, sample, source, *envs, title, pmid) for taxon, sample, source, *envs, title, pmid in taxonSampleSourcesEnv
                                     if envs[2] in self.environments['subtype']]
        if len(self.environments[args.env_partition]) < args.min_ubiquity:
            raise Exception('--min_ubiquity is set to {:d}. However, with your current settings, only {:d} different environments were selected.'.format(args.min_ubiquity, len(self.environments))) 
        ###Filter out taxa appearing in less than --min_ubiquity environments.
        taxaUbiquity = dict()
        for taxon, sample, source, envSupertype, envType, envSubtype, title, pmid in taxonSampleSourcesEnv:
            if args.env_partition == 'supertype': env = envSupertype
            elif args.env_partition == 'type': env = envType
            else: env = envSubtype
            if taxon not in taxaUbiquity:
                taxaUbiquity[taxon] = {env}
            else:
                taxaUbiquity[taxon].add(env)
        taxonSampleSourcesEnv = [(taxon, sample, source, *envs, title, pmid) for taxon, sample, source, *envs, title, pmid in taxonSampleSourcesEnv
                                 if len(taxaUbiquity[taxon]) >= args.min_ubiquity]
        ###Which samples are we finally including?
        self.samples = {taxonSampleSourceEnv[1]: {'source': taxonSampleSourceEnv[2], 'envs': taxonSampleSourceEnv[3:6], 'title': taxonSampleSourceEnv[6], 'pmid': taxonSampleSourceEnv[7]}
                        for taxonSampleSourceEnv in taxonSampleSourcesEnv}
        for sample in self.samples:
            if not self.samples[sample]: #If AMBIENTE.isource was NULL.
                self.samples[sample] = 'unknown'
        ###Populate graph with taxon nodes.
        for taxon, sample, source, *envs, title, pmid in taxonSampleSourcesEnv:
            if taxon not in self.baseNodes:
                self.baseNodes[taxon] = Node(taxon, self, None)
            self.baseNodes[taxon].addSample(sample)
        for node in self.baseNodes.values():
            node.addSourcesAndEnvs()


    def loadTaxaFromInputTable(self, fileName):
        """
        Parse taxa and sample information from a tab-formatted input table.
        TAXA AND RICHNESS FILTERS NOT IMPLEMENTED.
        """
        self.samples = dict()
        self.baseNodes = dict()
        #Rows are taxons, columns are samples.
        with open(fileName) as infile:
            self.samples = infile.readline().strip().split('\t')
            taxonSamplePairs = list()
            for line in infile:
                line = line.strip().split('\t')
                taxon = line[0]
                abundances = map(int, map(float, line[1:]))
                for sample, abundance in zip(self.samples, abundances):
                    if abundance:
                        if taxon not in self.baseNodes:
                            self.baseNodes[taxon] = Node(taxon, self, None)
                        self.baseNodes[taxon].addSample(sample)


    def runModel(self):
        """
        Build a presence/absence matrix from the data stored in self.baseNodes. Taxa in rows, samples in columns.
        Return a Model object containing the aggregation Score matrix. 
        """
        paMatrix = DataFrame(np.zeros([len(self.baseNodes), len(self.samples)]), index = sorted(self.baseNodes), columns = sorted(self.samples), dtype = np.int)
        for node in self.baseNodes.values():
            for sample in node.inSamples:
                paMatrix.ix[node.name, sample] = 1

        if not self.args.input:
            return Model(paMatrix, self.environments[self.args.env_partition], self.args.pval_cutoff, self.args.bootstrap, self.args.processors)
        else:
            return Model(paMatrix, None, self.args.pval_cutoff, self.args.bootstrap, self.args.processors)
