class Node:
    """Class representing either individual taxa (starting nodes in the directed graph) or combinations of co-occurring taxa."""
    def __init__(self, name, parentGraph, parentNodes, modelCoocSamples = None):
        self.name = name
        self.parentGraph = parentGraph #the graph object that spawned this node.
        self.args = self.parentGraph.args
        self.parentNodes = parentNodes if parentNodes else list() #the name of the nodes pointing to this node.
        self.components = self.getComponents()
        self.inSamples = set() #the samples this node is present in.
        self.envs = dict()
        self.sources = dict()
        self.titles = dict()
        self.pmids = dict()
        self.envs = dict()
        if parentNodes:
            if modelCoocSamples:
                self.inSamples = modelCoocSamples
            else:
                self.inSamples = self.parentNodes[0].inSamples.intersection(self.parentNodes[1].inSamples)
            if not self.parentGraph.args.input: #We are using microdb.
                self.addSourcesAndEnvs()


    def addSample(self, sample):
        """Add an additional sample to the self.inSamples list."""
        self.inSamples.add(sample)


    def addSourcesAndEnvs(self):
        """
        Add source and env tags, based on the samples on which this node appears. Must be called after self.inSamples is full.
        This method is called by the parent graph for the initial nodes, and by __init__ for the nodes generated during the clustering process.
        """
        self.sources = dict()
        self.envSupertypes = dict()
        self.envTypes = dict()
        self.envSubtypes = dict()
        self.consensusEnvSupertype = str()
        self.consensusEnvType = str()
        self.consensusEnvSubtype = str()
        self.titles = dict()
        self.pmids = dict()

        for sample in self.inSamples:
            source = self.parentGraph.samples[sample]['source']
            envSupertype = self.parentGraph.samples[sample]['envs'][0]
            envType = self.parentGraph.samples[sample]['envs'][1]
            envSubtype = self.parentGraph.samples[sample]['envs'][2]
            title = self.parentGraph.samples[sample]['title']
            pmid = self.parentGraph.samples[sample]['pmid']
            if not envSupertype: envSupertype = 'Undetermined'
            if not envType: envType = 'Undetermined'
            if not envSubtype: envSubtype = 'Undetermined'
            if source not in self.sources:
                self.sources[source] = 1
            else:
                self.sources[source] += 1
            if envSupertype not in self.envSupertypes:
                self.envSupertypes[envSupertype] = 1
            else:
                self.envSupertypes[envSupertype] += 1
            if envType not in self.envTypes:
                self.envTypes[envType] = 1
            else:
                self.envTypes[envType] += 1
            if envSubtype not in self.envSubtypes:
                self.envSubtypes[envSubtype] = 1
            else:
                self.envSubtypes[envSubtype] += 1
            if title not in self.titles:
                self.titles[title] = 1
            else:
                self.titles[title] += 1
            if pmid not in self.pmids:
                self.pmids[pmid] = 1
            else:
                self.pmids[pmid] += 1


        if len(self.envSupertypes) == 1:
            self.consensusEnvSupertype = list(self.envSupertypes.keys())[0]
        else:
            mainEnv, secondEnv = sorted(self.envSupertypes.keys(), key = lambda x: self.envSupertypes[x], reverse = True)[:2]
            if self.envSupertypes[mainEnv] / self.envSupertypes[secondEnv] > 1.5:
                self.consensusEnvSupertype = mainEnv
            else:
                self.consensusEnvSupertype = 'Undetermined'

        if len(self.envTypes) == 1:
            self.consensusEnvType = list(self.envTypes.keys())[0]
        else:
            mainEnv, secondEnv = sorted(self.envTypes.keys(), key = lambda x: self.envTypes[x], reverse = True)[:2]
            if self.envTypes[mainEnv] / self.envTypes[secondEnv] > 1.5:
                self.consensusEnvType = mainEnv
            else:
                self.consensusEnvType = 'Undetermined'

        if len(self.envSubtypes) == 1:
            self.consensusEnvSubtype = list(self.envSubtypes.keys())[0]
        else:
            mainEnv, secondEnv = sorted(self.envSubtypes.keys(), key = lambda x: self.envSubtypes[x], reverse = True)[:2]
            if self.envSubtypes[mainEnv] / self.envSubtypes[secondEnv] > 1.5:
                self.consensusEnvSubtype = mainEnv
            else:
                self.consensusEnvSubtype = 'Undetermined'
            

    def getComponents(self):
        components = set()
        for parentNode in self.parentNodes:
            components = components.union(parentNode.components)
        components = components if components else {self} #Return set(self) if the node has no children.
        return components
    

    def isSubset(self):
        """Check if this node is a subset of another node (all its component taxa are also included in another node."""
        for node in self.parentGraph.nodes.values():
            if node.name == self.name:
                continue
            else:
                if self.components.issubset(node.components):
                    return True  
