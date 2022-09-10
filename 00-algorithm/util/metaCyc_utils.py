import numpy as np
from pandas import read_csv, DataFrame
from copy import copy

import networkx as nx

MAX_ONTOLOGY_LEVELS = 7


def loadPathways(pathway_matrix = '/home/fer/Documentos/Projects/corNetworks/metNodes/merged6.2_goodPathways/goodMetaCyc/pathways_oneModelPerStrain.tsv',
                 pathway_list = '/home/fer/Documentos/Projects/corNetworks/metNodes/merged6.2_goodPathways/goodMetaCyc/metaCyc.names.tsv',
                 pathway_ontology = '/home/fer/Documentos/Projects/corNetworks/metNodes/merged6.2_goodPathways/goodMetaCyc/metaCyc.parsed.tsv',
                 threshold = 0.25):

    pathwayNames = {ID: name for ID, name in (line.strip().split('\t') for line in open(pathway_list))}
    pathwayTable = read_csv(pathway_matrix, sep='\t', index_col=0)

    ontology = {}
    for line in open(pathway_ontology):
        line = line.strip().split('\t')
        if len(line) > 4 and line[2] == 'Pathways':
            *ont, pwy = line[3:]
            ont = ont[:MAX_ONTOLOGY_LEVELS]
            for i in range(MAX_ONTOLOGY_LEVELS - len(ont)):
                ont.append(ont[-1]) #So that all the ontology lists have the same size.
            ont = [pathwayNames[x] if x in pathwayNames else x for x in ont] #Some ontology names are not in pathwayNames (bc they don't have a common name).

            if 'Reductants' in ont:
                # Now this is funny. There are two categories involved here, with their following IDs and names:
                # 1) Reductants -> Reductants biosynthesis
                # 2) REDUCTANTS -> Reductants
                # This will mess up the code below. REDUCTANTS actually is a category name for reductants degradation, and has only one pathway (PWYQT-4432).
                # I can't find this pathway in MetaCyc, nor in my data. So I will just ignore it exists.
                continue
            
            if pwy not in ontology: 
                ontology[pwy] = [ont]
            else:
                ontology[pwy].append(ont)


    ###Add ontology for metapathways.
    for pwy in list(ontology.keys()):
        for ont in ontology[pwy]:
            metapwy = ont[-1]
            ont = ont[:-1] #Metapathway ontologies are one level shorter than the rest.
            pathwayNames[metapwy] = metapwy #Ontology IDs were already translated into names before.
            if metapwy not in ontology:
                ontology[metapwy] = [ont]
            elif ont not in ontology[metapwy]:
                ontology[metapwy].append(ont)


    #Change pathway names in ontology dict:
    ontology = {pathwayNames[ID]: onts for (ID, onts) in ontology.items()}

    baseOntology = {o[0] for x in ontology.values() for o in x}

    ###Binarize pathway matrix.
    pathwayTable[pathwayTable > threshold] = 1
    pathwayTable[pathwayTable < 1] = 0

    ###Remove pathways with no ontology.
##    for pwy in pathwayTable.columns:
##        if pwy not in pathwayNames:
##            del pathwayTable[pwy]

    pathwayTable = pathwayTable.loc[:,np.array([pathwayNames[pwy] in ontology for pwy in pathwayTable.columns])]

    return pathwayTable, pathwayNames, ontology


def summarizePathways(pwySet, ontology, lvl, returnDict = False):
    """Summarize a list of pathways into broader categories"""
    pwySummary = {}
    for pwy in pwySet:
        base = {ont[0] for ont in ontology[pwy]}
        for ont in ontology[pwy]:
            cat = ont[lvl-1]
            if ont[0] == 'Superpathways' and len(base)>1: #Pwys belonging to the "Superpathways" category usually also belong to a more informative category, so we try to skip them.
                continue
            else:
                if cat not in pwySummary:
                    pwySummary[cat] = 1
                else:
                    pwySummary[cat] += 1
    if not returnDict:
        return ['{} | {}'.format(pwy, abund) for pwy, abund in sorted(pwySummary.items(), key = lambda x: x[1], reverse = True)]
    else:
        return pwySummary


def addBaseCatName(pwy, ontology, lvl): #Meow
    baseCats = [ont[lvl-1] for ont in ontology[pwy]]
    baseCats = '|'.join([simplifyName(cat) for cat in baseCats])
    return baseCats + ' - ' + pwy


def simplifyName(string):
    exclude = ('and')
    goodWords = [w for w in string.replace('-', ' ').split(' ') if w.lower() not in exclude]
    return ''.join([w[:4].capitalize() for w in goodWords[:3]])


def addPathways(G, pathwayTable, pathwayNames, ontology, min_support, addBaseCatStr = True):
    """
    Recursively assign pathways to graph nodes (in-place).
    - Node pathways will be the union of the sets pathways of their constituent genera.
    - Edge pathways will be the difference between the set of pathways of the parent node,
    and the union of the pathways of the rest of the parents of the child node.
    I.e. the pathways that the parent node (and only it) is bringing in to the community.
    """
    lvl2ontology = {o[1] for x in ontology.values() for o in x}
    #Define the recursion function.
    def addPathwaysToNode(node):
        if not G.predecessors(node):
            if 'pathways' not in G.node[node]:
                if len(G.node[node]['taxa']) != 1: print(node)
                assert len(G.node[node]['taxa']) == 1
                taxa = G.node[node]['taxa'][0]
                if taxa in pathwayTable.index:
                    pwyIDs = pathwayTable.columns[pathwayTable.loc[taxa] > 0]
                    pathways = {pathwayNames[ID] for ID in pwyIDs}
                else:
                    pathways = set()
                G.node[node]['pathways'] = pathways
                G.node[node]['numPathways'] = len(pathways)
                G.node[node]['lvl2pwys'] = summarizePathways(pathways, ontology, 2)
            else:
                pass #We already added the pathways for this node.
        else:
            pathways = set()
            for child in G.predecessors(node):
                if 'pathways' not in G.node[node]:
                    addPathwaysToNode(child)
                pathways = pathways.union(G.node[child]['pathways'])
            G.node[node]['pathways'] = pathways
            G.node[node]['numPathways'] = len(pathways)
            G.node[node]['lvl2pwys'] = summarizePathways(pathways, ontology, 2)

    #Run it.
    tipNodes = [node for node in G.nodes() if not G.successors(node)]
    for node in tipNodes:
        addPathwaysToNode(node)

    ###Assign pathways to graph edges.
    for node in G.nodes():
        parents = G.predecessors(node)
        havePathways = True
        if not parents:
            continue
        for parent in parents:
            if G.node[parent]['support'] < min_support: #Ignore parents with low support.
                G[parent][node]['pathways'] = set()
                G[parent][node]['numPathways'] = -2
                G[parent][node]['percentPathways'] = -2
                G[parent][node]['lvl2pwys'] = set()
                G[parent][node]['lvl3pwys'] = set()
                for cat in lvl2ontology:
                    G[parent][node][cat] = -2
                G[parent][node]['mainLvl1'] = 'no_data'
                G[parent][node]['mainLvl2'] = 'no_data'
                G[parent][node]['mainLvl3'] = 'no_data'

            else:
                exclusivePwys = copy(G.node[parent]['pathways'])
                if not exclusivePwys:
                    havePathways = False #Do not write info if any parent has no available genome.
                for partner in parents:
                    #Substract the pathways that appear in the rest of the parent nodes.
                    if partner == parent or G.node[parent]['support'] < min_support: #Ignore parents with low support.
                        continue
                    partnerPwys = copy(G.node[partner]['pathways'])
                    if not partnerPwys:
                        havePathways = False #Do not write info if any parent has no available genome.
                    exclusivePwys -= partnerPwys
                #Add edge information.
                if not havePathways:
                    G[parent][node]['pathways'] = set()
                    G[parent][node]['numPathways'] = -1
                    G[parent][node]['percentPathways'] = -1
                    G[parent][node]['lvl2pwys'] = set()
                    G[parent][node]['lvl3pwys'] = set()
                    for cat in lvl2ontology:
                        G[parent][node][cat] = -1
                    G[parent][node]['mainLvl1'] = 'no_data'
                    G[parent][node]['mainLvl2'] = 'no_data'
                    G[parent][node]['mainLvl3'] = 'no_data'

                else:
                    G[parent][node]['pathways'] = exclusivePwys
                    G[parent][node]['numPathways'] = len(exclusivePwys)
                    G[parent][node]['percentPathways'] = round(100*len(exclusivePwys)/len(G.node[parent]['pathways']))
                    G[parent][node]['lvl2pwys'] = summarizePathways(exclusivePwys, ontology, 2)
                    G[parent][node]['lvl3pwys'] = summarizePathways(exclusivePwys, ontology, 3)
                    lvl1 = summarizePathways(exclusivePwys, ontology, 1, returnDict = True)
                    lvl2 = summarizePathways(exclusivePwys, ontology, 2, returnDict = True)
                    lvl3 = summarizePathways(exclusivePwys, ontology, 3, returnDict = True)
                    for cat in lvl2ontology:
                        if cat in lvl2:
                            G[parent][node][cat] = lvl2[cat]
                        else:
                            G[parent][node][cat] = 0
                    if exclusivePwys:
                        try:
                            mainLvl1 = sorted(lvl1, key = lambda x: lvl1[x], reverse = True)[0]
                        except:
                            print(exclusivePwys)
                            print(lvl1)
                            raise
                        mainLvl2 = sorted(lvl2, key = lambda x: lvl2[x], reverse = True)[0]
                        mainLvl3 = sorted(lvl3, key = lambda x: lvl3[x], reverse = True)[0]
                    else:
                        mainLvl1 = 'no_pwys_added'
                        mainLvl2 = 'no_pwys_added'
                        mainLvl3 = 'no_pwys_added'
                    G[parent][node]['mainLvl1'] = mainLvl1
                    G[parent][node]['mainLvl2'] = mainLvl2
                    G[parent][node]['mainLvl3'] = mainLvl3

    ###If required, add the base functional category to the pathway names.
    if addBaseCatStr:
        for node in G.nodes():
            if G.node[node]['pathways']:
                G.node[node]['pathways'] = list(sorted([addBaseCatName(pwy, ontology, 2) for pwy in G.node[node]['pathways']]))
        for parent, child in G.edges():
            if G[parent][child]['pathways']:
                G[parent][child]['pathways'] = list(sorted([addBaseCatName(pwy, ontology, 2) for pwy in G[parent][child]['pathways']]))


def addEnvToEdges(G):
    """Add source network information to edges in place."""
    for parent, child in G.edges():
        G[parent][child]['sourceNetwork'] = G.node[parent]['sourceNetwork']
        G[parent][child]['numSourceNetworks'] = G.node[parent]['numSourceNetworks']


def countComplementingPathwayPairs(G, min_support):
    """
    Where G has already been annotated with addPathways.
    A pair of pathways will be considered to be complementing if they appear in different edges leading to the same node.
    """
    complementingPathwayPairs = {}
    for node in G.nodes():
        parents = G.predecessors(node)
        if len(parents) == 2 and G[parents[0]][node]['support'] >= min_support and G[parents[1]][node]['support'] >= min_support: #Avoid cases with several paths leading to the same node, to keep it simple.
            parent1pwys = G[parents[0]][node]['pathways']
            parent2pwys = G[parents[1]][node]['pathways']
            for p1 in parent1pwys:
                for p2 in parent2pwys:
                    mergedName = '*'.join(sorted([p1,p2]))
                    if mergedName in complementingPathwayPairs:
                        complementingPathwayPairs[mergedName] += 1
                    else:
                        complementingPathwayPairs[mergedName] = 1
    return complementingPathwayPairs


def addFunAnnotTaxa(G, pathwayTable):
    """Add the number of functionally annotated taxa to the nodes in place."""
    for node in G.nodes():
        taxa = G.node[node]['taxa']
        funAnnotTaxa = 0
        for t in taxa:
            if t in pathwayTable.index:
                funAnnotTaxa += 1
        G.node[node]['funAnnotTaxa'] = funAnnotTaxa


def addAvgPwysPerTaxa(G, pathwayTable):
    """Add the average number of pathways in the member taxa to the nodes in place."""
    for node in G.nodes():
        taxa = G.node[node]['taxa']
        numPathways = []
        for t in taxa:
            if t in pathwayTable.index:
                numPathways.append(sum(pathwayTable.loc[t,:]))
        if numPathways:
            G.node[node]['avgPwysPerTaxa'] = sum(numPathways)/len(numPathways)
        else:
            G.node[node]['avgPwysPerTaxa'] = 0
                
