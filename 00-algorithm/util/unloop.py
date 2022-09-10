#!/usr/bin/python3.5

"""
Remove loops from the graph.
When the same node can be reached by more than one path, pick the one with the most overall support, as long as it does not contain bifurcating nodes, or nodes coming from a different source network than the tip node.
"""

import networkx as nx
import sys
from copy import copy
from itertools import combinations, chain


def main():
    G = nx.read_gpickle(sys.argv[1])
    ###First pruning.
    G = prune(G)
    ###Unlooping.
    tipNodes = [node for node in G.nodes_iter() if not G.successors(node)]
    for node in tipNodes:
        sourceNetwork = G.node[node]['sourceNetwork']
        unloop(G, node, set(node.split('-')), sourceNetwork)
    nodesToRemove = [node for node in G.nodes_iter() if not G.successors(node) and not G.predecessors(node)]
    for node in nodesToRemove:
        G.remove_node(node)
    ###Second pruning.
    G = prune(G)
    ###Re-count the number of outgoing edges.
    for node in G.nodes():
        children = G.successors(node)
        for child in children:
            G[node][child]['parentOutEdges'] = len(children)
        G.node[node]['outEdges'] = len(children)
    ###Write output.
    outName = sys.argv[1].rsplit('.', 1)[0]
    write(G, '{}.unlooped.xml'.format(outName))
    nx.write_gpickle(G, '{}.unlooped.gpickle'.format(outName))


def prune(G):
    """Iteratively remove terminal nodes that share all their component taxa with another, bigger node, provided they appear in the same environments."""

    print('Pruning...')
    prunedGraph = copy(G)

    while True:
        nodeList = [set(node.split('-')) for node in prunedGraph.nodes()]
        termNodes = [node for node in prunedGraph.nodes_iter() if not prunedGraph.successors(node) and isSubset(node, nodeList, prunedGraph)]
        print(len(prunedGraph), len(termNodes))
        if not termNodes:
            break
        else:
            prunedGraph.remove_nodes_from(termNodes)

    return prunedGraph


def isSubset(node, nodeList, prunedGraph):
    nodeSet = set(node.split('-'))
    for partner in nodeList:
        if nodeSet==partner:
            continue
        if not nodeSet - partner:
            partner = '-'.join(sorted(partner))
            if 'sourceNetwork' in prunedGraph.node[node]:
                nodeSourceNets = set(prunedGraph.node[node]['sourceNetwork'].split('|'))
                partnerSourceNets = set(prunedGraph.node[partner]['sourceNetwork'].split('|'))
            else:
                nodeSourceNets, partnerSourceNets = set(), set()
            if not nodeSourceNets - partnerSourceNets:
                return True
    else:
        return False
    

def pairParents(parents, child):
    """Given a list of parents, group the ones that result in the child node when combined"""
    goodCombinations = []
    child = set(child.split('-'))
    for p1, p2 in combinations(parents,2):
        if set(p1.split('-')) ^ set(p2.split('-')) == child and len(p1.split('-')) + len(p2.split('-')) == len(child):
            goodCombinations.append(tuple([p1,p2]))
    try:
        assert len(set((chain(*goodCombinations)))) == len(set(parents))
    except:
        print(child)
        print(parents)
        print(goodCombinations)
        raise
    return goodCombinations


def unloop(G, node, tipNode, tipSourceNetwork):
    """Recursively remove redundant paths reaching a node, as long as they don't contain bifurcating nodes, or nodes coming from a different source network than the tip node."""
    ###Can we remove this node?
    forbidden = False
    if G.node[node]['sourceNetwork'] != tipSourceNetwork and len(G.node[node]['taxa']) > 1:
        forbidden = True #We can't if it comes from a different environment
        #print('{} forbidden due to different env'.format(node))
    children = G.successors(node)
    for child in children:
        if set(child.split('-')) - tipNode and len(G.node[node]['taxa']) > 1:
            forbidden = True #We can't if it leads to nodes with different taxa.
            #print('{} forbidden due to different children'.format(node))
    
    pairedParents = pairParents(G.predecessors(node), node)

    candidates = {}
    ###Now try to remove redundant predecessors.
    for pair in pairedParents:
        #Recursively check predecessors.
        forbidden1 = unloop(G, pair[0], tipNode, tipSourceNetwork)
        forbidden2 = unloop(G, pair[1], tipNode, tipSourceNetwork)
        support1 = G[pair[0]][node]['support']
        support2 = G[pair[1]][node]['support']
        candidates[pair] = {'support': support1 + support2}
        candidates[pair]['forbidden'] = forbidden1 or forbidden2
        #forbidden = forbidden or forbidden1 or forbidden2 #As long as one parent is forbidden to be removed, the child will be too.

    if len(candidates) > 1:
        for cand, values in sorted(candidates.items(), key = lambda x: x[1]['support'], reverse = True)[1:]: #We will keep the pair with the best support and try to remove the others.
            if not values['forbidden']:
                G.remove_edge(cand[0], node)
                G.remove_edge(cand[1], node)
                print(*cand)

    return forbidden
    
    
    
        
        





###XGMMLparser:
##    """
##    Modified from https://gist.github.com/samadlotia/8c987acb4de07ed86fba.
##    Copyright Samad Lotia, 2014.
##    Edge attribute processing added by Fernando Puente-Sanchez, 2016.
##    """
from xml.etree.ElementTree import ElementTree
import lxml.etree
from lxml.builder import E
import codecs

def _serialize_attrs(d):
    """
    Modified from https://gist.github.com/samadlotia/8c987acb4de07ed86fba.
    Samad Lotia, 2014.
    Edge and graph attribute processing added by Fernando Puente-Sanchez, 2016.
    """

    xattrs = list()
    for k, v in d.items():
        if type(v) == str:
            xattr = E.att({'name': k, 'value': v, 'type': 'string'})
        elif type(v) == float:
            xattr = E.att({'name': k, 'value': str(v), 'type': 'real'})
        elif type(v) == int:
            xattr = E.att({'name': k, 'value': str(v), 'type': 'integer'})
        elif type(v) == list or type(v) == set:
            xattr = E.att({'name': k, 'type': 'list'})
            for val in v:
                if type(val) == str:
                    xval = E.att({'name': k, 'value': val, 'type': 'string'})
                elif type(val) == int:
                    xval = E.att({'name': k, 'value': str(val), 'type': 'integer'})
                else:
                    raise Exception('Cannot serialize attribute of type %s' % type(val))
                xattr.append(xval)
        elif type(v) == type(None):
            xattr = E.att({'name': k, 'value': 'None', 'type': 'string'})
        elif type(v) == bool:
            if v:
                xattr = E.att({'name': k, 'value': 'True', 'type': 'string'})
            else:
                xattr = E.att({'name': k, 'value': 'False', 'type': 'string'})
        else:
            raise Exception('Cannot serialize attribute of type %s' % type(v))
        xattrs.append(xattr)
    return xattrs


def write(G, path):
    """
    Modified from https://gist.github.com/samadlotia/8c987acb4de07ed86fba.
    Copyright Samad Lotia, 2014.
    Edge and graph attribute processing added by Fernando Puente-Sanchez, 2016.
    """

    xG = E.graph({'label': 'Network', 'directed': '1'}, namespace='http://www.cs.rpi.edu/XGMML')
    for k, v in G.graph.items():
        xattrs = _serialize_attrs({k:v})
        xG.extend(xattrs)
    nodes = G.nodes()
    for node_id in range(len(nodes)):
        node = nodes[node_id]
        xnode = E.node({'id': str(node_id), 'label': node})
        xattrs = _serialize_attrs(G.node[node])
        xnode.extend(xattrs)
        xG.append(xnode)

    for src in G.edge:
        src_id = nodes.index(src)
        for trg in G.edge[src]:
            trg_id = nodes.index(trg)
            xedge = E.edge({'source': str(src_id), 'target': str(trg_id)})
            xattrs = _serialize_attrs(G.edge[src][trg])
            xedge.extend(xattrs)
            xG.append(xedge)

    t = ElementTree(xG)
    t.write(path, xml_declaration=True, encoding='UTF-8')


#########

main()
