#!/usr/bin/python3.5

"""
Merge several core-linkage graphs in the gpickle format.
Usage: merge-networks.py <G1> <G2>...<GN> <OUTPUT_PREFIX>
"""


import networkx as nx
import sys

MIN_SUPPORT = 0

def main():

    files = sys.argv[1:-1]
    outname = sys.argv[-1]
    mergedGraph = nx.DiGraph()

    for graphFile in files:
        graph = nx.read_gpickle(graphFile)
        for node, data in graph.node.items():
            del data['consensusEnvSupertype']
            del data['consensusEnvType']
            del data['consensusEnvSubtype']
            del data['pieSubtypeNames']
            del data['pieSubtype']
            if node in mergedGraph.node:
                mergedGraph.node[node]['support'] += graph.node[node]['support']
                mergedGraph.node[node]['sources'] = condenseAttributeLists(mergedGraph.node[node]['sources'], graph.node[node]['sources'])
                mergedGraph.node[node]['titles'] = condenseAttributeLists(mergedGraph.node[node]['titles'], graph.node[node]['titles'])
                mergedGraph.node[node]['pmids'] = condenseAttributeLists(mergedGraph.node[node]['pmids'], graph.node[node]['pmids'])
                mergedGraph.node[node]['envSupertypes'] = condenseAttributeLists(mergedGraph.node[node]['envSupertypes'], graph.node[node]['envSupertypes'])
                mergedGraph.node[node]['envTypes'] = condenseAttributeLists(mergedGraph.node[node]['envTypes'], graph.node[node]['envTypes'])
                mergedGraph.node[node]['envSubtypes'] = condenseAttributeLists(mergedGraph.node[node]['envSubtypes'], graph.node[node]['envSubtypes'])                
                mergedGraph.node[node]['samples'].extend(graph.node[node]['samples'])
                mergedGraph.node[node]['numSamples'] += graph.node[node]['numSamples']
                mergedGraph.node[node]['sourceNetwork'] += '|{}'.format(graphFile.split('.gpickle')[0].split('/')[-1])
                mergedGraph.node[node]['numSourceNetworks'] += 1
            else:
                mergedGraph.add_node(node, attr_dict = data)
                mergedGraph.node[node]['sourceNetwork'] = graphFile.split('.gpickle')[0].split('/')[-1]
                mergedGraph.node[node]['numSourceNetworks'] = 1
        for source, target in graph.edges():
            if (source, target) in mergedGraph.edges():
                mergedGraph[source][target]['support'] += graph[source][target]['support']
            else:
                mergedGraph.add_edge(source, target, attr_dict = graph[source][target])

    ###Delete nodes with low support.
    badNodes = [node for node in mergedGraph.nodes() if mergedGraph.node[node]['support'] < MIN_SUPPORT]
    for node in badNodes:
        mergedGraph.remove_node(node)
    ###Delete edges with low support.
    badEdges = [(s,t) for (s,t) in mergedGraph.edges() if mergedGraph[s][t]['support'] < MIN_SUPPORT]
    for s,t in badEdges:
        mergedGraph.remove_edge(s,t)
    ###Re-count the number of outgoing edges.
    for node in mergedGraph.nodes():
        children = mergedGraph.successors(node)
        for child in children:
            mergedGraph[node][child]['parentOutEdges'] = len(children)
        mergedGraph.node[node]['outEdges'] = len(children)


    write(mergedGraph, '{}.xml'.format(outname))
    nx.write_gpickle(mergedGraph, '{}.gpickle'.format(outname))



def condenseAttributeLists(list1, list2):
    """Get two list of strs in the format 'NAME | COUNTS' and add each element"""
    d1 = {x.split(' | ')[0]: int(x.split(' | ')[1]) for x in list1}
    d2 = {x.split(' | ')[0]: int(x.split(' | ')[1]) for x in list2}
    d3 = d1
    for k,v in d2.items():
        if k in d1:
            d3[k] += v
        else:
            d3[k] = v
    return sorted(['%s | %d'%(source, numSamples) for source, numSamples in d3.items()],
                            key = lambda x: int(x.split('| ')[1]), reverse = True)



#######WRITE THE NETWORKX GRAPH INTO AN XGMML FILE 

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
        elif type(v) == list:
            xattr = E.att({'name': k, 'type': 'list'})
            for val in v:
                if type(val) == str:
                    xval = E.att({'name': k, 'value': val, 'type': 'string'})
                elif type(val) == int:
                    xval = E.att({'name': k, 'value': str(val), 'type': 'integer'})
                elif type(val) == float:
                    xval = E.att({'name': k, 'value': str(val), 'type': 'float'})
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
