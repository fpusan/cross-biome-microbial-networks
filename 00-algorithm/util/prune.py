"""
Iteratively remove terminal nodes that share all their component taxa with another, bigger node, provided they appear in the same environments.
"""

import networkx as nx
import sys
from copy import copy


def main():
    G = nx.read_gpickle(sys.argv[1])

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
            
    ###Re-count the number of outgoing edges.
    for node in prunedGraph.nodes():
        children = prunedGraph.successors(node)
        for child in children:
            prunedGraph[node][child]['parentOutEdges'] = len(children)
        prunedGraph.node[node]['outEdges'] = len(children)
    
    outName = sys.argv[1].rsplit('.', 1)[0]
    write(prunedGraph, '{}.pruned.xml'.format(outName))
    nx.write_gpickle(prunedGraph, '{}.pruned.gpickle'.format(outName))


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
