"""
Add pathway information to the nodes and edges of a graph.
- Requires a table with fraction of genomes from each genera that have each pathway according to PathwayTools.
- This table will be binarized with a user-supplied presence/absence cuttof.
- Node pathways will be the union of the sets pathways of their constituent genera.
- Edge pathways will be the difference between the set of pathways of the parent node,
  and the union of the pathways of the rest of the parents of the child node.
  I.e. the pathways that the parent node (and only it) is bringing in to the community.
"""

import sys
import argparse
import networkx as nx

from metaCyc_utils import loadPathways, summarizePathways, addBaseCatName, addPathways




def main():

    parser = argparse.ArgumentParser(description='Add pathway information to the nodes and edges of a graph.')
    parser.add_argument('-i', '--input', type = str, required = True,
                        help = 'Input graph file in gpickle format')
    parser.add_argument('-p', '--pathway_matrix', type = str, default = '/home/fer/Software/opt/core-linkage/util/pathwaytable_genus.2.0.formatted.tsv',
                        help = 'Tab-separated matrix containing the pathway presence fraction for each genera')
    parser.add_argument('-n', '--pathway_list', type = str, default = '/home/fer/Documentos/Projects/corNetworks/metNodes/merged6.2_goodPathways/goodMetaCyc/metaCyc.names.tsv',
                        help = 'Tab-separated file with mappings from MetaCyc pathway IDs to pathway names')
    parser.add_argument('-o', '--pathway_onthology', type = str, default = '/home/fer/Software/opt/core-linkage/util/metaCyc.parsed.tsv',
                        help = 'File with onthology information for MetaCyc pathways')
    parser.add_argument('-m', '--min_support', type = int, default = 10,
                        help = 'Only consider nodes with more than min_support when assigning pathways to edges.')
    parser.add_argument('-t', '--threshold', type = float, default = 0.25,
                        help = 'Presence/absence threshold for binarizing the pathway matrix')

    args = parser.parse_args()
    
    G = nx.read_gpickle(args.input)

    ###Load MetaCyc info.
    pathwayTable, pathwayNames, onthology = loadPathways(args.pathway_matrix, args.pathway_list, args.pathway_onthology)
    ###Add pathways to the Graph.
    addPathways(G, pathwayTable, pathwayNames, onthology, args.min_support)
    ###Write results.
    G.graph['pathways'] = 'threshold = {:.2f}, min_support={}'.format(args.threshold, args.min_support)
    outname = '.'.join(args.input.split('.')[:-1]) + '.pathways'
    write(G, '{}.xml'.format(outname))
    nx.write_gpickle(G, '{}.gpickle'.format(outname))

                                           

            

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
