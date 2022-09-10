import argparse
from statistics import mean
from itertools import combinations
import networkx as nx
from pandas import read_csv


def main():

    parser = argparse.ArgumentParser(description='Add the average functional and phylogenetic relationships \
                                                  between node members as node attributes.')
    parser.add_argument('-i', '--input', type = str, required = True,
                         help = 'Input graph file in gpickle format')
    parser.add_argument('-f', '--funcorr_matrix', type = str, default = 'pathwaytable_genus.2.0.cors.tsv',
                         help = 'Tab-separated functional correlation matrix between genera')
    parser.add_argument('-p', '--phyldist_matrix', type = str, default = 'phylodist_clean.table.tsv',
                         help = 'Tab-separated phylogenetic distance matrix between genera')


    args = parser.parse_args()


    G = nx.read_gpickle(args.input)

    funCorrs = read_csv(args.funcorr_matrix, sep = '\t', index_col = 0)
    phyloDists = read_csv(args.phyldist_matrix, sep = '\t', index_col = 0)

    for node, data in G.node.items():
        taxa = data['taxa']
        #Functional correlations
        annotatedTaxa = len([x for x in taxa if x in funCorrs.index])
        if len(taxa) == 1 and taxa[0] in funCorrs.index:
            meanCorr = 1.0
        else:
            nodeCorrs = []
            for tax1, tax2 in combinations(taxa, 2):
                if tax1 in funCorrs.index and tax2 in funCorrs.index:
                    nodeCorrs.append(funCorrs.loc[tax1, tax2])
            if nodeCorrs:
                meanCorr = float(mean(nodeCorrs))
            else:
                meanCorr = 0.0 #No available correlations (only 1 taxa, or genomes not available).
        G.node[node]['meanFunCorr'] = meanCorr
        G.node[node]['funAnnotTaxa'] = annotatedTaxa
        #Phylogenetic correlations
        annotatedTaxa = len([x for x in taxa if x in phyloDists.index])
        if len(taxa) == 1 and taxa[0] in phyloDists.index:
            meanDist = 0.0
        else:
            nodeDists = []
            for tax1, tax2 in combinations(taxa, 2):
                if tax1 in phyloDists.index and tax2 in phyloDists.index:
                    nodeDists.append(phyloDists.loc[tax1, tax2])
            if nodeDists:
                meanDist = float(mean(nodeDists))
            else:
                meanDist = 1000.0 #No available correlations (only 1 taxa, or genomes not available).
        G.node[node]['meanPhyloDist'] = meanDist
        G.node[node]['phyloAnnotTaxa'] = annotatedTaxa

        
    G.graph['avgFunPhyl'] = True
    outname = '.'.join(args.input.split('.')[:-1]) + '.avgFunPhyl'
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

