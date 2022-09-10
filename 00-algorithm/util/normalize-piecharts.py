"""
Divide the values in the pieSubtype vector by the total number of samples in each environment.
This way, the fraction of the piecharts occupied by each environment will be normalized by
how well that environment is represented in our database.
"""

import sys
import argparse
import networkx as nx
import numpy as np




def main():

    parser = argparse.ArgumentParser(description='Divide the values in the pieSubtype vector \
                                                  by the total number of samples in each environment.')
    parser.add_argument('-i', '--input', type = str, required = True,
                         help = 'Input graph file in gpickle format')
    parser.add_argument('--remove_NAs', action = 'store_true',
                         help = 'Remove samples classified as "NA-NA-NA" from the pie charts')

    args = parser.parse_args()
    
    G = nx.read_gpickle(args.input)
    envs = {}
    for sample in G.graph['sampleEnvs']:
        sample, env = sample.split(' | ')
        if env not in envs:
            envs[env] = set([sample])
        else:
            envs[env].add(sample)

    pieNames = list(G.node.values())[0]['pieSubtypeNames'] #We draw them for one node, since they are always the same.
    total = np.array([len(envs[name]) for name in pieNames])
    print('{} total samples'.format(sum(total)))

    for node, data in G.node.items():
        normPie = np.array(data['pieSubtype'] / total)
        #Cytoscape wants integers for drawing piecharts.
        #Lets normalize so that the minimum non zero value becomes one and round the rest.
        minPie = min([x for x in normPie if x])
        G.node[node]['pieSubtype'] = [*map(lambda x: int(round(x)), normPie / minPie)]

    #Remove unclassified samples from the pie charts if needed.
    if args.remove_NAs:
        try:
            removeIndex = pieNames.index('NA-NA-NA')
            del(list(G.node.values())[0]['pieSubtypeNames'][removeIndex]) #Needed only one time, since all the nodes contain a reference to the same object.
            for node, data in G.node.items():
                del(data['pieSubtype'][removeIndex])
        except ValueError:
            print('No samples were classified as "NA-NA-NA"')
            args.remove_NAs = False

    outname = '.'.join(args.input.split('.')[:-1]) + '.normPie'
    G.graph['normPie'] = True
    if args.remove_NAs:
        outname = outname + '.noNAenvs'
        G.graph['noNAenvs'] = True
    write(G, '{}.xml'.format(outname))
    nx.write_gpickle(G, '{}.gpickle'.format(outname))
    with open('{}.envTotals.tsv'.format(outname), 'w') as outfile:
        for name, samples in envs.items():
            outfile.write('{}\t{}\n'.format(name, len(samples)))


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
