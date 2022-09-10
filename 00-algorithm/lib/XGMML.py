"""
Modified from https://gist.github.com/samadlotia/8c987acb4de07ed86fba.
Copyright Samad Lotia, 2014.
Edge and graph attribute processing added by Fernando Puente-Sanchez, 2016.
"""

from xml.etree.ElementTree import ElementTree
import lxml.etree
from lxml.builder import E
import codecs

def _serialize_attrs(d):
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
