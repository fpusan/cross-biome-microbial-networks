#!/usr/bin/python3.5

import networkx as nx
import copy
import sys

def main(classes_file, pathways_file):
    graph = nx.DiGraph()
    names = {}
    parseDatFile(classes_file, graph, names)
    parseDatFile(pathways_file, graph, names)
    parentNodes = [node for node in graph.nodes() if not graph.predecessors(node)]
    with open('metaCyc.parsed.tsv', 'w') as o:
        for parent in parentNodes:
            recurse(parent, graph, list(), o)
    with open('metaCyc.names.tsv', 'w') as o:
        for ID, name in names.items():
            o.write('{}\t{}\n'.format(ID, name))


def parseDatFile(inputFile, graph, names):
    links = [item.strip().split('\n') for item in open(inputFile, encoding='cp1252').read().strip().split('\n//\n')]
    for item in links:
        if item[0].startswith('#'):
            continue
        child = [line.split('UNIQUE-ID - ')[1] for line in item if line.startswith('UNIQUE-ID - ')]
        assert len(child) == 1
        child = child[0]
        name = [line.split('COMMON-NAME - ')[1] for line in item if line.startswith('COMMON-NAME - ')]
        if name:
            assert len(name) == 1
            names[child] = name[0]
        parents = [line.split('TYPES - ')[1] for line in item if line.startswith('TYPES - ')]
        for parent in parents:
            graph.add_edge(parent, child)

    
def recurse(node, graph, path, output):
    path.append(node)
    if not graph.successors(node):
        output.write('\t'.join(path) + '\n')
    else:
        for child in graph.successors(node):
            recurse(child, graph, copy.copy(path), output) #If we do not copy it all the children would be writing to the same list.


if __name__ == '__main__':
    main(*sys.argv[1:])
