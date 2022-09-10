#!/usr/bin/python3.5

"""
Generate a matrix with distribution of pathways from a given cathegory across genera.
"""

from metaCyc_utils import loadPathways, simplifyName

ONTHOLOGY_LEVEL = 2
OUTPUT_FOLDER = 'pwysPerCategory/'

pathwayTable, pathwayNames, onthology = loadPathways()
pathwayCodes = {v:k for k,v in pathwayNames.items()}

categories = {o[ONTHOLOGY_LEVEL-1] for ontList in onthology.values() for o in ontList} #The same pathway can have more than one onthology, values in the dict are lists of lists.

for cat in categories: #meow
    pwys = [pwy for pwy in onthology if any([o[ONTHOLOGY_LEVEL-1]==cat for o in onthology[pwy]])]
    pwys = [pathwayCodes[pwy] for pwy in pwys if pathwayCodes[pwy] in pathwayTable.columns]
    #outname = simplifyName(cat) + '.tsv'
    outname = cat.replace('/','-') + '.tsv'
    pathwayTable.loc[:,pwys].to_csv(OUTPUT_FOLDER+outname, sep='\t')
    
