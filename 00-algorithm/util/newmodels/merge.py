from pandas import read_csv, DataFrame, concat
from os import listdir

df = read_csv('../pathwaytable_genus.formatted.tsv', sep='\t')
print(df.shape)

pathways = {}

for f in listdir('.'):
    if f.endswith('.pathways'):
        name = f.replace('.pathways','')
        genus = name.split('_')[0]
        if genus not in pathways:
            pathways[genus] = {'numGenomes': 1, 'pwys': {}}
        else:
            pathways[genus]['numGenomes'] += 1
        pwys = [line.split('\t')[0] for line in open(f) if not line.startswith(';') and not line.startswith('FRAME')]
        for pwy in pwys:
            if pwy not in pathways[genus]['pwys']:
                pathways[genus]['pwys'][pwy] = 1
            else:
                pathways[genus]['pwys'][pwy] += 1
        

for genus in pathways:
    print(len(pathways[genus]['pwys']))
    pathways[genus]['pwys'] = {pwy: count/pathways[genus]['numGenomes'] for pwy, count in pathways[genus]['pwys'].items()}
    df = concat([df, DataFrame.from_dict({genus: pathways[genus]['pwys']}, orient = 'index')])

print(df.shape)

df.sort_index(inplace = True)
df.sort_index(1, inplace = True)
df.fillna(0).to_csv('../pathwaytable_genus.2.0.formatted.tsv', sep='\t')
