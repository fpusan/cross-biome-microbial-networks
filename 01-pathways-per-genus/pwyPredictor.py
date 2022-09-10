from copy import copy
from pandas import DataFrame
import networkx as nx
from networkx.algorithms import bipartite

MANUAL_KEY_RXNS = {
##                   'PWY-3781': set(['CYTOCHROME-C-OXIDASE-RXN'])
                  }

class Model:

    
    pgdbPaths = dict()
    pathways = dict()
    classes = dict()
    taxonomy = dict()
    ontology = dict()
    reactions = dict()


    @classmethod
    def parse_pgdbPaths(cls, pgdbFile):
        pgdbPaths = {}
        with open(pgdbFile) as infile:
            for line in infile:
                if 'MetaCyc' in line or 'Multi-Organism-Groupings' in line:
                    continue
                pgdb, tax, path = line.strip().split('\t')
                if not tax.startswith('TAX'):
                    continue
                if tax not in pgdbPaths:
                    pgdbPaths[tax] = {'name': pgdb, 'paths': [path]}
                else:
                    pgdbPaths[tax]['paths'].append(path) #Several pgdbs from the same strain.
                        
        cls.pgdbPaths = pgdbPaths
        

    @classmethod
    def parse_pathways(cls, datfile):
        pathways = cls.parse_dat(datfile)
        for pwy, data in pathways.items():
            data['REACTION-LIST'] = set(data['REACTION-LIST'])
            graph = nx.DiGraph()
            for rxn in data['REACTION-LAYOUT']:
                rxn, leftPrim, direction, rightPrim = [x.strip('() ') for x in rxn.split('(:')]
                assert 'LEFT' in leftPrim
                assert 'RIGHT' in rightPrim
                leftPrim = leftPrim.replace('LEFT-PRIMARIES', '').strip().split(' ')
                direction = direction.replace('DIRECTION', '').strip(' :')
                rightPrim = rightPrim.replace('RIGHT-PRIMARIES', '').strip().split(' ')
                if direction == 'NIL': #Happens in a few cases, for reactions that are not in MetaCyc.
                    continue
                if direction == 'L2R':
                    for c in leftPrim:
                        graph.add_edge(c, rxn)
                        graph.nodes[c]['reaction'] = 0
                    for c in rightPrim:
                        graph.add_edge(rxn, c)
                        graph.nodes[c]['reaction'] = 0
                else: #direction == 'R2L':
                    for c in rightPrim:
                        graph.add_edge(c, rxn)
                        graph.nodes[c]['reaction'] = 0
                    for c in leftPrim:
                        graph.add_edge(rxn, c)
                        graph.nodes[c]['reaction'] = 0

                graph.nodes[rxn]['reaction'] = 1
            data['REACTION-GRAPH'] = graph

        #Get pathways that are proper subsets of other pathways.
        for pwy, data in pathways.items():
            data['SUBSET-OF'] = []
            data['SUPERSET-OF'] = []
            for pwy2, data2 in pathways.items():
                if data['REACTION-LIST'] == data2['REACTION-LIST']:
                    continue
                if data['REACTION-LIST'].issubset(data2['REACTION-LIST']): #We don't care if the pathway is actually predicted or not, following https://academic.oup.com/bib/article/11/1/40/193600
                    data['SUBSET-OF'].append(pwy2)                      # 19-Oct-2020: pretty sure the ref I meant to include above was http://standardsingenomics.org/content/5/3/424/ instead
                if data['REACTION-LIST'].issuperset(data2['REACTION-LIST']):
                    data['SUPERSET-OF'].append(pwy2)

        cls.pathways = pathways


    @classmethod
    def parse_classes(cls, datfile):
        cls.classes = cls.parse_dat(datfile)


    @classmethod
    def parse_reactions(cls, datfile):
        reactions = cls.parse_dat(datfile)
        for rxn in reactions:
            reactions[rxn]['IN-PATHWAY'] = []
            reactions[rxn]['IN-PATHWAY-TYPES'] = set()

        for pwy, data in cls.pathways.items():  
            for rxn in data['REACTION-LIST']:
                if rxn not in reactions:
                    if rxn in cls.pathways:
                        reactions[rxn] = {'IN-PATHWAY': [pwy], 'IN-PATHWAY-TYPES': set(data['TYPES']), 'IS-PATHWAY': True}
                    else:
                        reactions[rxn] = {'IN-PATHWAY': [pwy], 'IN-PATHWAY-TYPES': set(data['TYPES']), 'NO-METACYC': True}
                else:
                    reactions[rxn]['IN-PATHWAY'].append(pwy)
                    reactions[rxn]['IN-PATHWAY-TYPES'] = reactions[rxn]['IN-PATHWAY-TYPES'] | set(data['TYPES'])
        cls.reactions = reactions


    @classmethod
    def parse_taxonomy(cls):

        def get_lineage(uid, lineage):
            parent = taxonomy[uid]['parent']
            if parent != 'Organisms':
                lineage.append(parent)
                get_lineage(parent, lineage)
                
            
        taxonomy = {}
        for uid, data in cls.classes.items():
            if uid.startswith('TAX'):
                name = data['COMMON-NAME'][0]
                synonyms = data['SYNONYMS'] if 'SYNONYMS' in data else []
                parent = data['TYPES'][0]
                taxonomy[uid] = {'name': name, 'synonyms': synonyms, 'parent': parent, 'lineage': [uid]}
        for uid, data in taxonomy.items():
            get_lineage(uid, data['lineage'])

        cls.taxonomy = taxonomy


    @classmethod
    def parse_ontology(cls):

        def get_lineage(pwy, lineage):
            lineage.append(pwy)
            if pwy in cls.pathways:
                parents = cls.pathways[pwy]['TYPES']
            else:
                parents = cls.classes[pwy]['TYPES']
            for parent in parents:
                if parent != 'FRAMES':
                    get_lineage(parent, copy(lineage))
                else:
                    ontology[lineage[0]].append(lineage)

        ontology = {}            
        for pwy in cls.pathways:
            ontology[pwy] = []
            get_lineage(pwy, [])

        metapathways = {}
        for pwy, data in cls.pathways.items():
            for parent in data['TYPES']:
                if parent not in metapathways:
                    name = cls.classes[parent]['COMMON-NAME'][0] if 'COMMON-NAME' in cls.classes[parent] else parent
                    metapathways[parent] = {'COMMON-NAME': name, 'PATHWAYS': [pwy]}
                else:
                    metapathways[parent]['PATHWAYS'].append(pwy)

        cls.ontology = ontology
        cls.metapathways = metapathways


    @staticmethod
    def parse_dat(datfile):
        results = {}
        with open(datfile, encoding='iso-8859-1') as infile:
            entries = ''.join([line for line in infile if not line.startswith('#')]) #Remove comments.
            entries = entries.strip().strip('/').split('\n//')
            for entry in entries:
                data = {}
                for line in entry.strip().split('\n'):
                    if not line.startswith('/'):
                        line = line.split(' - ', 1)
                        if len(line) == 2:
                            field, value = line
                            if field not in 'COMMENT':
                                if field not in data:
                                    data[field] = [value]
                                else:
                                    data[field].append(value)
                uid = data['UNIQUE-ID'][0]
                results[uid] = data
        return results
   

    def __init__(self, taxa):
        print(taxa, self.taxonomy[taxa]['name'])
        self.taxa = taxa
        self.enzymes = dict()
        self.enzrxns = {rxn: [] for rxn in self.reactions}
        self.inTaxRange = {}
        self.pwyPrediction = dict()
        self.metapwyPrediction = dict()
        for data in self.reactions.values():
            data['ENZYMES'] = []
        self.taxonomic_filter(self.taxa)
        self.get_pgdb_reactions()
        while True:
            pwyList = [pwy for pwy in self.pathways if self.is_predictable(pwy)]
            if not pwyList:
                break
            for pwy in pwyList:
                self.predict(pwy)
        self.refine_biosynthetic_pwys()
        self.filter_superset_pwys()
        self.predict_metapathways()
        print()


    def taxonomic_filter(self, taxa):
        
        validUIDs = {t for uid in self.taxonomy if taxa in self.taxonomy[uid]['lineage'] for t in self.taxonomy[uid]['lineage']}

        for rxn in self.reactions:
            self.inTaxRange[rxn] = False
        
        for pwy, data in self.pathways.items():
            if 'TAXONOMIC-RANGE' in data:
                pwyTaxa = set(data['TAXONOMIC-RANGE'])
                inRange = True if pwyTaxa & validUIDs else False
            else:
                inRange = True #If it has no declared explicit range assume it can be everywhere.
            self.inTaxRange[pwy] = inRange

        ### WE COMMENT THE LINES BELOW BECAUSE SOME PATHWAYS ARE REACTIONS OF OTHER PATHWAYS
        ### THE LINES BELOW WOULD MAKE A PATHWAY OUTSIDE THE TAXONOMIC RANGE BE OK IF IT ALSO COUNTS AS A REACTION OF A PATHWAY INSIDE THE TAXONOMIC RANGE
        ### WE PREFER TO TRUST THE METACYC GUYS, IF THEY SAY IT WAS OUTSIDE THEN IT IS
##        while True:
##            added = 0
##            for pwy, data in self.pathways.items(): # some reactions might not get properly assigned if they belonged to a pathway that was marked as positive by being itself a reaction of another positive pathway
##                if self.inTaxRange[pwy]:            # take care of them in this final pass
##                    for rxn in data['REACTION-LIST']:
##                        if not self.inTaxRange[rxn]:
##                            self.inTaxRange[rxn] = True
##                            added += 1
##            if not added:
##                break

        self.pgdbPaths =  {tax: data for tax, data in self.pgdbPaths.items() if taxa in self.taxonomy[tax]['lineage']}


    def get_pgdb_reactions(self):
        self.enzymes = {}
        paths = [path for data in self.pgdbPaths.values() for path in data['paths']]
        for path in paths:
            print('\t{}'.format(path))
            with open('{}/enzrxns.dat'.format(path), encoding='iso-8859-1') as infile:
                rxns = ''.join([line for line in infile if not line.startswith('#')]) #Remove comments.
                rxns = [rxn for rxn in rxns.strip().strip('/').split('//\n') if rxn] #Sometimes there are empty reactions.
                for entry in rxns:
                    data = {}
                    for line in entry.strip().split('\n'):
                        if not line.startswith('/'):
                            field, value = line.split(' - ', 1)
                            data[field] = value
                    if 'REACTION' in data and 'ENZYME' in data: #Sometimes the whole entry hasn't an assigned reaction or enzyme, bc of reasons.
                        rxn = data['REACTION']
                        enz = data['ENZYME']
                        data = {x: data[x] for x in data if x in ('BASIS-FOR-ASSIGNMENT', 'COMMON-NAME', 'ENZYME')}
                        if rxn in self.reactions: #Exclude reactions not assigned to any pathway.
                            self.enzrxns[rxn].append(enz) #Use enzrxns so as not to write in the class attributes.
                        if enz not in self.enzymes:
                            self.enzymes[enz] = data
                            self.enzymes[enz]['REACTION-LIST'] = [rxn]
                        else:
                            self.enzymes[enz]['REACTION-LIST'].append(rxn)


    def is_predictable(self, pwy):
        rxns = self.pathways[pwy]['REACTION-LIST']
        if pwy in self.pwyPrediction:
            return False
        for rxn in rxns:
            if rxn in self.pathways and rxn not in self.pwyPrediction:
                return False
        else:
            return True


    def predict(self, pwy):
        data = self.pathways[pwy]
        totalReactions = set(data['REACTION-LIST'])
        totalSpontaneous = set()
        totalEnzymatic = set()
        totalDiagnostic = set()
        totalKey = set(data['KEY-REACTIONS']) if 'KEY-REACTIONS' in data else set()
        totalKey = totalKey & totalReactions #Sometimes the key reaction is from a predecessor pathway because whatever.
        manualKey = MANUAL_KEY_RXNS[pwy] if pwy in MANUAL_KEY_RXNS else set()
        totalKey = totalKey | manualKey
        presentEnzymatic = set()
        presentEC = set()
        presentUnivocal = set()
        presentAmbiguous  = set()

        
        for rxn in totalReactions:
            if rxn in self.pathways:
                if self.pwyPrediction[rxn]['predicted']:
                    presentEC.add(rxn) #######
                    presentEnzymatic.add(rxn)
                    presentUnivocal.add(rxn)
            else:
                if 'SPONTANEOUS?' in self.reactions[rxn]:
                    totalSpontaneous.add(rxn)
                else:
                    if len(self.reactions[rxn]['IN-PATHWAY']) == 1 and 'EC-NUMBER' in self.reactions[rxn]:
                        totalDiagnostic.add(rxn)
                    for enz in self.enzrxns[rxn]:
                        presentEnzymatic.add(rxn)
                        if 'BASIS-FOR-ASSIGNMENT' in self.enzymes[enz] and self.enzymes[enz]['BASIS-FOR-ASSIGNMENT'] == ':EC-NUMBER': #In a few cases they don't have a 'basis for assignment' listed.
                            presentEC.add(rxn)
                        if len(self.enzymes[enz]['REACTION-LIST']) == 1:
                            presentUnivocal.add(rxn)


        totalEnzymatic = totalReactions - totalSpontaneous
        presentAmbiguous = presentEnzymatic - presentUnivocal - presentEC
        
        result = {'name': self.pathways[pwy]['COMMON-NAME'][0],
                  'totalReactions': totalReactions, 'totalSpontaneous': totalSpontaneous, 'totalEnzymatic': totalEnzymatic, 'totalKey': totalKey, 'totalDiagnostic': totalDiagnostic,
                  'presentEnzymatic': presentEnzymatic, 'presentEC': presentEC, 'presentUnivocal': presentUnivocal, 'presentAmbiguous': presentAmbiguous, 'predicted': False, 'basis': 'Not rejected'}

        presentEnzymaticUnambiguous = presentEnzymatic - presentAmbiguous
        

        ###Predict as positive unless explicit reason to reject.
        goodPwy = True
                
        if 'Electron-Transfer' in data['TYPES']:   #Is it an electron transport pathway?
            if presentEnzymatic < totalEnzymatic:
                goodPwy = False
                result['basis'] = 'Incomplete electron-transfer'

        elif not self.inTaxRange[pwy]:             #Is it not expected in this taxonomic range?
            if presentEnzymatic < totalEnzymatic:
                goodPwy = False
                result['basis'] = 'Incomplete outside of tax range'

        elif totalKey:                             #Does it have key reactions?
            if not totalKey.issubset(presentEnzymatic):
                goodPwy = False
                result['basis'] = 'Missing key reactions'
            else:
                result['basis'] = 'All key reactions present'
##
##        elif any(['Energy-Metabolism' in o for o in self.ontology[pwy]]) and not totalDiagnostic:
##            if len(presentEnzymatic) / len(totalEnzymatic) < 0.5:
##                goodPwy = False
##                result['basis'] = 'Energy-Metabolism, less than half-complete'
##            else:
##                result['basis'] = 'Energy-Metabolism, half-complete or better'

        else:
            if totalDiagnostic:                  #Does it have diagnostic reactions?
                if not (totalDiagnostic & presentEnzymatic): #At least one diagnostic (unique) reaction.
                    goodPwy = False
                    result['basis'] = 'No diagnostic reactions present'
            
            else:                                #If it doesn't, accept if complete, or if it has more than 1 reactions and is missing at most 1.
                if presentEnzymatic != totalEnzymatic:
                    if len(presentEnzymatic) <= 1 or len(totalEnzymatic - presentEnzymatic) > 1:
                        goodPwy = False
                        result['basis'] = 'Not nearly-complete'

        if totalSpontaneous == totalReactions:   #In any case accept if completely spontaneous.
            goodPwy = True
            result['basis'] = 'Spontaneous'


        result['predicted'] = goodPwy

        self.pwyPrediction[pwy] = result


    def refine_biosynthetic_pwys(self):
        """
        A biosynthetic pathway that is a proper subset of some other pathway and is missing its final two reactions is not inferred.
        Similarly, a degradative pathway that is a proper subset of some other pathway and is missing its initial step is not inferred.
        The intuition behind the preceding rules is that multiple biosynthetic pathways sometimes share common initial steps,
        and it is the final steps that are most indicative of the presence of the pathway (and conversely for catabolic pathways).

        Ignore if a pathway has key reactions.
        """
        
        for pwy, data in self.pathways.items():

            biosynthetic, degradative = False, False
            
            if any(['Biosynthesis' in o for o in self.ontology[pwy]]):
                biosynthetic = True
            elif any(['Degradation' in o for o in self.ontology[pwy]]):
                degradative = True
            else:
                continue
            
            assert not (biosynthetic and degradative)

            if data['SUBSET-OF'] and 'Super-Pathways' not in data['TYPES'] and 'KEY-REACTIONS' not in data: ##Check taxonomic range here?
                reactionGraph = data['REACTION-GRAPH']
                rxns = [n for n in reactionGraph.nodes() if reactionGraph.nodes[n]['reaction']]
                rxns = [n for n in rxns if 'SPONTANEOUS?' not in self.reactions[n]] #Don't consider spontaneous reactions for this.
                projectedReactionGraph = bipartite.projected_graph(reactionGraph, rxns)
                if biosynthetic:
                    lastRxns = [rxn for rxn in projectedReactionGraph.nodes() if not projectedReactionGraph.successors(rxn)]
                    if not lastRxns: #This happens for superpathways (which we are avoiding anyways), but also for pathways with reversible steps in the beginning/end.
                        continue
                    twoLastRxns = set(lastRxns)
                    for rxn in lastRxns:
                        twoLastRxns.update(projectedReactionGraph.predecessors(rxn))
                        
                    if all([len(self.enzrxns[rxn]) > 0 for rxn in twoLastRxns]):
                        self.pwyPrediction[pwy]['predicted'] = True
                        self.pwyPrediction[pwy]['basis'] = 'Biosynthetic proper subset with the two final reactions'
                    else:
                        self.pwyPrediction[pwy]['predicted'] = False
                        self.pwyPrediction[pwy]['basis'] = 'Biosynthetic proper subset without the two final reactions'
                        
                if degradative:
                    firstRxns = [rxn for rxn in projectedReactionGraph.nodes() if not projectedReactionGraph.predecessors(rxn)]
                    if not firstRxns: #This happens for superpathways (which we are avoiding anyways), but also for pathways with reversible steps in the beginning/end.
                        continue
                    if all([len(self.enzrxns[rxn]) > 0 for rxn in firstRxns]):
                        self.pwyPrediction[pwy]['predicted'] = True
                        self.pwyPrediction[pwy]['basis'] = 'Biosynthetic proper subset with the last reaction'
                    else:
                        self.pwyPrediction[pwy]['predicted'] = False
                        self.pwyPrediction[pwy]['basis'] = 'Biosynthetic proper subset without the last reaction'


    def filter_superset_pwys(self):
        """
        If the reactions in a pathway P1 are a superset of the reactions of an already inferred pathway P2, and P1 has more missing reactions than P2, P1 is not inferred.
        
        Ignore if a pathway has key reactions.
        """
        for pwy, data in self.pathways.items():
            if self.pwyPrediction[pwy]['predicted'] and data['SUPERSET-OF'] and 'KEY-REACTIONS' not in data:
                missing1 = self.pwyPrediction[pwy]['totalEnzymatic'] - self.pwyPrediction[pwy]['presentEnzymatic']
                for pwy2 in data['SUPERSET-OF']:
                    if self.pwyPrediction[pwy2]['predicted']:
                        missing2 = self.pwyPrediction[pwy2]['totalEnzymatic'] - self.pwyPrediction[pwy2]['presentEnzymatic']
                        if len(missing1) > len(missing2):
                            self.pwyPrediction[pwy]['predicted'] = False
                            self.pwyPrediction[pwy]['basis'] = 'Superset of better-predicted pathway'


    def predict_metapathways(self):
        for mpwy, data in self.metapathways.items():
            self.metapwyPrediction[mpwy] = {'predicted': False}
            if any([self.pwyPrediction[pwy]['predicted'] for pwy in data['PATHWAYS']]):
                self.metapwyPrediction[mpwy]['predicted'] = True


    def format_pwy_string(self, pwy, lvl):
        baseCats = set()
        for l in self.ontology[pwy]:
            try:
                baseCats.add(l[lvl])
            except IndexError:
                pass
        shortStrings = sorted([''.join([word[:4] for word in cat.replace('-',' ').split(' ')]) for cat in baseCats])
        return '{} - {} - {}'.format('|'.join(shortStrings), pwy, self.pathways[pwy]['COMMON-NAME'][0])
 
    
    def pwy_summary(self, pwy):
        assert pwy in self.pwyPrediction
        data = self.pwyPrediction[pwy]
        print(pwy)
        print('Name:', data['name'])
        print('Types:', self.pathways[pwy]['TYPES'])
        print('Predicted:', data['predicted'])
        print('Basis for prediction:', data['basis'])
        print('Reactions:')
        G = Model.pathways[pwy]['REACTION-GRAPH']
        rxns = [n for n in G.nodes() if G.nodes[n]['reaction']]
        G2 = bipartite.projected_graph(G, rxns)
        if nx.is_directed_acyclic_graph(G2):
            rxns = nx.topological_sort(G2)
        else:
            print('(non-acyclic graph, reactions are not in order)')  
        for rxn in rxns:
            rstring = ''
            if rxn in self.pathways: rstring += 'P'
            if rxn in data['totalSpontaneous']: rstring += 'S'
            if rxn in data['totalKey']: rstring += 'K'
            if rxn in data['totalDiagnostic']: rstring += 'D'
            if rstring: rstring += ' - '
            print('\t{}{}'.format(rstring, rxn))
            if rxn in self.pathways and self.pwyPrediction[rxn]['predicted']:
                print('\t\tPredicted')
            elif rxn in data['totalSpontaneous']:
                print('\t\tSpontaneous')
            else:
                for enz in self.enzrxns[rxn]:
                    enzd = self.enzymes[enz]
                    name = enzd['COMMON-NAME']
                    assignment = 'EC' if enzd['BASIS-FOR-ASSIGNMENT'] == ':EC-NUMBER' else 'name'
                    promiscuity = len(enzd['REACTION-LIST'])
                    promiscuity = 'Univocal' if promiscuity == 1 else '{} reactions'.format(promiscuity)
                    print('\t\t{}\t{}\t{}\t{}'.format(enz, name, assignment, promiscuity))


    def print_good_pwys(self, printAll = True):
        goodPwys = {p:v for p,v in self.pwyPrediction.items() if v['predicted']}
        if not printAll:
            print(len(goodPwys))
        else:
            for niceString in sorted([self.format_pwy_string(pwy, -4) for pwy in goodPwys]):
                print(niceString)


    def print_good_meta(self):
        goodMPwys = {data['COMMON-NAME'] for mpwy, data in self.metapathways.values() if self.metapwyPrediction[mpwy]['predicted']}
        for mpwy in sorted(goodMPwys):
            print(mpwy)
                



###Main2
Model.parse_pgdbPaths('pgdb_paths_SILVANI.tsv')
Model.parse_pathways('data/metacyc/pathways.dat')
Model.parse_classes('classes.parsed.dat')
Model.parse_reactions('data/metacyc/reactions.dat')
Model.parse_taxonomy()
Model.parse_ontology()

###Get all genera.
allGenera = {' '.join(data['name'].split(' ')[:2]) if data['name'].startswith('Candidatus') else data['name'].split(' ')[0] for data in Model.pgdbPaths.values() if data['name'] != 'MetaCyc'}

namesToUIDs = {md['name']: uid for uid, md in Model.taxonomy.items()}

#Detect "Bacillus <bacterium>" with name.split(' <')[0]
goodGenera = {name.split(' <')[0]: uid for name, uid in namesToUIDs.items() if name.split(' <')[0] in allGenera and ('TAX-2157' in Model.taxonomy[uid]['lineage'] or 'TAX-2' in Model.taxonomy[uid]['lineage'])}

goodGenera = {name: uid for name, uid in goodGenera.items() if 'ales' not in name and 'aceae' not in name}

models = {}
predictions = {}
metaPredictions = {}

for genus, g_uid in goodGenera.items():
##    if genus != 'Polaribacter':
##        continue
    print(genus, g_uid)
    pgdbs = [uid for uid, data in Model.pgdbPaths.items() if g_uid in Model.taxonomy[uid]['lineage']]
    models[genus] = [Model(uid) for uid in pgdbs]
    predictions[genus] = {}
    metaPredictions[genus] = {}
    for pwy in Model.pathways:
        prev = sum([mod.pwyPrediction[pwy]['predicted'] for mod in models[genus]])
        predictions[genus][pwy] = prev / len(models[genus])
    for mpwy in Model.metapathways:
        metaPrev = sum([mod.metapwyPrediction[mpwy]['predicted'] for mod in models[genus]])
        metaPredictions[genus][mpwy] = metaPrev / len(models[genus])


#Make a combined table, with some pathways summarized into metapathways.
ontologyToInclude = ['TCA-VARIANTS', 'Nucleotide-Biosynthesis', 'NUCLEO-DEG', 'Cofactor-Biosynthesis', 'COFACTOR-DEGRADATION',
                     'AROMATIC-COMPOUNDS-BIOSYN', 'AROMATIC-COMPOUNDS-DEGRADATION', 'Amino-Acid-Biosynthesis', 'Amino-Acid-Degradation',
                     'Polyamine-Biosynthesis', 'AMINE-DEG']


def getOntologySet(pwy):
    return {o for ont in Model.ontology[pwy] for o in ont}

metaPathwaysToInclude = [mpwy for mpwy in Model.metapathways if set(ontologyToInclude) & getOntologySet(Model.metapathways[mpwy]['PATHWAYS'][0])]
metaPathwaysToInclude = [mpwy for mpwy in metaPathwaysToInclude if not mpwy in ontologyToInclude] #Sometimes the high-rank ontology is also the direct parent of a pathway, we will avoid such cases.

combinedPredictions = {}

for genus, mpwys in metaPredictions.items():
    combinedPredictions[genus] = {}
    for mpwy in mpwys:
        if mpwy in metaPathwaysToInclude:
            combinedPredictions[genus][mpwy] = mpwys[mpwy]
        else:
            for pwy in Model.metapathways[mpwy]['PATHWAYS']:
                if 'Super-Pathways' not in Model.pathways[pwy]['TYPES']: #Avoid superpathways.
                    combinedPredictions[genus][pwy] = predictions[genus][pwy]


###Write output.
OUTNAME = 'oct2020_escherichia'

pwyTable = DataFrame.from_dict(predictions, orient='index').fillna(0) ######Nitrospira is missing?
pwyTable.loc[:,[pwy for pwy in pwyTable.columns if 'Super-Pathways' not in Model.pathways[pwy]['TYPES']]].to_csv(OUTNAME + '.tsv', sep='\t')

metaPwyTable = DataFrame.from_dict(metaPredictions, orient='index').fillna(0)
metaPwyTable.to_csv(OUTNAME + '.meta.tsv', sep='\t')

combinedTable = DataFrame.from_dict(combinedPredictions, orient='index').fillna(0)
combinedTable.to_csv(OUTNAME + '.combined.tsv', sep='\t')



check = ['GLYCOLYSIS', 'OXIDATIVEPENT-PWY', 'NONOXIPENT-PWY', 'TCA',
         'GLYOXYLATE-BYPASS', 'PWYQT-4429', 'PWY66-399', 'PWY0-1312',
         'PWY0-1313', 'PWY-5939', 'PWY-4261', 'GLYCOGENSYNTH-PWY', 'PWY-5675',
         'PWY-381', 'SO4ASSIM-PWY', 'DISSULFRED-PWY', 'PWY-5279', 'PWY-5276',
         'P381-PWY', 'PWY-5507', 'BIOTIN-BIOSYNTHESIS-PWY', 'ENTNER-DOUDOROFF-PWY']

##Model('TAX-313598').pwy_summary('GLYCOLYSIS')
##Model('TAX-313598').print_good_pwys(False)
##Model('TAX-313594').print_good_pwys(False)
#Model('TAX-716928')
