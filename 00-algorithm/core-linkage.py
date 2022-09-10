#!/usr/bin/python3.5

"""
--------------------------------------------------------------------------------------------------------------------------------------------

Build a community assemblage graph from a presence/absence matrix of different taxa in different samples.

REQUIREMENTS:
- Requires MySQLdb (if using the microdb database to build presence/absence matrices).
- Requires scipy.
- Requires pandas.
- Requires networkx.

USAGE:

OUTPUT:

PARAMETERS:

Distributed under the Modified BSD license.

--------------------------------------------------------------------------------------------------------------------------------------------

"""




__author__ = 'Fernando Puente-Sánchez'
__email__ = 'fpusan@gmail.com'
__version__ = '0.0.1'
__date__ = '11-Jan-2016'
__license__ = 'BSD-3'
__copyright__ = 'Copyright 2016 Fernando Puente-Sánchez'

BSD3_LICENSE = """
    Copyright (c) 2016, Fernando Puente Sánchez
    All rights reserved.
    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright notice, this
      list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.
    * Neither the name of the Centro Nacional de Biotecnologia, nor the names of
      its contributors may be used to endorse or promote products derived from this
      software without specific prior written permission.
    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
    FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

CITATION = """
Unpublished
"""

import argparse
from lib.metaGraph import MetaGraph
TAXON_HIERARCHY = ('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')


def parseArguments():
    """Parse the command line arguments and return an object containing them."""   
    def str2bool(value):
        return value.lower() in ("yes", "true", "t", "1")

    parser = argparse.ArgumentParser(description = 'Build a community assemblage graph from a presence/absence matrix of different taxa in different samples')
    parser.add_argument('-i', '--input', type = str,
                        help = 'Input tab-formatted table. If not provided microdb will be use as a source of data')
    parser.add_argument('-q', '--env_supertype', type = str, nargs='*',
                        help = 'Restrict to samples classified into the selected environmental supertypes (in the microdb hierarchy, requires microdb)')
    parser.add_argument('-w', '--env_type', type = str, nargs='*',
                        help = 'Restrict to samples classified into the selected environmental type (in the microdb hierarchy, requires microdb)')
    parser.add_argument('-e', '--env_subtype', type = str, nargs='*',
                        help = 'Restrict to samples classified into the selected environmental subtype (in the microdb hierarchy, requires microdb)')
    parser.add_argument('-p', '--env_partition', type = str, default = 'subtype',
                        choices = ('supertype', 'type', 'subtype'),
                        help = 'Calculate probability matrices independently for each environment at the provided environmental level')
    parser.add_argument('-l', '--tax_level', type = str, default = 'genus',
                        choices = TAXON_HIERARCHY,
                        help = 'Build presence/absence matrices using the selected taxonomic level (requires microdb)')
    parser.add_argument('-t', '--use_taxon', type = str,
                        help = 'Build co-occurrence tables using only sequences from this taxon. Its taxonomic level must be provided with --use_taxon_tax_level (requires microdb)')
    parser.add_argument('-x', '--use_taxon_tax_level', type = str,
                        choices = TAXON_HIERARCHY,
                        help = 'Taxonomic level for the --use_taxon command. Will be ignored if --use_taxon is not provided (requires microdb)')
    parser.add_argument('-r', '--min_richness', type = int, default = 10,
                        help = 'Avoid samples with less than --min_richness taxa')
    parser.add_argument('-m', '--min_cosmopolitanism', type = int, default = 10,
                        help = 'Avoid taxa in less than --min_cosmopolitanism samples')
    parser.add_argument('-s', '--min_samples', type = int, default = 10,
                        help = 'Avoid environments with less than --min_samples samples')
    parser.add_argument('-u', '--min_ubiquity', type = int, default = 1,
                        help = 'Avoid taxa in less than --min_ubiquity different environments')
    parser.add_argument('-c', '--pval_cutoff', type = float, default = 0.001,
                        help = 'pValue cutoff')
    parser.add_argument('-z', '--size_cutoff', type = int, default = 0,
                        help = 'Node size (samples in which two taxa co-occur) cutoff')
    parser.add_argument('-b', '--bootstrap', type = int, default = 5,
                        help = 'Number of random matrices used for score cutoff calculation')
    parser.add_argument('-g', '--graphs', type = int, default = 10,
                        help = 'Number of random-path graphs used for node support calculation')
    parser.add_argument('-k', '--best_scores', type = int, default = 10,
                        help = 'Randomly pick between the k most significantly aggregated pairs when generating random-path graphs.\
                                If zero, randomly select between all the significantly aggregated pairs, weighted by their aggregation score.')
    parser.add_argument('-o', '--output', type = str, default = 'metagraph',
                        help = 'Output file name.')
    parser.add_argument('-y', '--processors', type=int, default = 1,
                        help = 'Number of processesors to be used')
    parser.add_argument('--cluster_by_size', action = 'store_true',
                        help = 'Cluster significant pairs with a higher co-occurrence size first. By default, scores with a high aggregation score are clustered first')
    parser.add_argument('--profile', action = 'store_true',
                        help = 'Profile graph building with cProfile')
    
    args = parser.parse_args()
    
    ###Additional checks.
    if not args.input and args.use_taxon and not args.use_taxon_tax_level:
        parser.error('--use_taxon was provided. You also need to select a taxonomic level for this selection via --use_taxon_tax_level')
    if not args.input and args.use_taxon_tax_level:
        if TAXON_HIERARCHY.index(args.use_taxon_tax_level) >=  TAXON_HIERARCHY.index(args.tax_level):
            parser.error('The taxonomic level selected for building presence/absence matrices (--tax_level) is higher or equal than the taxonomic level we are restricted to (--use_taxon_tax_level)')
    return args



if __name__ == '__main__':
    #Graph(parseArguments())
    MetaGraph(parseArguments())
    

    
