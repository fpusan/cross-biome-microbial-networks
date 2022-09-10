PATHWAYS='goodMetaCyc/oct2020.combined_noPWY0-1324.tsv' # PWY0-1324 is duplicated with NAN-MANNACs-degradation
echo """
library(ade4)
ta = read.table("\'$PATHWAYS\'", header=T, row.names=1, check.names=F, sep='\t')
d = 1 - as.matrix(dist.binary(ta, method=1)) #Jaccard dist
write.table(d, 'goodMetaCyc/pwycors.tsv', col.names=NA, quote=F, sep='\t')
""" | R --slave

cd goodMetaCyc
python3.5 00-algorithm/util/parse-cpathways.py classes.parsed.dat pathways.parsed.dat
cd ..

python3.5 00-algorithm/util/merge-networks.py 01-output/*gpickle merged6.2
python3.5 00-algorithm/util/unloop.py merged6.2.gpickle
python3.5 00-algorithm/util/avg-fun-phylo-distance-in-node.py -p goodMetaCyc/phylodist_clean.table.tsv -f goodMetaCyc/pwycors.tsv -i merged6.2.unlooped.gpickle
python3.5 00-algorithm/util/add-pathways.py -i merged6.2.unlooped.avgFunPhyl.gpickle --threshold 0.25 -n goodMetaCyc/metaCyc.names.tsv -o goodMetaCyc/metaCyc.parsed.tsv -p $PATHWAYS

