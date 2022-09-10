library(ade4)

###Function definitions.

isSubset = function(node1,node2)
    {
    node1taxa = getNodeTaxa(node1)
    node2taxa = getNodeTaxa(node2)
    return(length(setdiff(node1taxa,node2taxa))==0)
    }

branchLength = function(c1, c2)
    {
    dists = c()
    for(i in c1)
        {
        for(j in c2)
            {
            dists = c(dists, allPhyloDists[i,j])
            }
        }
    return(mean(dists))
    }


plotClustPhylo = function(node)
    {
    nodeTaxa = getNodeTaxa(node)
    plot(hclust(as.dist(allPhyloDists[nodeTaxa,nodeTaxa])), main='Phylogeny', sub='', xlab='')
    }


plotClustFun = function(node)
    {
    nodeTaxa = getNodeTaxa(node)
    plot(hclust(as.dist(1-graphFunCorrs[nodeTaxa,nodeTaxa])), main='Function', sub='', xlab='')
    }




plotClust = function(node)
    {
    par(mfrow=c(1,2))
    plotClustPhylo(node)
    plotClustFun(node)
    leg = unlist(strsplit(node, '-'))
    splitIdx = floor((1:length(leg)-1)/3) #Three names per line.
    leg = split(leg,splitIdx)
    leg = paste(sapply(leg, paste, collapse=' '), collapse='\n')
    mtext(leg, outer=T, line=-4)
    dev.copy(png, paste('largestAnnotatedAssemblies/',node,'.png',sep=''))
    dev.off()
    }


grab_grob <- function()
    {
    library(gridGraphics)
    library(grid)
    grid.echo()
    grid.grab()
    }


plotUpset = function(node)
    {
    nodeTaxa = getNodeTaxa(node)
    l = list()
    for(taxon in nodeTaxa) { pwys = colnames(pathways)[pathways[taxon,]>0]; l[[taxon]] = pwys[!is.na(pwys)] }
    upset(fromList(l), nsets=length(nodeTaxa), nintersects = NA)
    }    

    


plotHeatmaps = function(node, minPhyloDist=0, maxPhyloDist=1, minFunDist=0, maxFunDist=1, save=F)
    {
    library(gplots)
    library(gridExtra)
    try(dev.off(), silent=T)
    nodeTaxa = getNodeTaxa(node)
    gl = list()
    pal = colorRampPalette(c('blue','red'))
    y = allPhyloDists[nodeTaxa,nodeTaxa]
    diag(y) = NA
    par(mar=c(7,4,4,2)+0.1)
    breaks = seq(minPhyloDist, maxPhyloDist, length.out = 10)
    heatmap.2(as.matrix(y), symm=T, na.col='black', trace='none', keysize=1, margins=c(15,15), main='Phylogenetic distance', col=pal, breaks=breaks)
    gl[[1]] = grab_grob()
    y = 1 - graphFunCorrs[nodeTaxa,nodeTaxa]
    diag(y) = NA
    par(mar=c(7,4,4,2)+0.1)
    breaks = seq(minFunDist, maxFunDist, length.out = 10)
    heatmap.2(as.matrix(y), symm=T, na.col='black', trace='none', keysize=1, margins=c(15,15), main='Functional distance', col=pal, breaks=breaks)
    gl[[2]] = grab_grob()
    if(save==T)
        {
        png(paste('largestAnnotatedAssemblies/',node,'.png'), width = 1000, height = 480)
        grid.arrange(grobs=gl, ncol=2, clip=TRUE)
        dev.off()
    }else
        {
        grid.arrange(grobs=gl, ncol=2, clip=TRUE)
        }
    }


testResults = function(index, allowDifferentEnv=F)
    {
    if(allowDifferentEnv==T)
        {
        res = res2
        }else
        {
        res = res
        }
    par(mfrow=c(1,2))
    randomDists = sapply(res[[index]], avgPhyloDist)
    realDist = avgPhyloDist(bigTipNodes[index])
    randomCorrs = sapply(res[[index]], avgFunCorr)
    realCorr = avgFunCorr(bigTipNodes[index])
    if(realCorr < min(randomCorrs)) {ylim = c(realCorr-0.01, max(randomCorrs)+0.01)
    }else if(realCorr > max(randomCorrs)) {ylim = c(min(randomCorrs)-0.01, realCorr+0.01)
    }else {ylim=NULL}
    boxplot(randomDists, ylab='Phylo dist')
    abline(h=realDist, col='red')
    boxplot(randomCorrs, ylab='Fun corr', ylim=ylim)
    abline(h=realCorr, col='red')
    leg = unlist(strsplit(bigTipNodes[index], '-'))
    splitIdx = floor((1:length(leg)-1)/3) #Three names per line.
    leg = split(leg,splitIdx)
    leg = paste(sapply(leg, paste, collapse=' '), collapse='\n')
    mtext(leg, outer=T, line=-4)
    }


avgFunCorr = function(node)
    {
    nodeTaxa = getNodeTaxa(node)
    ta = allFunCorrs[nodeTaxa, nodeTaxa]
    return(mean(ta[upper.tri(ta)]))
    }


avgPhyloDist = function(node)
    {
    nodeTaxa = getNodeTaxa(node)
    ta = allPhyloDists[nodeTaxa, nodeTaxa]
    return(mean(ta[upper.tri(ta)]))
    }

cvFunCorr = function(node)
    {
    nodeTaxa = getNodeTaxa(node)
    ta = graphFunCorrs[nodeTaxa, nodeTaxa]
    return( sd(ta[upper.tri(ta)]) / mean(ta[upper.tri(ta)]) )
    }
    
cvPhyloDist = function(node)
    {
    nodeTaxa = getNodeTaxa(node)
    ta = allPhyloDists[nodeTaxa, nodeTaxa]
    return( sd(ta[upper.tri(ta)]) / mean(ta[upper.tri(ta)]) )
    }

sharesEnvironment = function(taxa, nodeEnvs)
    {
    return(sapply(taxa, FUN=function(x) length(intersect(environments[[x]], nodeEnvs))>0))
    }

sharesSubtype = function(taxa, nodeEnvs)
    {
    return(sapply(taxa, FUN=function(x) length(intersect(subTypes[[x]], nodeEnvs))>0))
    }

generatePhyloEquiv = function(node, tol=PHYLO_SIM_TOLERANCE, allowDifferentEnv=F)
    {
    ###Generate random nodes in which the member taxa have the same distribution of pairwise phylogenetic distances than the original node.
    print(node)
    nodeTaxa = getNodeTaxa(node)
    if(allowDifferentEnv==T)
        {otherTaxa = goodTaxa[!goodTaxa %in% nodeTaxa]
    }else
        {
        nodeEnvs = getNodeEnvs(node)   
        otherTaxa = goodTaxa[!goodTaxa %in% nodeTaxa & sharesEnvironment(goodTaxa, nodeEnvs)]
        }
    clust = hclust(as.dist(allPhyloDists[nodeTaxa,nodeTaxa]), method='average')
    res = list()
    for(s in 1:(length(nodeTaxa)-1))
        {
        res[[s]] = character()
        merge = clust$merge[s,]
        targetDist = clust$height[s]
        if(merge[1]<0){cands1 = otherTaxa
        }else{cands1 = res[[merge[1]]]}
        if(merge[2]<0){cands2 = otherTaxa
        }else{cands2 = res[[merge[2]]]}
        print(paste(s, length(cands1), length(cands2)))
        #Check all possible pairings.
        for(i in 1:length(cands1))
            {
            for(j in 1:length(cands2))
                {
                c1 = unlist(strsplit(cands1[i], '-'))
                c2 = unlist(strsplit(cands2[j], '-'))
                if(length(intersect(c1, c2))==0) #No taxa in common.
                    {
                    if(abs(branchLength(c1, c2)-targetDist) < tol)
                        {
                        merged = paste(sort(c(c1, c2)), sep='-', collapse='-')
                        res[[s]] = union(res[[s]], merged)
                        }
                    }
                }
            }
        }
    return(res[[length(res)]])
    }


generateAvgPhyloEquiv = function(node, tol=PHYLO_SIM_TOLERANCE, replicates=AVG_PHY_EQ_REPLICATES, allowDifferentEnv=F, subtype=F, maxTries=10000, maxDist=0.05)
    ###Generate random nodes in which the member taxa have the same average pairwise phylogenetic distances than the original node.
    #For each node in our graph, will make 10 "phylogenetically equivalent" random nodes.
    {
    print(node)
    nodeTaxa = getNodeTaxa(node)
    if(allowDifferentEnv==T)
        {otherTaxa = goodTaxa[!goodTaxa %in% nodeTaxa]
    }else
        {
       	if(subtype) { nodeEnvs = getNodeSubtypes(node); otherTaxa = goodTaxa[!goodTaxa %in% nodeTaxa & sharesSubtype(goodTaxa, nodeEnvs)]
       	}else { nodeEnvs = getNodeEnvs(node); otherTaxa = goodTaxa[!goodTaxa %in% nodeTaxa & sharesEnvironment(goodTaxa, nodeEnvs)] }
        }

    meanPhyloDist = avgPhyloDist(node) #Loathing myself for using "mean" and "average" in the same line of code...
    res = c()
    tries = 0
    while(length(res) < replicates)
        {
        if(tol>=maxDist) 
            {
            break
            }
        if(length(otherTaxa) < length(nodeTaxa))
            {
            print('Not enough eligible taxa, skipping...')
            break
            }
        while(T)
            {
            if(tries == maxTries)
                {
                tol = tol + 0.001
                tries = 0
                print(paste('\tIncreasing tolerance to', tol))
                }
            if(tol>=maxDist)
                {
                print('Reached max dist. Skipping...')
                break
                }
            testComm = paste(sample(otherTaxa,length(nodeTaxa)), collapse='-')
            testMeanPhyloDist = avgPhyloDist(testComm)
            if(abs(testMeanPhyloDist-meanPhyloDist) < tol)
                {
                res = union(res, testComm)
                break
            }else
                {
                tries = tries + 1
                }
            }
        }
    return(res)
    }


generateAvgPwyEquiv = function(node, tol=0.05, replicates=100, allowDifferentEnv=F, maxTries=10000, maxDist=0.2, subtype=F)
    {
    ###Generate random nodes in which the member taxa have the same average number of pathways than the original node.
    print(node)
    nodeTaxa = getNodeTaxa(node)
    if(allowDifferentEnv==T)
        {otherTaxa = goodTaxa[!goodTaxa %in% nodeTaxa]
    }else
        {
       	if(subtype) { nodeEnvs = getNodeSubtypes(node); otherTaxa = goodTaxa[!goodTaxa %in% nodeTaxa & sharesSubtype(goodTaxa, nodeEnvs)]
       	}else { nodeEnvs = getNodeEnvs(node); otherTaxa = goodTaxa[!goodTaxa %in% nodeTaxa & sharesEnvironment(goodTaxa, nodeEnvs)] }
        }

    avgPathways = getAvgPwysPerNodeMember(node) #Loathing myself for using "mean" and "average" in the same line of code...
    res = c()
    tries = 0
    while(length(res) < replicates)
        {
        if(tol>=maxDist) 
            {
            break
            }
        if(length(otherTaxa) < length(nodeTaxa))
            {
            print('Not enough eligible taxa, skipping...')
            break
            }
        while(T)
            {
            if(tries == maxTries)
                {
                tol = tol + 0.01
                tries = 0
                print(paste('\tIncreasing tolerance to', tol))
                }
            if(tol>=maxDist)
                {
                print('Reached max dist. Skipping...')
                break
                }
            testComm = paste(sample(otherTaxa,length(nodeTaxa)), collapse='-')
            testAvgPathways = getAvgPwysPerNodeMember(testComm)
            if(abs(testAvgPathways-avgPathways)/avgPathways < tol)
                {
                res = union(res, testComm)
                break
            }else
                {
                tries = tries + 1
                }
            }
        }
    return(res)
    }

	


generateEnvEquiv = function(node, replicates=AVG_PHY_EQ_REPLICATES)
    {
    nodeTaxa = getNodeTaxa(node)
    size = length(nodeTaxa)
    nodeEnvs = getNodeEnvs(node) 
    envTaxa = setdiff(goodTaxa[sharesEnvironment(goodTaxa, nodeEnvs)], nodeTaxa)
    if(length(envTaxa)<size) { return(c()) }
    res = c()
    for(i in 1:replicates)
        {
        taxa = sample(envTaxa, size)
        res = c(res, paste(taxa, collapse='-'))
        }
    return(res)
    }


generateEnvEquivSubtype = function(node, replicates=AVG_PHY_EQ_REPLICATES)
    {
    nodeTaxa = getNodeTaxa(node)
    size = length(nodeTaxa)
    nodeSubtypes = getNodeSubtypes(node) 
    envTaxa = setdiff(goodTaxa[sharesSubtype(goodTaxa, nodeSubtypes)], nodeTaxa)
    if(length(envTaxa)<size) { return(c()) }
    res = c()
    for(i in 1:replicates)
        {
        taxa = sample(envTaxa, size)
        res = c(res, paste(taxa, collapse='-'))
        }
    return(res)
    }

generateAvgFunEquiv = function(node, tol=FUN_SIM_TOLERANCE, replicates=AVG_FUN_EQ_REPLICATES, allowDifferentEnv=F, maxTries=10000, maxDist=0.05)
    ###Generate random nodes in which the member taxa have the same average pairwise functional similarities than the original node.
    #For each node in our graph, will make 10 "functionally equivalent" random nodes.
    {
    print(node)
    nodeTaxa = getNodeTaxa(node)
    if(allowDifferentEnv==T)
        {otherTaxa = goodTaxa[!goodTaxa %in% nodeTaxa]
    }else
        {
        nodeEnvs = getNodeEnvs(node) 
        otherTaxa = goodTaxa[!goodTaxa %in% nodeTaxa & sharesEnvironment(goodTaxa, nodeEnvs)]
        }

    meanFunCorr = avgFunCorr(node) #Loathing myself for using "mean" and "average" in the same line of code...
    res = c()
    tries = 0
    while(length(res) < replicates)
        {
        if(tol>=maxDist) 
            {
            break
            }
        if(length(otherTaxa) < length(nodeTaxa))
            {
            print('Not enough eligible taxa, skipping...')
            break
            }
        while(T)
            {
            if(tries == maxTries)
                {
                tol = tol + 0.001
                tries = 0
                print(paste('\tIncreasing tolerance to', tol))
                }
            if(tol>=maxDist)
                {
                print('Reached max dist. Skipping...')
                break
                }
            testComm = paste(sample(otherTaxa,length(nodeTaxa)), collapse='-')
            testMeanFunCorr = avgFunCorr(testComm)
            if(abs(testMeanFunCorr-meanFunCorr) < tol)
                {
                res = union(res, testComm)
                break
            }else
                {
                tries = tries + 1
                }
            }
        }
    return(res)
    }



generateRandomComm = function(size, replicates=RAND_REPLICATES)
    {
    res = c()
    for(i in 1:replicates)
        {
        taxa = sample(goodTaxa, size)
        res = c(res, paste(taxa, collapse='-'))
        }
    return(res)
    }


getNodeTaxa = function(node) return(unlist(strsplit(as.character(node), '-', fixed=T)))
getNodeSize = function(node) return(length(getNodeTaxa(node)))
getNodeEnvs = function(node) return(unlist(strsplit(as.character(nodes[node, 'sourceNetwork']), '|', fixed=T)))
getNodeSubtypes = function(node) { st = unlist(strsplit(as.character(nodes[node, 'envSubtypes']), '|', fixed=T)); return(st[seq(1,length(st),2)])}
getNodeSamples = function(node) return(unlist(strsplit(as.character(nodes[node, 'samples']), '|', fixed=T)))
getNodeFunAnnotTaxa = function(node) { tx = getNodeTaxa(node); return(tx[tx %in% rownames(pathways)]) }
getNodePhyloAnnotTaxa = function(node) { tx = getNodeTaxa(node); return(tx[tx %in% rownames(allPhyloDists)]) }
getNodeAnnotTaxa = function(node) intersect(getNodeFunAnnotTaxa(node), getNodePhyloAnnotTaxa(node))
getAvgPwysPerNodeMember = function(node) return(mean(totalPathwaysPerGenus[getNodeTaxa(node)]))
getNodeUniquePathways = function(node) return(table(colSums(pathways[getNodeTaxa(node),])>0)[2])
getNodeTotalPathways = function(node) return(sum(totalPathwaysPerGenus[getNodeTaxa(node)]))
getNodePwyPrevalence = function(node) return(colSums(pathways[getNodeTaxa(node),]))



####DATA


###Static vars.

PATHWAY_PRESENCE_CUTOFF=0.25
MIN_NODE_SUPPORT = 10
RAND_REPLICATES=1000
AVG_PHY_EQ_REPLICATES = 10
PHYLO_SIM_TOLERANCE = 0.005
FUN_SIM_TOLERANCE = 0.005
AVG_FUN_EQ_REPLICATES = 10


###Preload data.
# The file `merged6.2.networkTable.csv` is a csv representation of the network.
# It can be obtained by loading the file `merged6.2.unlooped.avgFunPhyl.pathways.xml` into cytoscape and exporting the node table
nodes=read.table('../02-network-generation/02-output/merged6.2.networkTable.csv', header=T, sep=',', row.names=9)
goodNodes=subset(nodes, subset= support >= MIN_NODE_SUPPORT & funAnnotTaxa > 1 & phyloAnnotTaxa > 1 & phyloAnnotTaxa==funAnnotTaxa & funAnnotTaxa==numTaxa)
tipNodes=subset(goodNodes, subset=outEdges==0)
singleEnv = goodNodes[goodNodes$numSourceNetworks==1,]
multiEnv = goodNodes[goodNodes$numSourceNetworks>1,]
bifurcating = goodNodes[goodNodes$outEdges>1,]


goodTaxa = c()
for (i in 1:length(rownames(goodNodes))){
    taxa = unlist(strsplit(as.character(rownames(goodNodes)[i]), '-', fixed = TRUE))
    goodTaxa = c(goodTaxa, taxa)
    }

goodTaxa=unique(goodTaxa)

environments = list()
subTypes = list()
for(taxon in goodTaxa)
    {
    environments[[taxon]] = strsplit(as.character(nodes[taxon, 'sourceNetwork']), '|', fixed=T)[[1]]
    st = strsplit(as.character(nodes[taxon, 'envSubtypes']), '|', fixed=T)[[1]]
    subTypes[[taxon]] = st[seq(1,length(st),2)]
    }


envs = unlist(strsplit(as.character(goodNodes$sourceNetwork[i]), '|', fixed=T))

allPhyloDists = read.table('../02-network-generation/goodMetaCyc/phylodist_clean.table.tsv', header=T, row.names=1, sep='\t')

pathways = read.table('../02-network-generation/goodMetaCyc/oct2020.combined_noPWY0-1324.tsv', check.names=F, header=T, row.names=1, sep='\t')
binarize = function(table, cutoff)
	{
	ta2 = table
	ta2[ta2>cutoff]=1
        ta2[ta2<1]=0
        return(ta2)
	}
pathways = binarize(pathways, PATHWAY_PRESENCE_CUTOFF)
pathwayNames_ = read.table('../02-network-generation/goodMetaCyc/metaCyc.names.tsv', sep='\t', row.names=1, check.names=F, quote='')
pathwayNames = as.character(pathwayNames_[,1])
names(pathwayNames) = rownames(pathwayNames_)

noName = setdiff(colnames(pathways), names(pathwayNames)) #Pathways not in our names file.
names(noName) = noName
pathwayNames = c(pathwayNames, noName)

colnames(pathways) = pathwayNames[colnames(pathways)]

library(SQMtools)
exportTable(pathways, 'pathwaysPerGenus.binarized.tsv')
totalPathwaysPerGenus = rowSums(pathways)

allFunCorrs = 1-as.matrix(dist.binary(pathways, method=1)) #Jaccard
#allFunCorrs = as.matrix(dist.binary(pathways, method=1)) #Jaccard || USING DISTANCES INSTEAD OF CORRELATIONS!!!!!!!!!!!

graphPhyloDists = allPhyloDists[goodTaxa, goodTaxa]
graphFunCorrs = allFunCorrs[goodTaxa, goodTaxa]

maxFunAnnotTaxa = max(goodNodes$funAnnotTaxa)
maxPhyloAnnotTaxa = max(goodNodes$phyloAnnotTaxa)


### Genomes per genus in the genera included in the network
gpg=read.table('../01-pathways-per-genus/genomes_per_genus_raw.tsv', sep='\t')
gpg=gpg[gpg[,1] %in% row.names(pathways),]
gpg[,2] = as.numeric(as.character(gpg[,2]))
write.table(gpg[order(gpg[,2], decreasing=T),], file='genomes_per_genus.tsv')





###Generate fully random comms.
library(parallel)
set.seed(132)
randComms = unlist(lapply(goodNodes$numTaxa, generateRandomComm))
numTaxaRand = sapply(randComms, getNodeSize)
meanFunCorrsRand = sapply(randComms, avgFunCorr)
meanPhyloDistsRand = sapply(randComms, avgPhyloDist)
cvFunCorrsRand = sapply(randComms, cvFunCorr)              #Note NAs when numTaxa==2.
cvPhyloDistsRand = sapply(randComms, cvPhyloDist)          #Note NAs when numTaxa==2.
avgPwysRand = unlist(mclapply(randComms, getAvgPwysPerNodeMember, mc.cores=6))
totalPwysRand = unlist(mclapply(randComms, getNodeTotalPathways, mc.cores=6))


###Generate env-equiv comms
set.seed(139)
envEquivComms = unlist(lapply(rownames(goodNodes), generateEnvEquivSubtype)) # we could also use generateEnvEquivSubtype here
numTaxaEquivEnv = sapply(envEquivComms, getNodeSize)
meanFunCorrsEquivEnv = sapply(envEquivComms, avgFunCorr)
meanPhyloDistsEquivEnv = sapply(envEquivComms, avgPhyloDist)
cvFunCorrsEquivEnv = sapply(envEquivComms, cvFunCorr)     #Note NAs when numTaxa==2.
cvPhyloDistsEquivEnv = sapply(envEquivComms, cvPhyloDist) #Note NAs when numTaxa==2.
avgPwysEquivEnv = sapply(envEquivComms, getAvgPwysPerNodeMember)
totalPwysEquivEnv = sapply(envEquivComms, getNodeTotalPathways) 

###Generate avg-equiv comms.
set.seed(133)
avgEquivComms = unlist(lapply(rownames(goodNodes), generateAvgPhyloEquiv, allowDifferentEnv=T, maxTries=10000))
numTaxaEquivRand = sapply(avgEquivComms, getNodeSize)
meanFunCorrsEquivRand = sapply(avgEquivComms, avgFunCorr)
meanPhyloDistsEquivRand = sapply(avgEquivComms, avgPhyloDist)
cvFunCorrsEquivRand = sapply(avgEquivComms, cvFunCorr)     #Note NAs when numTaxa==2.
cvPhyloDistsEquivRand = sapply(avgEquivComms, cvPhyloDist) #Note NAs when numTaxa==2.
avgPwysEquivRand = sapply(avgEquivComms, getAvgPwysPerNodeMember)
totalPwysEquivRand = sapply(avgEquivComms, getNodeTotalPathways)

###Generate avg-fun equiv comms.
set.seed(134)
avgFunEquivComms = unlist(lapply(rownames(goodNodes), generateAvgFunEquiv, allowDifferentEnv=T, maxTries=10000))
numTaxaFunEquivRand = sapply(avgFunEquivComms, getNodeSize)
meanFunCorrsFunEquivRand = sapply(avgFunEquivComms, avgFunCorr)
meanPhyloDistsFunEquivRand = sapply(avgFunEquivComms, avgPhyloDist)
cvFunCorrsFunEquivRand = sapply(avgFunEquivComms, cvFunCorr)     #Note NAs when numTaxa==2.
cvPhyloDistsFunEquivRand = sapply(avgFunEquivComms, cvPhyloDist) #Note NAs when numTaxa==2.

###Generate avg-equiv comms, same environment.
set.seed(135)
avgEquivSameEnvComms = unlist(lapply(rownames(goodNodes), generateAvgPhyloEquiv, allowDifferentEnv=F, maxTries=10000, tol = 0.003, maxDist = 0.01, subtype=T))
numTaxaEquivSameEnvRand = sapply(avgEquivSameEnvComms, getNodeSize)
meanFunCorrsEquivSameEnvRand = sapply(avgEquivSameEnvComms, avgFunCorr)
meanPhyloDistsEquivSameEnvRand = sapply(avgEquivSameEnvComms, avgPhyloDist)
cvFunCorrsEquivSameEnvRand = sapply(avgEquivSameEnvComms, cvFunCorr)     #Note NAs when numTaxa==2.
cvPhyloDistsEquivSameEnvRand = sapply(avgEquivSameEnvComms, cvPhyloDist) #Note NAs when numTaxa==2.
totalPwysEquivSameEnvRand = sapply(avgEquivSameEnvComms, getNodeTotalPathways)
avgPwysEquivSameEnvRand = sapply(avgEquivSameEnvComms, getAvgPwysPerNodeMember)


### generate avg-pathway-equiv comms, same environments
set.seed(140)
avgPwyEquivSameEnvComms = unlist(lapply(rownames(goodNodes), generateAvgPwyEquiv, allowDifferentEnv=F, maxTries=10000, subtype=T))
numTaxaPwyEquivSameEnvRand = sapply(avgPwyEquivSameEnvComms, getNodeSize)
meanFunCorrsPwyEquivSameEnvRand = sapply(avgPwyEquivSameEnvComms, avgFunCorr)
meanPhyloDistsPwyEquivSameEnvRand = sapply(avgPwyEquivSameEnvComms, avgPhyloDist)
totalPwysPwyEquivSameEnvRand = sapply(avgPwyEquivSameEnvComms, getNodeTotalPathways)

###Single-environment communities.
set.seed(136)
numTaxaSE = sapply(rownames(singleEnv), getNodeSize)
meanFunCorrsSE = sapply(rownames(singleEnv), avgFunCorr)
meanPhyloDistsSE = sapply(rownames(singleEnv), avgPhyloDist)
cvFunCorrsSE = sapply(rownames(singleEnv), cvFunCorr)      #Note NAs when numTaxa==2.
cvPhyloDistsSE = sapply(rownames(singleEnv), cvPhyloDist)  #Note NAs when numTaxa==2.
avgPwysSE = sapply(rownames(singleEnv), getAvgPwysPerNodeMember)
totalPwysSE = sapply(rownames(singleEnv), getNodeTotalPathways)

###Multi-environment communities.
set.seed(137)
numTaxaME = sapply(rownames(multiEnv), getNodeSize)
meanFunCorrsME = sapply(rownames(multiEnv), avgFunCorr)
meanPhyloDistsME = sapply(rownames(multiEnv),avgPhyloDist) 
cvFunCorrsME = sapply(rownames(multiEnv), cvFunCorr)       #Note NAs when numTaxa==2.
cvPhyloDistsME = sapply(rownames(multiEnv),cvPhyloDist)    #Note NAs when numTaxa==2.
avgPwysME = sapply(rownames(multiEnv), getAvgPwysPerNodeMember)
totalPwysME = sapply(rownames(multiEnv), getNodeTotalPathways)






###Wilcoxon tests.
#Function: Single Env - Random
wilcox.test(meanFunCorrsSE[numTaxaSE==2], meanFunCorrsRand[numTaxaRand==2]) #2.2e-16
wilcox.test(meanFunCorrsSE[numTaxaSE==3], meanFunCorrsRand[numTaxaRand==3]) #2.2e-16
wilcox.test(meanFunCorrsSE[numTaxaSE==4], meanFunCorrsRand[numTaxaRand==4]) #3.747e-12
wilcox.test(meanFunCorrsSE[numTaxaSE==5], meanFunCorrsRand[numTaxaRand==5]) #1.739e-06
wilcox.test(meanFunCorrsSE[numTaxaSE==6], meanFunCorrsRand[numTaxaRand==6]) #0.0002515
wilcox.test(meanFunCorrsSE[numTaxaSE==7], meanFunCorrsRand[numTaxaRand==7]) #0.0007565
wilcox.test(meanFunCorrsSE[numTaxaSE==8], meanFunCorrsRand[numTaxaRand==8]) #0.1134
wilcox.test(meanFunCorrsSE[numTaxaSE==9], meanFunCorrsRand[numTaxaRand==9]) #0.08021
wilcox.test(meanFunCorrsSE[numTaxaSE==12], meanFunCorrsRand[numTaxaRand==12]) #0.2059

#Function: Single Env - phyloEquiv
wilcox.test(meanFunCorrsSE[numTaxaSE==2], meanFunCorrsEquivRand[numTaxaEquivRand==2]) #0.00173
wilcox.test(meanFunCorrsSE[numTaxaSE==3], meanFunCorrsEquivRand[numTaxaEquivRand==3]) #0.0001754
wilcox.test(meanFunCorrsSE[numTaxaSE==4], meanFunCorrsEquivRand[numTaxaEquivRand==4]) #0.004186
wilcox.test(meanFunCorrsSE[numTaxaSE==5], meanFunCorrsEquivRand[numTaxaEquivRand==5]) #0.1338
wilcox.test(meanFunCorrsSE[numTaxaSE==6], meanFunCorrsEquivRand[numTaxaEquivRand==6]) #0.2317
wilcox.test(meanFunCorrsSE[numTaxaSE==7], meanFunCorrsEquivRand[numTaxaEquivRand==7]) #0.1701
wilcox.test(meanFunCorrsSE[numTaxaSE==8], meanFunCorrsEquivRand[numTaxaEquivRand==8]) #0.1732
wilcox.test(meanFunCorrsSE[numTaxaSE==9], meanFunCorrsEquivRand[numTaxaEquivRand==9]) #0.7792
wilcox.test(meanFunCorrsSE[numTaxaSE==12], meanFunCorrsEquivRand[numTaxaEquivRand==12]) #0.5455

#Function: Single Env - phyloEquivSameEnv
wilcox.test(meanFunCorrsSE[numTaxaSE==2], meanFunCorrsEquivSameEnvRand[numTaxaEquivSameEnvRand==2]) #0.0331
wilcox.test(meanFunCorrsSE[numTaxaSE==3], meanFunCorrsEquivSameEnvRand[numTaxaEquivSameEnvRand==3]) #0.001646
wilcox.test(meanFunCorrsSE[numTaxaSE==4], meanFunCorrsEquivSameEnvRand[numTaxaEquivSameEnvRand==4]) #0.04957
wilcox.test(meanFunCorrsSE[numTaxaSE==5], meanFunCorrsEquivSameEnvRand[numTaxaEquivSameEnvRand==5]) #0.05345
wilcox.test(meanFunCorrsSE[numTaxaSE==6], meanFunCorrsEquivSameEnvRand[numTaxaEquivSameEnvRand==6]) #0.1135
wilcox.test(meanFunCorrsSE[numTaxaSE==7], meanFunCorrsEquivSameEnvRand[numTaxaEquivSameEnvRand==7]) #0.01697
wilcox.test(meanFunCorrsSE[numTaxaSE==8], meanFunCorrsEquivSameEnvRand[numTaxaEquivSameEnvRand==8]) #0.1039
wilcox.test(meanFunCorrsSE[numTaxaSE==9], meanFunCorrsEquivSameEnvRand[numTaxaEquivSameEnvRand==9]) #0.2597
wilcox.test(meanFunCorrsSE[numTaxaSE==12], meanFunCorrsEquivSameEnvRand[numTaxaEquivSameEnvRand==12]) #0.1818

#Function: Single Env - pwyEquivSameEnv
wilcox.test(meanFunCorrsSE[numTaxaSE==2], numTaxaPwyEquivSameEnvRand[numTaxaEquivSameEnvRand==2]) #2.2e-16
wilcox.test(meanFunCorrsSE[numTaxaSE==3], numTaxaPwyEquivSameEnvRand[numTaxaEquivSameEnvRand==3]) #2.2e-16
wilcox.test(meanFunCorrsSE[numTaxaSE==4], numTaxaPwyEquivSameEnvRand[numTaxaEquivSameEnvRand==4]) #2.2e-16
wilcox.test(meanFunCorrsSE[numTaxaSE==5], numTaxaPwyEquivSameEnvRand[numTaxaEquivSameEnvRand==5]) #2.034e-14
wilcox.test(meanFunCorrsSE[numTaxaSE==6], numTaxaPwyEquivSameEnvRand[numTaxaEquivSameEnvRand==6]) #3.436e-08
wilcox.test(meanFunCorrsSE[numTaxaSE==7], numTaxaPwyEquivSameEnvRand[numTaxaEquivSameEnvRand==7]) #2.236e-06
wilcox.test(meanFunCorrsSE[numTaxaSE==8], numTaxaPwyEquivSameEnvRand[numTaxaEquivSameEnvRand==8]) #0.002612
wilcox.test(meanFunCorrsSE[numTaxaSE==9], numTaxaPwyEquivSameEnvRand[numTaxaEquivSameEnvRand==9]) #0.006892
wilcox.test(meanFunCorrsSE[numTaxaSE==12], numTaxaPwyEquivSameEnvRand[numTaxaEquivSameEnvRand==12]) #0.03581


#Function: Single Env - MultiEnv
wilcox.test(meanFunCorrsSE[numTaxaSE==2], meanFunCorrsME[numTaxaME==2]) #0.0007302
wilcox.test(meanFunCorrsSE[numTaxaSE==3], meanFunCorrsME[numTaxaME==3]) #0.1211


#Phylogeny: Single Env - Random
wilcox.test(meanPhyloDistsSE[numTaxaSE==2], meanPhyloDistsRand[numTaxaRand==2]) #2.2e-16
wilcox.test(meanPhyloDistsSE[numTaxaSE==3], meanPhyloDistsRand[numTaxaRand==3]) #2.2e-16
wilcox.test(meanPhyloDistsSE[numTaxaSE==4], meanPhyloDistsRand[numTaxaRand==4]) #1.456e-10
wilcox.test(meanPhyloDistsSE[numTaxaSE==5], meanPhyloDistsRand[numTaxaRand==5]) #1.243e-06
wilcox.test(meanPhyloDistsSE[numTaxaSE==6], meanPhyloDistsRand[numTaxaRand==6]) #0.008084
wilcox.test(meanPhyloDistsSE[numTaxaSE==7], meanPhyloDistsRand[numTaxaRand==7]) #0.09589
wilcox.test(meanPhyloDistsSE[numTaxaSE==8], meanPhyloDistsRand[numTaxaRand==8]) #0.7341
wilcox.test(meanPhyloDistsSE[numTaxaSE==9], meanPhyloDistsRand[numTaxaRand==9]) #0.2167
wilcox.test(meanPhyloDistsSE[numTaxaSE==12], meanPhyloDistsRand[numTaxaRand==12]) #0.6317

#Phylogeny: Single Env - pwyEquivSameEnv
wilcox.test(meanPhyloDistsSE[numTaxaSE==2], meanPhyloDistsPwyEquivSameEnvRand[numTaxaPwyEquivSameEnvRand==2]) #2.2e-16
wilcox.test(meanPhyloDistsSE[numTaxaSE==3], meanPhyloDistsPwyEquivSameEnvRand[numTaxaPwyEquivSameEnvRand==3]) #2.2e-16
wilcox.test(meanPhyloDistsSE[numTaxaSE==4], meanPhyloDistsPwyEquivSameEnvRand[numTaxaPwyEquivSameEnvRand==4]) #1.042e-06
wilcox.test(meanPhyloDistsSE[numTaxaSE==5], meanPhyloDistsPwyEquivSameEnvRand[numTaxaPwyEquivSameEnvRand==5]) #2.795e-05
wilcox.test(meanPhyloDistsSE[numTaxaSE==6], meanPhyloDistsPwyEquivSameEnvRand[numTaxaPwyEquivSameEnvRand==6]) #0.0239
wilcox.test(meanPhyloDistsSE[numTaxaSE==7], meanPhyloDistsPwyEquivSameEnvRand[numTaxaPwyEquivSameEnvRand==7]) #0.3676
wilcox.test(meanPhyloDistsSE[numTaxaSE==8], meanPhyloDistsPwyEquivSameEnvRand[numTaxaPwyEquivSameEnvRand==8]) #0.841
wilcox.test(meanPhyloDistsSE[numTaxaSE==9], meanPhyloDistsPwyEquivSameEnvRand[numTaxaPwyEquivSameEnvRand==9]) #0.1397
wilcox.test(meanPhyloDistsSE[numTaxaSE==12], meanPhyloDistsPwyEquivSameEnvRand[numTaxaPwyEquivSameEnvRand==12]) #0.5257

#Phylogeny: Single Env - MultiEnv
wilcox.test(meanPhyloDistsSE[numTaxaSE==2], meanPhyloDistsME[numTaxaME==2]) #0.03546
wilcox.test(meanPhyloDistsSE[numTaxaSE==3], meanPhyloDistsME[numTaxaME==3]) #0.4743


#Avg pwys: Single Env - Random
wilcox.test(avgPwysSE[numTaxaSE==2], avgPwysRand[numTaxaRand==2]) #0.0005198
wilcox.test(avgPwysSE[numTaxaSE==3], avgPwysRand[numTaxaRand==3]) #0.2727
wilcox.test(avgPwysSE[numTaxaSE==4], avgPwysRand[numTaxaRand==4]) #0.0005599
wilcox.test(avgPwysSE[numTaxaSE==5], avgPwysRand[numTaxaRand==5]) #0.8894
wilcox.test(avgPwysSE[numTaxaSE==6], avgPwysRand[numTaxaRand==6]) #0.06559
wilcox.test(avgPwysSE[numTaxaSE==7], avgPwysRand[numTaxaRand==7]) #0.1333
wilcox.test(avgPwysSE[numTaxaSE==8], avgPwysRand[numTaxaRand==8]) #0.01836
wilcox.test(avgPwysSE[numTaxaSE==9], avgPwysRand[numTaxaRand==9]) #0.02012
wilcox.test(avgPwysSE[numTaxaSE==12], avgPwysRand[numTaxaRand==12]) #0.08576

#Avg pwys: Single Env - phyloEquivSameEnv
wilcox.test(avgPwysSE[numTaxaSE==2], avgPwysEquivSameEnvRand[numTaxaEquivSameEnvRand==2]) #2.968e-07
wilcox.test(avgPwysSE[numTaxaSE==3], avgPwysEquivSameEnvRand[numTaxaEquivSameEnvRand==3]) #2.628e-05
wilcox.test(avgPwysSE[numTaxaSE==4], avgPwysEquivSameEnvRand[numTaxaEquivSameEnvRand==4]) #0.8916
wilcox.test(avgPwysSE[numTaxaSE==5], avgPwysEquivSameEnvRand[numTaxaEquivSameEnvRand==5]) #0.1317
wilcox.test(avgPwysSE[numTaxaSE==6], avgPwysEquivSameEnvRand[numTaxaEquivSameEnvRand==6]) #0.03543
wilcox.test(avgPwysSE[numTaxaSE==7], avgPwysEquivSameEnvRand[numTaxaEquivSameEnvRand==7]) #0.06415
wilcox.test(avgPwysSE[numTaxaSE==8], avgPwysEquivSameEnvRand[numTaxaEquivSameEnvRand==8]) #0.008658
wilcox.test(avgPwysSE[numTaxaSE==9], avgPwysEquivSameEnvRand[numTaxaEquivSameEnvRand==9]) #0.0259
wilcox.test(avgPwysSE[numTaxaSE==12], avgPwysEquivSameEnvRand[numTaxaEquivSameEnvRand==12]) #0.1818


#Avg pwys: Single Env - MultiEnv
wilcox.test(avgPwysSE[numTaxaSE==2], avgPwysME[numTaxaME==2]) #0.4947
wilcox.test(avgPwysSE[numTaxaSE==3], avgPwysME[numTaxaME==3]) #0.3013



pwyAbund = function(node, pwy=NULL)
    {
    taxa = getNodeTaxa(node)
    if(is.null(pwy))
        {return(colSums(pathways[taxa,]))
    }else
        {return(sum(pathways[taxa,pwy]))}
    }



###PLOT STUFF.

library(ggplot2)

figureName = ''
outName = 'all.subtype'

df = data.frame(row.names=1:length(c(totalPwysSE, totalPwysEquivEnv, totalPwysRand[numTaxaRand<=maxFunAnnotTaxa], totalPwysME, totalPwysPwyEquivSameEnvRand, totalPwysEquivSameEnvRand)))
df['numTaxa'] = as.factor(c(numTaxaRand[numTaxaRand<=maxFunAnnotTaxa], numTaxaEquivEnv, numTaxaEquivSameEnvRand, numTaxaPwyEquivSameEnvRand, numTaxaSE, numTaxaME))
df['totalPwys'] = c(totalPwysRand[numTaxaRand<=maxFunAnnotTaxa], totalPwysEquivEnv, totalPwysEquivSameEnvRand, totalPwysPwyEquivSameEnvRand, totalPwysSE, totalPwysME)
df['avgPwys'] = df['totalPwys']/as.numeric(as.character(df[,'numTaxa']))
df['PhyloDists'] = c(meanPhyloDistsRand[numTaxaRand<=maxFunAnnotTaxa], meanPhyloDistsEquivEnv, meanPhyloDistsEquivSameEnvRand, meanPhyloDistsPwyEquivSameEnvRand, meanPhyloDistsSE, meanPhyloDistsME)
df['funDists'] = 1 - c(meanFunCorrsRand[numTaxaRand<=maxFunAnnotTaxa], meanFunCorrsEquivEnv, meanFunCorrsEquivSameEnvRand, meanFunCorrsPwyEquivSameEnvRand, meanFunCorrsSE, meanFunCorrsME)
df['type'] = factor(
                    c(
                      rep('Random', length(meanPhyloDistsRand[numTaxaRand<=maxFunAnnotTaxa])),
                      rep('Env-equiv random', length(meanFunCorrsEquivEnv)),
                      rep('Env-phylo-equiv random', length(meanFunCorrsEquivSameEnvRand)),
                      rep('Env-pathway-equiv random', length(meanFunCorrsPwyEquivSameEnvRand)),
                      rep('Single-Environment', length(meanPhyloDistsSE)),
                      rep('Multi-Environment', length(meanPhyloDistsME))
                     ), 
                    levels=c('Random', 'Env-equiv random', 'Env-phylo-equiv random', 'Env-pathway-equiv random', 'Single-Environment', 'Multi-Environment')
                   )

df = na.omit(df)


ggplot(df, aes(x = numTaxa, y = funDists, fill = type)) + 
       geom_boxplot() + 
       scale_fill_manual(values=c('grey30', 'grey60', 'gray80', 'gray90', 'royalblue1', 'darkolivegreen1')) +
       labs(x='Assemblage size', y='Intra-assemblage functional distance (Jaccard)', fill='Assemblage type') + 
       geom_hline(aes(yintercept=1-mean(as.matrix(graphFunCorrs), na.rm=T), col = 'Average pairwise\nfunctional distance')) + 
       labs(colour='') +
       guides(fill = guide_legend(order=1), colour = guide_legend(order=2)) +
       ylim(0,1) + 
       ggtitle(figureName)
ggsave(paste(outName, '.funDists.png', sep=''))



ggplot(df, aes(x = numTaxa, y = PhyloDists, fill = type)) + 
       geom_boxplot() + 
       scale_fill_manual(values=c('grey30', 'grey60', 'gray80', 'gray90', 'royalblue1', 'darkolivegreen1')) +
       labs(x='Assemblage size', y='Intra-assemblage phylogenetic distance (subs/site)', fill='Assemblage type') + 
       geom_hline(aes(yintercept=mean(as.matrix(graphPhyloDists)), col = 'Average pairwise\nphylogenetic distance')) + 
       labs(colour='') +
       guides(fill = guide_legend(order=1), colour = guide_legend(order=2)) +
       ggtitle(figureName)
ggsave(paste(outName, '.phyloDists.png', sep=''))


ggplot(df, aes(x = numTaxa, y = avgPwys, fill = type)) + 
       geom_boxplot() + 
       scale_fill_manual(values=c('grey30', 'grey60', 'gray80', 'gray90', 'royalblue1', 'darkolivegreen1')) +
       labs(x='Assembly size', y='Average number of pathways in the assemblage', fill='Assemblage type') + 
       guides(fill = guide_legend(order=1), colour = guide_legend(order=2)) +
       ggtitle(figureName)
ggsave(paste(outName, '.avgPathways.png', sep=''))




# Fun, distances, bifurcating
df = data.frame(row.names=1:length(c(meanFunCorrsSE, meanFunCorrsEquivEnv, meanFunCorrsRand[numTaxaRand<=maxFunAnnotTaxa], meanFunCorrsME, meanFunCorrsBF, meanFunCorrsEquivSameEnvRand)))
df['numTaxa'] = as.factor(c(numTaxaRand[numTaxaRand<=maxFunAnnotTaxa], numTaxaEquivEnv, numTaxaEquivSameEnvRand, numTaxaSE, numTaxaBF, numTaxaME))
df['funDists'] = 1 - c(meanFunCorrsRand[numTaxaRand<=maxFunAnnotTaxa], meanFunCorrsEquivEnv, meanFunCorrsEquivSameEnvRand, meanFunCorrsSE, meanFunCorrsBF, meanFunCorrsME)
df['type'] = factor(
                    c(
                      rep('Random', length(meanFunCorrsRand[numTaxaRand<=maxFunAnnotTaxa])),
                      rep('Env-equiv random', length(meanFunCorrsEquivEnv)),
                      rep('Env-phylo-equiv random', length(meanFunCorrsEquivSameEnvRand)),
                      rep('Single-Environment', length(meanFunCorrsSE)),
                      rep('Bifurcating', length(meanFunCorrsBF)),
                      rep('Multi-Environment', length(meanFunCorrsME))
                     ), 
                    levels=c('Random', 'Env-equiv random', 'Env-phylo-equiv random', 'Single-Environment', 'Bifurcating', 'Multi-Environment')
                   )
df = na.omit(df)

ggplot(df, aes(x = numTaxa, y = funDists, fill = type)) + 
       geom_boxplot() + 
       scale_fill_manual(values=c('grey30', 'grey60', 'gray80', 'royalblue1', 'yellow2', 'darkolivegreen1')) +
       labs(x='Assemblage size', y='Intra-assemblage functional distance (Jaccard)', fill='Assemblage type') + 
       geom_hline(aes(yintercept=1-mean(as.matrix(graphFunCorrs), na.rm=T), col = 'Average pairwise\nfunctional distance')) + 
       labs(colour='') +
       guides(fill = guide_legend(order=1), colour = guide_legend(order=2)) +
       ylim(0,1) + 
       ggtitle(figureName)
ggsave(paste(outName, '.funDists.bifurcating.png', sep=''))







####### INDIVIDUAL PATHWAYS

library(purrr)
allAnnotTaxa = intersect(rownames(nodes)[nodes$numTaxa==1], rownames(pathways))
goodPwys = colnames(pathways)[colSums(pathways[allAnnotTaxa,])/length(allAnnotTaxa)>0]
goodPwys = goodPwys[order(colSums(pathways[allAnnotTaxa, goodPwys]))]

consensusEdges = read.table('../02-network-generation/02-output/consensusNet.sif', sep='\t', as.is=T)[,-2] # the second column is empty
consensusNodes = nodes[unique(c(consensusEdges[,1], consensusEdges[,2])),]
tipConsensusNodes = rownames(consensusNodes)[!rownames(consensusNodes) %in% consensusEdges[,1]] # no successors
tipConsensusNodes = tipConsensusNodes[sapply(tipConsensusNodes, FUN=function(x) length(getNodeAnnotTaxa(x))) > 1]
failed = c()

NUMNULLS = 1000
allNulls = list()
set.seed(7)
for(real0 in tipConsensusNodes) {
    real = paste(getNodeAnnotTaxa(real0), collapse='-')
    realDist = avgPhyloDist(real)
    if(F) {
        nulls = generateAvgPwyEquiv(real, allowDifferentEnv=F, maxTries=100000, subtype=T, replicates=10000)
        if(length(nulls) < NUMNULLS) { failed = c(failed, real0); next }   
        nulls = nulls[abs(sapply(nulls, avgPhyloDist) - realDist) < 0.1]
    } else { nulls = generateAvgPhyloEquiv(real, allowDifferentEnv=F, maxTries=10000, replicates=1000, tol = 0.003, maxDist = 0.1, subtype=T) }
    if(length(nulls) < NUMNULLS) { failed = c(failed, real0); next }
    nulls = nulls[1:NUMNULLS]
    allNulls[[real0]] = nulls
    }

#allNullsSTORED = allNulls # allNullsSTORED has the nulls phyloEnvEquiv subtype=T for all the tipConsensus nodes with 2+ annotated taxa
tipConsensusNodes = names(allNulls)

res = matrix(NA, nrow=length(tipConsensusNodes), ncol=length(goodPwys), dimnames=list(tipConsensusNodes, goodPwys))
res_fdr = matrix(NA, nrow=length(tipConsensusNodes), ncol=length(goodPwys), dimnames=list(tipConsensusNodes, goodPwys))
resB = matrix(NA, nrow=length(tipConsensusNodes), ncol=9, dimnames=list(tipConsensusNodes, c('annot taxa', 'diff funDist', 'diff avgPwys', 'red sig', 'red ns', 'spe sig', 'spe ns', 'mis sig', 'mis ns')))
resC = matrix(NA, nrow=length(tipConsensusNodes), ncol=length(goodPwys), dimnames=list(tipConsensusNodes, goodPwys))
resD = matrix(NA, nrow=length(tipConsensusNodes), ncol=length(goodPwys), dimnames=list(tipConsensusNodes, goodPwys))
spe  = matrix(NA, nrow=length(tipConsensusNodes), ncol=length(goodPwys), dimnames=list(tipConsensusNodes, goodPwys))


for(real0 in tipConsensusNodes) {
    if(real0 %in% failed) { next }
    real = paste(getNodeAnnotTaxa(real0), collapse='-')
    nulls = allNulls[[real0]]

    pwyPrevRealA = getNodePwyPrevalence(real)
    pwyPrevNullA = simplify_all(transpose(lapply(nulls, getNodePwyPrevalence)))
    pwyPrevNullAvgA = sapply(pwyPrevNullA, mean)
    diffFromAvgA = pwyPrevRealA - pwyPrevNullAvgA

    pwyPrevReal = pwyPrevRealA[goodPwys]
    pwyPrevNull = pwyPrevNullA[goodPwys]
    pwyPrevNullAvg = sapply(pwyPrevNull, mean)
    diffFromAvg = pwyPrevReal - pwyPrevNullAvg

    diffFromAvg = pwyPrevReal - pwyPrevNullAvg
    # In null hypothesis significance testing, the p-value[note 1] is the probability of obtaining test results at least as extreme as the results actually observed, 
    pval_red = sapply(1:length(pwyPrevReal), FUN=function(x) sum(pwyPrevNull[[x]] >= pwyPrevReal[x])/length(pwyPrevNull[[x]]))
    pval_spe = sapply(1:length(pwyPrevReal), FUN=function(x) sum(pwyPrevNull[[x]] <= pwyPrevReal[x])/length(pwyPrevNull[[x]]))
    pval = pmin(pval_red, pval_spe)
    names(pval) = names(pwyPrevReal)
	
    nullAvgPwys = sapply(nulls, getAvgPwysPerNodeMember)
    realAvgPwys = getAvgPwysPerNodeMember(real)
    nullFunDist = 1 - sapply(nulls, avgFunCorr)
    realFunDist = 1 - avgFunCorr(real)

    print(c(realAvgPwys - mean(nullAvgPwys), sum(diffFromAvgA)/getNodeSize(real))) # note the use of diffFromAvgA, this is because getAvgPwysPerNodeMember uses all the pathways (not only goodPathways)
                                                                                   # shouldn't matter anymore since now I'm using all the pwys
    resC[real0,] = diffFromAvg/getNodeSize(real)
    resD[real0,] = pwyPrevReal

    CUTOFF = 0.05
    
    resB[real0, 'annot taxa'] = length(getNodeAnnotTaxa(real0))
    resB[real0, 'diff funDist'] = realFunDist - mean(nullFunDist)
    resB[real0, 'diff avgPwys'] = realAvgPwys - mean(nullAvgPwys) # esto no sumará exactamente lo mismo que rd sig + red ns + spe sig + spe ns + mis sig + mis ns (aquí usamos todos los pwys, en las otras goodPwys)
    resB[real0, 'red sig'] = sum(diffFromAvg[names(pwyPrevReal)[pwyPrevReal > pwyPrevNullAvg & pval < CUTOFF]])
    resB[real0, 'red ns' ] = sum(diffFromAvg[names(pwyPrevReal)[pwyPrevReal > pwyPrevNullAvg & pval >= CUTOFF]])
    resB[real0, 'spe sig'] = sum(diffFromAvg[names(pwyPrevReal)[pwyPrevReal < pwyPrevNullAvg & pwyPrevReal > 0 & pval < CUTOFF]])
    resB[real0, 'spe ns' ] = sum(diffFromAvg[names(pwyPrevReal)[pwyPrevReal < pwyPrevNullAvg & pwyPrevReal > 0 & pval >= CUTOFF]])
    resB[real0, 'mis sig'] = sum(diffFromAvg[names(pwyPrevReal)[pwyPrevReal == 0 & pval < CUTOFF]])
    resB[real0, 'mis ns' ] = sum(diffFromAvg[names(pwyPrevReal)[pwyPrevReal == 0 & pval >= CUTOFF]])

    for(pwy in names(pval)) {
        if(pval[pwy] < CUTOFF) {
            if(diffFromAvg[pwy]<0) {
                if(pwyPrevReal[pwy] > 0) { res[real0, pwy] = 1
                } else { res[real0, pwy] = 2 }
            } else { res[real0, pwy] = -1 }
        } else { res[real0, pwy] = 0 }
    }
}

res2 = res[!rownames(res) %in% failed,]

resB2 = resB[!rownames(resB) %in% failed,]

resC = resC[!rownames(resC) %in% failed,]
resD = resD[!rownames(resD) %in% failed,]
res2 = res2[,colSums(res2 != 0) > 0]

envcols = c(`host-associated` = '#990099', `soils` = '#33CC00', `marine-water` = '#0000FF', `oil` = '#666600', `hypothermal-nonmarine` = '#CCFFFF',
            `thermal` = '#FF0000', `freshwater` = '#00CCFF', `marine-sediment` = '#009999', `saline-hypersaline` = '#009966', `water-treatment` = '#CCCCFF')



# heatmap
res2F = res2[,colSums(res2 != 0)>=10]
res2F = res2F[rowSums(res2F != 0) > 0,]
# remove some pwys engulfed in others
res2F = res2F[,!colnames(res2F) %in% c('Guanosine Deoxyribonucleotides <i>De Novo</i> Biosynthesis',
                                       'Adenosine Deoxyribonucleotides <i>De Novo</i> Biosynthesis',
                                       'UDP-<i>N</i>-acetyl-D-glucosamine biosynthesis I')]

library(gplots)
library(dendextend)
library(svglite)

svglite('UgoHeatmap.phyloEnv.svg', height=20, width=50)
# custom dendrograms so we control line width
ddR = set(as.dendrogram(hclust(dist(res2F))), 'branches_lwd', 3)
ddC = set(as.dendrogram(hclust(dist(t(res2F)))), 'branches_lwd', 3)
rcols = envcols[as.character(nodes[rownames(res2F), 'sourceNetwork'])]
rcols[is.na(rcols)] = 'darkolivegreen1' # multi env
rownames(res2F) = sapply(rownames(res2F), FUN=function(x) paste(strwrap(gsub('-', ' ', x), width=100), collapse='\n'))

heatmap.2(res2F, trace='none', Rowv = ddR, Colv = ddC, na.color='white', col=colorRampPalette(c('#f9766e','#f9766e','#f9766e','#f9766e','white','#01ba38', '#01ba38','#01ba38','#619dff')),
          margins=c(80, 50), lwid=c(1,10), lhei=c(1,10), cexRow=0.5, cexCol=2, RowSideColors = rcols, colsep = 1:1000,  sepcolor = 'black')
dev.off()
###


resB2 = cbind(resB2,
              `red all` = resB2[,'red sig'] + resB2[,'red ns'],
              `spe all` = resB2[,'spe sig'] + resB2[,'spe ns'],
              `mis all` = resB2[,'mis sig'] + resB2[,'mis ns']
              )

x = c(resB2[,'diff avgPwys'], resB2[,'red all']/resB2[,'annot taxa'], resB2[,'spe all']/resB2[,'annot taxa'], resB2[,'mis all']/resB2[,'annot taxa'])
z = factor(c(rep('All', nrow(resB2)), rep('Redundant', nrow(resB2)), rep('Specific', nrow(resB2)), rep('Missing', nrow(resB2))), levels=c('All', 'Redundant', 'Specific', 'Missing'))
s = rep(resB2[,'annot taxa'], 4)

THRES = 5
xs = x[s<THRES]
zs = z[s<THRES]
ss = s[s<THRES]

xl = x[s>=THRES]
zl = z[s>=THRES]
sl = s[s>=THRES]

par(mfrow=c(2,2), mar=c(5.1, 7.1, 4.1, 2.1), cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
p = wilcox.test(xs[zs=='All'],xl[zl=='All'])$p.value
boxplot(xs[zs=='All'],xl[zl=='All'], names=c('Small (<5 members)', 'Large'),
        ylab='Average pathways per genome\n(difference with null model)', main='All pathways', col='grey40', sub=sprintf('p = %s', round(p, 4)))

p = wilcox.test(xs[zs=='Redundant'],xl[zl=='Redundant'])$p.value
boxplot(xs[zs=='Redundant'],xl[zl=='Redundant'],
        names=c('Small (<5 members)', 'Large'), ylab='Average pathways per genome\n(difference with null model)', main='Redundant pathways', col='#f9766e', sub=sprintf('p = %s', round(p, 4)))

p = wilcox.test(xs[zs=='Specific'],xl[zl=='Specific'])$p.value
boxplot(xs[zs=='Specific'],xl[zl=='Specific'],
        names=c('Small (<5 members)', 'Large'), ylab='Average pathways per genome\n(difference with null model)', main='Specific pathways', col='#01ba38', sub=sprintf('p = %s', round(p, 4)))

p = wilcox.test(xs[zs=='Missing'],xl[zl=='Missing'])$p.value
boxplot(xs[zs=='Missing'],xl[zl=='Missing'],
        names=c('Small (<5 members)', 'Large'), ylab='Average pathways per genome\n(difference with null model)', main='Missing pathways', col='#619dff', sub=sprintf('p = %s', round(p, 4)))

dev.copy(png, "smallVSlarge.png")
dev.off()




#### aminoacids

#Taken from http://www.pnas.org/content/99/6/3695.full
costs = c('L-alanine Biosynthesis' = 11.7,
          'L-arginine Biosynthesis' = 27.3,
          'L-asparagine Biosynthesis' = 14.7,
          'L-aspartate Biosynthesis' = 12.7,
          'L-cysteine Biosynthesis' = 24.7,
          'L-glutamate Biosynthesis' = 15.3,
          'L-glutamine Biosynthesis' = 16.3,
          'Glycine Biosynthesis' = 11.7,
          'L-histidine Biosynthesis' = 38.3,
          'L-isoleucine Biosynthesis' = 32.3,
          'L-leucine Biosynthesis' = 27.3,
          'L-lysine Biosynthesis' = 30.3,
          'L-methionine <i>De Novo</i> Biosynthesis' = 34.3,
          'L-phenylalanine Biosynthesis' = 52.0,
          'L-proline Biosynthesis' = 20.3,
          'L-serine Biosynthesis' = 11.7,
          'L-threonine Biosynthesis' = 18.7,
          'L-typtophan Biosynthesis' = 74.3, #This typo comes from metacyc... live with it.
          'L-tyrosine Biosynthesis' = 50.0,
          'L-valine Biosynthesis' = 23.3
         )




a = sapply(rownames(resD), FUN=function(x) length(getNodeAnnotTaxa(x)))

aux=1-colSums(resD)/sum(a) # SPECIFICITY/AUXOTROFY?, times we see it in a genome / genera involved

aux.aa = aux[names(costs)]

library(ggplot2)
gt = data.frame(row.names=names(costs))
gt[,'cost'] = costs
gt[,'aux'] = 100*aux.aa
gt[,'auxNoOutliers'] = 100*aux.aa
gt[c('L-typtophan Biosynthesis'), 'auxNoOutliers'] = NA

summary(lm(gt[,'aux']~gt[,'cost'])) # Adj R2=0.1528, p = 0.0497
summary(lm(gt[,'auxNoOutliers']~gt[,'cost'])) # Adj R2=0.485 p = 0.0005558

ggplot(gt, aes(x=cost, y=aux)) + 
  ggtitle('Aminoacid biosynthesis: average specificity vs biosynthetic cost') + xlab('Biosynthetic cost (high-energy phosphate bonds)') + ylab('Average specificity (%)') + 
  geom_point(shape=5, size=4) + 
  stat_smooth(geom='line', method='lm', col='red', linetype='dashed', alpha=1) +
  geom_smooth(method='lm', aes(x=cost, y=auxNoOutliers), col='royalblue1', fullrange=TRUE) +
  theme_light(base_size = 24)
ggsave('aaCost.png', height=4000, width=6000, unit='px')





# Summarize

annotTaxa = unique(unlist(lapply(rownames(res2), getNodeAnnotTaxa)))
res3 = data.frame(inGenomes = colSums(pathways[annotTaxa, colnames(res2)]),
                  redundant = sapply(colnames(res2), FUN=function(x) sum(res2[,x]==-1)),
                  specific  = sapply(colnames(res2), FUN=function(x) sum(res2[,x]== 1)),
                  missing   = sapply(colnames(res2), FUN=function(x) sum(res2[,x]== 2)),
                  nosig     = sapply(colnames(res2), FUN=function(x) sum(res2[,x]== 0)),
                  row.names = colnames(res2))


res3 = res3[order(res3$redundant),]

parseNodePathways = function(node) {
    prev = c()
    pwys = unlist(strsplit(as.character(nodes[node,'pathways']), split='|', fixed=T))
    goodPwys = c()
    for(pwy in pwys) {
        if(!grepl(' - ', pwy, fixed=T)) {
            prev = c(prev, pwy) # this is part of a larger annotation e.g. AcetCoaBios|CarbDegr - pyruvate decarboxylation to acetyl CoA
        }else{###############
            prev = c(prev, pwy) # yeah this is duplicated but I'd rather do it like this and explicit the special case above
            goodPwys = c(goodPwys, paste(prev, collapse='|'))
            prev = c()
        }
    }
    return(goodPwys)
}
            

longNames = unique(unlist(lapply(rownames(nodes), parseNodePathways)))
shortNames = sapply(strsplit(longNames, split=' - '), FUN = function(x) x[2])
names(shortNames) = longNames
names(longNames) = shortNames

res3$cat = sapply(strsplit(longNames[rownames(res3)], split=' - '), FUN = function(x) x[1])

# CarbDegr can mean Carbohydrate biosynthesis or Carboxylate biosynthesis. Manually correct here after having checked each pathway in MetaCyc
carboxylate = c('3-oxoadipate degradation', '2-oxopentenoate degradation', 'acetate formation from acetyl-CoA I', 'pyruvate decarboxylation to acetyl CoA',
                'L-malate degradation II', 'acetate conversion to acetyl-CoA', 'D-galactarate degradation I', '2-methylcitrate cycle I',
                'D-gluconate degradation', 'citrate degradation', '2-methylcitrate cycle II', 'oxalate degradation II',
                'glycolate and glyoxylate degradation II', 'D-fructuronate degradation', 'propionyl CoA degradation', 'itaconate degradation',
                'acetate formation from acetyl-CoA III (succinate)', 'cyclohexane-1-carboxylate degradation (anaerobic)', 'acrylate degradation',
                '2-amino-3-carboxymuconate semialdehyde degradation to 2-oxopentenoate', 'oxalate degradation V', '<i>D</i>-glucarate degradation I',
                'L-idonate degradation', 'glutaryl-CoA degradation', 'D-galactarate degradation II', 'D-malate degradation',
                'mevalonate degradation', 'propanoyl-CoA degradation II', 'D-galactonate degradation', 'D-glucosaminate degradation',
                '&beta;-D-glucuronide and D-glucuronate degradation', 'L-ascorbate degradation II (bacterial, aerobic)', 'L-ascorbate degradation I (bacterial, anaerobic)',
                'D-glucarate degradation II', 'D-galacturonate degradation I', 'nitrilotriacetate degradation', 'malonate degradation I (biotin-independent)',
                'acetate formation from acetyl-CoA II')
res3[carboxylate,'cat'] = gsub('CarbDegr', 'CarboxDegr', res3[carboxylate,'cat'], fixed=T)
                
# manually fix ferulate (no es ni carbohydrate ni carboxylate)
res3['ferulate degradation', 'cat'] = 'SecoMetaBios|AromCompDegr'
res3['Bifidobacterium shunt', 'cat'] = 'Ferm|CarbDegr|CarboxDegr' # es de las dos


categories = sort(unique(unlist(strsplit(res3$cat, '|', fixed=T))))

res3_summ = matrix(0, nrow=length(categories), ncol=4, dimnames=list(categories, c('Redundant', 'Specific', 'Missing', 'NoSig')))
for(i in 1:nrow(res3)) {
    cats = sort(unique(unlist(strsplit(res3[i,'cat'], '|', fixed=T)))) # meow
    for(cat in cats) { res3_summ[cat,] = res3_summ[cat,] + unlist(res3[i, c('redundant', 'specific', 'missing', 'nosig')]) }
}

res3_summ = res3_summ[,c('Redundant', 'Specific', 'Missing')]
top=names(sort(rowSums(res3_summ), decreasing=T))[1:15]
res3_summ = res3_summ[top,]
res3_summ = res3_summ[order(res3_summ[,'Specific'], decreasing=T),]

equivs = c('CarbDegr' = 'Carbohydrate Degradation',     'CarboxDegr' = 'Carboxylate Degradation', 'InorNutrMeta' = 'Inorganic Nutrient Metabolism', 'Resp' = 'Respiration',
           'CofaProsGrou' = 'Cofactors, Prosthetic Groups,\nElectron Carriers Biosynthesis', 'SecoMetaBios' = 'Secondary Metabolite Biosynthesis', 'FattAcidLipi' = 'Fatty Acid and Lipid Biosynthesis',
           'AromCompDegr' = 'Aromatic Compound Degradation', 'CarbBios' = 'Carbohydrate Biosynthesis', 'Ferm' = 'Fermentation', 'C1CompUtil' = 'C1 Compound Utilization and Assimilation',
           'ElecTran' = 'Electron Transfer',     'CellStruBios' = 'Cell Structures Biosynthesis', 'AminAcidDegr' = 'Amino Acids Degradation',
           'NuclNuclBios' = 'Nucleosides and Nucleotides Biosynthesis', 'AminAcidBios' = 'Amino Acids Biosynthesis')

rownames(res3_summ) = equivs[rownames(res3_summ)]

library(reshape2)
res3_summM = melt(res3_summ)



ggplot(data=res3_summM, aes(fill=Var2, x=Var1, y=value)) +
geom_bar(position='dodge', stat='identity') +
coord_flip() +
theme_light(base_size = 24) + 
ylab('Times redundant/specific/missing in the consensus network') + xlab('MetaCyc class') + labs(fill='Category')

ggsave('metabolicGroups.png', height=4000, width=4000, unit='px')


