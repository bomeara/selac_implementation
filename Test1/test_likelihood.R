
library(selac)
library(ape)


######################################################################################################################################
######################################################################################################################################
### Test 1 -- SelAC with GTR
######################################################################################################################################
######################################################################################################################################

set.seed(4)
tree <- read.tree("rokasYeast.tre")
phy <- drop.tip(tree, "Calb")
yeast.gene <- read.dna("gene1Yeast.fasta", format="fasta")
yeast.gene <- as.list(as.matrix(cbind(yeast.gene))[1:7,])
write.nexus.data(yeast.gene, file="dna.nex")
write.nexus(phy, file="tree.nex")

chars <- selac:::DNAbinToCodonNumeric(yeast.gene)
codon.data <- chars[phy$tip.label,]
aa.data <- selac:::ConvertCodonNumericDataToAAData(codon.data, numcode=1)
aa.optim <- apply(aa.data[, -1], 2, selac:::GetMaxName) #starting values for all, final values for majrule
aa.optim.full.list <- aa.optim
codon.freq.by.aa <- selac:::GetCodonFreqsByAA(codon.data[,-1], aa.optim, numcode=1)
codon.freq.by.gene <- selac:::GetCodonFreqsByGene(codon.data[,-1])
aa.optim.frame.to.add <- matrix(c("optimal", aa.optim), 1, dim(codon.data)[2])
colnames(aa.optim.frame.to.add) <- colnames(codon.data)
codon.data <- rbind(codon.data, aa.optim.frame.to.add)
codon.data <- selac:::SitePattern(codon.data, includes.optimal.aa=TRUE)
aa.optim = codon.data$optimal.aa
codon.index.matrix = selac:::CreateCodonMutationMatrixIndex()
selac.gtr <- selac:::GetLikelihoodSAC_CodonForManyCharGivenAllParams(log(c(4*4e-7*.5*5e6, 1.829272, 0.101799, .25, .25, .25, rep(1,5))), codon.data, phy, aa.optim_array=aa.optim, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, numcode=1, diploid=TRUE, aa.properties=NULL, volume.fixed.value=0.0003990333, nuc.model="GTR", codon.index.matrix, include.gamma=FALSE, ncats=4, k.levels=0, logspace=TRUE, verbose=FALSE, n.cores.by.gene.by.site=1)
comparison <- identical(round(selac.gtr, 3), -7063.585)
print(comparison)
######################################################################################################################################
######################################################################################################################################



######################################################################################################################################
######################################################################################################################################
### Test 2 -- SelAC with UNREST, but should lead to identical likelihood as Test 1
######################################################################################################################################
######################################################################################################################################
set.seed(4)
tree <- read.tree("rokasYeast.tre")
phy <- drop.tip(tree, "Calb")
yeast.gene <- read.dna("gene1Yeast.fasta", format="fasta")
yeast.gene <- as.list(as.matrix(cbind(yeast.gene))[1:7,])
chars <- selac:::DNAbinToCodonNumeric(yeast.gene)
codon.data <- chars[phy$tip.label,]
aa.data <- selac:::ConvertCodonNumericDataToAAData(codon.data, numcode=1)
aa.optim <- apply(aa.data[, -1], 2, selac:::GetMaxName) #starting values for all, final values for majrule
aa.optim.full.list <- aa.optim
codon.freq.by.aa <- selac:::GetCodonFreqsByAA(codon.data[,-1], aa.optim, numcode=1)
codon.freq.by.gene <- selac:::GetCodonFreqsByGene(codon.data[,-1])
aa.optim.frame.to.add <- matrix(c("optimal", aa.optim), 1, dim(codon.data)[2])
colnames(aa.optim.frame.to.add) <- colnames(codon.data)
codon.data <- rbind(codon.data, aa.optim.frame.to.add)
codon.data <- selac:::SitePattern(codon.data, includes.optimal.aa=TRUE)
aa.optim = codon.data$optimal.aa
codon.index.matrix = selac:::CreateCodonMutationMatrixIndex()
selac.unrest <- selac:::GetLikelihoodSAC_CodonForManyCharGivenAllParams(log(c(4*4e-7*.5*5e6, 1.829272, 0.101799, rep(1,11))), codon.data, phy, aa.optim_array=aa.optim, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, numcode=1, diploid=TRUE, aa.properties=NULL, volume.fixed.value=0.0003990333, nuc.model="UNREST", codon.index.matrix, include.gamma=FALSE, ncats=4, k.levels=0, logspace=TRUE, verbose=FALSE, n.cores.by.gene.by.site=1)
comparison <- identical(round(selac.unrest, 3), -7063.585)
print(comparison)
######################################################################################################################################
######################################################################################################################################



######################################################################################################################################
######################################################################################################################################
### Test 3 -- SelAC with UNREST, with among site variation set to "median"
######################################################################################################################################
######################################################################################################################################
set.seed(4)
tree <- read.tree("rokasYeast.tre")
phy <- drop.tip(tree, "Calb")
yeast.gene <- read.dna("gene1Yeast.fasta", format="fasta")
yeast.gene <- as.list(as.matrix(cbind(yeast.gene))[1:7,])
chars <- selac:::DNAbinToCodonNumeric(yeast.gene)
codon.data <- chars[phy$tip.label,]
aa.data <- selac:::ConvertCodonNumericDataToAAData(codon.data, numcode=1)
aa.optim <- apply(aa.data[, -1], 2, selac:::GetMaxName) #starting values for all, final values for majrule
aa.optim.full.list <- aa.optim
codon.freq.by.aa <- selac:::GetCodonFreqsByAA(codon.data[,-1], aa.optim, numcode=1)
codon.freq.by.gene <- selac:::GetCodonFreqsByGene(codon.data[,-1])
aa.optim.frame.to.add <- matrix(c("optimal", aa.optim), 1, dim(codon.data)[2])
colnames(aa.optim.frame.to.add) <- colnames(codon.data)
codon.data <- rbind(codon.data, aa.optim.frame.to.add)
codon.data <- selac:::SitePattern(codon.data, includes.optimal.aa=TRUE)
aa.optim = codon.data$optimal.aa
codon.index.matrix <- selac:::CreateCodonMutationMatrixIndex()
selac_gamma <- selac:::GetLikelihoodSAC_CodonForManyCharGivenAllParams(log(c(4*4e-7*.5*5e6, 1.829272, 0.101799, rep(1,11), 5)), codon.data, phy, aa.optim_array=aa.optim, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, numcode=1, diploid=TRUE, aa.properties=NULL, volume.fixed.value=0.0003990333, nuc.model="UNREST", codon.index.matrix, include.gamma=TRUE, gamma.type="median", ncats=4, k.levels=0, logspace=TRUE, verbose=FALSE)
comparison <- identical(round(selac_gamma, 3), -7006.621)
print(comparison)
######################################################################################################################################
######################################################################################################################################



######################################################################################################################################
######################################################################################################################################
### Test 4 -- SelAC with UNREST, with among site variation set to "quadrature"
######################################################################################################################################
######################################################################################################################################
set.seed(4)
tree <- read.tree("rokasYeast.tre")
phy <- drop.tip(tree, "Calb")
yeast.gene <- read.dna("gene1Yeast.fasta", format="fasta")
yeast.gene <- as.list(as.matrix(cbind(yeast.gene))[1:7,])
chars <- selac:::DNAbinToCodonNumeric(yeast.gene)
codon.data <- chars[phy$tip.label,]
aa.data <- selac:::ConvertCodonNumericDataToAAData(codon.data, numcode=1)
aa.optim <- apply(aa.data[, -1], 2, selac:::GetMaxName) #starting values for all, final values for majrule
aa.optim.full.list <- aa.optim
codon.freq.by.aa <- selac:::GetCodonFreqsByAA(codon.data[,-1], aa.optim, numcode=1)
codon.freq.by.gene <- selac:::GetCodonFreqsByGene(codon.data[,-1])
aa.optim.frame.to.add <- matrix(c("optimal", aa.optim), 1, dim(codon.data)[2])
colnames(aa.optim.frame.to.add) <- colnames(codon.data)
codon.data <- rbind(codon.data, aa.optim.frame.to.add)
codon.data <- selac:::SitePattern(codon.data, includes.optimal.aa=TRUE)
aa.optim = codon.data$optimal.aa
codon.index.matrix = selac:::CreateCodonMutationMatrixIndex()
selac_gamma_quad <- selac:::GetLikelihoodSAC_CodonForManyCharGivenAllParams(log(c(4*4e-7*.5*5e6, 1.829272, 0.101799, rep(1,11), 5)), codon.data, phy, aa.optim_array=aa.optim, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, numcode=1, diploid=TRUE, aa.properties=NULL, volume.fixed.value=0.0003990333, nuc.model="UNREST", codon.index.matrix, include.gamma=TRUE, gamma.type="quadrature", ncats=4, k.levels=0, logspace=TRUE, verbose=FALSE, n.cores.by.gene.by.site=1)
comparison <- identical(round(selac_gamma_quad, 3), -6998.186)
print(comparison)
######################################################################################################################################
######################################################################################################################################





######################################################################################################################################
######################################################################################################################################
### Test 5 -- Goldman Yang 94 implementation
######################################################################################################################################
######################################################################################################################################
tree <- read.tree("rokasYeast.tre")
phy <- drop.tip(tree, "Calb")
yeast.gene <- read.dna("gene1Yeast.fasta", format="fasta")
yeast.gene <- as.list(as.matrix(cbind(yeast.gene))[1:7,])
chars <- selac:::DNAbinToCodonNumeric(yeast.gene)
codon.data <- chars[phy$tip.label,]
codon.data <- selac:::SitePattern(codon.data)
gy94 <- selac:::GetLikelihoodGY94_YN98_CodonForManyCharGivenAllParams(log(c(1,1)), codon.data, phy, numcode=1, logspace=TRUE, verbose=FALSE, n.cores.by.gene.by.site=1)
comparison <- identical(round(gy94,3), -7826.123)
print(comparison)

######################################################################################################################################
######################################################################################################################################




results <- data.frame(model=c("SelAC with GTR", "SelAC with UNREST", "SelAC with UNREST, with among site variation set to median", "SelAC with UNREST, with among site variation set to quadrature", "GY94"), LL=c(selac.gtr, selac.unrest, selac_gamma, selac_gamma_quad, gy94), stringsAsFactors=FALSE)
write.csv(results, file="selac.csv", row.names=FALSE)
