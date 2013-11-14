library(picante)
library(ggplot2)
library(foreach)
library(reshape)
library(RJSONIO)
library(igraph0)
library(multtest)

library(phyloseq)
library(XML)
library(phyloch)
library(ape)

#setwd("/Users/pipgriffinold/Documents/Current Projects/454 Results/Gene Trees 6 runs 150M burnin removed resampled 100000 therefore 3000 final trees and mcc made")
setwd("/Applications/BEAST Suite/Runs Sept 2013/Combined Tree Files")



#Before reading in trees, tidy up the names in TextWrangler
#replace species codes with full names
#replace TAS with Tas
#don't remove underscores, this stuffs up the tip label detection
#and plot.phylo automatically plots them as spaces anyway

#################################
#         SPECIES TREE          #
#################################

sp_tree<-read.beast("All_except_Run_7_species.tree", digits=7)
lad_sp<-ladderize(sp_tree, right=FALSE)

#change branch lengths and HPD intervals from substitutions to time scale
factor<-13/0.0724261728

converted <- compute.brlen(lad_sp, lad_sp$edge.length*factor)
converted$"height_95%_HPD" <- converted$"height_95%_HPD"*factor
converted$"height_95%_HPD_MIN" <- converted$"height_95%_HPD_MIN"*factor
converted$"height_95%_HPD_MAX" <- converted$"height_95%_HPD_MAX"*factor

#Change tip labels from codes into names
converted$tip.label<-c("P. annua", "P. clivicola", "P. costiniana", "P. ensiformis", "P. fawcettiae", "P. gunnii", "P. helmsii", "P. hiemata", "P. hothamensis", "P. phillipsiana", "P. sieberiana")

pdf(file="StarBEAST Tree with HPD Intervals.pdf", width = 12, height = 7, useDingbats=FALSE)
plotdimensions<-plot(converted, edge.color=0, tip.color=0)
axisChrono(side=1, unit = "MY")
abline(v=13-0.1, col="red", lty=2)
abline(v=13-1, col="red", lty=2)
#abline(v=13-10, col="red", lty=2)
HPDbars(converted, col="lightgrey", lwd=4)
plot.phylo.upon(converted)
node.support(converted$posterior, digits=2)
node.support(converted$posterior, mode="dots", cutoff=0.8, col="black")
#node.support(converted$posterior, mode="dots", cutoff=0.95, col="blue")
node.support(converted$posterior, mode="dots", cutoff=0.99, col="red")
dev.off()


###########################################################

# Now want to do something similar for the bpp tree values

###########################################################

bpp_factor <- 13/0.04697

#The edge lengths (tau) here are only for one tree
bpp_tree<-read.tree(text="((((((((Pen: 0.000295, Phe: 0.000295): 0.000053, Pho: 0.000349): 0.000049, Pph: 0.000398): 0.000027, Phi: 0.000425): 0.000012, (Pcl: 0.000408, Psi: 0.000408): 0.000029): 0.000005, (Pco: 0.000042, Pfa: 0.000042): 0.000401): 0.000968, Pgu: 0.001411): 0.045782, Pan: 0.047193);")

#replace the edge lengths
setwd("/Applications/bpp2.2/Poa/Fixed Species Tree With Pan")
tau_values<-read.csv("bpp_tau_values.csv")

bpp_tree$edge.length<-tau_values$edge.length.mean[1:20]*bpp_factor
bpp_min<-tau_values$X95..HPD.lower[12:21]*bpp_factor
bpp_max<-tau_values$X95..HPD.upper[12:21]*bpp_factor
bpp_tree$"height_95%_HPD_MIN" <- bpp_min
bpp_tree$"height_95%_HPD_MAX" <- bpp_max
bpp_tree$tip.label<-c("P. ensiformis", "P. helmsii", "P. hothamensis", "P. phillipsiana", "P. hiemata", "P. clivicola", "P. sieberiana", "P. costiniana", "P. fawcettiae", "P. gunnii", "P. annua")
lad_bpp<-ladderize(bpp_tree, right=FALSE)

pdf(file="BPP Tree with HPD Intervals.pdf", width = 12, height = 7, useDingbats=FALSE)
plotdimensions<-plot(lad_bpp, edge.color=0, tip.color=0)
axisChrono(side=1, unit = "Ma")
abline(v=13-0.1, col="red", lty=2)
abline(v=13-1, col="red", lty=2)
HPDbars(lad_bpp, col="lightgrey", lwd=4)
plot.phylo.upon(lad_bpp)
dev.off()

#################################
#      CHLOROPLAST TREE         #
#################################

setwd("/Applications/BEAST Suite/Runs Sept 2013/Combined Tree Files")

cpdata<-read.beast.table("cp_9runs_downsampled_7000_mcc.tree", digits=2)

cptree<-read.beast("cp_9runs_downsampled_7000_mcc.tree")
cptree<-ladderize(cptree, right=FALSE)

plot.phylo(cptree, cex=0.3, root.edge=TRUE, x.lim=0.1, label.offset=0.003)

lab<-cptree$tip.label

pannsw<-grep("annua", lab)

pclnsw<-grep("clivicola_NSW", lab)
pclact<-grep("clivicola_ACT", lab)
pclvic<-grep("clivicola_Vic", lab)
pcltas<-grep("clivicola_Tas", lab)

pconsw<-grep("costiniana_NSW", lab)
pcotas<-grep("costiniana_Tas", lab)
pcovic<-grep("costiniana_Vic", lab)

pennsw<-grep("ensiformis_NSW", lab)
penact<-grep("ensiformis_ACT", lab)
penvic<-grep("ensiformis_Vic", lab)

pfansw<-grep("fawcettiae_NSW", lab)
pfaact<-grep("fawcettiae_ACT", lab)
pfavic<-grep("fawcettiae_Vic", lab)
pfatas<-grep("fawcettiae_Tas", lab)

pgutas<-grep("gunnii_Tas", lab)

phensw<-grep("helmsii_NSW", lab)
pheact<-grep("helmsii_ACT", lab)

phinsw<-grep("hiemata_NSW", lab)
phiact<-grep("hiemata_ACT", lab)
phivic<-grep("hiemata_Vic", lab)
phitas<-grep("hiemata_Tas", lab)

phovic<-grep("hothamensis_Vic", lab)

pphnsw<-grep("phillipsiana_NSW", lab)
pphact<-grep("phillipsiana_ACT", lab)
pphvic<-grep("phillipsiana_Vic", lab)
pphtas<-grep("phillipsiana_Tas", lab)

psinsw<-grep("sieberiana_NSW", lab)
psiact<-grep("sieberiana_ACT", lab)
psivic<-grep("sieberiana_Vic", lab)
psitas<-grep("sieberiana_Tas", lab)

pdf("CP MCC Tree With Support And Labels.pdf", width=10, height=8, useDingbats=FALSE)
plot.phylo(cptree, cex=0.3, label.offset=0.002)

add.scale.bar(x=0, y=15, length=0.05, cex=0.5)
node.support(as.data.frame(cpdata)$posterior, cutoff=0.8, pos="right", mode="dots", col="black")
node.support(as.data.frame(cpdata)$posterior, cutoff=0.95, pos="right", mode="dots", col="blue")
node.support(as.data.frame(cpdata)$posterior, cutoff=0.99, pos="right", mode="dots", col="red")

points(x=rep(0.001, times=11), y=seq(70, 25, length.out=11), pch=symbolsnsw, col=cols, cex=0.4)
points(x=rep(0.006, times=11), y=seq(70, 25, length.out=11), pch=symbolsvic, col=cols, cex=0.4)
points(x=rep(0.012, times=11), y=seq(70, 25, length.out=11), pch=symbolstas, col=cols, cex=0.4)
text(x=rep(0.015, times=11), y=seq(70, 25, length.out=11), labels=names, adj=0, cex=0.4, font=3)
text(x=c(0.001, 0.006, 0.012), y=rep(75, times=11), labels=c("NSW", "Vic", "Tas"), cex=0.4)

appendpannsw<-MYappend2tiplabel(cptree, tips=pannsw, col="black", pch=15, align=TRUE)
appendpclnsw<-MYappend2tiplabel(cptree, tips=pclnsw, col="mediumpurple", pch=15, align=TRUE)
appendpclact<-MYappend2tiplabel(cptree, tips=pclact, col="mediumpurple", pch=15, align=TRUE)
appendpclvic<-MYappend2tiplabel(cptree, tips=pclvic, col="mediumpurple", pch=19, align=TRUE)
appendpcoact<-MYappend2tiplabel(cptree, tips=pcoact, col="royalblue4", pch=15, align=TRUE)
appendpconsw<-MYappend2tiplabel(cptree, tips=pconsw, col="royalblue4", pch=15, align=TRUE)
appendpcovic<-MYappend2tiplabel(cptree, tips=pcovic, col="royalblue4", pch=19, align=TRUE)
appendpcotas<-MYappend2tiplabel(cptree, tips=pcotas, col="royalblue4", pch=17, align=TRUE)
appendpenact<-MYappend2tiplabel(cptree, tips=penact, col="black", bg="grey50", pch=22, align=TRUE)
appendpennsw<-MYappend2tiplabel(cptree, tips=pennsw, col="black", bg="grey50", pch=22, align=TRUE)
appendpenvic<-MYappend2tiplabel(cptree, tips=penvic, col="black", bg="grey50", pch=21, align=TRUE)
appendpfatas<-MYappend2tiplabel(cptree, tips=pfatas, col="springgreen4", pch=17, align=TRUE)
appendpfaact<-MYappend2tiplabel(cptree, tips=pfaact, col="springgreen4", pch=15, align=TRUE)
appendpfansw<-MYappend2tiplabel(cptree, tips=pfansw, col="springgreen4", pch=15, align=TRUE)
appendpfavic<-MYappend2tiplabel(cptree, tips=pfavic, col="springgreen4", pch=19, align=TRUE)
appendpgutas<-MYappend2tiplabel(cptree, tips=pgutas, col="black", pch=17, align=TRUE)
appendphensw<-MYappend2tiplabel(cptree, tips=phensw, col="turquoise3", pch=15, align=TRUE)
appendpheact<-MYappend2tiplabel(cptree, tips=pheact, col="turquoise3", pch=15)
appendphitas<-MYappend2tiplabel(cptree, tips=phitas, col="red3", pch=17, align=TRUE)
appendphiact<-MYappend2tiplabel(cptree, tips=phiact, col="red3", pch=15, align=TRUE)
appendphinsw<-MYappend2tiplabel(cptree, tips=phinsw, col="red3", pch=15, align=TRUE)
appendphivic<-MYappend2tiplabel(cptree, tips=phivic, col="red3", pch=19, align=TRUE)
appendphovic<-MYappend2tiplabel(cptree, tips=phovic, col="orangered2", pch=19, align=TRUE)
appendpphact<-MYappend2tiplabel(cptree, tips=pphact, col="goldenrod1", pch=15, align=TRUE)
appendpphnsw<-MYappend2tiplabel(cptree, tips=pphnsw, col="goldenrod1", pch=15, align=TRUE)
appendpphvic<-MYappend2tiplabel(cptree, tips=pphvic, col="goldenrod1", pch=19, align=TRUE)
appendpsiact<-MYappend2tiplabel(cptree, tips=psiact, col="gray50", pch=15, align=TRUE)
appendpsinsw<-MYappend2tiplabel(cptree, tips=psinsw, col="gray50", pch=15, align=TRUE)
appendpsivic<-MYappend2tiplabel(cptree, tips=psivic, col="gray50", pch=19, align=TRUE)

dev.off()

#################################
#         CDO504A TREE          #
#################################
setwd("/Applications/BEAST Suite/Runs Sept 2013/Combined Tree Files")
cdo504Adata<-read.beast.table("CDO504A_9runs_downsampled_7000_mcc.tree", digits=2)

cdo504Atree<-read.beast("CDO504A_9runs_downsampled_7000_mcc.tree")

cdo504Atree<-ladderize(cdo504Atree, right=FALSE)

plot.phylo(cdo504Atree, cex=0.3, label.offset=0.01)

lab<-cdo504Atree$tip.label

pannsw<-grep("annua", lab)

pclnsw<-grep("clivicola_NSW", lab)
pclact<-grep("clivicola_ACT", lab)
pclvic<-grep("clivicola_Vic", lab)
pcltas<-grep("clivicola_Tas", lab)

pconsw<-grep("costiniana_NSW", lab)
pcotas<-grep("costiniana_Tas", lab)
pcovic<-grep("costiniana_Vic", lab)

pennsw<-grep("ensiformis_NSW", lab)
penact<-grep("ensiformis_ACT", lab)
penvic<-grep("ensiformis_Vic", lab)

pfansw<-grep("fawcettiae_NSW", lab)
pfaact<-grep("fawcettiae_ACT", lab)
pfavic<-grep("fawcettiae_Vic", lab)
pfatas<-grep("fawcettiae_Tas", lab)

pgutas<-grep("gunnii_Tas", lab)

phensw<-grep("helmsii_NSW", lab)
pheact<-grep("helmsii_ACT", lab)

phinsw<-grep("hiemata_NSW", lab)
phiact<-grep("hiemata_ACT", lab)
phivic<-grep("hiemata_Vic", lab)
phitas<-grep("hiemata_Tas", lab)

phovic<-grep("hothamensis_Vic", lab)

pphnsw<-grep("phillipsiana_NSW", lab)
pphact<-grep("phillipsiana_ACT", lab)
pphvic<-grep("phillipsiana_Vic", lab)
pphtas<-grep("phillipsiana_Tas", lab)

psinsw<-grep("sieberiana_NSW", lab)
psiact<-grep("sieberiana_ACT", lab)
psivic<-grep("sieberiana_Vic", lab)
psitas<-grep("sieberiana_Tas", lab)

pdf("CDO504A MCC Tree With Support And Labels.pdf", width=10, height=8, useDingbats=FALSE)

plot.phylo(cdo504Atree, cex=0.3, label.offset=0.02)

add.scale.bar(x=0, y=30, length=0.05, cex=0.5)
node.support(as.data.frame(cdo504Adata)$posterior, cutoff=0.8, pos="right", mode="dots", col="black")
node.support(as.data.frame(cdo504Adata)$posterior, cutoff=0.95, pos="right", mode="dots", col="blue")
node.support(as.data.frame(cdo504Adata)$posterior, cutoff=0.99, pos="right", mode="dots", col="red")

points(x=rep(0.001, times=11), y=seq(80, 50, length.out=11), pch=symbolsnsw, col=cols, cex=0.3)
points(x=rep(0.01, times=11), y=seq(80, 50, length.out=11), pch=symbolsvic, col=cols, cex=0.3)
points(x=rep(0.019, times=11), y=seq(80, 50, length.out=11), pch=symbolstas, col=cols, cex=0.3)
text(x=rep(0.024, times=11), y=seq(80, 50, length.out=11), labels=names, adj=0, cex=0.3, font=3)
text(x=c(0.001, 0.01, 0.019), y=rep(85, times=11), labels=c("NSW", "Vic", "Tas"), cex=0.3)

appendpannsw<-MYappend2tiplabel(cdo504Atree, tips=pannsw, col="black", pch=15, align=TRUE)
appendpclnsw<-MYappend2tiplabel(cdo504Atree, tips=pclnsw, col="mediumpurple", pch=15, align=TRUE)
appendpclact<-MYappend2tiplabel(cdo504Atree, tips=pclact, col="mediumpurple", pch=15, align=TRUE)
appendpclvic<-MYappend2tiplabel(cdo504Atree, tips=pclvic, col="mediumpurple", pch=19, align=TRUE)
appendpcoact<-MYappend2tiplabel(cdo504Atree, tips=pcoact, col="royalblue4", pch=15, align=TRUE)
appendpconsw<-MYappend2tiplabel(cdo504Atree, tips=pconsw, col="royalblue4", pch=15, align=TRUE)
appendpcovic<-MYappend2tiplabel(cdo504Atree, tips=pcovic, col="royalblue4", pch=19, align=TRUE)
appendpcotas<-MYappend2tiplabel(cdo504Atree, tips=pcotas, col="royalblue4", pch=17, align=TRUE)
appendpenact<-MYappend2tiplabel(cdo504Atree, tips=penact, col="black", bg="grey50", pch=22, align=TRUE)
appendpennsw<-MYappend2tiplabel(cdo504Atree, tips=pennsw, col="black", bg="grey50", pch=22, align=TRUE)
appendpenvic<-MYappend2tiplabel(cdo504Atree, tips=penvic, col="black", bg="grey50", pch=21, align=TRUE)
appendpfatas<-MYappend2tiplabel(cdo504Atree, tips=pfatas, col="springgreen4", pch=17, align=TRUE)
appendpfaact<-MYappend2tiplabel(cdo504Atree, tips=pfaact, col="springgreen4", pch=15, align=TRUE)
appendpfansw<-MYappend2tiplabel(cdo504Atree, tips=pfansw, col="springgreen4", pch=15, align=TRUE)
appendpfavic<-MYappend2tiplabel(cdo504Atree, tips=pfavic, col="springgreen4", pch=19, align=TRUE)
appendpgutas<-MYappend2tiplabel(cdo504Atree, tips=pgutas, col="black", pch=17, align=TRUE)
appendphensw<-MYappend2tiplabel(cdo504Atree, tips=phensw, col="turquoise3", pch=15, align=TRUE)
appendpheact<-MYappend2tiplabel(cdo504Atree, tips=pheact, col="turquoise3", pch=15)
appendphitas<-MYappend2tiplabel(cdo504Atree, tips=phitas, col="red3", pch=17, align=TRUE)
appendphiact<-MYappend2tiplabel(cdo504Atree, tips=phiact, col="red3", pch=15, align=TRUE)
appendphinsw<-MYappend2tiplabel(cdo504Atree, tips=phinsw, col="red3", pch=15, align=TRUE)
appendphivic<-MYappend2tiplabel(cdo504Atree, tips=phivic, col="red3", pch=19, align=TRUE)
appendphovic<-MYappend2tiplabel(cdo504Atree, tips=phovic, col="orangered2", pch=19, align=TRUE)
appendpphact<-MYappend2tiplabel(cdo504Atree, tips=pphact, col="goldenrod1", pch=15, align=TRUE)
appendpphnsw<-MYappend2tiplabel(cdo504Atree, tips=pphnsw, col="goldenrod1", pch=15, align=TRUE)
appendpphvic<-MYappend2tiplabel(cdo504Atree, tips=pphvic, col="goldenrod1", pch=19, align=TRUE)
appendpsiact<-MYappend2tiplabel(cdo504Atree, tips=psiact, col="gray50", pch=15, align=TRUE)
appendpsinsw<-MYappend2tiplabel(cdo504Atree, tips=psinsw, col="gray50", pch=15, align=TRUE)
appendpsivic<-MYappend2tiplabel(cdo504Atree, tips=psivic, col="gray50", pch=19, align=TRUE)

dev.off()

#################################
#         CDO504B TREE          #
#################################
setwd("/Applications/BEAST Suite/Runs Sept 2013/Combined Tree Files")
cdo504Bdata<-read.beast.table("CDO504B_9runs_downsampled_7000_mcc.tree", digits=2)

cdo504Btree<-read.beast("CDO504B_9runs_downsampled_7000_mcc.tree")

cdo504Btree<-ladderize(cdo504Btree, right=FALSE)

plot.phylo(cdo504Btree, cex=0.3, label.offset=0.01)

lab<-cdo504Btree$tip.label

pannsw<-grep("annua", lab)

pclnsw<-grep("clivicola_NSW", lab)
pclact<-grep("clivicola_ACT", lab)
pclvic<-grep("clivicola_Vic", lab)
pcltas<-grep("clivicola_Tas", lab)

pconsw<-grep("costiniana_NSW", lab)
pcotas<-grep("costiniana_Tas", lab)
pcovic<-grep("costiniana_Vic", lab)

pennsw<-grep("ensiformis_NSW", lab)
penact<-grep("ensiformis_ACT", lab)
penvic<-grep("ensiformis_Vic", lab)

pfansw<-grep("fawcettiae_NSW", lab)
pfaact<-grep("fawcettiae_ACT", lab)
pfavic<-grep("fawcettiae_Vic", lab)
pfatas<-grep("fawcettiae_Tas", lab)

pgutas<-grep("gunnii_Tas", lab)

phensw<-grep("helmsii_NSW", lab)
pheact<-grep("helmsii_ACT", lab)

phinsw<-grep("hiemata_NSW", lab)
phiact<-grep("hiemata_ACT", lab)
phivic<-grep("hiemata_Vic", lab)
phitas<-grep("hiemata_Tas", lab)

phovic<-grep("hothamensis_Vic", lab)

pphnsw<-grep("phillipsiana_NSW", lab)
pphact<-grep("phillipsiana_ACT", lab)
pphvic<-grep("phillipsiana_Vic", lab)
pphtas<-grep("phillipsiana_Tas", lab)

psinsw<-grep("sieberiana_NSW", lab)
psiact<-grep("sieberiana_ACT", lab)
psivic<-grep("sieberiana_Vic", lab)
psitas<-grep("sieberiana_Tas", lab)

pdf("CDO504B MCC Tree With Support And Labels.pdf", width=10, height=10, useDingbats=FALSE)

plot.phylo(cdo504Btree, cex=0.3, label.offset=0.01)

add.scale.bar(x=0, y=30, length=0.05, cex=0.5)
node.support(as.data.frame(cdo504Bdata)$posterior, cutoff=0.8, pos="right", mode="dots", col="black")
node.support(as.data.frame(cdo504Bdata)$posterior, cutoff=0.95, pos="right", mode="dots", col="blue")
node.support(as.data.frame(cdo504Bdata)$posterior, cutoff=0.99, pos="right", mode="dots", col="red")

points(x=rep(0.001, times=11), y=seq(70, 40, length.out=11), pch=symbolsnsw, col=cols, cex=0.4)
points(x=rep(0.01, times=11), y=seq(70, 40, length.out=11), pch=symbolsvic, col=cols, cex=0.4)
points(x=rep(0.019, times=11), y=seq(70, 40, length.out=11), pch=symbolstas, col=cols, cex=0.4)
text(x=rep(0.024, times=11), y=seq(70, 40, length.out=11), labels=names, adj=0, cex=0.4, font=3)
text(x=c(0.001, 0.01, 0.019), y=rep(75, times=11), labels=c("NSW", "Vic", "Tas"), cex=0.4)

appendpannsw<-MYappend2tiplabel(cdo504Btree, tips=pannsw, col="black", pch=15, align=TRUE)
appendpclnsw<-MYappend2tiplabel(cdo504Btree, tips=pclnsw, col="mediumpurple", pch=15, align=TRUE)
appendpclact<-MYappend2tiplabel(cdo504Btree, tips=pclact, col="mediumpurple", pch=15, align=TRUE)
appendpclvic<-MYappend2tiplabel(cdo504Btree, tips=pclvic, col="mediumpurple", pch=19, align=TRUE)
appendpcoact<-MYappend2tiplabel(cdo504Btree, tips=pcoact, col="royalblue4", pch=15, align=TRUE)
appendpconsw<-MYappend2tiplabel(cdo504Btree, tips=pconsw, col="royalblue4", pch=15, align=TRUE)
appendpcovic<-MYappend2tiplabel(cdo504Btree, tips=pcovic, col="royalblue4", pch=19, align=TRUE)
appendpcotas<-MYappend2tiplabel(cdo504Btree, tips=pcotas, col="royalblue4", pch=17, align=TRUE)
appendpenact<-MYappend2tiplabel(cdo504Btree, tips=penact, col="black", bg="grey50", pch=22, align=TRUE)
appendpennsw<-MYappend2tiplabel(cdo504Btree, tips=pennsw, col="black", bg="grey50", pch=22, align=TRUE)
appendpenvic<-MYappend2tiplabel(cdo504Btree, tips=penvic, col="black", bg="grey50", pch=21, align=TRUE)
appendpfatas<-MYappend2tiplabel(cdo504Btree, tips=pfatas, col="springgreen4", pch=17, align=TRUE)
appendpfaact<-MYappend2tiplabel(cdo504Btree, tips=pfaact, col="springgreen4", pch=15, align=TRUE)
appendpfansw<-MYappend2tiplabel(cdo504Btree, tips=pfansw, col="springgreen4", pch=15, align=TRUE)
appendpfavic<-MYappend2tiplabel(cdo504Btree, tips=pfavic, col="springgreen4", pch=19, align=TRUE)
appendpgutas<-MYappend2tiplabel(cdo504Btree, tips=pgutas, col="black", pch=17, align=TRUE)
appendphensw<-MYappend2tiplabel(cdo504Btree, tips=phensw, col="turquoise3", pch=15, align=TRUE)
appendpheact<-MYappend2tiplabel(cdo504Btree, tips=pheact, col="turquoise3", pch=15)
appendphitas<-MYappend2tiplabel(cdo504Btree, tips=phitas, col="red3", pch=17, align=TRUE)
appendphiact<-MYappend2tiplabel(cdo504Btree, tips=phiact, col="red3", pch=15, align=TRUE)
appendphinsw<-MYappend2tiplabel(cdo504Btree, tips=phinsw, col="red3", pch=15, align=TRUE)
appendphivic<-MYappend2tiplabel(cdo504Btree, tips=phivic, col="red3", pch=19, align=TRUE)
appendphovic<-MYappend2tiplabel(cdo504Btree, tips=phovic, col="orangered2", pch=19, align=TRUE)
appendpphact<-MYappend2tiplabel(cdo504Btree, tips=pphact, col="goldenrod1", pch=15, align=TRUE)
appendpphnsw<-MYappend2tiplabel(cdo504Btree, tips=pphnsw, col="goldenrod1", pch=15, align=TRUE)
appendpphvic<-MYappend2tiplabel(cdo504Btree, tips=pphvic, col="goldenrod1", pch=19, align=TRUE)
appendpsiact<-MYappend2tiplabel(cdo504Btree, tips=psiact, col="gray50", pch=15, align=TRUE)
appendpsinsw<-MYappend2tiplabel(cdo504Btree, tips=psinsw, col="gray50", pch=15, align=TRUE)
appendpsivic<-MYappend2tiplabel(cdo504Btree, tips=psivic, col="gray50", pch=19, align=TRUE)

dev.off()

#################################
#         trx TREE          #
#################################

trxdata<-read.beast.table("trx_9runs_downsampled_7000_mcc.tree", digits=2)

trxtree<-read.beast("trx_9runs_downsampled_7000_mcc.tree")
trxtree<-ladderize(trxtree, right=FALSE)

plot.phylo(trxtree, cex=0.3, label.offset=0.01)

lab<-trxtree$tip.label

pannsw<-grep("annua", lab)

pclnsw<-grep("clivicola_NSW", lab)
pclact<-grep("clivicola_ACT", lab)
pclvic<-grep("clivicola_Vic", lab)
pcltas<-grep("clivicola_Tas", lab)

pconsw<-grep("costiniana_NSW", lab)
pcotas<-grep("costiniana_Tas", lab)
pcovic<-grep("costiniana_Vic", lab)

pennsw<-grep("ensiformis_NSW", lab)
penact<-grep("ensiformis_ACT", lab)
penvic<-grep("ensiformis_Vic", lab)

pfansw<-grep("fawcettiae_NSW", lab)
pfaact<-grep("fawcettiae_ACT", lab)
pfavic<-grep("fawcettiae_Vic", lab)
pfatas<-grep("fawcettiae_Tas", lab)

pgutas<-grep("gunnii_Tas", lab)

phensw<-grep("helmsii_NSW", lab)
pheact<-grep("helmsii_ACT", lab)

phinsw<-grep("hiemata_NSW", lab)
phiact<-grep("hiemata_ACT", lab)
phivic<-grep("hiemata_Vic", lab)
phitas<-grep("hiemata_Tas", lab)

phovic<-grep("hothamensis_Vic", lab)

pphnsw<-grep("phillipsiana_NSW", lab)
pphact<-grep("phillipsiana_ACT", lab)
pphvic<-grep("phillipsiana_Vic", lab)
pphtas<-grep("phillipsiana_Tas", lab)

psinsw<-grep("sieberiana_NSW", lab)
psiact<-grep("sieberiana_ACT", lab)
psivic<-grep("sieberiana_Vic", lab)
psitas<-grep("sieberiana_Tas", lab)

pdf("Trx MCC Tree With Support And Labels.pdf", width=10, height=10, useDingbats=FALSE)

plot.phylo(trxtree, cex=0.3, label.offset=0.01)

add.scale.bar(x=0, y=30, length=0.05, cex=0.5)
node.support(as.data.frame(trxdata)$posterior, cutoff=0.8, pos="right", mode="dots", col="black")
node.support(as.data.frame(trxdata)$posterior, cutoff=0.95, pos="right", mode="dots", col="blue")
node.support(as.data.frame(trxdata)$posterior, cutoff=0.99, pos="right", mode="dots", col="red")

points(x=rep(0.001, times=11), y=seq(70, 40, length.out=11), pch=symbolsnsw, col=cols, cex=0.4)
points(x=rep(0.01, times=11), y=seq(70, 40, length.out=11), pch=symbolsvic, col=cols, cex=0.4)
points(x=rep(0.019, times=11), y=seq(70, 40, length.out=11), pch=symbolstas, col=cols, cex=0.4)
text(x=rep(0.024, times=11), y=seq(70, 40, length.out=11), labels=names, adj=0, cex=0.4, font=3)
text(x=c(0.001, 0.01, 0.019), y=rep(75, times=11), labels=c("NSW", "Vic", "Tas"), cex=0.4)

appendpannsw<-MYappend2tiplabel(trxtree, tips=pannsw, col="black", pch=15, align=TRUE)
appendpclnsw<-MYappend2tiplabel(trxtree, tips=pclnsw, col="mediumpurple", pch=15, align=TRUE)
appendpclact<-MYappend2tiplabel(trxtree, tips=pclact, col="mediumpurple", pch=15, align=TRUE)
appendpclvic<-MYappend2tiplabel(trxtree, tips=pclvic, col="mediumpurple", pch=19, align=TRUE)
appendpcoact<-MYappend2tiplabel(trxtree, tips=pcoact, col="royalblue4", pch=15, align=TRUE)
appendpconsw<-MYappend2tiplabel(trxtree, tips=pconsw, col="royalblue4", pch=15, align=TRUE)
appendpcovic<-MYappend2tiplabel(trxtree, tips=pcovic, col="royalblue4", pch=19, align=TRUE)
appendpcotas<-MYappend2tiplabel(trxtree, tips=pcotas, col="royalblue4", pch=17, align=TRUE)
appendpenact<-MYappend2tiplabel(trxtree, tips=penact, col="black", bg="grey50", pch=22, align=TRUE)
appendpennsw<-MYappend2tiplabel(trxtree, tips=pennsw, col="black", bg="grey50", pch=22, align=TRUE)
appendpenvic<-MYappend2tiplabel(trxtree, tips=penvic, col="black", bg="grey50", pch=21, align=TRUE)
appendpfatas<-MYappend2tiplabel(trxtree, tips=pfatas, col="springgreen4", pch=17, align=TRUE)
appendpfaact<-MYappend2tiplabel(trxtree, tips=pfaact, col="springgreen4", pch=15, align=TRUE)
appendpfansw<-MYappend2tiplabel(trxtree, tips=pfansw, col="springgreen4", pch=15, align=TRUE)
appendpfavic<-MYappend2tiplabel(trxtree, tips=pfavic, col="springgreen4", pch=19, align=TRUE)
appendpgutas<-MYappend2tiplabel(trxtree, tips=pgutas, col="black", pch=17, align=TRUE)
appendphensw<-MYappend2tiplabel(trxtree, tips=phensw, col="turquoise3", pch=15, align=TRUE)
appendpheact<-MYappend2tiplabel(trxtree, tips=pheact, col="turquoise3", pch=15)
appendphitas<-MYappend2tiplabel(trxtree, tips=phitas, col="red3", pch=17, align=TRUE)
appendphiact<-MYappend2tiplabel(trxtree, tips=phiact, col="red3", pch=15, align=TRUE)
appendphinsw<-MYappend2tiplabel(trxtree, tips=phinsw, col="red3", pch=15, align=TRUE)
appendphivic<-MYappend2tiplabel(trxtree, tips=phivic, col="red3", pch=19, align=TRUE)
appendphovic<-MYappend2tiplabel(trxtree, tips=phovic, col="orangered2", pch=19, align=TRUE)
appendpphact<-MYappend2tiplabel(trxtree, tips=pphact, col="goldenrod1", pch=15, align=TRUE)
appendpphnsw<-MYappend2tiplabel(trxtree, tips=pphnsw, col="goldenrod1", pch=15, align=TRUE)
appendpphvic<-MYappend2tiplabel(trxtree, tips=pphvic, col="goldenrod1", pch=19, align=TRUE)
appendpsiact<-MYappend2tiplabel(trxtree, tips=psiact, col="gray50", pch=15, align=TRUE)
appendpsinsw<-MYappend2tiplabel(trxtree, tips=psinsw, col="gray50", pch=15, align=TRUE)
appendpsivic<-MYappend2tiplabel(trxtree, tips=psivic, col="gray50", pch=19, align=TRUE)

dev.off()

#################################
#         Waxy A TREE           #
#################################

waxyAdata<-read.beast.table("waxyA_9runs_downsampled_7000_mcc.tree", digits=2)

waxyAtree<-read.beast("waxyA_9runs_downsampled_7000_mcc.tree")
waxyAtree<-ladderize(waxyAtree, right=FALSE)

plot.phylo(waxyAtree, cex=0.3, label.offset=0.005)

lab<-waxyAtree$tip.label

pannsw<-grep("annua", lab)

pclnsw<-grep("clivicola_NSW", lab)
pclact<-grep("clivicola_ACT", lab)
pclvic<-grep("clivicola_Vic", lab)
pcltas<-grep("clivicola_Tas", lab)

pconsw<-grep("costiniana_NSW", lab)
pcotas<-grep("costiniana_Tas", lab)
pcovic<-grep("costiniana_Vic", lab)

pennsw<-grep("ensiformis_NSW", lab)
penact<-grep("ensiformis_ACT", lab)
penvic<-grep("ensiformis_Vic", lab)

pfansw<-grep("fawcettiae_NSW", lab)
pfaact<-grep("fawcettiae_ACT", lab)
pfavic<-grep("fawcettiae_Vic", lab)
pfatas<-grep("fawcettiae_Tas", lab)

pgutas<-grep("gunnii_Tas", lab)

phensw<-grep("helmsii_NSW", lab)
pheact<-grep("helmsii_ACT", lab)

phinsw<-grep("hiemata_NSW", lab)
phiact<-grep("hiemata_ACT", lab)
phivic<-grep("hiemata_Vic", lab)
phitas<-grep("hiemata_Tas", lab)

phovic<-grep("hothamensis_Vic", lab)

pphnsw<-grep("phillipsiana_NSW", lab)
pphact<-grep("phillipsiana_ACT", lab)
pphvic<-grep("phillipsiana_Vic", lab)
pphtas<-grep("phillipsiana_Tas", lab)

psinsw<-grep("sieberiana_NSW", lab)
psiact<-grep("sieberiana_ACT", lab)
psivic<-grep("sieberiana_Vic", lab)
psitas<-grep("sieberiana_Tas", lab)

pdf("WaxyA MCC Tree With Support And Labels.pdf", width=10, height=12, useDingbats=FALSE)

plot.phylo(waxyAtree, cex=0.3, label.offset=0.005)

add.scale.bar(x=0, y=35, length=0.05, cex=0.5)
node.support(as.data.frame(waxyAdata)$posterior, cutoff=0.8, pos="right", mode="dots", col="black")
node.support(as.data.frame(waxyAdata)$posterior, cutoff=0.95, pos="right", mode="dots", col="blue")
node.support(as.data.frame(waxyAdata)$posterior, cutoff=0.99, pos="right", mode="dots", col="red")

points(x=rep(0.001, times=11), y=seq(80, 50, length.out=11), pch=symbolsnsw, col=cols, cex=0.4)
points(x=rep(0.006, times=11), y=seq(80, 50, length.out=11), pch=symbolsvic, col=cols, cex=0.4)
points(x=rep(0.011, times=11), y=seq(80, 50, length.out=11), pch=symbolstas, col=cols, cex=0.4)
text(x=rep(0.015, times=11), y=seq(80, 50, length.out=11), labels=names, adj=0, cex=0.4, font=3)
text(x=c(0.001, 0.006, 0.011), y=rep(85, times=11), labels=c("NSW", "Vic", "Tas"), cex=0.4)

appendpannsw<-MYappend2tiplabel(waxyAtree, tips=pannsw, col="black", pch=15, align=TRUE)
appendpclnsw<-MYappend2tiplabel(waxyAtree, tips=pclnsw, col="mediumpurple", pch=15, align=TRUE)
appendpclact<-MYappend2tiplabel(waxyAtree, tips=pclact, col="mediumpurple", pch=15, align=TRUE)
appendpclvic<-MYappend2tiplabel(waxyAtree, tips=pclvic, col="mediumpurple", pch=19, align=TRUE)
appendpcoact<-MYappend2tiplabel(waxyAtree, tips=pcoact, col="royalblue4", pch=15, align=TRUE)
appendpconsw<-MYappend2tiplabel(waxyAtree, tips=pconsw, col="royalblue4", pch=15, align=TRUE)
appendpcovic<-MYappend2tiplabel(waxyAtree, tips=pcovic, col="royalblue4", pch=19, align=TRUE)
appendpcotas<-MYappend2tiplabel(waxyAtree, tips=pcotas, col="royalblue4", pch=17, align=TRUE)
appendpenact<-MYappend2tiplabel(waxyAtree, tips=penact, col="black", bg="grey50", pch=22, align=TRUE)
appendpennsw<-MYappend2tiplabel(waxyAtree, tips=pennsw, col="black", bg="grey50", pch=22, align=TRUE)
appendpenvic<-MYappend2tiplabel(waxyAtree, tips=penvic, col="black", bg="grey50", pch=21, align=TRUE)
appendpfatas<-MYappend2tiplabel(waxyAtree, tips=pfatas, col="springgreen4", pch=17, align=TRUE)
appendpfaact<-MYappend2tiplabel(waxyAtree, tips=pfaact, col="springgreen4", pch=15, align=TRUE)
appendpfansw<-MYappend2tiplabel(waxyAtree, tips=pfansw, col="springgreen4", pch=15, align=TRUE)
appendpfavic<-MYappend2tiplabel(waxyAtree, tips=pfavic, col="springgreen4", pch=19, align=TRUE)
appendpgutas<-MYappend2tiplabel(waxyAtree, tips=pgutas, col="black", pch=17, align=TRUE)
appendphensw<-MYappend2tiplabel(waxyAtree, tips=phensw, col="turquoise3", pch=15, align=TRUE)
appendpheact<-MYappend2tiplabel(waxyAtree, tips=pheact, col="turquoise3", pch=15)
appendphitas<-MYappend2tiplabel(waxyAtree, tips=phitas, col="red3", pch=17, align=TRUE)
appendphiact<-MYappend2tiplabel(waxyAtree, tips=phiact, col="red3", pch=15, align=TRUE)
appendphinsw<-MYappend2tiplabel(waxyAtree, tips=phinsw, col="red3", pch=15, align=TRUE)
appendphivic<-MYappend2tiplabel(waxyAtree, tips=phivic, col="red3", pch=19, align=TRUE)
appendphovic<-MYappend2tiplabel(waxyAtree, tips=phovic, col="orangered2", pch=19, align=TRUE)
appendpphact<-MYappend2tiplabel(waxyAtree, tips=pphact, col="goldenrod1", pch=15, align=TRUE)
appendpphnsw<-MYappend2tiplabel(waxyAtree, tips=pphnsw, col="goldenrod1", pch=15, align=TRUE)
appendpphvic<-MYappend2tiplabel(waxyAtree, tips=pphvic, col="goldenrod1", pch=19, align=TRUE)
appendpsiact<-MYappend2tiplabel(waxyAtree, tips=psiact, col="gray50", pch=15, align=TRUE)
appendpsinsw<-MYappend2tiplabel(waxyAtree, tips=psinsw, col="gray50", pch=15, align=TRUE)
appendpsivic<-MYappend2tiplabel(waxyAtree, tips=psivic, col="gray50", pch=19, align=TRUE)

dev.off()

#################################
#         Waxy B TREE           #
#################################

waxyBdata<-read.beast.table("waxyB_9runs_downsampled_7000_mcc.tree", digits=2)

waxyBtree<-read.beast("waxyB_9runs_downsampled_7000_mcc.tree")

plot.phylo(waxyBtree, cex=0.3, label.offset=0.005)
waxyBtree<-ladderize(waxyBtree, right=FALSE)

lab<-waxyBtree$tip.label

pannsw<-grep("annua", lab)

pclnsw<-grep("clivicola_NSW", lab)
pclact<-grep("clivicola_ACT", lab)
pclvic<-grep("clivicola_Vic", lab)
pcltas<-grep("clivicola_Tas", lab)

pconsw<-grep("costiniana_NSW", lab)
pcotas<-grep("costiniana_Tas", lab)
pcovic<-grep("costiniana_Vic", lab)

pennsw<-grep("ensiformis_NSW", lab)
penact<-grep("ensiformis_ACT", lab)
penvic<-grep("ensiformis_Vic", lab)

pfansw<-grep("fawcettiae_NSW", lab)
pfaact<-grep("fawcettiae_ACT", lab)
pfavic<-grep("fawcettiae_Vic", lab)
pfatas<-grep("fawcettiae_Tas", lab)

pgutas<-grep("gunnii_Tas", lab)

phensw<-grep("helmsii_NSW", lab)
pheact<-grep("helmsii_ACT", lab)

phinsw<-grep("hiemata_NSW", lab)
phiact<-grep("hiemata_ACT", lab)
phivic<-grep("hiemata_Vic", lab)
phitas<-grep("hiemata_Tas", lab)

phovic<-grep("hothamensis_Vic", lab)

pphnsw<-grep("phillipsiana_NSW", lab)
pphact<-grep("phillipsiana_ACT", lab)
pphvic<-grep("phillipsiana_Vic", lab)
pphtas<-grep("phillipsiana_Tas", lab)

psinsw<-grep("sieberiana_NSW", lab)
psiact<-grep("sieberiana_ACT", lab)
psivic<-grep("sieberiana_Vic", lab)
psitas<-grep("sieberiana_Tas", lab)

pdf("WaxyB MCC Tree With Support And Labels.pdf", width=10, height=12, useDingbats=FALSE)

plot.phylo(waxyBtree, cex=0.3, label.offset=0.005)

add.scale.bar(x=0, y=35, length=0.02, cex=0.5)
node.support(as.data.frame(waxyBdata)$posterior, cutoff=0.8, pos="right", mode="dots", col="black")
node.support(as.data.frame(waxyBdata)$posterior, cutoff=0.95, pos="right", mode="dots", col="blue")
node.support(as.data.frame(waxyBdata)$posterior, cutoff=0.99, pos="right", mode="dots", col="red")

points(x=rep(0.001, times=11), y=seq(80, 50, length.out=11), pch=symbolsnsw, col=cols, cex=0.4)
points(x=rep(0.006, times=11), y=seq(80, 50, length.out=11), pch=symbolsvic, col=cols, cex=0.4)
points(x=rep(0.011, times=11), y=seq(80, 50, length.out=11), pch=symbolstas, col=cols, cex=0.4)
text(x=rep(0.015, times=11), y=seq(80, 50, length.out=11), labels=names, adj=0, cex=0.4, font=3)
text(x=c(0.001, 0.006, 0.011), y=rep(85, times=11), labels=c("NSW", "Vic", "Tas"), cex=0.4)

appendpannsw<-MYappend2tiplabel(waxyBtree, tips=pannsw, col="black", pch=15, align=TRUE)
appendpclnsw<-MYappend2tiplabel(waxyBtree, tips=pclnsw, col="mediumpurple", pch=15, align=TRUE)
appendpclact<-MYappend2tiplabel(waxyBtree, tips=pclact, col="mediumpurple", pch=15, align=TRUE)
appendpclvic<-MYappend2tiplabel(waxyBtree, tips=pclvic, col="mediumpurple", pch=19, align=TRUE)
appendpcoact<-MYappend2tiplabel(waxyBtree, tips=pcoact, col="royalblue4", pch=15, align=TRUE)
appendpconsw<-MYappend2tiplabel(waxyBtree, tips=pconsw, col="royalblue4", pch=15, align=TRUE)
appendpcovic<-MYappend2tiplabel(waxyBtree, tips=pcovic, col="royalblue4", pch=19, align=TRUE)
appendpcotas<-MYappend2tiplabel(waxyBtree, tips=pcotas, col="royalblue4", pch=17, align=TRUE)
appendpenact<-MYappend2tiplabel(waxyBtree, tips=penact, col="black", bg="grey50", pch=22, align=TRUE)
appendpennsw<-MYappend2tiplabel(waxyBtree, tips=pennsw, col="black", bg="grey50", pch=22, align=TRUE)
appendpenvic<-MYappend2tiplabel(waxyBtree, tips=penvic, col="black", bg="grey50", pch=21, align=TRUE)
appendpfatas<-MYappend2tiplabel(waxyBtree, tips=pfatas, col="springgreen4", pch=17, align=TRUE)
appendpfaact<-MYappend2tiplabel(waxyBtree, tips=pfaact, col="springgreen4", pch=15, align=TRUE)
appendpfansw<-MYappend2tiplabel(waxyBtree, tips=pfansw, col="springgreen4", pch=15, align=TRUE)
appendpfavic<-MYappend2tiplabel(waxyBtree, tips=pfavic, col="springgreen4", pch=19, align=TRUE)
appendpgutas<-MYappend2tiplabel(waxyBtree, tips=pgutas, col="black", pch=17, align=TRUE)
appendphensw<-MYappend2tiplabel(waxyBtree, tips=phensw, col="turquoise3", pch=15, align=TRUE)
appendpheact<-MYappend2tiplabel(waxyBtree, tips=pheact, col="turquoise3", pch=15)
appendphitas<-MYappend2tiplabel(waxyBtree, tips=phitas, col="red3", pch=17, align=TRUE)
appendphiact<-MYappend2tiplabel(waxyBtree, tips=phiact, col="red3", pch=15, align=TRUE)
appendphinsw<-MYappend2tiplabel(waxyBtree, tips=phinsw, col="red3", pch=15, align=TRUE)
appendphivic<-MYappend2tiplabel(waxyBtree, tips=phivic, col="red3", pch=19, align=TRUE)
appendphovic<-MYappend2tiplabel(waxyBtree, tips=phovic, col="orangered2", pch=19, align=TRUE)
appendpphact<-MYappend2tiplabel(waxyBtree, tips=pphact, col="goldenrod1", pch=15, align=TRUE)
appendpphnsw<-MYappend2tiplabel(waxyBtree, tips=pphnsw, col="goldenrod1", pch=15, align=TRUE)
appendpphvic<-MYappend2tiplabel(waxyBtree, tips=pphvic, col="goldenrod1", pch=19, align=TRUE)
appendpsiact<-MYappend2tiplabel(waxyBtree, tips=psiact, col="gray50", pch=15, align=TRUE)
appendpsinsw<-MYappend2tiplabel(waxyBtree, tips=psinsw, col="gray50", pch=15, align=TRUE)
appendpsivic<-MYappend2tiplabel(waxyBtree, tips=psivic, col="gray50", pch=19, align=TRUE)

dev.off()

#################################
#         DMC A TREE           #
#################################

DMCAdata<-read.beast.table("DMCA_9runs_downsampled_7000_mcc.tree", digits=2)

DMCAtree<-read.beast("DMCA_9runs_downsampled_7000_mcc.tree")
DMCAtree<-ladderize(DMCAtree, right=FALSE)

plot.phylo(DMCAtree, cex=0.3, label.offset=0.005)

lab<-DMCAtree$tip.label

pannsw<-grep("annua", lab)

pclnsw<-grep("clivicola_NSW", lab)
pclact<-grep("clivicola_ACT", lab)
pclvic<-grep("clivicola_Vic", lab)
pcltas<-grep("clivicola_Tas", lab)

pconsw<-grep("costiniana_NSW", lab)
pcotas<-grep("costiniana_Tas", lab)
pcovic<-grep("costiniana_Vic", lab)

pennsw<-grep("ensiformis_NSW", lab)
penact<-grep("ensiformis_ACT", lab)
penvic<-grep("ensiformis_Vic", lab)

pfansw<-grep("fawcettiae_NSW", lab)
pfaact<-grep("fawcettiae_ACT", lab)
pfavic<-grep("fawcettiae_Vic", lab)
pfatas<-grep("fawcettiae_Tas", lab)

pgutas<-grep("gunnii_Tas", lab)

phensw<-grep("helmsii_NSW", lab)
pheact<-grep("helmsii_ACT", lab)

phinsw<-grep("hiemata_NSW", lab)
phiact<-grep("hiemata_ACT", lab)
phivic<-grep("hiemata_Vic", lab)
phitas<-grep("hiemata_Tas", lab)

phovic<-grep("hothamensis_Vic", lab)

pphnsw<-grep("phillipsiana_NSW", lab)
pphact<-grep("phillipsiana_ACT", lab)
pphvic<-grep("phillipsiana_Vic", lab)
pphtas<-grep("phillipsiana_Tas", lab)

psinsw<-grep("sieberiana_NSW", lab)
psiact<-grep("sieberiana_ACT", lab)
psivic<-grep("sieberiana_Vic", lab)
psitas<-grep("sieberiana_Tas", lab)

pdf("DMCA MCC Tree With Support And Labels.pdf", width=10, height=12, useDingbats=FALSE)

plot.phylo(DMCAtree, cex=0.3, label.offset=0.005)

add.scale.bar(x=0, y=30, length=0.05, cex=0.5)
node.support(as.data.frame(DMCAdata)$posterior, cutoff=0.8, pos="right", mode="dots", col="black")
node.support(as.data.frame(DMCAdata)$posterior, cutoff=0.95, pos="right", mode="dots", col="blue")
node.support(as.data.frame(DMCAdata)$posterior, cutoff=0.99, pos="right", mode="dots", col="red")

points(x=rep(0.001, times=11), y=seq(80, 50, length.out=11), pch=symbolsnsw, col=cols, cex=0.4)
points(x=rep(0.006, times=11), y=seq(80, 50, length.out=11), pch=symbolsvic, col=cols, cex=0.4)
points(x=rep(0.011, times=11), y=seq(80, 50, length.out=11), pch=symbolstas, col=cols, cex=0.4)
text(x=rep(0.015, times=11), y=seq(80, 50, length.out=11), labels=names, adj=0, cex=0.4, font=3)
text(x=c(0.001, 0.006, 0.011), y=rep(85, times=11), labels=c("NSW", "Vic", "Tas"), cex=0.4)

appendpannsw<-MYappend2tiplabel(DMCAtree, tips=pannsw, col="black", pch=15, align=TRUE)
appendpclnsw<-MYappend2tiplabel(DMCAtree, tips=pclnsw, col="mediumpurple", pch=15, align=TRUE)
appendpclact<-MYappend2tiplabel(DMCAtree, tips=pclact, col="mediumpurple", pch=15, align=TRUE)
appendpclvic<-MYappend2tiplabel(DMCAtree, tips=pclvic, col="mediumpurple", pch=19, align=TRUE)
appendpcoact<-MYappend2tiplabel(DMCAtree, tips=pcoact, col="royalblue4", pch=15, align=TRUE)
appendpconsw<-MYappend2tiplabel(DMCAtree, tips=pconsw, col="royalblue4", pch=15, align=TRUE)
appendpcovic<-MYappend2tiplabel(DMCAtree, tips=pcovic, col="royalblue4", pch=19, align=TRUE)
appendpcotas<-MYappend2tiplabel(DMCAtree, tips=pcotas, col="royalblue4", pch=17, align=TRUE)
appendpenact<-MYappend2tiplabel(DMCAtree, tips=penact, col="black", bg="grey50", pch=22, align=TRUE)
appendpennsw<-MYappend2tiplabel(DMCAtree, tips=pennsw, col="black", bg="grey50", pch=22, align=TRUE)
appendpenvic<-MYappend2tiplabel(DMCAtree, tips=penvic, col="black", bg="grey50", pch=21, align=TRUE)
appendpfatas<-MYappend2tiplabel(DMCAtree, tips=pfatas, col="springgreen4", pch=17, align=TRUE)
appendpfaact<-MYappend2tiplabel(DMCAtree, tips=pfaact, col="springgreen4", pch=15, align=TRUE)
appendpfansw<-MYappend2tiplabel(DMCAtree, tips=pfansw, col="springgreen4", pch=15, align=TRUE)
appendpfavic<-MYappend2tiplabel(DMCAtree, tips=pfavic, col="springgreen4", pch=19, align=TRUE)
appendpgutas<-MYappend2tiplabel(DMCAtree, tips=pgutas, col="black", pch=17, align=TRUE)
appendphensw<-MYappend2tiplabel(DMCAtree, tips=phensw, col="turquoise3", pch=15, align=TRUE)
appendpheact<-MYappend2tiplabel(DMCAtree, tips=pheact, col="turquoise3", pch=15)
appendphitas<-MYappend2tiplabel(DMCAtree, tips=phitas, col="red3", pch=17, align=TRUE)
appendphiact<-MYappend2tiplabel(DMCAtree, tips=phiact, col="red3", pch=15, align=TRUE)
appendphinsw<-MYappend2tiplabel(DMCAtree, tips=phinsw, col="red3", pch=15, align=TRUE)
appendphivic<-MYappend2tiplabel(DMCAtree, tips=phivic, col="red3", pch=19, align=TRUE)
appendphovic<-MYappend2tiplabel(DMCAtree, tips=phovic, col="orangered2", pch=19, align=TRUE)
appendpphact<-MYappend2tiplabel(DMCAtree, tips=pphact, col="goldenrod1", pch=15, align=TRUE)
appendpphnsw<-MYappend2tiplabel(DMCAtree, tips=pphnsw, col="goldenrod1", pch=15, align=TRUE)
appendpphvic<-MYappend2tiplabel(DMCAtree, tips=pphvic, col="goldenrod1", pch=19, align=TRUE)
appendpsiact<-MYappend2tiplabel(DMCAtree, tips=psiact, col="gray50", pch=15, align=TRUE)
appendpsinsw<-MYappend2tiplabel(DMCAtree, tips=psinsw, col="gray50", pch=15, align=TRUE)
appendpsivic<-MYappend2tiplabel(DMCAtree, tips=psivic, col="gray50", pch=19, align=TRUE)

dev.off()

#################################
#         DMC B TREE           #
#################################

DMCBdata<-read.beast.table("DMCB_9runs_downsampled_7000_mcc.tree", digits=2)

DMCBtree<-read.beast("DMCB_9runs_downsampled_7000_mcc.tree")
DMCBtree<-ladderize(DMCBtree, right=FALSE)

plot.phylo(DMCBtree, cex=0.3, label.offset=0.005)

lab<-DMCBtree$tip.label

pannsw<-grep("annua", lab)

pclnsw<-grep("clivicola_NSW", lab)
pclact<-grep("clivicola_ACT", lab)
pclvic<-grep("clivicola_Vic", lab)
pcltas<-grep("clivicola_Tas", lab)

pconsw<-grep("costiniana_NSW", lab)
pcotas<-grep("costiniana_Tas", lab)
pcovic<-grep("costiniana_Vic", lab)

pennsw<-grep("ensiformis_NSW", lab)
penact<-grep("ensiformis_ACT", lab)
penvic<-grep("ensiformis_Vic", lab)

pfansw<-grep("fawcettiae_NSW", lab)
pfaact<-grep("fawcettiae_ACT", lab)
pfavic<-grep("fawcettiae_Vic", lab)
pfatas<-grep("fawcettiae_Tas", lab)

pgutas<-grep("gunnii_Tas", lab)

phensw<-grep("helmsii_NSW", lab)
pheact<-grep("helmsii_ACT", lab)

phinsw<-grep("hiemata_NSW", lab)
phiact<-grep("hiemata_ACT", lab)
phivic<-grep("hiemata_Vic", lab)
phitas<-grep("hiemata_Tas", lab)

phovic<-grep("hothamensis_Vic", lab)

pphnsw<-grep("phillipsiana_NSW", lab)
pphact<-grep("phillipsiana_ACT", lab)
pphvic<-grep("phillipsiana_Vic", lab)
pphtas<-grep("phillipsiana_Tas", lab)

psinsw<-grep("sieberiana_NSW", lab)
psiact<-grep("sieberiana_ACT", lab)
psivic<-grep("sieberiana_Vic", lab)
psitas<-grep("sieberiana_Tas", lab)

pdf("DMCB MCC Tree With Support And Labels.pdf", width=10, height=12, useDingbats=FALSE)

plot.phylo(DMCBtree, cex=0.3, label.offset=0.005)

add.scale.bar(x=0, y=30, length=0.05, cex=0.5)
node.support(as.data.frame(DMCBdata)$posterior, cutoff=0.8, pos="right", mode="dots", col="black")
node.support(as.data.frame(DMCBdata)$posterior, cutoff=0.95, pos="right", mode="dots", col="blue")
node.support(as.data.frame(DMCBdata)$posterior, cutoff=0.99, pos="right", mode="dots", col="red")

points(x=rep(0.001, times=11), y=seq(80, 50, length.out=11), pch=symbolsnsw, col=cols, cex=0.4)
points(x=rep(0.006, times=11), y=seq(80, 50, length.out=11), pch=symbolsvic, col=cols, cex=0.4)
points(x=rep(0.011, times=11), y=seq(80, 50, length.out=11), pch=symbolstas, col=cols, cex=0.4)
text(x=rep(0.015, times=11), y=seq(80, 50, length.out=11), labels=names, adj=0, cex=0.4, font=3)
text(x=c(0.001, 0.006, 0.011), y=rep(85, times=11), labels=c("NSW", "Vic", "Tas"), cex=0.4)

appendpannsw<-MYappend2tiplabel(DMCBtree, tips=pannsw, col="black", pch=15, align=TRUE)
appendpclnsw<-MYappend2tiplabel(DMCBtree, tips=pclnsw, col="mediumpurple", pch=15, align=TRUE)
appendpclact<-MYappend2tiplabel(DMCBtree, tips=pclact, col="mediumpurple", pch=15, align=TRUE)
appendpclvic<-MYappend2tiplabel(DMCBtree, tips=pclvic, col="mediumpurple", pch=19, align=TRUE)
appendpcoact<-MYappend2tiplabel(DMCBtree, tips=pcoact, col="royalblue4", pch=15, align=TRUE)
appendpconsw<-MYappend2tiplabel(DMCBtree, tips=pconsw, col="royalblue4", pch=15, align=TRUE)
appendpcovic<-MYappend2tiplabel(DMCBtree, tips=pcovic, col="royalblue4", pch=19, align=TRUE)
appendpcotas<-MYappend2tiplabel(DMCBtree, tips=pcotas, col="royalblue4", pch=17, align=TRUE)
appendpenact<-MYappend2tiplabel(DMCBtree, tips=penact, col="black", bg="grey50", pch=22, align=TRUE)
appendpennsw<-MYappend2tiplabel(DMCBtree, tips=pennsw, col="black", bg="grey50", pch=22, align=TRUE)
appendpenvic<-MYappend2tiplabel(DMCBtree, tips=penvic, col="black", bg="grey50", pch=21, align=TRUE)
appendpfatas<-MYappend2tiplabel(DMCBtree, tips=pfatas, col="springgreen4", pch=17, align=TRUE)
appendpfaact<-MYappend2tiplabel(DMCBtree, tips=pfaact, col="springgreen4", pch=15, align=TRUE)
appendpfansw<-MYappend2tiplabel(DMCBtree, tips=pfansw, col="springgreen4", pch=15, align=TRUE)
appendpfavic<-MYappend2tiplabel(DMCBtree, tips=pfavic, col="springgreen4", pch=19, align=TRUE)
appendpgutas<-MYappend2tiplabel(DMCBtree, tips=pgutas, col="black", pch=17, align=TRUE)
appendphensw<-MYappend2tiplabel(DMCBtree, tips=phensw, col="turquoise3", pch=15, align=TRUE)
appendpheact<-MYappend2tiplabel(DMCBtree, tips=pheact, col="turquoise3", pch=15)
appendphitas<-MYappend2tiplabel(DMCBtree, tips=phitas, col="red3", pch=17, align=TRUE)
appendphiact<-MYappend2tiplabel(DMCBtree, tips=phiact, col="red3", pch=15, align=TRUE)
appendphinsw<-MYappend2tiplabel(DMCBtree, tips=phinsw, col="red3", pch=15, align=TRUE)
appendphivic<-MYappend2tiplabel(DMCBtree, tips=phivic, col="red3", pch=19, align=TRUE)
appendphovic<-MYappend2tiplabel(DMCBtree, tips=phovic, col="orangered2", pch=19, align=TRUE)
appendpphact<-MYappend2tiplabel(DMCBtree, tips=pphact, col="goldenrod1", pch=15, align=TRUE)
appendpphnsw<-MYappend2tiplabel(DMCBtree, tips=pphnsw, col="goldenrod1", pch=15, align=TRUE)
appendpphvic<-MYappend2tiplabel(DMCBtree, tips=pphvic, col="goldenrod1", pch=19, align=TRUE)
appendpsiact<-MYappend2tiplabel(DMCBtree, tips=psiact, col="gray50", pch=15, align=TRUE)
appendpsinsw<-MYappend2tiplabel(DMCBtree, tips=psinsw, col="gray50", pch=15, align=TRUE)
appendpsivic<-MYappend2tiplabel(DMCBtree, tips=psivic, col="gray50", pch=19, align=TRUE)

dev.off()

########################

#making a legend

cols<-c("black", "mediumpurple", "royalblue4", "grey50", "springgreen4", "black", "turquoise3", "red3", "orangered2", "goldenrod1", "gray50")
species<-c("Poa annua", "Poa clivicola", "Poa costiniana", "Poa ensiformis", "Poa fawcettiae", "Poa gunnii", "Poa helmsii", "Poa hiemata", "Poa hothamensis" ,"Poa phillipsiana", "Poa sieberiana")
yvalues<-c(11:1)
xvalues<-rep(1,11)

pdf(file="Legend for Poa trees.pdf", width=10, height=15, useDingbats=FALSE)
plot(yvalues~xvalues, col=cols, pch=17, xlim=c(-6, 20), ylim=c(0, 37))
points(yvalues~c(xvalues+1), col=cols, pch=19)
points(yvalues~c(xvalues+2), col=cols, pch=15)
text(x=-4, y=yvalues, labels=species, adj=c(0,0), vfont=c("sans serif", "italic"), cex=0.7)
text(x=c(1,2,3), y=12, labels=c("Tas", "Vic", "NSW"), vfont=c("sans serif", "plain"), cex=0.7)
dev.off()