setwd("~/project/data")

# LOAD NORMALIZED DATA, COVARIATES AND FUNCTIONS
load("indoRNA.logCPM.TMM.filtered.Rda")
load("covariates.Rda")
source("CILP_functions.R")


###########################
### PRE-PROCESSING DATA ###
###########################

library(dplyr)

# correcting samples' names to a specific format
samples <- colnames(lcpmFilt)
individuals <- sapply(strsplit(samples, "_"), "[[",1) # get sample's id without batch
individuals <- gsub("([A-Z]{3})(\\d{3})$", "\\1-\\2", individuals) # adjust format (snug a hyphen where it's needed)
colnames(lcpmFilt) <- individuals

# removing duplicated individuals from TMM and covariates data
dupl <- which(duplicated(colnames(lcpmFilt))) # identify columns with repeated id
lcpmFilt <- subset(lcpmFilt, select = -c(dupl)) # remove duplicated columns
covariates <- covariates[!duplicated(covariates$ID), ] # remove duplicated IDs

# changing island names to make further processing easier
covariates$Island[covariates$Island == "Sumba"] <- "SMB"
covariates$Island[covariates$Island == "Mentawai"] <- "MTW"
covariates$Island[covariates$Island == "West Papua"] <- "KOR"
colnames(lcpmFilt) <- gsub("MPI", "\\1KOR", colnames(lcpmFilt))

# correcting the individual from Sumba without an age record
sumba_ages <- (subset(covariates, Island == "SMB", select = Age))$Age # get all the age records for Sumba
no_age <- which(is.na(covariates$Age) == TRUE) # get the index for the unaged individual
covariates$Age[no_age] <- round(mean(sumba_ages, na.rm=TRUE)) # allocate the average age for the island to the individual with no age


#######################
### PROCESSING DATA ###
#######################

# specify the groups to be compared
pop1 <- "MTW"
pop2 <- "KOR"

# subset normalized data according to each population
pop1Sub <- subset(lcpmFilt, select = grep(paste0("^",pop1), colnames(lcpmFilt)))
pop2Sub <- subset(lcpmFilt, select = grep(paste0("^",pop2), colnames(lcpmFilt)))

# load expression data comparison between populations
comp <- as.data.frame(read.table(paste0("topTable.voomNoNorm.tmm.filtered.dup_corrected.",pop1,"-",pop2,".txt")))


#### CILP ####

library(ggplot2)

# create the expression file containing expression values for only genes with significant differences (p<0.01; logFC>0.5) between populations
compFilt1 <- filter(comp, abs(adj.P.Val) < 0.01)
compFilt <- filter(comp, abs(adj.P.Val) < 0.01, abs(logFC) > 0.5) # filter genes according to p-value and logFC
compGenes <- compFilt$genes # save significant genes' names
compCILP <- cbind( subset(pop1Sub, rownames(pop1Sub) %in% compGenes) , subset(pop2Sub, rownames(pop2Sub) %in% compGenes)) # get expression values for significant genes in both pops
write.table(compCILP, paste0(pop1,"_",pop2,"_exp_dataset.txt"), append= FALSE, sep="\t", dec=".", row.names=TRUE, col.names=TRUE) # write the expression file

# create the predictor dataset
predictor <- c() # columns for the predictor file
covariate1 <- c() # island
covariate2 <- c() # age
covariate3 <- c() # RIN
covariate4 <- c() # gran cells
covariate5 <- c() # b cells
covariate6 <- c() # CD4T
covariate7 <- c() # CD8T
covariate8 <- c() # nk
covariate9 <- c() # monocytes

for (sample in colnames(compCILP)) { # add info in the predictor file (class and covariate value) for each sample
  if (startsWith(sample, pop1)) { # in this case, class is equal to the original location (island) of the sample
    predictor <- c(predictor,0)
  } else {
    predictor <- c(predictor,1)
  }
  if (covariates$Island[covariates$Sample.ID == sample] == pop1) {
    covariate1 <- c(covariate1,0)
  } else {
    covariate1 <- c(covariate1,1)
  }
  covariate2 <- c(covariate2, covariates$Age[covariates$Sample.ID == sample])
  covariate3 <- c(covariate3, covariates$RIN[covariates$Sample.ID == sample])
  covariate4 <- c(covariate4, covariates$Gran[covariates$Sample.ID == sample])
  covariate5 <- c(covariate5, covariates$Bcell[covariates$Sample.ID == sample])
  covariate6 <- c(covariate6, covariates$CD4T[covariates$Sample.ID == sample])
  covariate7 <- c(covariate7, covariates$CD8T[covariates$Sample.ID == sample])
  covariate8 <- c(covariate8, covariates$NK[covariates$Sample.ID == sample])
  covariate9 <- c(covariate9, covariates$Mono[covariates$Sample.ID == sample])
}
pred_dataset <- data.frame(predictor, covariate1, covariate2, covariate3, covariate4, covariate5, covariate6, covariate7, covariate8, covariate9) # create the predictor file dataframe
write.table(pred_dataset, paste0(pop1,"_",pop2,"_pred_dataset.txt"), append= FALSE, sep="\t", dec=".", row.names=FALSE, col.names= TRUE )

# cilp method
expressiondata <- read.delim(paste0(pop1,"_",pop2,"_exp_dataset.txt"))
yvars=read.delim(paste0(pop1,"_",pop2,"_pred_dataset.txt"))
CILP_withLM(expressiondata,yvars,1,2,10,pop1,pop2)

# save RESULTS (correlation p-values) and gene pairs tested
results <- as.data.frame(read.table(paste0(pop1,"_",pop2,"_completemod_results.txt"), header=TRUE))
pairs <- as.data.frame(read.table(paste0(pop1,"_",pop2,"_completemod_pairs_tested.txt"), header=TRUE))

# plot
png(file=paste0("~/project/output/qqplot_completemod_",pop1,"_",pop2,".png"), width=600, height=350)
qqplot(-log10(runif(dim(pairs)[1])),-log10(results$p_value),pch=20,main=paste0(pop1," vs. ",pop2),xlab='p-value, uniform distribution',ylab='p-value, differential correlation test',bty='n',xlim=c(0,max(-log10(results$p_value))),ylim=c(0,max(-log10(results$p_value))))
x=c(0,1);y=c(0,1);abline(lm(y~x),lty=2)
dev.off()

# plot p value distribution
png(file=paste0("~/project/output/pvalue_completemod_",pop1,"_",pop2,".png"), width=600, height=350)
hist(results$p_value, main=paste0("p-value distribution\n(",pop1,"-",pop2,")"),
     col="red", xlab="p-value",
     bins=20)
dev.off()


#######################
### PROCESS RESULTS ###
#######################

# correct the p-values for multiple testing (BH correction) and create dataframe with significant genes
library(qvalue)

p <- results$p_value
results$qvalue <- qvalue(p, fdr.level=0.05, pi0=1)$qvalues  # p-value correction
results$V1 <- pairs$V1
results$V2 <- pairs$V2
results$gene1 <- compGenes[results$V1]
results$gene2 <- compGenes[results$V2]
allGenes <- (unique(c(results$gene1, results$gene2))) # list of unique genes

sigResults <- results[results$qvalue <= 0.05,] # get the significant correlations
sigPairs <- select(sigResults, gene1, gene2) # pairs of genes with significant correlations
sigGenes <- unique(c(sigPairs$gene1, sigPairs$gene2)) # list of unique genes in significant pairs






###########################################
### IF THERE ARE SIGNIFICANT GENE PAIRS ###
###########################################

# count how many times a gene appears in the DE pairings
Count <- c()
Gene <- sigGenes
for (gene in sigGenes) {
  occurrences <- length(subset(c(sigPairs$gene1, sigPairs$gene2),  c(sigPairs$gene1, sigPairs$gene2)== gene))
  Count <- c(Count, occurrences)
}
sigGeneCount <- data.frame(Gene, Count) # create a data frame for plotting

# plot a histogram (see distribution, horizontal axis is number of times a gene appears in all the DE pairings)
png(file=paste0("~/project/output/occurrences_",pop1,"_",pop2,".png"), width=600, height=350)
hist(Gene_count_all$Count, main=paste0("Occurrence Distribution of Genes in\n DE Pairwirse Comparisons (",pop1,"-",pop2,")"),
     xlab="Occurence",
     col="darkmagenta",
     breaks=(seq(0,5,1)), xaxp=c(0,5,5),
     xlim=c(0,5), ylim=c(0,12))
dev.off()


# gene interactions: plot a network of interactions between genes in DE pairings
library(igraph)
library(biomaRt)

edgeList <- c()
for (i in seq(1,dim(sigPairs)[1],by=1)) {edgeList <- c(edgeList, sigPairs$gene1[i], sigPairs$gene2[i])}

firstTest <- results[1:19,]
ensembl75 <- useMart(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
genenames <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), mart=ensembl75, filter='ensembl_gene_id', values=edgeList)
firstTest$hgnc1 <- sapply(firstTest$gene1, function(x) (genenames[genenames$ensembl_gene_id %in% x,2]))
firstTest$hgnc2 <- sapply(firstTest$gene2, function(x) (genenames[genenames$ensembl_gene_id %in% x,2]))
gsub("^$", "na", genenames$hgnc_symbol) # changes blank spaces for "na"s
gsub("^$", "na", firstTest$hgcn2) # changes blank spaces for "na"s
firstTest[grepl("^$", firstTest$hgnc2),] <- NA # change for NAs

edgeList2 <- genenames$hgnc_symbol[match(edgeList,genenames$ensembl_gene_id)]
G <- graph(edgeList2, directed = TRUE )

# visualization
png(file=paste0("~/project/output/gene_interactions_",pop1,"_",pop2,".png"), width=600, height=350)
plot(G, layout = layout.fruchterman.reingold,
     vertex.size = 25,vertex.color="yellow",vertex.frame.color= "black",vertex.label.color = "black", vertex.label.family = "sans",
     edge.width=2, edge.color="black")
dev.off()


# visualization
plot(G, layout = layout.fruchterman.reingold,
     vertex.size = 25,vertex.color="maroon",vertex.frame.color= "black",vertex.label.color = "black", vertex.label.family = "sans",
     edge.width=2, edge.color="black")


# GO  and KEGG analyses
library(clusterProfiler)
library(organism, character.only = TRUE)
library(wordcloud)
library(enrichplot)

# GO over representation analysis
genes <- sigGeneCount$Count
genes <- sort(genes, decreasing=TRUE)
names(genes) <- sigGeneCount$Gene
organism <-"org.Hs.eg.db"

go_enrich <- enrichGO(gene = names(genes),
                      universe = names(all_genes),
                      OrgDb= organism,
                      keyType= "ENSEMBL",
                      readable = T,
                      ont="MF",
                      pvalueCutoff=1,
                      qvalueCutoff=1)
x1 <- pairwise_termsim(go_enrich)
emapplot(x1)

# KEGG enrichment analysis
ids <- bitr(names(all_genes), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb= organism)
dedup_ids <- ids[!duplicated(ids[c("ENSEMBL")]),] # remove duplicate IDs
df2 <- Gene_count_all[Gene_count_all$Gene %in% dedup_ids$ENSEMBL,]
df2$EntrezIds <- dedup_ids$ENTREZID
kegg_gene_list <- df2$Count
names(kegg_gene_list) <- df2$EntrezIds
kegg_gene_list <- na.omit(kegg_gene_list)
kegg_gene_list <- sort(kegg_gene_list, decreasing=TRUE)

kegg_organism <- "hsa"
kk2 <- gseKEGG(geneList = kegg_gene_list,
               organism = kegg_organism,
               nPerm = 10000,
               minGSSize = 1,
               maxGSSize = 2000,
               pvalueCutoff = 1,
               pAdjustMethod = "none",
               keyType = "ncbi-geneid")

x2 <- pairwise_termsim(kk2)
emapplot(x2)



#############
### EXTRA ###
#############

# map figure
library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library(sf)

population <- c("Mentawai", "Sumba", "Korowai")
latitude <- c(-3,-9.8,-6)
longitude <- c(100.3,120,140)
popsdf <- data.frame(population,latitude,longitude)

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

cleanup <- 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white', colour = 'white'), 
        axis.line = element_line(colour = "black"), legend.position="none",
        axis.ticks=element_blank(), axis.text.x=element_blank(),
        axis.text.y=element_blank())

png(file="~/project/output/indonesia_pops.png", width=600, height=350)
ggplot(data = world) +cleanup + xlab("") + ylab("") +
  geom_sf(color="darkgray", fill=ifelse(world$sovereignt == "Indonesia", "light blue", "lightgray")) +
  annotate(geom="text", x= 99, y=-4, label="Mentawai", fontface="bold", color="grey22", size= 5) +
  annotate(geom="text", x= 120, y=-11, label="Sumba", fontface="bold", color="grey22", size= 5) +
  annotate(geom="text", x= 141, y=-6.8, label="Korowai", fontface="bold", color="grey22", size= 5) +
  coord_sf(ylim = c(-15,15), xlim=c(148,93), expand=FALSE) +
  geom_point(data=popsdf, aes(x=longitude, y=latitude), colour="Deep Pink", fill="Pink", pch=21, size=5, alpha=I(0.7)) 
dev.off()

