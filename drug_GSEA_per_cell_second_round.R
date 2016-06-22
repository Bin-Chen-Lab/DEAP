setwd("/Users/binchen1/Documents/stanford/hcc/code/DEAP/")
#in the second round, using the targets predicted in the first rounds.
#initital evaluation indicates it does not give better results.

#input:
#ranked_file
#output_file
#target_file

#DEAP package used to find enriched substructures, targets, target class, pathways for the predictions from LINCS


#may used to identify potential targets
GSEA.EnrichmentScore <- function(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) {  
  #
  # Computes the weighted GSEA score of gene.set in gene.list. 
  # The weighted score type is the exponent of the correlation 
  # weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). When the score type is 1 or 2 it is 
  # necessary to input the correlation vector with the values in the same order as in the gene list.
  #
  # Inputs:
  #   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)  
  #   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset) 
  #   weighted.score.type: Type of score: weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)  
  #  correl.vector: A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list 
  #
  # Outputs:
  #   ES: Enrichment score (real number between -1 and +1) 
  #   arg.ES: Location in gene.list where the peak running enrichment occurs (peak of the "mountain") 
  #   RES: Numerical vector containing the running enrichment score for all locations in the gene list 
  #   tag.indicator: Binary vector indicating the location of the gene sets (1's) in the gene list 
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.
  
  tag.indicator <- sign(match(gene.list, gene.set, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
  no.tag.indicator <- 1 - tag.indicator 
  N <- length(gene.list) 
  Nh <- length(gene.set) 
  Nm <-  N - Nh 
  if (weighted.score.type == 0) {
    correl.vector <- rep(1, N)
  }
  alpha <- weighted.score.type
  correl.vector <- abs(correl.vector**alpha)
  sum.correl.tag    <- sum(correl.vector[tag.indicator == 1])
  norm.tag    <- 1.0/sum.correl.tag
  norm.no.tag <- 1.0/Nm
  RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)      
  max.ES <- max(RES)
  min.ES <- min(RES)
  if (max.ES > - min.ES) {
    #      ES <- max.ES
    ES <- signif(max.ES, digits = 5)
    arg.ES <- which.max(RES)
  } else {
    #      ES <- min.ES
    ES <- signif(min.ES, digits=5)
    arg.ES <- which.min(RES)
  }
  return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))    
}

GSEA.EnrichmentScore2 <- function(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) {  
  #
  # Computes the weighted GSEA score of gene.set in gene.list. It is the same calculation as in 
  # GSEA.EnrichmentScore but faster (x8) without producing the RES, arg.RES and tag.indicator outputs.
  # This call is intended to be used to asses the enrichment of random permutations rather than the 
  # observed one.
  # The weighted score type is the exponent of the correlation 
  # weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). When the score type is 1 or 2 it is 
  # necessary to input the correlation vector with the values in the same order as in the gene list.
  #
  # Inputs:
  #   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)  
  #   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset) 
  #   weighted.score.type: Type of score: weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)  
  #  correl.vector: A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list 
  #
  # Outputs:
  #   ES: Enrichment score (real number between -1 and +1) 
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.
  
  N <- length(gene.list) 
  Nh <- length(gene.set) 
  Nm <-  N - Nh 
  
  loc.vector <- vector(length=N, mode="numeric")
  peak.res.vector <- vector(length=Nh, mode="numeric")
  valley.res.vector <- vector(length=Nh, mode="numeric")
  tag.correl.vector <- vector(length=Nh, mode="numeric")
  tag.diff.vector <- vector(length=Nh, mode="numeric")
  tag.loc.vector <- vector(length=Nh, mode="numeric")
  
  loc.vector[gene.list] <- seq(1, N)
  tag.loc.vector <- loc.vector[gene.set]
  
  tag.loc.vector <- sort(tag.loc.vector, decreasing = F)
  
  if (weighted.score.type == 0) {
    tag.correl.vector <- rep(1, Nh)
  } else if (weighted.score.type == 1) {
    tag.correl.vector <- correl.vector[tag.loc.vector]
    tag.correl.vector <- abs(tag.correl.vector)
  } else if (weighted.score.type == 2) {
    tag.correl.vector <- correl.vector[tag.loc.vector]*correl.vector[tag.loc.vector]
    tag.correl.vector <- abs(tag.correl.vector)
  } else {
    tag.correl.vector <- correl.vector[tag.loc.vector]**weighted.score.type
    tag.correl.vector <- abs(tag.correl.vector)
  }
  
  norm.tag <- 1.0/sum(tag.correl.vector)
  tag.correl.vector <- tag.correl.vector * norm.tag
  norm.no.tag <- 1.0/Nm
  tag.diff.vector[1] <- (tag.loc.vector[1] - 1) 
  tag.diff.vector[2:Nh] <- tag.loc.vector[2:Nh] - tag.loc.vector[1:(Nh - 1)] - 1
  tag.diff.vector <- tag.diff.vector * norm.no.tag
  peak.res.vector <- cumsum(tag.correl.vector - tag.diff.vector)
  valley.res.vector <- peak.res.vector - tag.correl.vector
  max.ES <- max(peak.res.vector)
  min.ES <- min(valley.res.vector)
  ES <- signif(ifelse(max.ES > - min.ES, max.ES, min.ES), digits=5)
  
  return(list(ES = ES))
  
}

library(gplots)
library(ggplot2)
library(RColorBrewer)

#predictions
#pred = read.csv("data/lincs_predictions_all_cmpd.csv")
files = list.files("~/Documents/stanford/hcc/data/LIHC/moa/drug_neighbor", pattern = "hepg2")
drugs = sapply(files, function(drug_name){
  unlist(strsplit(drug_name, "_"))[1]
})

#drugs = drugs[14:length(drugs)] #c("niclosamide")
drugs = c("niclosamide") #, "wortmannin", "sunitinib" , "tamoxifen", "strophanthidin", "sorafenib", "selamectin", "sirolimus","pyrvinium-pamoate","oligomycin-a","obatoclax","gefitinib")
for (drug in drugs){
#input drug list
#pred = read.csv(paste("~/Documents/stanford/hcc/data/LIHC/moa/drug_neighbor/", drug, "_sim_drugs_all_cells.csv", sep=""))
pred = read.csv(paste("~/Documents/stanford/hcc/data/LIHC/moa/drug_neighbor/", drug, "_sim_drugs_hepg2_cells.csv", sep=""))

if (nrow(pred)<1000){next}
#suppose the one with itself is highest correlated
pred = rbind(pred, data.frame(X=100000, pert_iname=drug, cor=1))
pred$pert_iname = tolower(pred$pert_iname)

#load drug targets
load("data/lincs_cmpd_info_mesh_target_sea_v1.RData") #
#lincs_cmpd_info_mesh_target = read.csv("data/lincs_cmpd_info_mesh_target_hepg2_manual.csv")

#remove pert_id
lincs_cmpd_info_mesh_target = unique(lincs_cmpd_info_mesh_target[, !colnames(lincs_cmpd_info_mesh_target) %in% c("meshes", "V2", "pubchem_cid", "pert_id")])
#lincs_cmpd_info_mesh_target$pert_iname = lincs_cmpd_info_mesh_target$pert_iname_lower

pred_terms = merge( pred, lincs_cmpd_info_mesh_target, by.x="pert_iname", by.y="pert_iname")

pred_targets = subset(pred_terms, !is.na(sea_targets_1) & sea_targets_1 != "")
pred_targets$targets = pred_targets$sea_targets_1

#construct target sets
target_cmpd = data.frame()
for(i in 1:nrow(pred_targets)){
  target_cmpd = rbind(target_cmpd, data.frame(pert_iname = pred_targets$pert_iname[i], target= unlist(strsplit(as.character(pred_targets$targets[i]), "\t"))))
}
target_cmpd = unique(target_cmpd)

#select target sets (remove rare targets)
targets = sort(table(unlist(strsplit(paste(pred_targets$targets, collapse="\t"), "\t"))))
targets = names(targets[targets>2])

#construct random list
pred_targets$rank = rank(-1 * pred_targets$cor, ties.method= "random") #high correlation ranks better

target_cmpd_rank = merge(target_cmpd, subset(pred_targets, select=c( "rank", "pert_iname"), by.x="pert_iname", by.y="pert_iname"))

drug.list = sort(pred_targets$rank)
nperm = 1000

Ng = length(targets)
gs.names = targets
Obs.ES <- vector(length = Ng, mode = "numeric")
Obs.ES.norm <- vector(length = Ng, mode = "numeric")

phi <- matrix(nrow = Ng, ncol = nperm)
phi.norm <- matrix(nrow = Ng, ncol = nperm)
obs.phi <- matrix(nrow = Ng, ncol = nperm)
obs.phi.norm <- matrix(nrow = Ng, ncol = nperm)

print("Computing ES...")


sizes = NULL
arg.ESs = NULL
for (i in 1:Ng){
  target = targets[i]
  print(i)
  drug.set = target_cmpd_rank$rank[target_cmpd_rank$target == target]
  gsea_o = GSEA.EnrichmentScore(drug.list, drug.set, 0) 
  #plot(x=c(1:length(gsea$RES)), y = gsea$RES, type = "l" , xlab = "", ylab = "score")
  
  sizes = c(sizes, length(drug.set))
  arg.ESs = c(arg.ESs, gsea_o$arg.ES)
  
  Obs.ES[i] = gsea_o$ES
  for (r in 1:nperm){
    #target_cmpd_rank_shuffle = target_cmpd_rank
    #target_cmpd_rank_shuffle$rank = sample(target_cmpd_rank_shuffle$rank, length(target_cmpd_rank_shuffle$rank))
    drug.set = sample(drug.list, length(drug.set))
    gsea = GSEA.EnrichmentScore(drug.list, drug.set, 0)
    if (is.na(gsea$ES)){
      gsea = GSEA.EnrichmentScore2(drug.list, drug.set, 0)
    }
    phi[i, r] = gsea$ES

  }
  
  #we are not going to resample the observerd label. In our case, our observation is fixed
  obs.phi[i, 1] <- gsea_o$ES
  for (r in 2:nperm) {
    obs.phi[i, r] <- obs.phi[i, 1]
  }
  
}

print("Computing nominal p-values...")

p.vals <- matrix(0, nrow = Ng, ncol = 2)

for (i in 1:Ng) {
  pos.phi <- NULL
  neg.phi <- NULL
  for (j in 1:nperm) {
    if (phi[i, j] >= 0) {
      pos.phi <- c(pos.phi, phi[i, j]) 
    } else {
      neg.phi <- c(neg.phi, phi[i, j]) 
    }
  }
  ES.value <- Obs.ES[i]
  if (ES.value >= 0) {
    p.vals[i, 1] <- signif(sum(pos.phi >= ES.value)/length(pos.phi), digits=5)
  } else {
    p.vals[i, 1] <- signif(sum(neg.phi <= ES.value)/length(neg.phi), digits=5)
  }
}

print("Computing rescaling normalization for each gene set null...")
for (i in 1:Ng) {
  pos.phi <- NULL
  neg.phi <- NULL
  for (j in 1:nperm) {
    if (phi[i, j] >= 0) {
      pos.phi <- c(pos.phi, phi[i, j]) 
    } else {
      neg.phi <- c(neg.phi, phi[i, j]) 
    }
  }
  pos.m <- mean(pos.phi)
  neg.m <- mean(abs(neg.phi))
  
  pos.phi <- pos.phi/pos.m
  neg.phi <- neg.phi/neg.m
  for (j in 1:nperm) {
    if (phi[i, j] >= 0) {
      phi.norm[i, j] <- phi[i, j]/pos.m
    } else {
      phi.norm[i, j] <- phi[i, j]/neg.m
    }
  }
  for (j in 1:nperm) {
    if (obs.phi[i, j] >= 0) {
      obs.phi.norm[i, j] <- obs.phi[i, j]/pos.m
    } else {
      obs.phi.norm[i, j] <- obs.phi[i, j]/neg.m
    }
  }
  
  if (Obs.ES[i] >= 0) {
    Obs.ES.norm[i] <- Obs.ES[i]/pos.m
  } else {
    Obs.ES.norm[i] <- Obs.ES[i]/neg.m
  }
}


print("Computing FDR q-values...")

NES <- vector(length=Ng, mode="numeric")
phi.norm.mean  <- vector(length=Ng, mode="numeric")
obs.phi.norm.mean  <- vector(length=Ng, mode="numeric")
phi.norm.median  <- vector(length=Ng, mode="numeric")
obs.phi.norm.median  <- vector(length=Ng, mode="numeric")
phi.norm.mean  <- vector(length=Ng, mode="numeric")
obs.phi.mean  <- vector(length=Ng, mode="numeric")
FDR.mean <- vector(length=Ng, mode="numeric")
FDR.median <- vector(length=Ng, mode="numeric")
phi.norm.median.d <- vector(length=Ng, mode="numeric")
obs.phi.norm.median.d <- vector(length=Ng, mode="numeric")

Obs.ES.index <- order(Obs.ES.norm, decreasing=T)
Orig.index <- seq(1, Ng)
Orig.index <- Orig.index[Obs.ES.index]
Orig.index <- order(Orig.index, decreasing=F)
Obs.ES.norm.sorted <- Obs.ES.norm[Obs.ES.index]
gs.names.sorted <- gs.names[Obs.ES.index]

for (k in 1:Ng) {
  print(k)
  NES[k] <- Obs.ES.norm.sorted[k]
  ES.value <- NES[k]
  count.col <- vector(length=nperm, mode="numeric")
  obs.count.col <- vector(length=nperm, mode="numeric")
  for (i in 1:nperm) {
    phi.vec <- phi.norm[,i]
    obs.phi.vec <- obs.phi.norm[,i]
    if (ES.value >= 0) {
      count.col.norm <- sum(phi.vec >= 0)
      obs.count.col.norm <- sum(obs.phi.vec >= 0)
      count.col[i] <- ifelse(count.col.norm > 0, sum(phi.vec >= ES.value)/count.col.norm, 0)
      obs.count.col[i] <- ifelse(obs.count.col.norm > 0, sum(obs.phi.vec >= ES.value)/obs.count.col.norm, 0)
    } else {
      count.col.norm <- sum(phi.vec < 0)
      obs.count.col.norm <- sum(obs.phi.vec < 0)
      count.col[i] <- ifelse(count.col.norm > 0, sum(phi.vec <= ES.value)/count.col.norm, 0)
      obs.count.col[i] <- ifelse(obs.count.col.norm > 0, sum(obs.phi.vec <= ES.value)/obs.count.col.norm, 0)
    }
  }
  phi.norm.mean[k] <- mean(count.col)
  obs.phi.norm.mean[k] <- mean(obs.count.col)
  phi.norm.median[k] <- median(count.col)
  obs.phi.norm.median[k] <- median(obs.count.col)
  FDR.mean[k] <- ifelse(phi.norm.mean[k]/obs.phi.norm.mean[k] < 1, phi.norm.mean[k]/obs.phi.norm.mean[k], 1)
  FDR.median[k] <- ifelse(phi.norm.median[k]/obs.phi.norm.median[k] < 1, phi.norm.median[k]/obs.phi.norm.median[k], 1)
}

# adjust q-values
adjust.FDR.q.val = T
if (adjust.FDR.q.val == T) {
  pos.nes <- length(NES[NES >= 0])
  min.FDR.mean <- FDR.mean[pos.nes]
  min.FDR.median <- FDR.median[pos.nes]
  for (k in seq(pos.nes - 1, 1, -1)) {
    if (FDR.mean[k] < min.FDR.mean) {
      min.FDR.mean <- FDR.mean[k]
    }
    if (min.FDR.mean < FDR.mean[k]) {
      FDR.mean[k] <- min.FDR.mean
    }
  }
  
  neg.nes <- pos.nes + 1
  min.FDR.mean <- FDR.mean[neg.nes]
  min.FDR.median <- FDR.median[neg.nes]
  for (k in seq(neg.nes + 1, Ng)) {
    if (FDR.mean[k] < min.FDR.mean) {
      min.FDR.mean <- FDR.mean[k]
    }
    if (min.FDR.mean < FDR.mean[k]) {
      FDR.mean[k] <- min.FDR.mean
    }
  }
}

obs.phi.norm.mean.sorted <- obs.phi.norm.mean[Orig.index]
phi.norm.mean.sorted <- phi.norm.mean[Orig.index]
FDR.mean.sorted <- FDR.mean[Orig.index]
FDR.median.sorted <- FDR.median[Orig.index]

Obs.ES <- signif(Obs.ES, digits=3)
Obs.ES.norm <- signif(Obs.ES.norm, digits=3)
p.vals <- signif(p.vals, digits=3)
FDR.mean.sorted <- signif(FDR.mean.sorted, digits=3)
FDR.median.sorted <- signif(FDR.median.sorted, digits=3)

#bias towards to the genes with bigger size
report <- data.frame(gs.names, Obs.ES, Obs.ES.norm, p.vals[,1], FDR.mean.sorted, FDR.median.sorted, sizes, arg.ESs)
names(report) <- c("GS", "ES", "NES", "NOM p-val", "FDR q-val",  "FDR (median)", "size", "arg.ES")
report = report[order(report[,"NOM p-val"], decreasing=F), ]
#report = subset(report, ES > 0)
#head(report, 20)
report$drug = drug

write.csv(report, paste("~/Documents/stanford/hcc/data/LIHC/moa/hepg2/", drug, "_sea_hepg2_t2.csv", sep=""))
}

#visualize one enrichment
target = "CHEMBL4128" #STAT3
target_score = subset(report, GS ==  target)
#target_cmpd_rank[target_cmpd_rank$target == target,]
drug.set = target_cmpd_rank$rank[target_cmpd_rank$target == target]
gsea_o = GSEA.EnrichmentScore(drug.list, drug.set, 0) 


#pdf(paste("~/Documents/stanford/hcc/data/LIHC/moa/drug_enrich/", drug, "_", target, ".pdf", sep="")) 
  layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE), heights= c(1,10))
  par(mar=c(0, 4, 2, 0.5))
  matched = rep(2, length(unique(target_cmpd_rank$pert_iname)))
  matched[drug.set] = 1

  image(as.matrix(matched, nrow=1, ncol=length(matched)), col=  redgreen(100), axes=F, srt=45)
  par(mar=c(12, 4, 2, 0.5))
  plot(x=c(1:length(gsea_o$RES)), y = gsea_o$RES, type = "l" , xlab = "", 
       ylab = "score", col = "red", main = paste(target, "enrichment"))
  abline(v=target_score$arg.ES[1])
  usr <- par( "usr" )
  #text( usr[ 2 ] , usr[ 4 ],  paste("ES =", target_score$ES[1]), adj = c( 3, 2 ), col = "blue" )
  #text( usr[ 2 ], usr[ 4 ], paste("NES =", target_score$NES[1]), adj = c( 3, 2), col = "blue" )
  #text( usr[ 2 ], usr[ 4 ], paste("q value =", target_score[1, "FDR q-val"]), adj = c( 3, 2 ), col = "blue" )
  legend("topright", legend=c(paste("ES =", target_score$ES[1]),
                              paste("NES =", target_score$NES[1]),
                              paste("q value =", target_score[1, "FDR q-val"])))
#dev.off()

#gefitinib: EGFR ranked 1st using chembl while ranked >100 using sticth