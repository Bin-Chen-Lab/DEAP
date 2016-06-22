#annotate each LINCS drugs
#using mesh terms, targets, etc
setwd("/Users/binchen1/Documents/stanford/hcc/code/DEAP/")

library(stringr)

lincs_cmpd_info = read.csv("data/lincs_cmpd_info.csv")
lincs_cmpd_info = unique(subset(lincs_cmpd_info, select=c("pert_iname", "pubchem_cid", "inchikey")))


#prepare targets from SEA
lincs_cmpd_pert_id_smiles = read.csv("data/lincs_cmpd_info.csv")
lincs_cmpd_pert_id_smiles = unique(subset(lincs_cmpd_pert_id_smiles, select=c("pert_id", "canonical_smiles")))

#since LINCS has so many isoforms, SEA only predicts one, need to match back through smiles
lincs_cmpd_sea_pred = read.csv("data/lincs.cond.sea_res.csv")
lincs_cmpd_sea_pred = merge(lincs_cmpd_sea_pred, lincs_cmpd_pert_id_smiles, by.x="cpd_id", by.y="pert_id")

lincs_cmpd = read.csv("data/lincs_cmpd_info.csv")
lincs_pert_iname_sea = merge(lincs_cmpd, lincs_cmpd_sea_pred, by.x="canonical_smiles", by.y="canonical_smiles")

lincs_pert_iname_sea_merged = unique(subset(lincs_pert_iname_sea, select=c("pert_iname", "chembl_id")))
lincs_pert_iname_sea_merged$pert_iname = tolower(lincs_pert_iname_sea_merged$pert_iname)
sea_pert_iname_targets = tapply(lincs_pert_iname_sea_merged$chembl_id, lincs_pert_iname_sea_merged$pert_iname, function(symbols){paste(symbols, collapse="\t")})
sea_pert_iname_targets = data.frame(pert_iname = tolower(names(sea_pert_iname_targets)), sea_targets = sea_pert_iname_targets)

#prepare targets from chembl
#map chembl cmpds to LINCS by inchikey
chembl_drug_targets = read.csv("data/drug_targets_chembl_v20.csv")
chembl_drug_targets = subset(chembl_drug_targets, !is.na(Symbol))
chembl_pert_iname_targets = tapply(chembl_drug_targets$Symbol, chembl_drug_targets$pert_iname, function(symbols){paste(symbols, collapse="\t")})
chembl_pert_iname_targets = data.frame(pert_iname = tolower(names(chembl_pert_iname_targets)), chembl_targets = chembl_pert_iname_targets)

#prepare targets from CTD
#map CTD cmpds to LINCS by name
ctd_drugs = read.csv("~/Downloads/CTD_chemicals (1).csv", skip = 26, stringsAsFactors=F)

'ctd_drug_synonyms = data.frame()
for(i in 1:nrow(ctd_drugs)){
  print(i)
  if (ctd_drugs$Synonyms[i] != ""){
    write.table(data.frame(chemicalID = ctd_drugs$ChemicalID[i], synonyms= unlist(strsplit(as.character(ctd_drugs$Synonyms[i]), "\\|"))), "~/Downloads/ctd_synoyms.txt", append=T, sep="\t", quote=F, col.names=F, row.names=F)
  }
}'
#mapping synonyms
ctd_drug_synonyms = read.delim("~/Downloads/ctd_synoyms.txt", sep="\t", header=F)
ctd_drugs_subset = subset(ctd_drugs, select = c("ChemicalID", "X..ChemicalName"))
names(ctd_drugs_subset) = c("V1", "V2")
ctd_drug_synonyms = unique(rbind(ctd_drugs_subset, ctd_drug_synonyms))
ctd_drug_synonyms$V2 = tolower(ctd_drug_synonyms$V2)
lincs_cmpd_info$pert_iname_lower = tolower(lincs_cmpd_info$pert_iname)
ctd_drug_synonyms_lincs = merge(ctd_drug_synonyms, lincs_cmpd_info, by.x="V2", by.y="pert_iname_lower")

ctd_drug_targets = read.csv("~/Downloads/CTD_chem_gene_ixns (1).csv", skip = 26 )
ctd_drug_targets = unique(subset(ctd_drug_targets, !is.na(GeneSymbol) , select=c("ChemicalID", "GeneSymbol")))
names(ctd_drug_targets) = c("ChemicalID", "Symbol")
ctd_drug_targets$mesh = paste("MESH:", ctd_drug_targets$ChemicalID, sep="")

ctd_drug_targets = merge(ctd_drug_targets, ctd_drug_synonyms_lincs, by.x="mesh", by.y="V1")
ctd_drug_targets = unique(subset(ctd_drug_targets, select=c("pert_iname", "Symbol")))

ctd_pert_iname_targets = tapply(ctd_drug_targets$Symbol, ctd_drug_targets$pert_iname, function(symbols){paste(symbols, collapse="\t")})
ctd_pert_iname_targets = data.frame(pert_iname = tolower(names(ctd_pert_iname_targets)), ctd_targets = ctd_pert_iname_targets)
ctd_pert_iname_targets = unique(subset(ctd_pert_iname_targets, !is.na(ctd_targets)))


#prepare targets from drugbank
#map drugbank to LINCS by inchikey
drugbank_drug_targets = read.csv("~/Desktop/drugbank_drug_targets_V2.csv")
drugbank_drug_targets$inchikey_new = sapply(as.character(drugbank_drug_targets$inchikey), function(inchikey){
  items = unlist(strsplit(inchikey, "")) 
  if (length(items)>0){
    as.character(paste(items[10:length(items)], collapse=""))
  }else{
    ""
  }
})
drugbank_drug_targets = unique(subset(drugbank_drug_targets, !is.na(gene_symbol) & inchikey_new %in% lincs_cmpd_info$inchikey, select=c("name", "inchikey_new", "gene_symbol")))
names(drugbank_drug_targets) = c("name","inchikey", "Symbol")
drugbank_drug_targets = merge(drugbank_drug_targets, lincs_cmpd_info, by="inchikey")
drugbank_pert_iname_targets = tapply(drugbank_drug_targets$Symbol, drugbank_drug_targets$pert_iname, function(symbols){paste(unique(symbols), collapse="\t")})
drugbank_pert_iname_targets = data.frame(pert_iname = tolower(names(drugbank_pert_iname_targets)), drugbank_targets = drugbank_pert_iname_targets)
drugbank_pert_iname_targets = unique(subset(drugbank_pert_iname_targets, !is.na(drugbank_targets)))

#from stitch
# stitch_lincs_cmpd = read.csv("~/Downloads/lincs_in_stitch.csv")
# stitch_cmpd_target = read.csv("~/Downloads/9606.protein_chemical.links.detailed.v4.0.tsv", sep="\t")
# stitch_cmpd_targets = merge(stitch_cmpd_target, stitch_lincs_cmpd, by.x="chemical", by.y="stereo_chemical_id")
# stitch_cmpd_targets = subset(stitch_cmpd_targets, combined_score>500)
# stitch_cmpd_targets$protein_human = sapply(as.character(stitch_cmpd_targets$protein), function(id){unlist(strsplit(id, "\\."))[2]})
# 
# ensembl_peptide_id = sapply(as.character(unique(stitch_cmpd_targets$protein)), function(id){unlist(strsplit(id, "\\."))[2]})
# 
# library(biomaRt)
# mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# ensembl_symbol <- getBM(attributes = c( "hgnc_symbol", "ensembl_peptide_id"),
#                  filters = "ensembl_peptide_id", values = ensembl_peptide_id,
#                  mart = mart)
# 
# stitch_cmpd_targets_symbol = merge(stitch_cmpd_targets, ensembl_symbol, by.x="protein_human", by.y="ensembl_peptide_id")
# stitch_drug_targets = unique(subset(stitch_cmpd_targets_symbol, select=c("pert_iname", "inchikey", "hgnc_symbol")))
# names(stitch_drug_targets) = c("pert_iname","inchikey", "Symbol")
# stitch_pert_iname_targets = tapply(stitch_drug_targets$Symbol, stitch_drug_targets$pert_iname, function(symbols){paste(unique(symbols), collapse="\t")})
# stitch_pert_iname_targets = data.frame(pert_iname = tolower(names(stitch_pert_iname_targets)), stitch_targets = stitch_pert_iname_targets)
# stitch_pert_iname_targets = unique(subset(stitch_pert_iname_targets, !is.na(stitch_targets)))
# stitch_pert_iname_targets$pert_iname = tolower(stitch_pert_iname_targets$pert_iname)

#prepare target class
cid_names = read.delim("data/CID-MeSH.txt", sep="\t", header=F)
name_mesh = read.delim("data/MeSH-Pharm.txt", sep="|", header=F, stringsAsFactors=F)
name_mesh_reformat = data.frame()
for (i in 1:nrow(name_mesh)){
  fields = unlist(strsplit(name_mesh$V1[i], "\t"))
  name_mesh_reformat = rbind(name_mesh_reformat, data.frame(name = fields[1], meshes = paste(fields[-c(1)], collapse="\t")))
}


lincs_cmpd_info_mesh = merge(lincs_cmpd_info, cid_names, by.x="pubchem_cid", by.y="V1", all.x=T)
lincs_cmpd_info_mesh = merge(lincs_cmpd_info_mesh, name_mesh_reformat, by.x= "V2", by.y="name", all.x=T)

lincs_cmpd_info_mesh_target = merge(lincs_cmpd_info_mesh, chembl_pert_iname_targets, by.x= "pert_iname_lower", by.y="pert_iname", all.x=T)
lincs_cmpd_info_mesh_target = merge(lincs_cmpd_info_mesh_target, ctd_pert_iname_targets, by.x= "pert_iname_lower", by.y="pert_iname", all.x=T)
lincs_cmpd_info_mesh_target = merge(lincs_cmpd_info_mesh_target, drugbank_pert_iname_targets, by.x= "pert_iname_lower", by.y="pert_iname", all.x=T)
#lincs_cmpd_info_mesh_target = merge(lincs_cmpd_info_mesh_target, stitch_pert_iname_targets, by.x= "pert_iname_lower", by.y="pert_iname", all.x=T)
lincs_cmpd_info_mesh_target = merge(lincs_cmpd_info_mesh_target, sea_pert_iname_targets, by.x= "pert_iname_lower", by.y="pert_iname", all.x=T)


# lincs_cmpd_info_mesh_target$combined = sapply(1:nrow(lincs_cmpd_info_mesh_target), function(id){
#   combined = paste(lincs_cmpd_info_mesh_target$chembl_targets[id], lincs_cmpd_info_mesh_target$drugbank_targets[id], lincs_cmpd_info_mesh_target$ctd_targets[id], lincs_cmpd_info_mesh_target$stitch_targets[id], sep="\t")
#   combined = unique(unlist(strsplit(combined, "\t")))
#   paste(unique(combined[!is.na(combined) & combined != "NA"]), collapse="\t")
# })


save(lincs_cmpd_info_mesh_target, file="data/lincs_cmpd_info_mesh_target_sea.RData")

#for cell line 
lincs_sig_info = read.csv("data/lincs_sig_info.csv")
pert_inames = unique(lincs_sig_info$pert_iname[lincs_sig_info$cell_id == "HEPG2" & lincs_sig_info$pert_type == "trt_cp" & lincs_sig_info$is_gold == 1])
lincs_cmpd_info_mesh_target_subset = subset(lincs_cmpd_info_mesh_target, pert_iname_lower %in% tolower(pert_inames))
write.csv(lincs_cmpd_info_mesh_target_subset, "data/lincs_cmpd_info_mesh_target_hepg2.csv")



