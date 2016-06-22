#build target_cell line for each drug

library(pheatmap)

drugs = "niclosamide"
files = list.files("drug_targets")
drugs = sapply(files, function(drug_name){
  unlist(strsplit(unlist(strsplit(drug_name, "_"))[1], "\\."))[1]
})
drugs = unique(drugs)

for (drug in drugs){
files = list.files("drug_targets", pattern = drug)
cell_lines = sapply(files, function(drug_name){
  unlist(strsplit(unlist(strsplit(drug_name, "_"))[3], "\\."))[1]
})

lincs_cmpd_sea_pred = read.csv("lincs.cond.sea_res.csv")
lincs_cmpd_sea_pred_subset = unique(lincs_cmpd_sea_pred[, c("uniprot_id", "targ_desc", "chembl_id")])

preds = data.frame()
for (cell_line in cell_lines){
  pred1 = read.csv(paste("drug_targets/", drug, "_sea_", cell_line, ".csv", sep=""), stringsAsFactors = F)
  pred1$cell_line = cell_line
  pred1$NES = scale(pred1$NES)
  preds = rbind(preds, pred1)
}

#targets 

preds = merge(preds, lincs_cmpd_sea_pred_subset, by.x="GS", by.y="chembl_id")

hits = subset(preds, ES>0 & FDR.q.val<0.05)
targets = table(hits$targ_desc)
targets = names(targets[targets>1])

target_cell_matrix = matrix(NA, nrow=length(targets), ncol=length(cell_lines))
rownames(target_cell_matrix) = targets
colnames(target_cell_matrix) = cell_lines

preds_subset = subset(preds, targ_desc %in% targets)

if (nrow(preds_subset)>1 & length(targets)>1 & length(cell_lines)>1){
  for (i in 1:nrow(preds_subset)){
    target_cell_matrix[as.character(preds_subset$targ_desc[i]), as.character(preds_subset$cell_line[i])] = preds_subset$NES[i]
  }
  
  pheatmap(target_cell_matrix,cellheight = 10, cellwidth = 10, file=paste("drug/", drug, ".pdf", sep=""))
}

######################
#limited targtes from sea
lincs_cmpd = read.csv("lincs_cmpd_info.csv", stringsAsFactors = F)
drug_pert_id = lincs_cmpd$pert_id[lincs_cmpd$pert_iname == drug][1]
targets = unique(as.character(lincs_cmpd_sea_pred$targ_desc[lincs_cmpd_sea_pred$cpd_id == drug_pert_id]))

targets = targets[targets %in% preds$targ_desc]

target_cell_matrix = matrix(NA, nrow=length(targets), ncol=length(cell_lines))
rownames(target_cell_matrix) = targets
colnames(target_cell_matrix) = cell_lines

preds_subset = subset(preds, targ_desc %in% targets)
if (nrow(preds_subset)>1 & length(targets)>1 & length(cell_lines)>1){
  for (i in 1:nrow(preds_subset)){
    target_cell_matrix[as.character(preds_subset$targ_desc[i]), as.character(preds_subset$cell_line[i])] = preds_subset$NES[i]
  }
  
  pheatmap(target_cell_matrix,cellheight = 10, cellwidth = 10, file=paste("drug/", drug, "_sea_targets.pdf", sep=""))
}
}
