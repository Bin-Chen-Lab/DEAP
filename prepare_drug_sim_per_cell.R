
#this code is used to build compound similarity matrix for a given cell line. 
#it will generate a file for each compound
#make sure the compounds in one pair should be tested under the same condition

cell_lines = c("HEPG2")
drugs = c("gefitinib")

#load LINCS compound profiles
load("lincs_signatures_cmpd_landmark.RData")
drug_lincs_signatures = lincs_signatures

#LINCS compound meta data
lincs_sig_info = read.csv("lincs_sig_info.csv", stringsAsFactors=F)
#comment out if looking single cell lines
cell_lines = table(lincs_sig_info$cell_id[lincs_sig_info$pert_type == "trt_cp"])
cell_lines = names(tail(sort(cell_lines), 15))

#comment out if looking single cell lines
drugs = unique(lincs_sig_info$pert_iname[lincs_sig_info$cell_id %in% (cell_lines) & lincs_sig_info$pert_type == "trt_cp" & lincs_sig_info$is_gold==1]) 
#

for (cell_line in cell_lines){
for (drug in drugs ){
  print(drug)
  #find drug profiles
  pert_ids = unique(lincs_sig_info$pert_id[lincs_sig_info$pert_iname == drug])
  
  drug_sigs = unique(subset(lincs_sig_info, pert_iname == drug & is_gold == 1 & cell_id %in% c(cell_line)))
  
  if (nrow(drug_sigs)==0){next}
  
  rownames(drug_sigs) = drug_sigs$id
  
  #compute similarity between target profile and all other profiles
  drug_drug_cor = list()
  
  for (i in 1:nrow(drug_sigs)){
    lincs_sig_info_subset = lincs_sig_info[lincs_sig_info$cell_id == cell_line & lincs_sig_info$is_gold == 1 & lincs_sig_info$pert_type == "trt_cp"  ,]
    
    drug_lincs_signatures_subset = drug_lincs_signatures[, colnames(drug_lincs_signatures) %in%  lincs_sig_info_subset$id]
    cors = cor(drug_lincs_signatures_subset[, colnames(drug_lincs_signatures_subset) == drug_sigs$id[i]], drug_lincs_signatures_subset[, colnames(drug_lincs_signatures_subset) != drug_sigs$id[i]], method="pearson")
    drug_cor = cors[1,]
    drug_cor = data.frame(id = names(drug_cor), cor = drug_cor)
    drug_cor = merge(drug_cor, subset(lincs_sig_info_subset, select=c("id", "sig_id","pert_iname")), by="id")
    drug_cor = cbind(drug_cor, cell_line, sig_id.1 = drug_sigs$sig_id[i])
    drug_drug_cor[[i]] = drug_cor
  }
  drug_drug_cor_list=list(cor = drug_drug_cor, cell = drug_sigs$cell_id)

  #find the most correlated drug
  drug_drug_cor_frame = data.frame()
  for (i in 1:length(drug_drug_cor)){
    drug_drug_cor_frame = rbind(drug_drug_cor_frame, drug_drug_cor[[i]])
  }
  
  if (nrow(drug_drug_cor_frame) == 0 ){next}
  
  drug_drug_cor_frame$sig_id_time = sapply(drug_drug_cor_frame$sig_id, function(id){
    unlist(strsplit(unlist(strsplit(id, ":"))[1], "_"))[3]
  })
  drug_drug_cor_frame$sig_id1_time = sapply(as.character(drug_drug_cor_frame$sig_id.1), function(id){
    unlist(strsplit(unlist(strsplit(id, ":"))[1], "_"))[3]
  })
  
  drug_drug_cor_frame$sig_id_dose = sapply(drug_drug_cor_frame$sig_id, function(id){
    round(as.numeric(unlist(strsplit(id, ":"))[3]),1)
  })
  drug_drug_cor_frame$sig_id1_dose = sapply(as.character(drug_drug_cor_frame$sig_id.1), function(id){
    round(as.numeric(unlist(strsplit(id, ":"))[3]),1)
  })
  
  #compare similarity under the same treatment condition (time, duration)
  drug_drug_cor_frame_subset = subset(drug_drug_cor_frame, pert_iname != drug & sig_id_time == sig_id1_time & sig_id_dose == sig_id1_dose)
  if (nrow(drug_drug_cor_frame_subset) ==0 ) {next}
  
  #take median, if one drug have two profiles
  pert_iname_cor = aggregate(cor ~ pert_iname, drug_drug_cor_frame_subset, median)
  
  pert_iname_cor = pert_iname_cor[order(pert_iname_cor$cor, decreasing=T), ]

  write.csv(pert_iname_cor, paste("drug_sim/", drug, "_sim_drugs_", cell_line, "_cells.csv", sep=""))
}
}

