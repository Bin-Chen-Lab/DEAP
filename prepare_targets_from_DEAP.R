setwd("/Users/binchen1/Documents/stanford/hcc/code/DEAP/")

files = list.files("~/Documents/stanford/hcc/data/LIHC/moa/hepg2/", pattern = "hepg2")
drugs = sapply(files, function(drug_name){
  unlist(strsplit(drug_name, "_"))[1]
})

drug_targets = data.frame()
drug_targets_frame = data.frame()


for (drug in drugs){
  pred = read.csv(paste("~/Documents/stanford/hcc/data/LIHC/moa/hepg2/", drug, "_sea_hepg2_t1.csv", sep=""))
  pred_sig = pred[pred$NES>0 & pred$FDR.q.val< 0.05, ]
  if (nrow(pred_sig)>0){
    drug_targets = rbind(drug_targets, data.frame(pert_iname = tolower(drug), sea_targets_1 = paste(pred_sig$GS, collapse="\t")))
    drug_targets_frame = rbind(drug_targets_frame, data.frame(pert_iname = drug, pred_sig ))
  }else{
    drug_targets = rbind(drug_targets, data.frame(pert_iname = tolower(drug), sea_targets_1 = NA))
    #drug_targets_frame = pred_sig
  }
}
lincs_cmpd_info_mesh_target = drug_targets

#filter out SEA targets using DEAP
lincs_cmpd_sea_pred = read.csv("data/lincs.cond.sea_res.csv")
lincs_cmpd_sea_pred =  subset(lincs_cmpd_sea_pred, max_tc > 0.999)


save(lincs_cmpd_info_mesh_target, file="data/lincs_cmpd_info_mesh_target_sea_v1.RData")


lincs_cmpd_sea_pred = read.csv("data/lincs.cond.sea_res.csv")
lincs_cmpd_sea_pred = unique(lincs_cmpd_sea_pred[, c("uniprot_id", "targ_desc", "chembl_id")])
drug_targets_frame_anno = merge(drug_targets_frame, lincs_cmpd_sea_pred, by.x="GS", by.y="chembl_id")

tail(sort(table(drug_targets_frame_anno$targ_desc)), 30)
tail(sort(table(drug_targets_frame_anno$pert_iname)), 30)
drug_targets_frame_anno = drug_targets_frame_anno[order(drug_targets_frame_anno$GS), ]
tail(drug_targets_frame_anno, 30)

write.csv(drug_targets_frame_anno, "data/drug_targets_frame_anno_hepg2.csv")
