drug = "vorinostat"
cell_line = "HEPG2"
pred1 = read.csv(paste("drug_targets/", drug, "_sea_", cell_line, ".csv", sep=""))

lincs_cmpd_sea_pred = read.csv("lincs.cond.sea_res.csv")
lincs_cmpd_sea_pred = unique(lincs_cmpd_sea_pred[, c("uniprot_id", "targ_desc", "chembl_id")])

pred1 = merge(pred1, lincs_cmpd_sea_pred, by.x="GS", by.y="chembl_id")
pred1 = subset(pred1, ES>0)
pred1 = pred1[order(pred1$FDR.q.val ), ]
head(pred1)

