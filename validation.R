#validate DEAP using gold standard (not really gold standard, as the targets are barely known for many drugs in a particular cell line)

library("ROCR")

files = list.files( "drug_targets")

drug_targets = data.frame()
for (file in files){
  drug_target = read.csv(paste("drug_targets/", file, sep=""))
  drug = unlist(strsplit(file, "_"))[1]
  cell_line = unlist(strsplit(file, "_"))[3]
  target_rank = rank(-1 * drug_target$NES)
  drug_targets = rbind(drug_targets, data.frame(drug_target, drug, cell_line, target_rank))
}

#need to match drug targets....
chembl_targets = read.csv("chembl_targets.csv")

drug_targets = merge(drug_targets, chembl_targets, by.x="GS", by.y="chembl_id")
#
drugbank_targets = read.csv("drugbank_drug_targets_V2.csv")

drugbank_targets = subset(drugbank_targets, tolower(name) %in% drug_targets$drug)
drugbank_targets = subset(drugbank_targets, type == "target", select=c("name", "gene_symbol", "action"))
drugbank_targets$name = tolower(drugbank_targets$name)
drugbank_targets$target = 1

drug_targets_validate = merge(drug_targets,drugbank_targets, by.x=c("drug", "Symbol"), by.y=c("name", "gene_symbol"), all.x=T )
drug_targets_validate$target[is.na(drug_targets_validate$target)] = 0

pred = prediction(drug_targets_validate$NES , drug_targets_validate$target)
performance(pred, "auc")
plot( performance( pred, "tpr", "fpr" ) , col="red", add=T)

#merage base cell lines
drug_targets_validate_aggre = aggregate(NES ~ drug + Symbol + target, drug_targets_validate, mean )
pred = prediction(drug_targets_validate_aggre$NES , drug_targets_validate_aggre$target)
performance(pred, "auc")


drug_targets_validate_subset = subset(drug_targets_validate, target == 1)
write.csv(drug_targets_validate_subset, "drug_targets_validate_subset.csv")

targets_for_one_drug = subset(drug_targets_validate, drug == "niclosamide")
targets_for_one_drug = targets_for_one_drug[order(targets_for_one_drug$FDR.q.val), ]
head(targets_for_one_drug, 20)

drug_targets_validate = drug_targets_validate[order(drug_targets_validate$arg.ES), ]
head(drug_targets_validate, 100)

tail(sort(table(drug_targets_validate$GS[drug_targets_validate$arg.ES <0.05 & drug_targets_validate$NES])), 30)

drugs = unique(drug_targets_validate$drug[drug_targets_validate$target == 1 & drug_targets_validate$NES>0])
ps = data.frame()
for (drug in drugs){
  drug_targets_validate_single = drug_targets_validate[drug_targets_validate$drug == drug ,]
  cell_lines = unique(drug_targets_validate_single$cell_line)
  for (cell_line in cell_lines){
    drug_targets_validate_single_cell_line = drug_targets_validate_single[drug_targets_validate_single$cell_line == cell_line, ]
    if (sum(drug_targets_validate_single_cell_line$target == 1) > 0){
    pred = prediction(drug_targets_validate_single_cell_line$NES , drug_targets_validate_single_cell_line$target)
    perf = performance(pred, "auc")
    ps = rbind(ps, data.frame(drug, cell_line, p=perf@y.values[[1]]))
    
    }
  }
}
#confirmed: NES lead to the best prediction

sort(by(ps$p, ps$drug, mean))
sort(by(ps$p, ps$cell_line, mean))

NEN = ps
ES = ps

