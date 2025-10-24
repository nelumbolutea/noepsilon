
################r dosage sensitivity########################

#1. Infer dosage sensitivity of an OG by Pearson correlation (R)#
#R https://www.r-project.org/#

##1.1 Input observed copy numbers (e.g. column 1~28) and expected copy numbers (e.g. column 29~56) of each OG (row)#
setwd("/YourWorkFolder/")
data=read.delim("dosage.obs.exp.txt",header=T, row.names = 1)
  
##1.2 Initialize a matrix to store the results  
results <- matrix(NA, nrow = nrow(data), ncol = 2)  
colnames(results) <- c("Pearson_r", "p_value")  
  
##1.3 Loop through each row  
for (i in 1:nrow(data)) {    # Extract the relevant columns for this row  
  x <- as.numeric(data[i, 1:28]) # 
  y <- as.numeric(data[i, 29:56]) # 
     
##1.4 Use cor.test() for Pearson correlation##

pearson_test <- cor.test(x, y, method = "pearson")  

pearson_cor <- pearson_test$estimate  

pearson_pvalue <- pearson_test$p.value  
    
##1.5 Store the Pearson r and p-value  
  results[i, "Pearson_r"] <- pearson_cor # Extract the Pearson r value  
  results[i, "p_value"] <- pearson_pvalue # Extract the p-value  
}  
  
##1.6 Convert the results matrix to a data frame and add row names  
results_df <- data.frame(results, row.names = rownames(data))  

##1.7 Write the results to a CSV file  
write.csv(results_df, "dosage.obs.exp,correlation_results.csv")

#########################end###########################################




################WGD by Gene tree and Species tree Reconcilation########################

#2. Phylogenetic placement of duplication events by gene tree and species tree reconciliation (Linux)#
#Notung https://www.cs.cmu.edu/~durand/notung/ #

##2.1. save all gene tree file names to a file##
##Note that gene IDs should have species prefix, e.g. nnu_geneID, ata_geneID, etc. Please see Notung guide###
ls OG*.nwk > z_use_OG_trees

##2.2. root the gene tree by Notung##
##speciesTreeBig.nwk is species file##
##Please use prefix to represent the species names, such as nnu, ata etc.##
java -jar Notung-3.0-beta.jar  -s speciesTreeBig.nwk -b z_use_OG_trees --speciestag prefix  --treeoutput newick --events  --root --outputdir z_OGrooted_iq_nobs

##2.3. reconcile the gene tree given node support threshold 0.9##
ls z_OGrooted_iq_nobs/*rooting.0>z_rooted

java -jar Notung-3.0-beta.jar  -s speciesTreeBig.nwk -b z_rooted --speciestag prefix  --treeoutput newick --events  --reconcile --outputdir z_OGreconciled_iq_nobs --threshold 0.9 --edgeweights nhx --saveweakedgespng

##2.4. remove node supports as they occupied node IDs##
mkdir  z_OGreconciled_iq_nobs_final

cp  z_OGreconciled_iq_nobs/*.nwk.rooting.0.reconciled  z_OGreconciled_iq_nobs_final

find z_OGreconciled_iq_nobs_final/ -name "*.nwk.rooting.0.reconciled" -exec sed -i 's/)[0-9]*\.[0-9]*/)/g' {} +

##2.5. name the nodes of gene trees by Notung ##
cd z_OGreconciled_iq_nobs_final/

ls *.fasttree.nwk.rooting.0.reconciled>z_nobs

cp ../speciesTreeBig.nwk .

java -jar Notung-3.0-beta.jar -s speciesTreeBig.nwk  -b z_nobs --speciestag prefix  --treeoutput newick --events --root --outputdir ../z_OGfinal_iq_nobs

##2.6. obtain phylogenetic placement information of gene duplication nodes####
cd ..

find z_OGfinal_iq_nobs/ -name "*.events.txt" -exec grep -P "^n\d" {} + > z_OGfinal_dup_placement_iq


###2.7. obtain child nodes obtain###############################

find z_OGfinal_iq_nobs/ -name "*.rooting.0"  > z_finallist_iq_nobs

perl z_bulk_childnodes.pl z_finallist_iq_nobs

##2.8 put all tee nodes into one file###
cd z_OGfinal_iq_nobs

grep "$" *0.nodes>z_all.nodes

sed -i 's/\.fasttree\.nwk\.rooting\.0\.reconciled\.rooting\.0\.nodes//g' z_all.nodes

cd ..

##2.9. finally get all duplication nodes with phylogenetic placement of a target species, e.g. atrSC###
grep -P 'atrSC_\S*\t\S*atrSC_' z_OGfinal_iq_nobs/z_all.nodes >z_atrSC.nodes
##########################end#########################################




################WGD by WGDgc estimation########################

#3. Gene-count-based estimation of WGD retention rate (q), likelihood of a WGD model (R language)#
#WGDgc, manual in https://github.com/cecileane/WGDgc/tree/master/man#
 
##3.1. install WGDgc##
setwd("/YourWorkFolder/")

if (!requireNamespace("WGDgc", quietly = TRUE)) {
  devtools::install_github("cecileane/WGDgc")  
}

library(WGDgc)

##3.2. input gene count data of OGs##
gene_count_data <- read.delim("OGcount.txt", sep = "\t", check.names = FALSE,header=TRUE)

##3.3. input simmap format species tree##

tre_with_wgd.string <- "Simmap_format_speciestree"

tree_with_wgd <- read.simmap(text = tre_with_wgd.string)  


##3.4. MLE estimation##
result_with_wgd <- MLEGeneCount(
  tr = tree_with_wgd,
  geneCountData = gene_count_data ,
  geomMean = 1.240687, #Input geometric mean of ancestral copy numbers of OGs, e.g. 1.240687
  conditioning = "oneOrMore",
  fixedRetentionRates = FALSE 
)



##3.5. Output likelihood estimation##
cat(
  "Birth Rate (λ):", result_with_wgd$birthrate, "\n",
  "Death Rate (μ):", result_with_wgd$deathrate, "\n",
  "Log-Likelihood:", result_with_wgd$loglikelihood, "\n",
  "Convergence:", result_with_wgd$convergence, "\n",
  "mMax:", result_with_wgd$mMax, "\n",
  file = "result_with_wgd_summary.txt"
)

write.csv(result_with_wgd$WGDtable, "result_with_wgd_WGD_retention_rates.csv", row.names = FALSE)

write.csv(result_with_wgd$phyloMat, "result_with_wgd_phylogenetic_matrix.csv", row.names = FALSE)
###########################end########################################







