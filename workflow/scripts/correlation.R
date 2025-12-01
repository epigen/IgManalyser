snakemake@source("init.R")
snakemake@source("QN_functions.R")

#save.image("correlation.RData")


array_df <- read.csv(snakemake@input[['array_subset']], row.names = 1)
colnames(array_df) <- gsub("X", "", colnames(array_df))

metadata_df <- read.csv(snakemake@input[['metadata_file']], row.names = 1)

outpath <- snakemake@params[['project_dir']]
dir.create(outpath, recursive = TRUE)

params_to_plot <- snakemake@params[["param_columns"]]
params_to_plot <- as.numeric(strsplit(params_to_plot, ",")[[1]])
ldl_col <- as.numeric(snakemake@config[['ldl_col']])
gensini_col <- as.numeric(snakemake@config[['gensini_score']])

params <- colnames(metadata_df)[6:ncol(metadata_df)]
params <- params[params != snakemake@params[["split_group"]]]


row.names(array_df) <- array_df$coordinate_unique
peptide_annot <- array_df[,c("coordinate_unique", "sequence_label")]
peptide_annot$sequence_label_short <- sapply(peptide_annot$sequence_label, 
                                    function(x) strsplit(x, "-")[[1]][1])
array_df <- array_df[,1:(ncol(array_df)-2)]
#if(is.numeric(metadata_df$Sample_full_ID[[1]])){
#    metadata_df$Sample_full_ID <- paste0("X", metadata_df$Sample_full_ID)   
#}

metadata_df$Sample_full_ID <- as.character(metadata_df$Sample_full_ID)                                             
### calculating correlations
corcoef_df <- data.table()
for(peptide_id in row.names(array_df)){
    v1 <- as.numeric(t(array_df[peptide_id, metadata_df$Sample_full_ID])) ### making sure it's in the same order
    for(cname in params){
        if(is.numeric(metadata_df[,cname])){
            v2 <- as.numeric(metadata_df[,cname])
            corcoef <- cor.test(v2, v1, use="complete.obs")
            corcoef_df <- rbind(corcoef_df, 
                    data.frame("peptide" = peptide_id, "column" = cname, 
                               "corr" = corcoef$estimate[["cor"]], 
                               "pval" = corcoef$p.value), fill = TRUE)
        }
        
    }
}


## plotting full heatmap with all params:
corr_df <- corcoef_df %>% select(corr, peptide, column)  %>% pivot_wider(values_from = corr, names_from = peptide)
corr_df <- as.data.frame(corr_df)
row.names(corr_df) <- corr_df$column
corr_df <- corr_df[-1]
corr_df[is.na(corr_df)] <- 0
pheatmap(corr_df, 
         fontsize_col = 3,
         width = 25,
         height = 5,
         filename = file.path(outpath,  "corr_heatmap_full.png")
        )
write.csv(corr_df, snakemake@output[["corr_df"]])

## now subselecting for only "important" parameters
corcoef_df_imp <- corcoef_df %>% filter(column %in% colnames(metadata_df)[params_to_plot]) %>% 
    mutate(p_adj = p.adjust(pval, method = "BH")) %>% mutate(adjpval_label = ifelse(p_adj < 0.05, "*", ""))
corcoef_df_imp <- left_join(corcoef_df_imp, peptide_annot, by = c("peptide" = "coordinate_unique")) #so we have also labels

write.csv(corcoef_df_imp, snakemake@output[["corr_df_imp"]])

corr_df <- corcoef_df_imp %>%  select(corr, peptide, column)  %>% pivot_wider(values_from = corr, names_from = peptide)
corr_df <- as.data.frame(corr_df)
row.names(corr_df) <- corr_df$column
corr_df <- corr_df[-1]
pheatmap(corr_df, 
         fontsize_col = 3,
         width = 25,
         height = 5,
         filename = file.path(outpath,  "corr_heatmap_imp_params.png")
        )

### show only highly-correlated
peptides_higly_corr <- corcoef_df_imp %>% filter(column!= colnames(metadata_df)[ldl_col]) %>% filter(abs(corr) > 0.5) %>% pull(peptide)
peptides_higly_corr <- unique(peptides_higly_corr)
print(paste0("detected highly-correlated peptides: ", length(peptides_higly_corr)))

if(length(peptides_higly_corr) > 2){
pheatmap(corr_df[, peptides_higly_corr], 
         labels_col = as.character(peptide_annot[colnames(corr_df[, peptides_higly_corr]),]$sequence_label),
         fontsize_row = 15,
         fontsize_col = 15,
         cutree_row = 1,
         cutree_col = 1,
         width = 30,
         height = 12,
         filename = file.path(outpath,  "corr_heatmap_imp_params_highly_corr.pdf")
        )

pheatmap(corr_df[, peptides_higly_corr], 
         labels_col = as.character(peptide_annot[colnames(corr_df[, peptides_higly_corr]),]$sequence_label_short),
         fontsize_row = 15,
         fontsize_col = 15,
         cutree_row = 1,
         cutree_col = 1,
         width = 30,
         height = 12,
         filename = file.path(outpath,  "corr_heatmap_imp_params_highly_corr_sequence.pdf")
        )


pheatmap(corr_df[, peptides_higly_corr], 
         fontsize_row = 15,
         fontsize_col = 15,
         cutree_row = 1,
         cutree_col = 1,
         width = 30,
         height = 12,
         filename = file.path(outpath,  "corr_heatmap_imp_params_highly_corr_nolabel.pdf")
        )
}
                                             
### show only highly-correlated without the gensini score
peptides_higly_corr <- corcoef_df_imp %>% filter((column!= colnames(metadata_df)[ldl_col]) & (column!= colnames(metadata_df)[gensini_col])) %>% filter(abs(corr) > 0.5) %>% pull(peptide)

peptides_higly_corr <- unique(peptides_higly_corr)
print(paste0("detected highly-correlated peptides: ", length(peptides_higly_corr)))
if(length(peptides_higly_corr) > 2){
pheatmap(corr_df[row.names(corr_df)!=colnames(metadata_df)[gensini_col], peptides_higly_corr], 
         labels_col = as.character(peptide_annot[colnames(corr_df[row.names(corr_df)!=colnames(metadata_df)[gensini_col], peptides_higly_corr]),]$sequence_label),
         fontsize_row = 15,
         fontsize_col = 15,
         cutree_row = 1,
         cutree_col = 1,
         width = 30,
         height = 12,
         filename = file.path(outpath,  "corr_heatmap_imp_params_highly_corr_noGensini.pdf")
        )

pheatmap(corr_df[row.names(corr_df)!=colnames(metadata_df)[gensini_col], peptides_higly_corr], 
         labels_col = as.character(peptide_annot[colnames(corr_df[row.names(corr_df)!=colnames(metadata_df)[gensini_col], peptides_higly_corr]),]$sequence_label_short),
         fontsize_row = 15,
         fontsize_col = 15,
         cutree_row = 1,
         cutree_col = 1,
         width = 30,
         height = 12,
         filename = file.path(outpath,  "corr_heatmap_imp_params_highly_corr_sequence_noGensini.pdf")
        )


pheatmap(corr_df[row.names(corr_df)!=colnames(metadata_df)[gensini_col], peptides_higly_corr], 
         fontsize_row = 15,
         fontsize_col = 15,
         cutree_row = 1,
         cutree_col = 1,
         width = 30,
         height = 12,
         filename = file.path(outpath,  "corr_heatmap_imp_params_highly_corr_noGensini_nolabel.pdf")
        )
}
