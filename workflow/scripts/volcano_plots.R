snakemake@source("init.R")
#save.image("volcano.RData")


filepath <- snakemake@input[["array_subset"]]
annotpath <- snakemake@input[["annot_df"]]

outdir <- snakemake@params[["project_dir"]]
dir.create(outdir, recursive = TRUE)

gensini_thr <- snakemake@params[["gensini_thr"]]
print(gensini_thr)

if(is.na(gensini_thr) == FALSE){
    gensini_thr <- as.numeric(strsplit(gensini_thr, ",")[[1]])
}
print(gensini_thr)
## reading in data for comparison:
#df <- fread(file.path("../processed_data", id, paste0(id, "_pool"), "summary", "full_avrg_paired.tsv"))
#df <- fread(file.path("../processed_data", species, id, "summary", "full_avrg.tsv"))
df <- fread(filepath, header = TRUE)
head(df)


##annotation to define phenotypes:
annot <- fread(annotpath)

annot_full <- fread(snakemake@params[['metadata']])
gensini_col_id <- snakemake@params[['gensini_col']]
genisin_col <- colnames(annot_full)[gensini_col_id]
genisin_col <- gsub(" ", "_", genisin_col)
print(genisin_col)

sub_col <- annot %>% pull(.data[[genisin_col]])

if(is.na(gensini_thr)[1] == FALSE){
    sub_col[sub_col < gensini_thr[2] & sub_col > gensini_thr[1]] <- NA
    sub_col[sub_col <= gensini_thr[1]] <- 0
    sub_col[sub_col >= gensini_thr[2]] <- 1
}

#sub_col[is.na(sub_col)] <- -1
annot <- annot[!is.na(sub_col)]
sub_col <- sub_col[!is.na(sub_col)]


group1_ids <- as.character(annot$Sample_full_ID[sub_col == sort(unique(sub_col))[1]])

group2_ids <- as.character(annot$Sample_full_ID[sub_col == sort(unique(sub_col))[2]])
ids <- c(group1_ids, group2_ids)
### 
if(length(group1_ids) <2 |length(group2_ids) <2 ){
    cat("not enough data for comparison", file = snakemake@output[["volcano_df"]])
}else{


#df_m <- df[, -c(24,25)]# %>% select(name, avrg,unique_peptide_id) %>%
#pivot_wider(names_from = name, values_from = avrg)
df_m <- as.data.frame(df[,..ids])
row.names(df_m) <- df$coordinate_unique


##exluding full zeros
df_m <- df_m[rowSums(df_m[,ids]) != 0,]
#print(head(df_m[, group1_ids]))
## simple t-test comparing rows
pval_list <- list()
for(i in row.names(df_m)){
    pval_list[[i]] <- t.test(as.numeric(df_m[i,group1_ids]),
     as.numeric(df_m[i,group2_ids]))$p.value
}

df_m$pval <- as.numeric(pval_list)
df_m$pval_adj <- p.adjust(df_m$pval, method = "fdr")

log2FC_list <- list()
for(i in seq(1, NROW(df_m))){
    log2FC_list[[i]] <- log2(mean(as.numeric(df_m[i,group2_ids]))) - log2(mean(as.numeric(df_m[i,group1_ids])))
}
df_m$log2foldchange <- as.numeric(log2FC_list)

df_m$coordinate_unique <- row.names(df_m)
add_cols <- c("sequence_label", "coordinate_unique")
df_m <- left_join(df_m, df[,..add_cols])

df_m$sequence_label_short <- sapply(df_m$sequence_label, 
                                    function(x) strsplit(x, "-")[[1]][1])
                                    
###pvalue
ggplot(df_m, aes(y = -log10(pval), x = log2foldchange)) + geom_point() + geom_hline(yintercept = -log10(0.05), linetype = "dashed") + geom_text_repel(data = df_m[df_m$pval < 0.05,], 
          aes(x = log2foldchange, y =-log10(pval), 
              label = coordinate_unique ))

ggsave(file.path(outdir, "volcano_ttest_not_adj.pdf"), width = 5, height = 5)

ggplot(df_m, aes(y = -log10(pval), x = log2foldchange)) + geom_point() + geom_hline(yintercept = -log10(0.05), linetype = "dashed") + geom_text_repel(data = df_m[df_m$pval < 0.05,], 
          aes(x = log2foldchange, y =-log10(pval), 
              label = sequence_label_short ), max.overlaps = 10)

ggsave(file.path(outdir, "volcano_ttest_not_adj_label.pdf"), width = 5, height = 5)


#####adjusted
ggplot(df_m, aes(y = -log10(pval_adj), x = log2foldchange)) +  geom_point() + 
geom_text_repel(data = df_m[df_m$pval_adj < 1,], 
          aes(x = log2foldchange, y =-log10(pval_adj), 
              label = coordinate_unique	 )) + geom_hline(yintercept = -log10(0.05), linetype = "dashed")

ggsave(file.path(outdir, "volcano_ttest_adjusted.pdf"), width = 5, height = 5)

ggplot(df_m, aes(y = -log10(pval_adj), x = log2foldchange)) +  geom_point() + 
geom_text_repel(data = df_m[df_m$pval_adj < 1,], 
          aes(x = log2foldchange, y =-log10(pval_adj), 
              label = sequence_label_short )) + geom_hline(yintercept = -log10(0.05), linetype = "dashed")

ggsave(file.path(outdir, "volcano_ttest_adjusted_labels.pdf"), width = 5, height = 5)


write.csv(df_m[, c("pval_adj","pval", "log2foldchange", "coordinate_unique", "sequence_label", "sequence_label_short")], snakemake@output[["volcano_df"]])

}
