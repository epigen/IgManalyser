snakemake@source("init.R")
snakemake@source("plot_functions.R")
library(ggpubr)
#save.image("plots.RData")

filepath <- snakemake@input[["array_subset"]]
annotpath <- snakemake@input[["annot_df"]]
corr_path <- snakemake@input[["corr_df"]]
volcano_results <- snakemake@input[["volcano_df"]]

outdir <- snakemake@params[["project_dir"]]
dir.create(outdir)
corr_thr <- as.numeric(snakemake@params[["CORR_THR"]])
lfc_thr <- as.numeric(snakemake@params[["LFC_THR"]])

gensini_thr <- snakemake@params[["gensini_thr"]]
print(gensini_thr)

if(is.na(gensini_thr) == FALSE){
    gensini_thr <- as.numeric(strsplit(gensini_thr, ",")[[1]])
}
print(gensini_thr)



volcano_df <- fread(volcano_results)
corr_df <- fread(corr_path)
array_df <- fread(filepath, header = TRUE)


##annotation to define phenotypes:
annot <- fread(annotpath)
annot$Sample_full_ID <- as.character(annot$Sample_full_ID)
annot_full <- fread(snakemake@params[['metadata']])
gensini_col_id <- snakemake@params[['gensini_col']]
gensini_col <- colnames(annot_full)[gensini_col_id]
gensini_col <- gsub(" ", "_", gensini_col)
cols_to_add <- c("Sample_full_ID", gensini_col)
annot_subs <- annot[,c()]
sub_col <- annot %>% pull(.data[[gensini_col]])
if(is.na(gensini_thr)[1] == FALSE){
    sub_col[sub_col < gensini_thr[2] & sub_col > gensini_thr[1]] <- NA
    sub_col[sub_col <= gensini_thr[1]] <- 0
    sub_col[sub_col >= gensini_thr[2]] <- 1
}

sub_col[is.na(sub_col)] <- -1

group1_ids <- annot$Sample_full_ID[sub_col == unique(annot[[gensini_col]])[0]]
group2_ids <- annot$Sample_full_ID[sub_col == unique(annot[[gensini_col]])[1]]

ids <- c(group1_ids, group2_ids)


array_df_long <- array_df[,-1] %>% 
    pivot_longer(!coordinate_unique	 & !sequence_label ) %>%
    left_join(annot[,..cols_to_add], by = c("name" = "Sample_full_ID")) %>% 
    filter(!is.na(.data[[gensini_col]]))
array_df_long$sequence_label_short <- sapply(array_df_long$sequence_label, 
      function(x) strsplit(x, "-")[[1]][1])

top_corr_df <- corr_df %>% filter(corr > corr_thr)

for(group in unique(top_corr_df$column)){
    print(group)
    pept_list <- top_corr_df[column == group,]$peptide
    plot_peptides_boxplots(array_df_long, pept_list, 
                           outdir, paste0("highly_corr_", group), gensini_col)
}


top_anticorr_df <- corr_df %>% filter(corr < -corr_thr)

for(group in unique(top_anticorr_df$column)){
    print(group)
    pept_list <- top_anticorr_df[column == group,]$peptide
    if(length(pept_list) > 0){
    plot_peptides_boxplots(array_df_long, pept_list, 
                           outdir, paste0("highly_anticorr_", group), gensini_col)
    }
}

if(NROW(volcano_df) > 1){                                             
diff_low <- volcano_df %>% filter(pval < 0.05 & log2foldchange < -lfc_thr)
                                             
plot_peptides_boxplots(array_df_long, diff_low$coordinate_unique, 
                           outdir, paste0("volcano_down"), gensini_col)

diff_high <- volcano_df %>% filter(pval < 0.05 & log2foldchange > lfc_thr)
plot_peptides_boxplots(array_df_long, diff_high$coordinate_unique, 
                           outdir, paste0("volcano_up"), gensini_col)
###combine top hits
diff_high$group <- "volcano_up"
diff_low$group <- "volcano_down"
                                            
}else{
     diff_low <- data.frame("coordinate_unique" = numeric(), "sequence_label" = numeric(),
                        "group" = numeric(), "pval" =numeric())
     diff_high <- data.table("coordinate_unique" = numeric(), "sequence_label" = numeric(),
                        "group" = numeric(), "pval" =numeric())
}
                                             
colnames(top_corr_df)[2] <- "coordinate_unique"
colnames(top_anticorr_df)[2] <- "coordinate_unique"                                           
                                             
top_corr_df$group <- paste0("high_corr_", top_corr_df$column)
top_anticorr_df$group <- paste0("high_anticorr_", top_anticorr_df$column)                                              
top_peptide_df <- rbind(rbind(diff_high[, c("coordinate_unique", "sequence_label", "group", "pval")], 
      diff_low[, c("coordinate_unique", "sequence_label", "group", "pval")]), 
      rbind(top_corr_df[, c("coordinate_unique", "sequence_label", "group", "pval")], 
      top_anticorr_df[, c("coordinate_unique", "sequence_label", "group", "pval")]))
write.csv(top_peptide_df, snakemake@output[["pept_hits"]])                                             
                                             
