snakemake@source("init.R")
snakemake@source("QN_functions.R")

#save.image("correlation_per_peptide.RData")

split_subsets <- snakemake@params[['split_subsets']]
split_subsets <- strsplit(split_subsets, ",")[[1]]

params_to_plot <- snakemake@params[["param_columns"]]
params_to_plot <- as.numeric(strsplit(params_to_plot, ",")[[1]])

outpath <- snakemake@params[["project_dir"]]

dir.create(file.path(outpath, "paramOnAssay"))
dir.create(file.path(outpath, "corr_per_peptide"))



corr_df_list <- lapply(seq_along(snakemake@input[["corr_df_imp"]]), function(x) {df <- fread(snakemake@input[["corr_df_imp"]][x], header = TRUE); df$group <- split_subsets[x]; return(df)})
corr_df <- rbindlist(corr_df_list)

mtd_df_list <- lapply(seq_along(snakemake@input[["mtd_df"]]), function(x) {df <- read.csv(snakemake@input[["mtd_df"]][x], row.names = 1); df$group <- split_subsets[x]; return(df)})
mtd_df <- rbindlist(mtd_df_list)
mtd_df$Sample_full_ID <- as.character(mtd_df$Sample_full_ID)

array_df_list <- lapply(seq_along(snakemake@input[["array_subsets"]]), function(x) {df <- fread(snakemake@input[["array_subsets"]][x], header = TRUE); df <- df[,2:(ncol(df))]; return(df)})
array_df <- Reduce( left_join, array_df_list)
peptide_annot <- array_df[, c("coordinate_unique", "sequence_label")]

for(param_val in unique(corr_df$column)){
    corr_sub <- corr_df[column == param_val]
    
    corr_sub$library <- sapply(corr_sub$peptide, function(x) strsplit(x, "_")[[1]][3])
    corr_sub$row <- sapply(corr_sub$peptide, function(x) strsplit(x, "_")[[1]][1])
    corr_sub$col <- sapply(corr_sub$peptide, function(x) strsplit(x, "_")[[1]][2])

    corr_sub$col <- factor(corr_sub$col, levels = paste0("C", 1:24))
    corr_sub$row <- factor(corr_sub$row, levels = rev(paste0("L", 1:16)))

    ggplot(corr_sub, aes(y = row, x = col, color = corr)) + 
        geom_point(size = 5) + facet_grid(group~library) + 
        scale_color_distiller(palette = "RdBu", limits = c(-0.9, 0.9))
        ggtitle(param_val)
    ggsave(file.path(outpath, "paramOnAssay", paste0(param_val, ".pdf")), width = 14, height = 10)
}

#### make param-expression correlation plots
array_df <- as.data.frame(array_df)
row.names(array_df) <- array_df$coordinate_unique
array_df$coordinate_unique <- NULL       
array_df$sequence_label <- NULL      
 array_df <- as.data.frame(t(array_df))  
array_df$Sample_full_ID <- row.names(array_df) 
                           
select_col <- c(colnames(mtd_df)[params_to_plot], c("Sample_full_ID", "group"))                           
array_df_annot <- left_join(array_df, mtd_df[,..select_col])



for(peptide_id in colnames(array_df)[1:(ncol(array_df)-1)]){
    sel_cols <- c(colnames(mtd_df)[params_to_plot],"group", peptide_id)
    df_plot <- array_df_annot[,sel_cols]
    colnames(df_plot)[ncol(df_plot)] <- 'peptide'

    df_plot <- df_plot %>% pivot_longer(!(peptide|group))

    sub_cor_df <- corr_df[corr_df$peptide == peptide_id]

    sub_cor_df$label = paste0(round(sub_cor_df$corr, 2), "(", round(sub_cor_df$pval, 2), "/", round(sub_cor_df$p_adj, 2), ")")
    
    sub_cor_df$x = 0.8*max(df_plot$peptide)

    colnames(sub_cor_df)[3] <- "name"

    
    sub_cor_df <- left_join(sub_cor_df, df_plot %>% group_by(name) %>% summarise(y = 1.1*max(value, na.rm = TRUE)))
    
    sub_cor_df[sub_cor_df$group == split_subsets[[1]],]$y <- 1.1*sub_cor_df[sub_cor_df$group == split_subsets[[2]],]$y
    
    options(repr.plot.width = 10, repr.plot.height = 9)
    ggplot(df_plot, aes(x = peptide, y = value, color = as.factor(group))) + geom_point()  +
        geom_text(data = sub_cor_df, aes(x = x, y = y, label = label)) + facet_wrap(~name, scale = "free_y") + 
        scale_color_manual(values = setNames( c("#67a9cf","#ef8a62"), split_subsets)) + ggtitle(sub_cor_df$sequence_label[1])
    ggsave(file.path(outpath,"corr_per_peptide", paste0(strsplit(sub_cor_df$sequence_label[1], "-")[[1]][1], ".pdf")), width = 20, height = 20)
    }


##distribution of correlations
ggplot(corr_df, aes(x = corr)) + geom_histogram(color = "black", fill = "lightblue")+
    facet_grid(column~group)    
ggsave(snakemake@output[["corr_dist"]], width =6, height = 30)
