if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")

  BiocManager::install("preprocessCore", 
                     configure.args = "--disable-threading", 
                     force = TRUE,
                     type = "source")
}

#save.image("QN.RData")

snakemake@source("init.R")
snakemake@source("QN_functions.R")

library(preprocessCore)
dataPath <- snakemake@input[[1]]
metadataPath <- snakemake@params[["metadataPath"]]
annot <- fread(metadataPath)

group_by <- snakemake@params[["group_by"]]

outpath <- snakemake@params[['normPath']]
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)

if(group_by == 'NA'){
    annot$group_col = "group"} ##TODO: grouping by phenotype, etc
##TODO-make the color choice fixed
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
for(val in unique(annot$group_col)){
    if(val %in% names(group_colors)) print(paste0("color defined for ", val, ": ", group_colors[val]))
    else group_colors[val] <- color[sample(length(color),1)]
    }
group_colors <- group_colors[unique(annot$group_col)]
### defining group splits for the QN (i.e. if we have heathy and disease we would want to do the normalisation in each of the groups independantly)

group1_id <- unique(annot$group_col)[[1]]

if(length(unique(annot$group_col)) > 1){
group2_id <- unique(annot$group_col)[[2]]}else group2_id <- NA  ## TODO: add an option to have more than two groups?


full_df <- fread(dataPath)
full_df$coordinate_unique <- paste0(full_df$Coordinate, "_", full_df$library)

##removing the samples that failed
annot <- annot[sample_name %in% unique(full_df$sample_name),]
full_array <- full_df %>% select(coordinate_unique, avrg, sample_name ) %>%
    pivot_wider(names_from = sample_name, values_from = avrg)
full_array <- as.data.frame(full_array)
row.names(full_array) <- full_array$coordinate_unique


full_array <-full_array[as.character(annot$sample_name)]
group1 <- which(annot$group_col == group1_id)
group2 <- which(annot$group_col == group2_id) 
if(!is.na(group2_id)){
    full_array_norm <- get_QN_norm(full_array_L1, group1, group2)
}else{
    full_array_norm <- get_QN_norm_onegroup(full_array, group1)
}


#visualising ALL peptides
##TO DO add annotations
full_array_norm$coordinate_unique <- row.names(full_array_norm)
full_array_norm <- left_join(full_array_norm,
          unique(full_df[,c("coordinate_unique", "sequence_label", "library", "color_col")]))

col_annot <- as.data.frame(annot[,c(2,4)])
row.names(col_annot) <- as.character(annot$sample_name)
row.names(full_array_norm) <- full_array_norm$sequence_label
p <- pheatmap(full_array_norm[c(group1, group2)], cluster_cols = FALSE, 
              annotation_col = col_annot,
                color = colorRampPalette(c("white","lightblue", "blue"))(50),
                fontsize_col = 5,
             fontsize_row = 3)
save_pheatmap_pdf(p,file.path(outpath, "peptide_QN.pdf"), 
                  width = 3 + round(NCOL(full_array_norm)/8), height = 30)

##visualising distributions:
p1 <- full_array_norm %>% pivot_longer(!(sequence_label|coordinate_unique|library|color_col)) %>%
            ggplot(aes(x = value, color = color_col)) +
            geom_density() + 
            scale_color_manual(values = group_colors)


p2 <- full_array_norm %>% pivot_longer(!(sequence_label|coordinate_unique|library|color_col)) %>%
            ggplot(aes(x = value, color = color_col)) +
            geom_density() + 
            scale_color_manual(values = group_colors)+ facet_wrap(~library)
g <- p1/p2
ggsave(file.path(outpath, "average_dist.pdf"),g,  width = 8, height = 6)


### plot top peptides
top50 <- full_array_norm %>% pivot_longer(!(sequence_label|coordinate_unique|library|color_col)) %>% 
group_by(color_col, coordinate_unique, sequence_label) %>% summarize(mean_QN = mean(value)) %>% 
ungroup() %>% group_by(color_col) %>% slice_max(order_by = mean_QN, n = 50)
my_wt(top50, file.path(outpath, "top50.csv"))

p <- pheatmap(full_array_norm[unique(top50$sequence_label),c(group1, group2)], cluster_cols = FALSE, 
              annotation_col = col_annot,
                color = colorRampPalette(c("white","lightblue", "blue"))(50),
                fontsize_col = 5,
             fontsize_row = 5)
save_pheatmap_pdf(p,file.path(outpath, "peptide_QN_top50.pdf"), 
                  width = 4 + round(NCOL(full_array_norm)/8), height = 10)



top15 <- full_array_norm %>% pivot_longer(!(sequence_label|coordinate_unique|library|color_col)) %>% 
group_by(color_col, coordinate_unique, sequence_label) %>% summarize(mean_QN = mean(value)) %>% 
ungroup() %>% group_by(color_col) %>% slice_max(order_by = mean_QN, n = 15)

p <- pheatmap(full_array_norm[unique(top15$sequence_label),c(group1, group2)], cluster_cols = FALSE, 
              annotation_col = col_annot,
                color = colorRampPalette(c("white","lightblue", "blue"))(50),
                fontsize_col = 5,
             fontsize_row = 5)
save_pheatmap_pdf(p,file.path(outpath, "peptide_QN_top15.pdf"), 
                  width = 4 + round(NCOL(full_array_norm)/8), height = 5)
my_wt(top50, file.path(outpath, "top15.csv"))

my_wt(full_array_norm, file.path(outpath, "full_samples_normalised.csv"))

### making a PCA 
df_rot <- run_pca(full_array_norm[,c(group1, group2)], outpath)
df_rot$sample_name <- as.character(row.names(df_rot))
annot$sample_name <- as.character(annot$sample_name)

df_rot <- left_join(df_rot, annot)
ggplot(df_rot, aes(x = PC1, y = PC2, color = group_col)) + geom_point() + 
    geom_text_repel(aes(label = sample_name), size = 2) + scale_color_manual(values = group_colors)
ggsave(file.path(outpath,  "PCA_plot.pdf"), width = 5, height = 5)
