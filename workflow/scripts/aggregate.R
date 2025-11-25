snakemake@source("init.R")
outpath <- snakemake@params[["aggregatePath"]]
dir.create(outpath)

metadataPath <- snakemake@params[["metadataPath"]]
CORTHR <- snakemake@params[["CORTHR"]]
group_by <- snakemake@params[["group_by"]]
## collecting summary stats:
summary_df <- data.table()
for(dataPath in snakemake@input[["corr_files"]]){
    df <- fread(dataPath)
    summary_df <- rbind(summary_df, df)
}

if(group_by == 'NA'){
    summary_df$color_col = "group"} ##TODO: grouping by phenotype, etc

color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
for(val in unique(summary_df$color_col)){
    if(val %in% names(group_colors)) print(paste0("color defined for ", val, ": ", group_colors[val]))
    else group_colors[val] <- color[sample(length(color),1)]
    }
group_colors <- group_colors[unique(summary_df$color_col)]

p <- ggplot(summary_df, aes(x = cor_coef, fill = color_col)) + 
    geom_histogram(binwidth = 0.005,position = "identity", alpha = 0.5) + scale_fill_manual(values = group_colors) + facet_wrap(~library) + geom_vline(xintercept = CORTHR, linetype = "dashed", color = "red")
ggsave(paste0(outpath, "/cv_comparison.pdf"),p, height = 5, width = 10)
my_wt(summary_df, paste0(outpath, "/cv_comparison.csv"))

summary_df[, status:="FAIL",]
summary_df[summary_df$cor_coef > CORTHR, status:="PASS",]

## if one library failed, other also:
summary_df <- unique(summary_df[, c("sample_name","color_col", "status")])
one_lib_failed <- summary_df[duplicated(summary_df$sample_name),]$sample_name
summary_df[summary_df$sample_name %in% one_lib_failed,]$status <- "FAIL"
summary_df <- unique(summary_df) 

annot <- fread(metadataPath)
annot <- left_join(annot, summary_df)

print(paste0("passed cv threshold ", round(100*NROW(annot[status=="PASS"])/NROW(annot)), "% samples"))

## collecting data:
full_df <- data.table()
for(dataPath in snakemake@input[["array_files"]]){
    df <- fread(dataPath)
    full_df <- rbind(full_df, df)
}
full_df <- left_join(full_df, annot)

full_df <- full_df[status == "PASS"]
print(NROW(full_df))

##adding the annotation to the peptides:
library_annot <- list()

library_annot <- fread(snakemake@params[["array_anotation_path"]])
if(!("sequence_label" %in% colnames(library_annot))|!("is_duplicate" %in% colnames(library_annot))){
print("Critical feature annotation missing. Please design your metadata according to the requirements")
}
save.image("aggregate.RData")
full_df$Coordinate <- paste0(full_df$peptide, "_", full_df$variable)
full_df <- left_join(full_df, library_annot[,c("Coordinate", "library","sequence", "sequence_label", "is_duplicate")])

full_df <- full_df[!is.na(sequence),] ###removing empty ones

##Visulising duplicates
duplicate_peptides <- unique(full_df[is_duplicate=="yes"]$sequence)

if(length(duplicate_peptides) > 0){
duplicate_plot_list <- list()
for(peptide_seq in duplicate_peptides){
    
    print(peptide_seq)
    
    df <- full_df[sequence == peptide_seq]
    df$coordinate_unique <- paste0(df$Coordinate, "_", df$library) ## in case duplicates are spreaded across libraries
    
    df <- unique(df) %>% select(coordinate_unique, sample_name, avrg) %>%
            pivot_wider(names_from = sample_name, values_from = avrg)

    if(NROW(df) > 1){
    df <- as.data.frame(df)
    row.names(df) <- df$coordinate_unique
    
    p <- pheatmap(df[, c(-1)], cluster_cols = FALSE, 
                  main = peptide_seq, color = colorRampPalette(c("white", "blue"))(50),
                 breaks = seq(0,0.5, length.out = 51), fontsize_col = 4)
    duplicate_plot_list[[peptide_seq]] = p[[4]]
    }    
}

nCol <- 2
g <- do.call(grid.arrange,c(duplicate_plot_list, ncol=nCol))
ggsave(file.path(outpath, "duplicates.pdf"),g, width = 20, height = round(0.5*length(duplicate_peptides))) 
}

##dropping duplicates:
full_df <- full_df[is_duplicate == "no"]

full_df_m <- full_df %>% select(sequence_label, sample_name, avrg) %>%
            pivot_wider(names_from = sample_name, values_from = avrg)
full_df_m <- as.data.frame(full_df_m)
row.names(full_df_m) <- full_df_m$sequence_label

##visualising ALL peptides
##TO DO add annotations
col_annot <- as.data.frame(annot[,c(2,4)])
row.names(col_annot) <- as.character(annot$sample_name)
p <- pheatmap(full_df_m[, c(-1)], cluster_cols = FALSE, 
              annotation_col = col_annot,
                color = colorRampPalette(c("white","lightblue", "blue"))(50),
                  fontsize_col = 4,
             fontsize_row = 3)
save_pheatmap_pdf(p,file.path(outpath, "peptide_avrg.pdf"), 
                  width = 2 + round(NCOL(full_df_m)/4), height = 30)

##visualising distributions:
p1 <- ggplot(full_df, aes(x = avrg, color = color_col)) +
            geom_density() + 
            scale_color_manual(values = group_colors)
p2 <- ggplot(full_df, aes(x = avrg, color = color_col)) +
            geom_density() + 
            scale_color_manual(values = group_colors) + facet_wrap(~library)
g <- p1/p2
ggsave(file.path(outpath, "average_dist.pdf"),g,  width = 6, height = 6)

p3 <- ggplot(full_df, aes(x = avrg, y = sdev, color = color_col)) +
            geom_point(size = 0.1, alpha = 0.4) + 
            scale_color_manual(values = group_colors) + facet_wrap(~library)
ggsave(file.path(outpath, "avrg_sd.pdf"), p3, width = 8, height = 4)

my_wt(full_df, file.path(outpath,"full_samples_combined.csv"))



