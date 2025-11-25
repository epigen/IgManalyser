## run from the project directory! Rscript src/array_preprecessing.R
### script for preprocessing the data from peptide array screening
### for each array two replicas are provided 
### we compare the values for each well in each replica, calulate the "error rate" and correlation between them

print("executing QC script!")

snakemake@source("init.R")
library(patchwork)

##function that normalizes the array by the max value 
scale_array <- function(df){
  M <- max(df)
  print(M)
  df_norm <- df/M
  return(df_norm)
}

### inputs

dataPath <- snakemake@input[[1]]
print(dataPath)
outpath <- snakemake@params[["QC_path"]]
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)

replicas_full <- data.table()
for(dataPath in unlist(snakemake@input)){
    
    print(dataPath)
    array_df <- fread(dataPath)
    colnames(array_df)[1] <- "indx"

    repl1 <- as.data.frame(array_df[, c(2:25)])
    row.names(repl1) <- array_df$indx

    repl2 <- as.data.frame(array_df[, c(27:50)])
    row.names(repl2) <- array_df$indx

    ##normalising for the visualisation and mean-cv plots and calculating the average (otherwise won't work)
    repl1_n <- scale_array(repl1)
    ##saving the visualisations
    p1 <- pheatmap(repl1_n, cluster_rows = F, cluster_cols = F, 
             main = "Replica 1", angle_col = 45)
    repl1_n$peptide <- row.names(repl1_n) ## reformat for downstream reshaping

    repl2_n <- scale_array(repl2)
    p2 <- pheatmap(repl2_n, cluster_rows = F, cluster_cols = F, 
             main = "Replica 2", angle_col = 45)
    repl2_n$peptide <- row.names(repl2_n)

    ##reshaping and joining
    repl1_m <- reshape2::melt(repl1_n, id.vars = "peptide")
    colnames(repl1_m)[3] <- "replica1"

    colnames(repl2_n) <- colnames(repl1_n)
    repl2_m <- reshape2::melt(repl2_n, id.vars = "peptide")
    colnames(repl2_m)[3] <- "replica2"

    repl_merged_m <- full_join(repl1_m, repl2_m)
    repl_merged_m$delta <- abs(repl_merged_m$replica2 - repl_merged_m$replica1)
    repl_merged_m$sdev <- apply(repl_merged_m[, c("replica1", "replica2")], 1, sd)
    repl_merged_m$avrg <- apply(repl_merged_m[, c("replica1", "replica2")], 1, mean)
    repl_merged_m$cov_var <- abs(repl_merged_m$sdev / repl_merged_m$avrg)
    ##if the CV sd/mean ratio is more than 10 it most likely indicates mean too close to zero -> for the vizualization we should filter them out from coloring, 
      ##by assigning the cov_var = 0
    if(NROW(repl_merged_m[repl_merged_m$avrg < 0.01, ]) > 0) repl_merged_m[repl_merged_m$avrg < 0.01, ]$cov_var <- 0
    ##average should also be plotted here, so creating a matrix for it
    replicas_avrg <- repl_merged_m[,c("peptide", "variable", "avrg")] %>% spread(variable, avrg)
    replicas_avrg_m <- as.matrix(replicas_avrg[,c(-1)])
    row.names(replicas_avrg_m) <- gsub(" ", "", replicas_avrg$peptide)   
    p3 <- pheatmap(replicas_avrg_m[row_levels,col_levels], cluster_rows = F, cluster_cols = F, 
             main = "avrg", angle_col = 45)

    ## also visualizing the CV: 
    replicas_CV <- repl_merged_m[,c("peptide", "variable", "cov_var")] %>% spread(variable, cov_var)
    replicas_CV_m <- as.matrix(replicas_CV[,c(-1)])
    row.names(replicas_CV_m) <- gsub(" ", "", replicas_avrg$peptide)  
    cols = colorRampPalette(c("white", "red"))(30) ## color palette (red is bad)
    p4 <- pheatmap(replicas_CV_m[row_levels,col_levels], cluster_rows = F, cluster_cols = F, 
             main = "cov_var", angle_col = 45, color = cols)

    ##and raw 
    p5 <- pheatmap(repl1, cluster_rows = F, cluster_cols = F, 
             main = "Replica 1 raw", angle_col = 45)
    p6 <- pheatmap(repl2, cluster_rows = F, cluster_cols = F, 
             main = "Replica 2 raw", angle_col = 45)
    b<-list(p5[[4]], p6[[4]], p1[[4]], p2[[4]], p3[[4]], p4[[4]])
    z <-  grid.arrange(arrangeGrob(grobs=b, ncol=2))

    libid <- ifelse(grepl("_L1", dataPath), "L1", "L2")
    try(dev.off())
    pdf(paste0(outpath, "/replica_heatmaps_", libid, ".pdf"), width = 10, height = 15)
    plot(z)
    dev.off()


    repl_merged_m$library <- libid
    replicas_full <- rbind(replicas_full, repl_merged_m)
}

#defining plot scale
lim_max <- round(max(c(replicas_full$replica1, replicas_full$replica2)), 1) + 0.1
lim_min <- min(round(min(c(replicas_full$replica1, replicas_full$replica2)), 1) - 0.1,0)
  
    

replicas_full <- replicas_full %>% group_by(library) %>% mutate(cor_coef = round(cor(replica1, replica2),3))

replicas_full$lib_and_cor <- paste0(replicas_full$library, ", corr = ", replicas_full$cor_coef)
g1 <- ggplot(replicas_full, aes(x = replica1, y = replica2 )) + 
    geom_point(aes(color = cov_var)) +
    xlim(c(lim_min, lim_max)) + ylim(c(lim_min,lim_max)) +
    coord_fixed() + 
    geom_abline(slope = 1, linetype = "dashed", color = "grey", alpha = 0.5) + 
    theme_bw() + 
    #ggtitle(paste0("corr = ", round(corr_coef, 3))) +
    scale_color_continuous(high = "red", low = "blue") + 
    labs( color = "CV = sd/mean") + facet_wrap(~lib_and_cor)
  
g2 <- ggplot(replicas_full, aes(x = replica1, y = replica2 )) + geom_point(aes(color = cov_var)) +
    xlim(c(lim_min, lim_max)) + ylim(c(lim_min,lim_max)) +
    coord_fixed() + geom_abline(slope = 1, linetype = "dashed", color = "grey", alpha = 0.5) + 
    theme_bw() + 
    #ggtitle(paste0("corr = ", round(corr_coef, 3))) +
    scale_color_continuous(high = "red", low = "blue") + 
    geom_text_repel(data = replicas_full[replicas_full$cov_var > 0.9 | replicas_full$avrg > 0.8, ], 
             aes(x = replica1, y = replica2, label = paste(peptide, variable, sep = "_")), size = 3) + 
    labs( color = "CV = sd/mean") + facet_wrap(~lib_and_cor)
g <- g1/g2
ggsave(paste0(outpath, "/replica_correlation.pdf"),g,width = 8, height = 6)
  
  #saving error distribution:
p_hist <- ggplot(replicas_full, aes(x = cov_var)) + geom_histogram(fill = "grey", color = "black") +
    geom_vline(xintercept = 0.5, color = "red", linetype = "dashed", alpha = 0.5) + 
    geom_vline(xintercept = 1, color = "red", linetype = "dashed", alpha = 0.5) + theme_bw() + 
    xlab ("CV = sd/mean") +  facet_wrap(~library)

ggsave(file.path(outpath, "replica_error.pdf"), 
        plot = p_hist, device = "pdf", width = 10, height = 5)
    
##saving the table
replicas_full$sample_name <- snakemake@params[["sample_id"]]
write.table(replicas_full, file.path(outpath, "replica_compared.csv"), 
              quote = F, row.names = F, sep = ";")
#summarizing
    
summary_data <- unique(replicas_full[c("library", "cor_coef")])
summary_data$sample_name <- snakemake@params[["sample_id"]]
    
write.table(summary_data, file.path(outpath, "/corr_summary.csv"), 
            quote = F, sep = ";", row.names = F)


