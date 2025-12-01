snakemake@source("init.R")
snakemake@source("QN_functions.R")

#save.image("split.RData")

array_df <- fread(snakemake@input[['array_file_QN']], header = TRUE)
#colnames(array_df) <- gsub("V", "", colnames(array_df))
mtd <- fread(snakemake@params[['metadataPath']])
outdir <- snakemake@params[["project_dir"]]
group_id <- snakemake@params[['split_group']]
dir.create(outdir)
###splits
split_subsets <- snakemake@params[['split_subsets']]
split_subsets <- strsplit(split_subsets, ";")[[1]]


#### plot_pca_categorical and qualitive
plot_params_pca <- function(pca_df, annotation_df, params_qual = NA, params_quant = NA){
     pca_df$Sample_full_ID <- gsub("X", "", pca_df$Sample_full_ID )                          
    full_df <- left_join(pca_df, annotation_df)

    if(!is.na(params_qual)[1]){
    plot_list <- c()
    for(param_value in params_qual){
        plot_list[[param_value]] <- ggplot(full_df, aes(x = PC1, y = PC2, color = as.factor(.data[[param_value]]))) + 
            geom_point(size = 2) + ggtitle(param_value) + labs(color = param_value)
    }
    g1 <- wrap_plots(plot_list, ncol = 4)
    }else{
        g1 = NULL
    }

    if(!is.na(params_qual)[1]){
    plot_list_quant <- c()
        for(param_value in params_quant){
            plot_list_quant[[param_value]] <- ggplot(full_df, aes(x = PC1, y = PC2, color = .data[[param_value]])) + 
                geom_point(size = 2) +  ggtitle(param_value) + labs(color = param_value)
            }
        if(length(plot_list_quant) > 1){
              g2 <- wrap_plots(plot_list_quant, ncol = 4)  
        }else{
            g2 <- plot_list_quant[[1]]
        }
    }else{
        g2 = NULL
    }
    return(list("qual" = g1, "quant" = g2))
}


### mtd cleanup
mtd <- mtd[mtd$Sample_full_ID %in% colnames(array_df),]
mtd <- as.data.frame(mtd)

mtd$Sample_full_ID <- as.character(mtd$Sample_full_ID)

## select params
param_quant <- c()
param_qual <- c()
for(param in colnames(mtd)[6:(ncol(mtd))]){
    if(is.numeric(mtd[[param]])){
        param_quant <- c(param_quant, param)
    }else{
        if(sum(is.na(mtd[[param]])) == NROW(mtd)){
            print("only NAs in the column")
        }else{
           param_qual <- c(param_qual, param) 
        }        
    }   
}

###prepare subsets and run PCA
for(spl_group in split_subsets){
    dir.create(file.path(outdir, spl_group))
    patients_subgr <- mtd[mtd[group_id] == spl_group,]$Sample_full_ID
    idx <- length(patients_subgr)
    patients_subgr <-c( patients_subgr, 'coordinate_unique', 'sequence_label')
    array_sub <- array_df[,..patients_subgr]
    
    pca_df <- run_pca(array_sub[,1:idx], file.path(outdir, spl_group))
    write.csv(array_sub, file.path(outdir, spl_group, paste0(spl_group, '_samples_subsetted.csv')))
    write.csv(pca_df, file.path(outdir, spl_group, paste0(spl_group, '_pca.csv')))
    pca_df$Sample_full_ID <- row.names(pca_df)
    plots <- plot_params_pca(pca_df, mtd, params_qual =  param_qual, params_quant = param_quant)
    ggsave(file.path(outdir, spl_group, "params_quant_pca.pdf"), plots$quant, width = 20, height = 3*round(length(param_quant)/4))
    ggsave(file.path(outdir, spl_group, "params_qual_pca.pdf"), plots$qual, width = 20, height = 3*round(length(param_qual)/4))

    write.csv(mtd[mtd[group_id] == spl_group,], file.path(outdir, spl_group, paste0(spl_group, '_metadata_subsetted.csv')))
}
