snakemake@source("init.R")
library(RColorBrewer)
#save.image("CFA_PCA.RData")


annotpath <- snakemake@input[["metadata_file"]]
outdir <- snakemake@params[["project_dir"]]

metadata_df <- fread(annotpath)
head(metadata_df)

params <- colnames(metadata_df)[6:ncol(metadata_df)]
params <- params[params != snakemake@params[["split_group"]]]

pca_df <-as.data.frame(fread(snakemake@input[['pca_file']]))
row.names(pca_df) <- pca_df$V1
pca_df <- pca_df[,-1]


params <- colnames(metadata_df)[6:ncol(metadata_df)]
params <- params[params != snakemake@params[["split_group"]]]

pca_ids <- colnames(pca_df)

## calculating pvalues
pval_data <- data.table()
for(i in 1:(length(params))){
    for(j in 1:length(pca_ids)){
        x <- metadata_df[[params[i]]]
        y <- pca_df[as.character(metadata_df$Sample_full_ID), pca_ids[j]]

        type_x <- if(is.numeric(x)) "numeric" else check.categorical(x, y)
        type_y <- "numeric"

            
        print(paste0(params[i], " and ", params[j])) 
        
        if(type_x=="numeric" & type_y=="numeric" ) p.val <- сor.test.numeric(x,y)
        else if(type_y=="numeric" & type_x=="categorical") p.val <- cor.test.cat_num(y,x)
        else{
            print("at least one column is only NAs")
            p.val <- NA
        }
            
        if(!is.na(p.val)){
                pval_data <- rbind(pval_data, data.frame("var1" = params[i], 
                                                        "var2" =  pca_ids[j],
                                                        "pval" = p.val))
            }
    }
    }

p_values <- pval_data
# adjust p-values for multiple testing
p_values$p_values_adjusted <- p.adjust(as.vector(p_values$pval), method = "BH")
# handle adjusted p-value 0 exception
min_nonzero <- min(p_values$p_values_adjusted[p_values$p_values_adjusted > 0])
p_values$p_values_adjusted[p_values$p_values_adjusted == 0] <- min_nonzero


# Transform p-values to -log10(p-values)
p_values$log_p_values <- -log10(p_values$p_values_adjusted)
write.csv(p_values, snakemake@output[['cfa_df']])


mat <- p_values[,c(1,2,5)] %>% pivot_wider(values_from = log_p_values, names_from = var2)
mat <- as.data.frame(mat)
row.names(mat) <- mat$var1
mat <- mat[,-1]

breaksList = seq(0, max(-log10(0.01), max(p_values$log_p_values)), by = 0.1)


pheatmap(mat,
        filename = file.path(outdir, "metadata_params_PCA.pdf"), 
        width = round(NROW(mat)/2) + 1, 
        height = round(NROW(mat)/2),
        color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList)
