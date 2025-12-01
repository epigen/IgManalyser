snakemake@source("init.R")
#save.image("CFA_annot.RData")

annotpath <- snakemake@input[["metadata_file"]]
outdir <- snakemake@params[["project_dir"]]

metadata_df <- fread(annotpath)
head(metadata_df)

params <- colnames(metadata_df)[6:ncol(metadata_df)]
params <- params[params != snakemake@params[["split_group"]]]

## calculating pvalues
pval_data <- data.table()
for(i in 1:(length(params)-1)){
    for(j in (i+1):length(params)){
        x <- metadata_df[[params[i]]]
        y <- metadata_df[[params[j]]]

        type_x <- if(is.numeric(x)) "numeric" else check.categorical(x, y)
        type_y <- if(is.numeric(y)) "numeric" else check.categorical(y, x)

            
        print(paste0(params[i], " and ", params[j])) 
        
        if(type_x=="numeric" & type_y=="numeric" ) p.val <- сor.test.numeric(x,y)
        else if(type_x=="categorical" & type_y=="categorical") p.val <- cor.test.categorical(x,y)
        else if(type_x=="numeric" & type_y=="categorical") p.val <- cor.test.cat_num(x,y)
        else if(type_y=="numeric" & type_x=="categorical") p.val <- cor.test.cat_num(y,x)
        else{
            print("at least one column is only NAs")
            p.val <- NA
        }
            
        if(!is.na(p.val)){
                pval_data <- rbind(pval_data, data.frame("var1" = params[i], 
                                                        "var2" = params[j],
                                                        "pval" = p.val))
            }
    }
    }


p_values <- pval_data
### visualising test

# adjust p-values for multiple testing
p_values$p_values_adjusted <- p.adjust(as.vector(p_values$pval), method = "BH")

# handle adjusted p-value 0 exception
min_nonzero <- min(p_values$p_values_adjusted[p_values$p_values_adjusted > 0])
p_values$p_values_adjusted[p_values$p_values_adjusted == 0] <- min_nonzero

# Transform p-values to -log10(p-values)
p_values$log_p_values <- -log10(p_values$p_values_adjusted)
write.csv(p_values, snakemake@output[['cfa_df']])
    
# create symmetric matrix of -log10(adjusted p values)
var_names <- unique(c(p_values$var1, p_values$var2))
mat <- matrix(NA, nrow=length(var_names), ncol=length(var_names), dimnames=list(var_names, var_names))
for(i in seq_len(nrow(p_values))){
  mat[p_values$var1[i], p_values$var2[i]] <- p_values$log_p_values[i]
  mat[p_values$var2[i], p_values$var1[i]] <- p_values$log_p_values[i]
}

# set diagonal to the maximum of -log10(adjusted p values) in the data for clustering
diag(mat) <- max(p_values$log_p_values, na.rm=TRUE)

pheatmap(mat,
        filename = file.path(outdir, "metadata_params.pdf"), 
        width = round(NROW(mat)/2) + 1, 
        height = round(NROW(mat)/2))

