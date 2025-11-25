## always needed libraries
### working with data
library(data.table)
library(dplyr)
library(tidyr)

##data visualisation
library(ggplot2)
library(ggrepel)
library(patchwork)
library(pheatmap)
library(gridExtra)
##other
library(stringr)


##ggplot2 setup
theme_set(theme_classic())

## color setup
group_colors = c("HCD" = "#d73027", "CHOW" = "#4575b4", "Unstable_plaque" = "#d73027", "Stable_plaque" = "#4575b4", "group" = "darkblue")

timepoint_levels <- c("DAY0", "WEEK2", "WEEK8", "WEEK13")

##levels to sort 
col_levels <- paste0("C", c(1:24))
row_levels <- paste0("L", c(1:16))



my_wt <-function(df, path){
    write.table(df, path, sep = "\t", quote = F, row.names = F)
}


save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}


check.categorical <- function(c, s){
    type_c <- "categorical"
    if(sum(is.na(c[!is.na(s)])) == length(c[!is.na(s)])) type_c <- "no_values"
    if(length(unique(c[!is.na(s)])) < 2) type_c <- "no_values"
    return(type_c)
}

сor.test.numeric <- function(x,y){
    return(cor.test(x, y, method="kendall")$p.value)
}

cor.test.categorical <- function(x,y){
    return(fisher.test(table(x, y), simulate.p.value=TRUE, B=10000)$p.value)
}

cor.test.cat_num <- function(x,y){
    return( kruskal.test(x ~ as.factor(y))$p.value)
}
