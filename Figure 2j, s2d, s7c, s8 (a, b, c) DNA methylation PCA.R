### For Figure 2j, supplementary figure 2d, supplementary figure 7c, supplementary figure 8 (a, b, c). For subset analyses, substitute filenames below
### DNA methylation, PCA analysis
library(factoextra)
library(data.table)
library(writexl)
library(readxl)

N_sheets = 6
setwd('')
pr_store = list()
big_mat = data.frame()
for (i in seq(1,6)) {
    
    pd = read_excel('sample sheet here with beta values',col_names = TRUE, sheet = i)
    genes = pd[,1]
    pd = as.data.frame (pd[,-1], row.names = genes)
    pd2 = transpose(pd)
    colnames(pd2) <- as.matrix(genes)
    rownames(pd2) <- colnames(pd)
    pr = prcomp(pd2)
    pr_store[[i]] = pr
    big_mat = rbind(big_mat,pd2)
    fviz_pca_ind(pr)
    
    write_xlsx(as.data.frame(cbind(rownames(pd2),pr$x)),paste('PCA_sheet_',i,'.xlsx'))

}

pr = prcomp(big_mat)
fviz_pca_ind(pr)
write_xlsx(as.data.frame(cbind(rownames(big_mat),pr$x)),paste('PCA_sheet_','combined','.xlsx'))

fviz_pca_ind(pr_store[[1]])
fviz_pca_ind(pr_store[[2]])
fviz_pca_ind(pr_store[[3]])
fviz_pca_ind(pr_store[[4]])
fviz_pca_ind(pr_store[[5]])
fviz_pca_ind(pr_store[[6]])