library(Matrix)
preprocess <- function(expressionmatrix, min_cell_ratio = 0.05, min_gene_ratio = 0.05, normalization = TRUE
                       outdir, outname, outtype){
  in_x <- expressionmatrix
  in_ngene <- nrow(in_x)
  in_ncell <- ncol(in_x)
  cell_chosen <- colSums(in_x) >= min_gene_ratio * in_ngene
  gene_chosen <- rowSums(in_x) >= min_cell_ratio * in_ncell
  X <- in_x[gene_chosen, cell_chosen]
  if(normalization){
  cell_Sums <- colSums(X)
  # x*10^6/cell_Sums: transform count into TPM
  tpm <-  sweep(X, MARGIN = 2, 10^6/cell_Sums, FUN = "*")
  out_x <- log10(tpm + 1)
  return(out_x)}else{
    out_x <- X
    return(out_x)
  }
  
  if(outtype == "RDS"){
    saveRDS(out_x , file = paste0(out_directory,"/", outname,".RDS"))
  }
  else if(intype == "RData"){
    save(out_x , file = paste0(out_directory,"/", outname,".RData"))
  }
  else{
    write.csv(out_x , file = paste0(input_directory, filename))
  }
}

data <- usoskin@assays[["data"]]@listData[["count"]]
X <- data
normalize.by.size.effect <- FALSE
M <- ncol(X)
N <- nrow(X)
min.expressed.gene <- 5
max.expressed.ratio <- 0.5
min.expressed.cell <- 5
m <- Matrix::colSums(X > 1) >= min.expressed.gene	# cells that have at least min.expressed.gene expreseed genes
n <- Matrix::rowSums(X > 1) <= max.expressed.ratio * M & Matrix::rowSums(X > 1) >= min.expressed.cell	# genes that are detected in at least min.expressed.cell or at most max.expressed.ratio cells
if (normalize.by.size.effect){
  sf <- apply((X[n, m] + 1) / exp(Matrix::rowMeans(log(X[n, m] + 1))), 2, median)
  X <- t(t(X[n, m]) / sf)
}else
  X <- X[n, m]
data <- X
rm(X,usoskin,m,n,M,N)
rm(normalize.by.size.effect)
rm(min.expressed.cell,min.expressed.gene,max.expressed.ratio)

