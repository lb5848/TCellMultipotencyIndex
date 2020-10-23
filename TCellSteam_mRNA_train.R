library(gelnet)
library(dplyr)
library(biomaRt)
library(synapser)
synLogin("*******", "*******")

# Maps ENSEMBL IDs to HUGO
# Use srcType = "ensembl_gene_id" for Ensembl IDs
# Use srcType = "entrezgene" for Entrez IDs
genes2hugo <- function( v, srcType = "ensembl_gene_id" )
{
  ## Retrieve the EMSEMBL -> HUGO mapping
  ensembl <- biomaRt::useMart( "ENSEMBL_MART_ENSEMBL", host="www.ensembl.org", dataset="hsapiens_gene_ensembl" )
  ID <- biomaRt::getBM( attributes=c(srcType, "hgnc_symbol"), filters=srcType, values=v, mart=ensembl )
  ## Make sure there was at least one mapping
  if( nrow(ID) < 1 ) top( "No IDs mapped successfully" )
  
  ## Drop empty duds
  j <- which( ID[,2] == "" )
  if( length(j) > 0 ) ID <- ID[-j,]
  stopifnot( all( ID[,1] %in% v ) )
  ID
}

main.train <- function( fnOut = "pcbc-stemsig.tsv", fnGenes = NULL )
{
  #Train the model
  # Load RNAseq data
  synRNA <- synGet("syn2701943")
  
  X <- read.delim(synRNA$path) %>%
    tibble::column_to_rownames("tracking_id") %>% as.matrix
  
  # Retrieve metadata
  request <- synTableQuery('select * from syn3156503')
  synMeta <- request$asDataFrame()
  
  Y <- synMeta %>%
    mutate( UID = gsub("-", ".", UID) ) %>%
    tibble::column_to_rownames( "UID" )
  
  # Retrieve the labels from the metadata
  y <- Y[colnames(X),]
  rownames(y) <- colnames(X)
  
  # Fix the missing labels by hand
  y["SC11.014BEB.133.5.6.11"] <- "EB"
  y["SC12.039ECTO.420.436.92.16"] <- "ECTO"
  
  ## Drop the splice form ID from the gene names
  v <- strsplit( rownames(X), "\\." ) %>% lapply( "[[", 1 ) %>% unlist()
  rownames(X) <- v
  
  # Map Ensembl IDs to HUGO
  V <- genes2hugo(rownames(X))
  
  #Change the row names of `X` from the gene name to the hgnc 
  #(HUGO Gene Nomenclature Committee) symbol
  X <- X[V[,1],]
  rownames(X) <- V[,2]
  
  #Reduce gene set to the provides list (if any)
  if(is.null(fnGenes ) == FALSE ){
    vGenes <- read.delim( fnGenes, header=FALSE ) %>% as.matrix() %>% drop()
    VE <- genes2hugo( vGenes, "entrezgene" )
    X <- X[intersect( rownames(X), VE[,2] ),]
  }
  
  #Find the mean center by subtracting the mean of each gene (`m`) from the RNA-seq data (`X`)
  m <- apply( X, 1, mean )
  X <- X - m
  
  #Identify stem cells and break up all samples into 2 groups:
  j <- which(y$Diffname_short == "SC")
  X.tr <- X[,j]
  
  X.bk <- X[,-j]
  
  ## Train a one-class model 
  #(transpose the matrix so that the genes are listed as rows and samples as columns)
  mm <- gelnet(t(X.tr), NULL, 0, 1)
  
  ## Store the signature to a file
  write.table(mm$w, file = fnOut, sep = "\t", quote = FALSE, col.names = FALSE)
  
  auc <- c()
  for(i in 1:ncol(X.tr)){
    ## Train a model on non-left-out data
    X1 <- X.tr[,-i]
    m1 <- gelnet(t(X1), NULL, 0, 1 )
    
    ## Score the left-out sample against the background
    s.bk <- apply( X.bk, 2, function(z) {cor( m1$w, z, method="sp" )} )
    s1 <- cor( m1$w, X.tr[,i], method="sp" )
    
    ## AUC = P( left-out sample is scored above the background )
    auc[i] <- sum( s1 > s.bk ) / length(s.bk)
    cat( "Current AUC: ", auc[i], "\n" )
    cat( "Average AUC: ", mean(auc), "\n" )
  }
  return(auc)
}

