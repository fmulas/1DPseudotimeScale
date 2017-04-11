
require("sva")
# select paths
path <-  "data/"
path_res <- "res/"

# variance-based filtering. Returns data filtered for top genes according 
# to the selected scores if ngenes is chosen; data with score > threshold 
# if min.score or min.perc are specified
var.filter <- function(dat,
                       score = "mad",
                       score.fun = "mad",
                       dir= "top",
                       ngenes=NULL,
                       min.score=NULL, 
                       min.perc=NULL)
{
  if (is.null(ngenes) && is.null(min.score)  && is.null(min.perc) )
    stop( "must specify ngenes or min.score or or min.perc" )
  if (length(which(c(!is.null(ngenes), !is.null(min.score), !is.null(min.perc)))) > 1 )
    stop( "cannot specify more than one among ngenes, min.score and min.perc" )
  if (!is.null(ngenes) && ngenes>nrow(dat) )
    stop( "ngenes is too large" )
  
  ctr <- if (score=="mad") apply( dat, 1, median ) else rowMeans( dat )
  SC1 <- apply( dat, 1, score.fun )
  
  if ( is.null(ngenes) ) {
    if (! is.null(min.perc)) {
      scorev <- quantile(SC1, min.perc)
    } else {
      scorev <- min.score
    }
    
    idx <- if(dir=="top") SC1>=scorev else SC1<=scorev
    if (sum(idx)==0)
      stop( "no genes passed the filtering criteria" )
  }
  else {
    if (dir=="top") SC1 <- -SC1
    idx <- order(SC1)[1:ngenes]
  }
  
  datff <- dat[idx,,drop=FALSE]
  
  return(datff)
  
}



#---------------------------------- 

hkg <- read.table(paste(path, 'human_HK_genes.txt', sep = ""), header=TRUE, sep='\t')
hkl <- as.character(hkg$Gene.Symbol)

data_all <- read.table(paste(path, "ALL_RPKM_sc.txt", sep = ""))
genes <- rownames(data_all)

data_bulk <- read.table(paste(path, "bulk_RPKM_nognames.txt", sep = ""), header=TRUE, row.names=genes)
# log!
data <- log2(data_bulk+1)

# select expressed genes
rowm = rowMeans(data) 
keepIndex=which(rowm > 5)
dataok = as.matrix(data[keepIndex,])

# select genes based on low variance
h <- var.filter(dataok, min.perc=0.05, dir="bottom")
SC1 <- apply( h, 1, "mad" )
ord <- order(SC1)
mad_ord <- rownames(h)[ord]

dataok_m <- as.matrix(data_all)
dataok_all <- dataok_m[which(rowSums(dataok_m)>0),]

# in bulk
data_flat <- dataok[which(rownames(dataok) %in% mad_ord),]
genes_flat <- rownames(data_flat)

# in single cells
data_flat_sc <- dataok_all[which(rownames(dataok_all) %in% mad_ord),]

# expressed in sc
SC3 <- apply( data_flat_sc, 1, function(x) length(x[log2(x+1)>5]))
keepIndex=which(SC3 > 2)
dataokf_sc = as.matrix(data_flat_sc[keepIndex,])

#ord3 <- order(SC3, decreasing = TRUE)
SC12 <- apply( dataokf_sc, 1, "mad" )

ord2 <- order(SC12)
mad_ord2 <- rownames(dataokf_sc)[ord2]

batches_f <- read.table(paste(path, "ALL_plates_libraries.txt", sep = ""), header=TRUE, row.name=1)

mod1 = model.matrix(~1, data=batches_f)
mod0 = cbind(mod1[,1])

# from bmc bioinfo paper
# this function regresses SVs out of genomic data
cleaningP = function(y, mod, svaobj,  P=ncol(mod)) {
  X=cbind(mod,svaobj$sv)
  Hat=solve(t(X)%*%X)%*%t(X)
  beta=(Hat%*%t(y))
  cleany=y-t(as.matrix(X[,-c(1:P)])%*%beta[-c(1:P),])
  return(cleany)
}

pp <- 0.9 
h2 <- var.filter(dataokf_sc, min.perc=pp, dir="bottom")
genes_flat <- rownames(h2) 
inters_hk <- length(which(toupper(genes_flat) %in% hkl))

# in order to retain only genes for which low variance was detected in the single cells data and to exclude markers of 
# beta cells function, we applied an additional filter on the control genes, by excluding the upper 10th percentile of 
# the genes ranked by MAD in the sc data. The resulting list contains a high number of known housekeeping genes (23 out 
# of 45) as well as other genes showing the least variation across time in our bulk data

controls <- rownames(dataok_all) %in% genes_flat

sup_svseq = svaseq(dataok_all,mod1,mod0,controls=controls,n.sv=1)

data_clean <- cleaningP(log2(dataok_all+1), mod1, sup_svseq)
datag <- var.filter(data_clean, min.perc=0.75)

write.table(data_clean, file = paste(path_res, 'data_ALL_stages_batchcorr_svaseq.txt', sep=''), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".",
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "", row.name=TRUE)

write.table(datag, file = paste(path_res, 'data_ALL_stages_batchcorr_svaseq_upper_quart_variant.txt', sep=''), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".",
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "", row.name=TRUE)




