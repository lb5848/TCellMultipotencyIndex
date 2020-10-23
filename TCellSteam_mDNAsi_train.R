library(gelnet)
library(dplyr)
library(gdata)
library(DT)

steam_dir <- getwd() 
work_dir <- paste(steam_dir, "TCellMultipotencyIndex", sep = "/")
setwd(work_dir)

replace.NA <-function(data,type.info,by = "mean"){
  if(!"group" %in% colnames(type.info)) stop("type.info must have group column")
  if(!"sample" %in% colnames(type.info)) stop("type.info must have a sample column")
  
  # Do we have NAs?
  if(is.na(table(is.na(data))["TRUE"])){
    message("No NAs were found")
    return(data)
  }
  # get NAs index 
  idx <- which(is.na(data) == TRUE,arr.ind=TRUE)
  count <- table(rownames(idx))
  message("======= Status Number of NA in probes ========")
  message("--------------------- Summary------------------")
  print(summary(as.numeric(count)))
  message("\n----------- Probes with more nb of NAs -----------")
  print(head(sort(count,decreasing = T)))
  message("===============================================")
  
  idx <- cbind(idx, mean = NA, median = NA)
  
  # For each NA value calculate the mean for the same probe for the samples
  # where it belongs
  for(line in 1:nrow(idx)){
    row <- idx[line,1]
    col <- idx[line,2]
    probe <- rownames(idx)[line]
    sample <- colnames(data)[col]
    group <- type.info[type.info$sample == sample,"group"]
    samples.in.group <- type.info[type.info$group == group,]$sample
    
    # get the probe value for all samples in the same group 
    aux <- data[rownames(data) %in% probe, colnames(data) %in% samples.in.group] 
    
    idx[line,3] <- mean(as.numeric(aux),na.rm = TRUE)
    idx[line,4] <- median(as.numeric(aux),na.rm = TRUE)
  }
  # Step 2 replace
  for(line in 1:nrow(idx)){
    row <- idx[line,1]
    col <- idx[line,2]
    if(by == "mean"){
      data[idx[line,1],idx[line,2]] <- idx[line,3]  
    } else if(by == "median") { 
      data[idx[line,1],idx[line,2]] <- idx[line,4]
    }
  }
  return(data)
}

# Data (99 samples as columns and 219 met probes as rows)
load("pcbc.data.Rda")
datatable(pcbc.data[1:3,1:4],
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 40), 
          rownames = FALSE)

load("pcbc.pd.f.Rda")
datatable(pcbc.pd.f[1:3,1:4],
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 40), 
          rownames = FALSE)



#Find the mean center by subtracting the mean of each probe from the entire pcbc data
#The mean of each probe just be in a numeric vector the same size as the number of probes (in this case 219)

m <- apply(pcbc.data, 1, mean)
m[1:5]

pcbc.data.2 <- pcbc.data - m
datatable(head(pcbc.data.2),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 40), 
          rownames = FALSE)

#Identify stem cells and break up all samples into 2 groups:
#* Stem cell (X.tr object)
#* not stem cell (X.bk object)

# Define PCBC groups (SC and non.SC)
M1_smp <- pcbc.pd.f[pcbc.pd.f$Diffname_short %in% "SC",] #SC
M2_smp <- pcbc.pd.f[!(pcbc.pd.f$Diffname_short %in% "SC"),] #non-SC

# Select PCBC data
X.tr <- pcbc.data.2[, as.character(M1_smp$UID)] # 44 samples
X.bk <- pcbc.data.2[, as.character(M2_smp$UID)] # 55 samples

## Train 

## Train a one-class model (transpose the matrix!)
mm <- gelnet(t(X.tr), NULL, 0, 1) #NULL for a one-class task 
## Store the signature to a file
save( mm, file = "pcbc-stemsig.p219.Rda")

## Leave One Out Cross-Validation 

# Cross-validation with linear model:
# Perform leave-one-out cross-validation
auc <- c()
for(i in 1:ncol(X.tr)) {
  ## Train a model on non-left-out data
  X1 <- X.tr[,-i]
  X1 <- as.matrix(X1)
  K <- t(X1) %*% X1 / nrow(X1)
  m1 <- gelnet.ker(K, NULL, lambda = 1)
  w1 <- X1 %*% m1$v
  
  ## Score the left-out sample against the background
  X.bk <- X.bk[rownames(X.tr),]
  X.bk <- as.matrix(X.bk)
  s.bk <- t(w1) %*% X.bk
  s.bk <- unmatrix(s.bk)
  
  s1 <- t(w1) %*% X.tr[,i]
  s1 <- unmatrix(s1)
  
  ## AUC = P( left-out sample is scored above the background )
  auc[i] <- sum(s1 > s.bk) / length(s.bk)
  cat( "Current AUC: ", auc[i], "\n" )
  cat( "Average AUC: ", mean(auc), "\n" )
}

head(auc)
all(auc == 1)
