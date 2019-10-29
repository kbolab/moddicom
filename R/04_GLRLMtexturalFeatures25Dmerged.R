######################################
### BEGIN FUNCION ####################
######################################


#' @import radiomics data.table

glrlmTexturalFeatures25Dmerged <- function(imgObj, n_grey){

  # compute number of non-NA voxels

  #nVoxel <- dim(imgObj)[1]*dim(imgObj)[2]*dim(imgObj)[3] - sum(is.na(imgObj))
  Nv <- length(imgObj[which(!is.na(imgObj))])*4

  ### compute Gray Levels Cooccurrence Matrices

  R_list_0 <- list()
  R_list_45 <- list()
  R_list_90 <- list()
  R_list_135 <- list()
  #Nv <- numeric()
  #Compute grey level cooccurrence matrices for 4 different directions within each slice
  for (i in 1:dim(imgObj)[3]){
    if(length(imgObj[,,i])*.005 >= sum(!is.na(imgObj[,,i]))) next
    imgObj[,,i] <- round(imgObj[,,i])
    #if(length(table(unique(imgObj[,,i])))<n_grey) n_grey <- length(table(unique(imgObj[,,i])))-1
    if (min(imgObj[,,i],na.rm = T) < 0) {imgObj[,,i] <- imgObj[,,i] + abs(min(imgObj[,,i],na.rm = T))}
    R_list_0[[i]] <- as.matrix(glrlm(imgObj[,,i], angle = 0, verbose=F,truncate = F, n_grey = n_grey))
    R_list_45[[i]] <- as.matrix(glrlm(imgObj[,,i], angle = 45, verbose=F,truncate = F, n_grey = n_grey))
    R_list_90[[i]] <- as.matrix(glrlm(imgObj[,,i], angle = 90, verbose=F,truncate = F, n_grey = n_grey))
    R_list_135[[i]] <- as.matrix(glrlm(imgObj[,,i], angle = 135, verbose=F,truncate = F, n_grey = n_grey))
  }

  
  #elimino gli elementi NULL della lista
  if(length(R_list_0)!=0){
    if(length(which(sapply(R_list_0, is.null)))!=0){
      R_list_0 = R_list_0[-which(sapply(R_list_0, is.null))]
    }
  }
  
  #elimino gli elementi NULL della lista
  if(length(R_list_45)!=0){
    if(length(which(sapply(R_list_45, is.null)))!=0){
      R_list_45 = R_list_45[-which(sapply(R_list_45, is.null))]
    }
  }
  
  #elimino gli elementi NULL della lista
  if(length(R_list_90)!=0){
    if(length(which(sapply(R_list_90, is.null)))!=0){
      R_list_90 = R_list_90[-which(sapply(R_list_90, is.null))]
    }
  }
  
  #elimino gli elementi NULL della lista
  if(length(R_list_135)!=0){
    if(length(which(sapply(R_list_135, is.null)))!=0){
      R_list_135 = R_list_135[-which(sapply(R_list_135, is.null))]
    }
  }
  
  
  sumtot <- list()
  ### DIRECTION 0
  matrix.df <- ldply(R_list_0, data.table::melt, varnames=c("row", "col"))
  sumtot[[1]] <- acast(matrix.df, row ~ col, sum)
  ### DIRECTION 45
  matrix.df <- ldply(R_list_45, data.table::melt, varnames=c("row", "col"))
  sumtot[[2]] <- acast(matrix.df, row ~ col, sum)
  ### DIRECTION 90
  matrix.df <- ldply(R_list_90, data.table::melt, varnames=c("row", "col"))
  sumtot[[3]] <- acast(matrix.df, row ~ col, sum)
  ### DIRECTION 135
  matrix.df <- ldply(R_list_135, data.table::melt, varnames=c("row", "col"))
  sumtot[[4]] <- acast(matrix.df, row ~ col, sum)
  
  matrix.df <- ldply(sumtot, data.table::melt, varnames=c("row", "col"))
  sumtottutte <- acast(matrix.df, row ~ col, sum)

  #Initialise data table for storing GLCM features; I have added a few
  featNames <- c("F_rlm.2.5Dmerged.sre","F_rlm.2.5Dmerged.lre","F_rlm.2.5Dmerged.lgre","F_rlm.2.5Dmerged.hgre","F_rlm.2.5Dmerged.srlge",
                 "F_rlm.2.5Dmerged.srhge","F_rlm.2.5Dmerged.lrlge","F_rlm.2.5Dmerged.lrhge","F_rlm.2.5Dmerged.glnu","F_rlm.2.5Dmerged.glnu.norm","F_rlm.2.5Dmerged.rlnu","F_rlm.2.5Dmerged.rlnu.norm", "F_rlm.2.5Dmerged.r.perc",
                 "F_rlm.2.5Dmerged.gl.var","F_rlm.2.5Dmerged.rl.var","F_rlm.2.5Dmerged.rl.entr")
  F_rlm <- data.table(matrix(NA, nrow=1, ncol=length(featNames)))
  colnames(F_rlm) <- featNames


    Ns <- sum(unlist(sumtot),na.rm=T)

    #Convert matrix to data table
    df.R    <- data.table(sumtottutte)

    #Add row grey level intensity
    df.R$i <- as.numeric(row.names(sumtottutte))

    #Convert from wide to long table. This is the preferred format for data tables and data frames
    df.R   <- data.table::melt(df.R, id.vars="i", variable.name="j", value.name="n", variable.factor=FALSE)

    #Convert j from string to numeric
    df.R$j <- as.numeric(df.R$j)

    #Remove combinations with 0 counts
    df.R   <- df.R[n>0,]

    df.R <- df.R[order(rank(i))]

    #Convert Grey level coccurrence matrix to joint probability
    #df.r_ij <- df.R[,.(r_ij=n/sum(df.R$n)), by=.(i,j)] #joint probability

    df.r_i  <- df.R[,.(r_i=sum(n)), by=i]        #marginal probability over columns
    df.r_j  <- df.R[,.(r_j=sum(n)), by=j]        #marginal probability over rows

    #Diagonal probabilities (p(i-j))
    #First, we create a new column k which contains the absolute value of i-j.
    #Second, we sum the joint probability where k is the same.
    #This can written as one line by chaining the operations.
    # df.p_imj <- copy(df.p_ij)
    # df.p_imj <- df.p_imj[,"k":=abs(i-j)][,.(p_imj=sum(p_ij)), by=k]

    #Cross-diagonal probabilities (p(i+j))
    #Again, we first create a new column k which contains i+j
    #Second, we sum the probability where k is the same.
    #This is written in one line by chaining the operations.
    # df.p_ipj <- copy(df.p_ij)
    # df.p_ipj <- df.p_ipj[,"k":=i+j][,.(p_ipj=sum(p_ij)), by=k]

    #Merger of df.p_ij, df.p_i and df.p_j
    df.R <- merge(x=df.R, y=df.r_i, by="i")
    df.R <- merge(x=df.R, y=df.r_j, by="j")

    #Thus we have five probability matrices
    #Joint probability:          df.p_ij with probabilities p_ij, p_i and p_j, and indices i, j
    #Marginal probability:       df.p_i with probability p_i, and index i
    #Marginal probability:       df.p_j with probability p_j, and index j
    #Diagonal probability:       df.p_imj with probability p_imj and index k
    #Cross-diagonal probability: df.p_ipj with probability p_ipj and index k

    #Calculate features
    #Short runs emphasis
    F_rlm$F_rlm.2.5Dmerged.sre         <- (1/Ns) * sum(df.r_j$r_j/(df.r_j$j^2))

    #Long runs emphasis
    F_rlm$F_rlm.2.5Dmerged.lre         <- (1/Ns) * sum(df.r_j$r_j*(df.r_j$j^2))

    #Low grey level run emphasis
    F_rlm$F_rlm.2.5Dmerged.lgre         <- (1/Ns) * sum(df.r_i$r_i/(df.r_i$i^2))

    #High grey level run emphasis
    F_rlm$F_rlm.2.5Dmerged.hgre         <- (1/Ns) * sum(df.r_i$r_i*(df.r_i$i^2))

    #Short run low grey level emphasis
    F_rlm$F_rlm.2.5Dmerged.srlge <- (1/Ns) * sum(df.R$n/((df.R$i^2)*(df.R$j^2)))

    #Short run high grey level emphasis
    F_rlm$F_rlm.2.5Dmerged.srhge <- (1/Ns) * sum((df.R$n)*(df.R$i^2)/(df.R$j^2))

    #Long run low grey level emphasis
    F_rlm$F_rlm.2.5Dmerged.lrlge <- (1/Ns) * sum((df.R$n)*(df.R$j^2)/(df.R$i^2))

    #Long run high grey level emphasis
    F_rlm$F_rlm.2.5Dmerged.lrhge <- (1/Ns) * sum((df.R$n)*(df.R$j^2)*(df.R$i^2))

    #Grey level non-uniformity
    F_rlm$F_rlm.2.5Dmerged.glnu <- (1/Ns) * sum(df.r_i$r_i^2)

    #Grey level non-uniformity normalized
    F_rlm$F_rlm.2.5Dmerged.glnu.norm <- (1/Ns^2) * sum(df.r_i$r_i^2)

    #Run length non-uniformity
    F_rlm$F_rlm.2.5Dmerged.rlnu <- (1/Ns) * sum(df.r_j$r_j^2)

    #Run length non-uniformity normalized
    F_rlm$F_rlm.2.5Dmerged.rlnu.norm <- (1/Ns^2) * sum(df.r_j$r_j^2)

    #Run percentage

    F_rlm$F_rlm.2.5Dmerged.r.perc <- Ns/Nv

    #Grey level variance
    p_ij <- df.R$n/sum(df.R$n)
    mu_i <- sum(df.R$i * p_ij)
    F_rlm$F_rlm.2.5Dmerged.gl.var <- sum((df.R$i - mu_i)^2 * p_ij)

    #Run length variance
    mu_j <- sum(df.R$j * p_ij)
    F_rlm$F_rlm.2.5Dmerged.rl.var <- sum((df.R$j - mu_j)^2 * p_ij)

    #Run entropy
    F_rlm$F_rlm.2.5Dmerged.rl.entr <- - sum(p_ij * log2(p_ij))

  return(F_rlm)
}


