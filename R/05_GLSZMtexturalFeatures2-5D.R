##############################################################################################
########################## GLSZM ############################################################
############################################################################################

glszmTexturalFeatures25D <- function(imgObj,px=2,py=2,pz=2,n_grey){

  # compute number of non-NA voxels

  Nv <- dim(imgObj)[1]*dim(imgObj)[2]*dim(imgObj)[3] - sum(is.na(imgObj))
  #Nv <- length(imgObj[which(!is.na(imgObj))])

  ### compute Gray Levels Cooccurrence Matrices

  S_list <- list()
  #Compute grey level cooccurrence matrices for 4 different directions within each slice
  for (i in 1:dim(imgObj)[3]){
    if(length(imgObj[,,i])*.005 >= sum(!is.na(imgObj[,,i]))) next
    #if(length(table(unique(imgObj[,,i])))<n_grey) n_grey <- length(table(unique(imgObj[,,i])))-1
    if (min(imgObj[,,i],na.rm = T) < 0) {imgObj[,,i] <- imgObj[,,i] + abs(min(imgObj[,,i],na.rm = T))}
    S_list[[i]] <- as.matrix(glszm(imgObj[,,i], verbose=F,n_grey = n_grey))
  }

  #elimino gli elementi NULL della lista
  if(length(S_list)!=0){
    if(length(which(sapply(S_list, is.null)))!=0){
      S_list = S_list[-which(sapply(S_list, is.null))]
    }
  }
  
  matrix.df <- ldply(S_list, data.table::melt, varnames=c("row", "col"))
  sumtot <- acast(matrix.df, row ~ col, sum)


  #Initialise data table for storing GLCM features; I have added a few
  featNames <- c("F_szm_2.5D.sze","F_szm_2.5D.lze","F_szm_2.5D.lgze","F_szm_2.5D.hgze","F_szm_2.5D.szlge","F_szm_2.5D.szhge","F_szm_2.5D.lzlge",
                 "F_szm_2.5D.lzhge","F_szm_2.5D.glnu","F_szm_2.5D.glnu.norm", "F_szm_2.5D.zsnu","F_szm_2.5D.zsnu.norm", "F_zsm_2.5D.z.perc","F_szm_2.5D.gl.var","F_szm_2.5D.zs.var",
                 "F_szm_2.5D.z.entr")

  F_szm <- data.table(matrix(NA, nrow=1, ncol=length(featNames)))
  colnames(F_szm) <- featNames

  #Iterate over grey level cooccurrence matrices
  #The idea is that basically every GLCM is the same, i.e. we can just perform the same operations on every glcm.

    Ns <- sum(sumtot,na.rm=T)
    #Convert matrix to data table
    df.S    <- data.table(sumtot)

    #Add row grey level intensity
    df.S$i <- as.numeric(row.names(sumtot))


    #Convert from wide to long table. This is the preferred format for data tables and data frames
    df.S   <- data.table::melt(df.S, id.vars="i", variable.name="j", value.name="n", variable.factor=FALSE)

    #Convert j from string to numeric
    df.S$j <- as.numeric(df.S$j)

    #Remove combinations with 0 counts
    df.S   <- df.S[n>0,]

    #Convert Grey level coccurrence matrix to joint probability
    #df.r_ij <- df.R[,.(r_ij=n/sum(df.R$n)), by=.(i,j)] #joint probability

    df.s_i  <- df.S[,.(s_i=sum(n)), by=i]        #marginal probability over columns
    df.s_j  <- df.S[,.(s_j=sum(n)), by=j]        #marginal probability over rows

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
    df.S <- merge(x=df.S, y=df.s_i, by="i")
    df.S <- merge(x=df.S, y=df.s_j, by="j")

    #Thus we have five probability matrices
    #Joint probability:          df.p_ij with probabilities p_ij, p_i and p_j, and indices i, j
    #Marginal probability:       df.p_i with probability p_i, and index i
    #Marginal probability:       df.p_j with probability p_j, and index j
    #Diagonal probability:       df.p_imj with probability p_imj and index k
    #Cross-diagonal probability: df.p_ipj with probability p_ipj and index k

    #Calculate features

    #Small zone emphasis
    F_szm$F_szm_2.5D.sze        <- (1/Ns) * sum(df.s_j$s_j/(df.s_j$j^2))

    #Large zone emphasis
    F_szm$F_szm_2.5D.lze        <- (1/Ns) * sum(df.s_j$s_j * df.s_j$j^2)

    #Low grey level zone emphasis
    F_szm$F_szm_2.5D.lgze       <- (1/Ns) * sum(df.s_i$s_i/(df.s_i$i^2))

    #High grey level zone emphasis
    F_szm$F_szm_2.5D.hgze        <- (1/Ns) * sum(df.s_i$s_i * df.s_i$i^2)

    #Small zone low grey level emphasis
    F_szm$F_szm_2.5D.szlge <- (1/Ns) * sum(df.S$n/((df.S$i^2)*(df.S$j^2)))

    #Small zone high grey level emphasis
    F_szm$F_szm_2.5D.szhge <- (1/Ns) * sum((df.S$n)*(df.S$i^2)/(df.S$j^2))

    #Large zone low grey level emphasis
    F_szm$F_szm_2.5D.lzlge <- (1/Ns) * sum((df.S$n)*(df.S$j^2)/(df.S$i^2))

    #Large zone high grey level emphasis
    F_szm$F_szm_2.5D.lzhge <- (1/Ns) * sum((df.S$n)*(df.S$j^2)*(df.S$i^2))

    #Grey level non-uniformity
    F_szm$F_szm_2.5D.glnu <- (1/Ns) * sum(df.s_i$s_i^2)

    #Grey level non-uniformity normalised
    F_szm$ F_szm_2.5D.glnu.norm <- (1/Ns^2) * sum(df.s_i$s_i^2)

    #Zone size non-uniformity
    F_szm$F_szm_2.5D.zsnu <- (1/Ns) * sum(df.s_j$s_j^2)

    #Zone size non-uniformity normalized
    F_szm$F_szm_2.5D.zsnu.norm <- (1/Ns^2) * sum(df.s_j$s_j^2)

    #Zone percentage
    F_szm$F_zsm_2.5D.z.perc <- Ns/Nv

    #Grey level variance
    p_ij <- df.S$n/sum(df.S$n)
    mu_i <- sum(df.S$i * p_ij)
    F_szm$F_szm_2.5D.gl.var <- sum((df.S$i - mu_i)^2 * p_ij)

    #Zone size variance
    mu_j <- sum(df.S$j * p_ij)
    F_szm$F_szm_2.5D.zs.var <- sum((df.S$j - mu_j)^2 * p_ij)

    #Zone size entropy
    F_szm$F_szm_2.5D.z.entr <- - sum(p_ij * log2(p_ij))

  return(F_szm)
}

