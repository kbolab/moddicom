######################################
### BEGIN FUNCION ####################
######################################


#' @import radiomics data.table

glrlmTexturalFeaturesMerged <- function(imgObj,n_grey){

  # compute number of non-NA voxels
  Nv <- numeric()
  sumtot <- list()
  for (i in 1:dim(imgObj)[3]){
    R_list <- list()
    if(length(imgObj[,,i])*.005 >= sum(!is.na(imgObj[,,i]))) next
    imgObj[,,i] <- round(imgObj[,,i])
    #if(length(table(unique(imgObj[,,i])))<n_grey) n_grey <- length(table(unique(imgObj[,,i])))-1
    if (min(imgObj[,,i],na.rm = T) < 0) {imgObj[,,i] <- imgObj[,,i] + abs(min(imgObj[,,i],na.rm = T))}
    R_list[[1]] <- as.matrix(glrlm(imgObj[,,i], angle = 0, verbose=F,truncate = F, n_grey = n_grey))
    R_list[[2]] <- as.matrix(glrlm(imgObj[,,i], angle = 45, verbose=F,truncate = F, n_grey = n_grey))
    R_list[[3]] <- as.matrix(glrlm(imgObj[,,i], angle = 90, verbose=F,truncate = F, n_grey = n_grey))
    R_list[[4]] <- as.matrix(glrlm(imgObj[,,i], angle = 135, verbose=F,truncate = F, n_grey = n_grey))
    Nv[i] <- length(imgObj[,,i][which(!is.na(imgObj[,,i]))])*4
    matrix.df <- ldply(R_list, data.table::melt, varnames=c("row", "col"))
    sumtot[[i]] <- acast(matrix.df, row ~ col, sum)
  }

  #elimino gli elementi NULL della lista
  if(length(sumtot)!=0){
    if(length(which(sapply(sumtot, is.null)))!=0){
      sumtot = sumtot[-which(sapply(sumtot, is.null))]
    }
  }

  #Initialise data table for storing GLCM features; I have added a few
  featNames <- c("F_rlm_merged.sre","F_rlm_merged.lre","F_rlm_merged.lgre","F_rlm_merged.hgre","F_rlm_merged.srlge",
                 "F_rlm_merged.srhge","F_rlm_merged.lrlge","F_rlm_merged.lrhge","F_rlm_merged.glnu","F_rlm_merged.glnu.norm","F_rlm_merged.rlnu","F_rlm_merged.rlnu.norm", "F_rlm_merged.r.perc",
                 "F_rlm_merged.gl.var","F_rlm_merged.rl.var","F_rlm_merged.rl.entr")
  F_rlm <- data.table(matrix(NA, nrow=length(sumtot), ncol=length(featNames)))
  colnames(F_rlm) <- featNames

  #Iterate over grey level cooccurrence matrices
  #The idea is that basically every GLCM is the same, i.e. we can just perform the same operations on every glcm.
  for (iter in seq_len(length(sumtot))){

    # if(iter==1 | iter==5 | iter==9 | iter==13) {    Nv <- length(imgObj[,,1][which(!is.na(imgObj[,,1]))]) }
    # else if(iter==2 | iter==6 | iter==10 | iter==14) {    Nv <- length(imgObj[,,2][which(!is.na(imgObj[,,2]))]) }
    # else if(iter==3 | iter==7 | iter==11 | iter==15) {    Nv <- length(imgObj[,,3][which(!is.na(imgObj[,,3]))]) }
    # else if(iter==4 | iter==8 | iter==12 | iter==16) {    Nv <- length(imgObj[,,4][which(!is.na(imgObj[,,4]))]) }

    Ns <- sum(sumtot[[iter]],na.rm=T)

    #Convert matrix to data table
    df.R    <- data.table(sumtot[[iter]])

    #Add row grey level intensity
    df.R$i <- as.numeric(row.names(sumtot[[iter]]))

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
    F_rlm$F_rlm_merged.sre[iter]         <- (1/Ns) * sum(df.r_j$r_j/(df.r_j$j^2))

    #Long runs emphasis
    F_rlm$F_rlm_merged.lre[iter]         <- (1/Ns) * sum(df.r_j$r_j*(df.r_j$j^2))

    #Low grey level run emphasis
    F_rlm$F_rlm_merged.lgre[iter]         <- (1/Ns) * sum(df.r_i$r_i/(df.r_i$i^2))

    #High grey level run emphasis
    F_rlm$F_rlm_merged.hgre[iter]         <- (1/Ns) * sum(df.r_i$r_i*(df.r_i$i^2))

    #Short run low grey level emphasis
    F_rlm$F_rlm_merged.srlge[iter] <- (1/Ns) * sum(df.R$n/((df.R$i^2)*(df.R$j^2)))

    #Short run high grey level emphasis
    F_rlm$F_rlm_merged.srhge[iter] <- (1/Ns) * sum((df.R$n)*(df.R$i^2)/(df.R$j^2))

    #Long run low grey level emphasis
    F_rlm$F_rlm_merged.lrlge[iter] <- (1/Ns) * sum((df.R$n)*(df.R$j^2)/(df.R$i^2))

    #Long run high grey level emphasis
    F_rlm$F_rlm_merged.lrhge[iter] <- (1/Ns) * sum((df.R$n)*(df.R$j^2)*(df.R$i^2))

    #Grey level non-uniformity
    F_rlm$F_rlm_merged.glnu[iter] <- (1/Ns) * sum(df.r_i$r_i^2)

    #Grey level non-uniformity normalized
    F_rlm$F_rlm_merged.glnu.norm[iter] <- (1/Ns^2) * sum(df.r_i$r_i^2)

    #Run length non-uniformity
    F_rlm$F_rlm_merged.rlnu[iter] <- (1/Ns) * sum(df.r_j$r_j^2)

    #Run length non-uniformity normalized
    F_rlm$F_rlm_merged.rlnu.norm[iter] <- (1/Ns^2) * sum(df.r_j$r_j^2)

    #Run percentage

    F_rlm$F_rlm_merged.r.perc[iter] <- Ns/Nv[iter]

    #Grey level variance
    p_ij <- df.R$n/sum(df.R$n)
    mu_i <- sum(df.R$i * p_ij)
    F_rlm$F_rlm_merged.gl.var[iter] <- sum((df.R$i - mu_i)^2 * p_ij)

    #Run length variance
    mu_j <- sum(df.R$j * p_ij)
    F_rlm$F_rlm_merged.rl.var[iter] <- sum((df.R$j - mu_j)^2 * p_ij)

    #Run entropy
    F_rlm$F_rlm_merged.rl.entr[iter] <- - sum(p_ij * log2(p_ij))
  }

  return(F_rlm)
}
