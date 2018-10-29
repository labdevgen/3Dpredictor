MatToVec <- function(dat){
  
  mat = as.matrix(dat)
  
  nc = ncol(mat)
  rc = nrow(mat)
  
  test = matrix(0, nc*rc, 3)
  test[,3] = as.vector(mat)
  test[,2] = as.double(rep(rownames(mat), nc))
  
  tmp = rep(as.double(colnames(mat)), each=rc)
  test[,1] = tmp
  return(test)
}


vstran  <- function(d){
  
  x1r = rank(d[,1], ties.method = "random")
  x2r = rank(d[,2], ties.method = "random")
  x1.cdf.func = ecdf(x1r); x2.cdf.func = ecdf(x2r)
  x1.cdf = x1.cdf.func(x1r)
  x2.cdf = x2.cdf.func(x2r)
  new_d = cbind(x1.cdf, x2.cdf)
  
  return(new_d)
}



depth.adj = function(d, size, resol, out = 0){
  
  cd = d[,-c(1,2,3)]
  rownames(cd) = colnames(cd) = d[,3]-resol/2
  
  temp = MatToVec1(cd)
  p1 = temp[,3]/sum(temp[,3])+.Machine$double.eps
  
  subrd = sample(nrow(temp), size, prob=p1, replace=TRUE)
  freq = table(subrd)
  idx = as.double(names(freq))
  vec = as.vector(freq)
  temp[,3] = 0
  temp[idx,3] = vec
  
  ##turn it back to matrix
  
  ntemp = temp[which(temp[,3]!=0),]
  ntemp[,1] = (ntemp[,1]+resol/2)/resol
  ntemp[,2] = (ntemp[,2]+resol/2)/resol
  cd[cd>0] = 0
  cd[ntemp[,c(1,2)]] = ntemp[,3]
  
  cdm = cbind(d[,c(1,2,3)], cd)
  colnames(cdm) = rownames(cdm) = NULL
  if(out == 1){return(temp)}
  else return(cdm)
}




htrain1 <- function(R1, R2, resol, max, range){
  
  corr = matrix(0, max(range)+1, 2)
  corr[,1] = range
  for (i in range){
    message(c("smoothing:", i))
    pre = prep1(R1, R2, resol, i)
    s_cor = array()
    for (j in 1:10){
      idx = sample(1:nrow(pre), floor(nrow(pre)*0.1), replace=FALSE)
      sub = pre[idx,]
      s_cor[j] = get.scc1(sub, resol, max)[[3]]
    }
    corr[i+1, 2] = round(mean(s_cor),4)
    if (i > 0){
      if ((corr[i+1,2] - corr[i,2])<0.01){
        break
      }
    }
  }
  if (i == max(range)){
    warning("Note: It's likely that your searching range is too narrow. 
                    Try to expand the range and rerun it")
  }
  return(corr[i,1])
}


get.scc1 <- function (dat, resol, max){
  
  ub <- floor(max/resol)
  corr <- array(ub)
  cov <- array(ub)
  wei <- array(ub)
  n <- array(ub)
  
  gdist = abs(dat[,2]-dat[,1])
  #print("hire1")
  est.scc = function(idx){
    
    if (length(idx) != 0){
      
      n = length(idx)
      ffd = dat[idx,c(3,4)]
      nd = vstran(ffd)
      
      if (length(unique(ffd[,1])) != 1 & length(unique(ffd[,2])) != 1) {
        corr = cor(ffd[,1], ffd[,2])
        wei = sqrt(var(nd[,1])*var(nd[,2]))*n
      } else {
        corr = NA
        wei = NA
      }
    } else {
      corr = NA 
      wei = NA
    }
    
    return(list(corr = corr, wei = wei))
  }
  grp <- match(gdist, seq_len(ub) * resol)
  idx <- split(seq_len(length(gdist)), grp)

  st = sapply(idx, est.scc)
  corr0 = unlist(st[1,])
  wei0 = unlist(st[2,])

  corr = corr0[!is.na(corr0)]
  wei = wei0[!is.na(wei0)]

  scc = corr %*% wei/sum(wei)

  std = sqrt(sum(wei^2*(1-corr^2)^2/(n-3))/(sum(wei))^2)
  
  return(list(corr = corr, wei = wei, scc = scc, std = std))
}

smoothMat1 <- function(dat, h){

    matr = as.matrix(dat)

    c = ncol(matr)
    r = nrow(matr)

    smd_matr = matrix(0,r,c)

    i <- seq_len(r)
    rlb <- ifelse(i - h > 0, i - h, 1)
    rrb <- ifelse(i + h < r, i + h, r)

    j <- seq_len(c)
    clb <- ifelse(j - h > 0, j - h, 1)
    crb <- ifelse(j + h < c, j +  h, c)

    for (i in seq_len(r)){
		#message(c("i:", i))
        for (j in seq_len(c)){
            if(abs(r-c)<60 - h){
            smd_matr[i,j] = mean(matr[rlb[i]:rrb[i], clb[j]:crb[j]])}
        }
    }

    colnames(smd_matr)=colnames(dat)
    rownames(smd_matr)=rownames(dat)

    return(smd_matr)
}

prep1 <- function(R1, R2, resol, h, max){

  pro_Rep1 = R1[,-c(1,2,3)]
  rownames(pro_Rep1) = colnames(pro_Rep1) = R1[,3]-resol/2

  pro_Rep2 = R2[,-c(1,2,3)]
  rownames(pro_Rep2)=colnames(pro_Rep2)=R2[,3]-resol/2

  if(h == 0){
    vec_Rep1=MatToVec(pro_Rep1)

    vec_Rep2=MatToVec(pro_Rep2)
  } else {

    smt_Rep1 = smoothMat1(pro_Rep1, h)
    smt_Rep2 = smoothMat1(pro_Rep2, h)
    vec_Rep1 = MatToVec(smt_Rep1)
    vec_Rep2 = MatToVec(smt_Rep2)
  }

  comb = data.frame(vec_Rep1,vec_Rep2[,3])
  colnames(comb) = c("V1", "V2", "V3", "V4")
  eidx = which(comb[,3] == 0 & comb[,4] == 0)

  if (length(eidx) == 0) {
    filt = comb
  } else {
    filt = comb[-eidx,]
  }
  
  return(filt)
}

prep1_half_smoof <- function(R1, R2, resol, h, max){

  pro_Rep1 = R1[,-c(1,2,3)]
  rownames(pro_Rep1) = colnames(pro_Rep1) = R1[,3]-resol/2

  pro_Rep2 = R2[,-c(1,2,3)]
  rownames(pro_Rep2)=colnames(pro_Rep2)=R2[,3]-resol/2

  if(h == 0){
    vec_Rep1=MatToVec(pro_Rep1)

    vec_Rep2=MatToVec(pro_Rep2)
  } else {

    smt_Rep1 = smoothMat1(pro_Rep1, h)
    smt_Rep2 = pro_Rep2
    vec_Rep1 = MatToVec(smt_Rep1)
    vec_Rep2 = MatToVec(smt_Rep2)
  }

  comb = data.frame(vec_Rep1,vec_Rep2[,3])
  colnames(comb) = c("V1", "V2", "V3", "V4")
  eidx = which(comb[,3] == 0 & comb[,4] == 0)

  if (length(eidx) == 0) {
    filt = comb
  } else {
    filt = comb[-eidx,]
  }
  
  return(filt)
}

loop_calc <- function(M1){
  max1 <- max(M1[,2])
  min1 <- min(M1[ ,1])
  min2 <- min(M1[ ,2])
  wide <- max(M1[,2] - M1[,1])
  maxdist = max1 - min1 - 2*binsize
  size <- max1 - min1
  Z2 <- seq(0,size - binsize,by = binsize)
  Z1 <- cbind(1:length(Z2))
  Z3 <- seq(binsize ,size,by = binsize)
  Z <- cbind(Z1,Z2,Z3)
  Z4 <- matrix(0, nrow = length(Z2), ncol = length(Z2))
  for (t in 1:length(M1[,3])) { 
    r = (M1[t,1]-min1)/25000
    c = (M1[t,2]-min2)/25000
    if(abs(r-c)<60 && M1[t,5] == 1){
      Z4[r,c] = M1[t,3]
      Z4[c,r] = M1[t,3]}}
  H1_loop	<- cbind(Z,Z4)
  
	Z5 <- matrix(0, nrow = length(Z2), ncol = length(Z2))
	for (t in 1:length(M1[,4])) {
	  r = (M1[t,1]-min1)/25000
	  c = (M1[t,2]-min2)/25000
	  if(abs(r-c)<60 && M1[t,5] == 1){
	  Z5[r,c] = M1[t,4]
	  Z5[c,r] = M1[t,4]}
	}
	H2_loop <- cbind(Z,Z5)
  H <- prep1(H1_loop, H2_loop, binsize, 0, maxdist)
  corr_loop <- cor(as.double(H[,3]),as.double(H[,4]))
  j_loop = get.scc1(H, binsize, maxdist)
  return(paste(" scc_loop_only = ",as.numeric(j_loop$scc)," corr_loop_only = ",corr_loop,sep = " "))
}



args = commandArgs(trailingOnly=TRUE)
in_fname = args[1]
h = as.numeric(args[2])
message(c("in_fname =:", in_fname))
binsize <- 25000
M1 <- read.table(in_fname,head=T)
M1 <- as.matrix.data.frame(M1)
M1 <- as.matrix.data.frame(M1)

max1 <- max(M1[,2])
min1 <- min(M1[,1])
wide <- max(M1[,2] - M1[,1])
maxdist = max1 - min1 - 2*binsize

size <- max1 - min1
Z2 <- seq(0,size - binsize,by = binsize)
Z1 <- cbind(1:length(Z2))
Z3 <- seq(binsize ,size,by = binsize)
Z <- cbind(Z1,Z2,Z3)
if (length(M1[1,])>4){loop <- loop_calc(M1)} else loop = ""
message(c("loop = ",loop))
Z4 <- matrix(0, nrow = length(Z2), ncol = length(Z2))
Z5 <- matrix(0, nrow = length(Z2), ncol = length(Z2))
min1 <- min(M1[ ,1])
min2 <- min(M1[ ,2])
k = 0
euc_sq = 0
euc_mod_avr = 0
for (t in 1:length(M1[,3])) { 
  r = (M1[t,1]-min1)/25000
  c = (M1[t,2]-min2)/25000
  if(abs(r-c)<60){
  Z4[r,c] = M1[t,3]
  Z4[c,r] = M1[t,3]
  Z5[r,c] = M1[t,4]
  Z5[c,r] = M1[t,4]
  }
}

#norm
s= 0
k= 0
for (t in 1:length(Z4[,1])) { 
  s= s + sum(as.numeric(Z4[t,]))
  k = k + 1
}
s = s/k/10000
Z4 <- Z4[,]/s
#norm
s= 0
k= 0
for (t in 1:length(Z5[,1])) { 
  s= s + sum(as.numeric(Z5[t,]))
  k = k + 1
}
s = s/k/10000
Z5 <- Z5[,]/s


for (t in 1:length(M1[,3])) { 
  r = (M1[t,1]-min1)/25000
  c = (M1[t,2]-min2)/25000
  if(abs(r-c)<60){
    euc_sq = euc_sq + (M1[t,3] - M1[t,4])*(M1[t,3] - M1[t,4])
    euc_mod_avr = euc_mod_avr + abs(M1[t,3] - M1[t,4])
    k = k + 1
  }
}

euc_sq = sqrt(euc_sq)
euc_sq_avr = sqrt(euc_sq)/k
euc_mod_avr = euc_mod_avr/k
H1 <- cbind(Z,Z4)
H2 <- cbind(Z,Z5)
k <- 0
c <- 0
for (t in 1:length(Z2)) {

  m = cor(as.double(Z5[t,]),as.double(Z4[t,]))

  if (!is.na(m)) {
    c <- c + as.double(cor(as.double(Z5[t,]),as.double(Z4[t,])))

    k <- k + 1
  }
}

c = c/k



#h_hat <- htrain1(H1, H2, binsize, maxdist, 0:10)
processed <- prep1(H1, H2, binsize, h, maxdist)
j = get.scc1(processed, binsize, maxdist)

message(c("scc =:", j$scc))
message(c("cor =:", c))
out_fname = paste(in_fname,"out",sep=".")
message(c("loop = ", loop))
write.table(paste("corr = ",c,"scc = ",as.numeric(j$scc),loop, " euc_sq = ", euc_sq," euc_sq_avr = ",euc_sq_avr, " euc_mod_avr = ", euc_mod_avr, sep=" "),out_fname)
