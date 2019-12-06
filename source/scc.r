#message(c("step:", 0))
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
    #message(c("smoothing:", i))
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
  #print("hire2")
  grp <- match(gdist, seq_len(ub) * resol)
  idx <- split(seq_len(length(gdist)), grp)
  #print("hire3")
  st = sapply(idx, est.scc)
  corr0 = unlist(st[1,])
  wei0 = unlist(st[2,])
  #print("hire4")
  corr = corr0[!is.na(corr0)]
  wei = wei0[!is.na(wei0)]
  #print("hire5")
  scc = corr %*% wei/sum(wei)
  #print("hire6")
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
            if((abs(i-j)<1500000/binsize-h)&&abs(i-j)>2 + h){
            smd_matr[i,j] = mean(matr[rlb[i]:rrb[i], clb[j]:crb[j]])
            }
        }
    }

    colnames(smd_matr)=colnames(dat)
    rownames(smd_matr)=rownames(dat)

    return(smd_matr)
}




binsize_def <- function(M0){
  buf = unique(M0[order(M0[,1]),1])
  for (i in 1:(length(buf)-1)){
    buf[i] = buf[i + 1] - buf[i]
  }
  buf = buf[1:(length(buf)-1)]
  binsize = min(buf)
  
  return(binsize)
}

prep1 <- function(R1, R2, resol, h, max){
  #print("start_prep")
  #message(c("start_prep:", 0))
  pro_Rep1 = R1[,-c(1,2,3)]
  rownames(pro_Rep1) = colnames(pro_Rep1) = R1[,3]-resol/2
  #print("start_prep_1")
  #message(c("start_prep_1:", 0))
  pro_Rep2 = R2[,-c(1,2,3)]
  rownames(pro_Rep2)=colnames(pro_Rep2)=R2[,3]-resol/2
  #print("start_prep_2")
  if(h == 0){
    vec_Rep1=MatToVec(pro_Rep1)
    ##print("vec_Rep1")
    ##print(vec_Rep1)
    vec_Rep2=MatToVec(pro_Rep2)
  } else {
    #print("start_prep_3")
    #message(c("start_prep_3:", 0))
    smt_Rep1 = smoothMat1(pro_Rep1, h)
    smt_Rep2 = smoothMat1(pro_Rep2, h)
    vec_Rep1 = MatToVec(smt_Rep1)
    vec_Rep2 = MatToVec(smt_Rep2)
  }
  #message(c("start_prep_4:", 0))
  #print("start_prep_4")
  comb = data.frame(vec_Rep1,vec_Rep2[,3])
  colnames(comb) = c("V1", "V2", "V3", "V4")
  eidx = which(comb[,3] == 0 & comb[,4] == 0)
  #print("start_prep_5")
  if (length(eidx) == 0) {
    filt = comb
  } else {
    filt = comb[-eidx,]
  }
  
  return(filt)
}


chrs_length = c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566)
#message(c("step:", 0))
#print("0")



args = commandArgs(trailingOnly=TRUE)
in_fname = args[1]
out_file = args[2]
h = as.numeric(args[3])
message(c("h=",h))
chr = args[4]
promoters = args[4]
enhancers = args[5]
message(c("in_fname =:", in_fname))
#chr <- "1"
#cell_type1 <- "NHEK"
#cell_type2 <- "IMR90"
#m1 <- paste(chr,".5MB.",cell_type1,".contacts", sep = "")
#m2 <- paste(chr,".5MB.",cell_type2,".contacts", sep = "")
#m1 <- paste(chr,".5MB.",cell_type1,".contacts.gz", sep = "")
#m2 <- paste(chr,".5MB.",cell_type2,".contacts.gz", sep = "")
M1 <- read.table(in_fname,head=T)
M2 <- read.table(in_fname,head=T)
M2[,3] = M1[,4]
M1[is.na(M1)] <- 0
M2[is.na(M2)] <- 0
binsize = binsize_def(M1)
chrs_length = chrs_length%/%binsize
maxx = max(max(M1[,2]),max(M2[,2]))
minn = min(min(M1[,1]),min(M2[,1]))
mx = max(max(M1[,2]),max(M2[,2]))/binsize
mn = min(min(M1[,1]),min(M2[,1]))/binsize
#print(M1[1,3])
wide = binsize*1500000/binsize
maxdist = maxx - minn - 2*binsize
#print(maxdist)
#message(c("maxM1:", max(M1[,2])))
#message(c("maxM2:", max(M2[,2])))
#message(c(":", length(M1[,3])))
#message(c(":", length(M2[,3])))
size <- maxx
Z2 <- seq(0,size - binsize + binsize,by = binsize)
Z1 <- cbind(1:length(Z2))
Z3 <- seq(binsize ,size + binsize,by = binsize)
Z <- cbind(Z1,Z2,Z3)
Z4 <- matrix(0, nrow = length(Z2), ncol = length(Z2))
Z5 <- matrix(0, nrow = length(Z2), ncol = length(Z2))
Z6 <- matrix(0, nrow = length(Z2), ncol = length(Z2))
#print("1")
	

for (t in 1:length(M1[,3])) {
  r = (M1[t,1])/binsize + 1
  c = (M1[t,2])/binsize + 1
  if(abs(r-c)<1500000/binsize&&abs(r-c)>2){
    Z4[r,c] = M1[t,3]
    Z4[c,r] = M1[t,3]}
}

#message(c("step:", 0))
for (t in 1:length(M2[,3])) {
  r = (M2[t,1])/binsize + 1
  c = (M2[t,2])/binsize + 1
  if(abs(r-c)<1500000/binsize&&abs(r-c)>2){
    Z5[r,c] = M2[t,3]
    Z5[c,r] = M2[t,3]}
}

#pr_en
mse_pr_en = 0
mod_avr_pr_en = 0
mod_ps_pr_en = 0
k_pr_en = 0
k_mod_pr_en = 0
k_mod_ps_pr_en = 0

mse = 0
mod_avr = 0
mod_ps = 0
k = 0
k_mod = 0
k_mod_ps = 0
if (!is.na(promoters)&&!is.na(enhancers)){
Pr <- read.table(promoters,head=T)
En <- read.table(enhancers,head=T)
message(c("scc =:", 0))
for(t in which(En[,1]==chr)) {
  for(l in which(Pr[,1]==chr)) {
    r = Pr[l,2]%/%binsize + 1
    c = En[t,2]%/%binsize + 1
    if(r<mx&&c<mx&&r>mn&&c>mn&&abs(r-c)<1500000/binsize&&abs(r-c)>2){
      Z6[r,c] = 1
      Z6[c,r] = 1}
  }}
message(c("scc =:", 1))

#evaluate easy metrics


for (r in 1:length(Z4[,1])){
  for (c in 1:length(Z4[,1])){
    #pr_en
    if(abs(r-c)<=1500000/binsize&&abs(r-c)>2&&Z6[r,c]==1){
      mse_pr_en = mse_pr_en + (Z4[r,c] - Z5[r,c])*(Z4[r,c] - Z5[r,c])
      mod_avr_pr_en = mod_avr_pr_en + abs(Z4[r,c] - Z5[r,c])
      if(Z4[r,c]!=0){
        mod_ps_pr_en = mod_ps_pr_en + abs((Z4[r,c] - Z5[r,c])/Z4[r,c])
        k_mod_ps_pr_en = k_mod_ps_pr_en + 1}
      if(r>mn&&r<mx&&c>mn&&c<mx){k_pr_en = k_pr_en + 1}
      if(abs(Z4[r,c])+abs(Z5[r,c]!=0)){k_mod_pr_en = k_mod_pr_en + 1}
    }
  }}

#pr_en
mse_pr_en = sqrt(mse_pr_en/k_pr_en)
mod_avr_pr_en = mod_avr_pr_en/k_mod_pr_en
mod_ps_pr_en = mod_ps_pr_en/k_mod_pr_en
message(c("mse_pr_en",mse_pr_en))
message(c("mod_avr_pr_en",mod_avr_pr_en))
message(c("mod_ps_pr_en",mod_ps_pr_en))
}
if (1==1){
#evaluate easy metrics


for (r in 1:length(Z4[,1])){
  for (c in 1:length(Z4[,1])){
    #all
    if(abs(r-c)<=1500000/binsize&&abs(r-c)>2){
      mse = mse + (Z4[r,c] - Z5[r,c])*(Z4[r,c] - Z5[r,c])
      mod_avr = mod_avr + abs(Z4[r,c] - Z5[r,c])
      if(Z4[r,c]!=0){mod_ps = mod_ps + abs((Z4[r,c] - Z5[r,c])/Z4[r,c])
      k_mod_ps = k_mod_ps + 1}
      if(r>mn&&r<mx&&c>mn&&c<mx){k = k + 1}
      if(abs(Z4[r,c])+abs(Z5[r,c]!=0)){k_mod = k_mod + 1}
    }
  }}


#all
mse = sqrt(mse/k)
mod_avr = mod_avr/k
mod_ps = mod_ps/k_mod_ps
message(c("mse",mse))
message(c("mod_avr",mod_avr))
message(c("mod_ps",mod_ps))
}
#correlation
H1 <- cbind(Z,Z4)
H2 <- cbind(Z,Z5)

#pearson
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


#scc  
#h_hat <- htrain1(H1, H2, binsize, maxdist, 0:10)
processed <- prep1(H1, H2, binsize, h, maxdist)
j = get.scc1(processed, binsize, maxdist)
message(c("scc =:", j$scc))
message(c("cor =:", c))
out_fname = paste(in_fname,"out",sep=".")
#if (length(M1[1,])>4){
#  out1 = paste("corr","scc","MSE","mod_avr","Pr_en_MSE","Pr_en_mod_avr","mod_ps","mod_ps_pr_en",sep = "	")
#  out2 = paste(c,as.numeric(j$scc), mse, mod_avr,mse_pr_en,mod_avr_pr_en,mod_ps,mod_ps_pr_en, sep =  "	")
#}
out1 = paste("corr","scc","MSE","mod_avr","Pr_en_MSE","Pr_en_mod_avr","mod_ps","mod_ps_pr_en", sep = "	")
out2 = paste(c,as.numeric(j$scc), mse, mod_avr,mse_pr_en,mod_avr_pr_en,mod_ps,mod_ps_pr_en, sep =  "	")
out = paste(out1,"\n",out2,sep = "")
write(out,out_file,sep = " ")