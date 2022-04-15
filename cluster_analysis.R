purity<-function(standard, in_cluster){
  label_s <- unique(standard)
  label_c <- unique(in_cluster)
  sum <- 0
  bs <- matrix(ncol=length(label_c),nrow=length(label_s))
  for(i in 1:length(label_s)){
    a <- which(standard == label_s[i])
    for(j in 1:length(label_c)){
      b <- which(in_cluster == label_c[j])
      bs[i,j] <- length(intersect(a,b))
    }
  }
  for(k in 1:length(label_s)){
    sum <- sum+max(bs[k,])
  }
  purity <- sum/length(standard)
  return(purity)
}

ari<-function(standard, in_cluster){
  label_s <- unique(standard)
  label_c <- unique(in_cluster)
  nlen_s <- length(label_s)
  nlen_c <- length(label_c)
  contingency_table <- matrix(ncol = nlen_s, nrow = nlen_c)
  for(i in 1:nlen_c){
    a <- label_c[i]
    v2 <- which(a == cluster_c)
    for (j in 1:nlen_s){
      b<-label_s[j]
      v1<-which(b == cluster_s)
      contingency_table[i,j] <- length(intersect(v1,v2))
    }
  }
  n<-sum(contingency_table)
  a<-colSums(contingency_table)
  b<-rowSums(contingency_table)
  suma<-sum(choose(a,2))
  sumb<-sum(choose(b,2))
  sumn<-sum(choose(contingency_table,2))
  part2<-(suma*sumb)/choose(n,2)
  ari<-(sumn-part2)/(0.5*suma+0.5*sumb-part2)
  return(ari)
}

jaccard<-function(standard, in_cluster){
  label_s <- unique(standard)
  label_c <- unique(in_cluster)
  A <- matrix(nrow=length(standard),ncol=length(standard))
  for(i in 1:length(label_s)){
    a <- which(standard == label_s[i])
    A[a,a]=1
  }
  A[is.na(A)] <- 0
  B <- matrix(nrow = length(cluster_in), ncol = length(cluster_in))
  for(i in 1:length(label_c)){
    a <- which(cluster2 == label_c[i])
    B[a,a]=1
  }
  B[is.na(B)] <- 0
  A_1 <- which(A==1)
  A_0 <- which(A==0)
  B_1 <- which(B==1)
  B_0 <- which(B==0)
  
  a<-length(intersect(A_1,B_1))
  b<-length(intersect(A_1,B_0))
  c<-length(intersect(A_0,B_1))
  jaccard=a/(a+b+c)
  return(jaccard)
}

H<-function(prob_count){
  n<-sum(prob_count)
  H<-0
  prob_info<-(-1)*(prob_count/n)*log(prob_count/n)
  H<-sum(prob_info)
  return(H)
}

I<-function(contingency_table){
  prob_count1<-colSums(contingency_table)
  prob_count2<-rowSums(contingency_table)
  n<-sum(prob_count1)
  cal_table<-matrix(0,ncol=ncol(contingency_table),nrow = nrow(contingency_table))
  contingency_table[contingency_table ==0] <- 0.000001
  for(i in 1:nrow(cal_table)){
    for (j in 1:ncol(cal_table)){
      cal_table[i,j]<-(contingency_table[i,j]/n)*log(n*contingency_table[i,j]/(prob_count2[i]*prob_count1[j]))
    }
  }
  I<-sum(cal_table)
  return(I)
}

nmi<-function(standard, in_cluster){
  num<-length(standard)
  label_s<-unique(standard)
  label_c<-unique(in_cluster)
  row_contigency<-length(label_s)
  col_contigency<-length(label_c)
  # row for standard types
  contingency_table<-matrix(0,nrow = row_contigency, ncol = col_contigency)
  for (i in 1:num){
    row<-as.integer(standard[i])+1
    col<-as.integer(in_cluster[i])+1
    contingency_table[row,col]<-contingency_table[row,col]+1
  }
  #H(U)
  prob_count1<-colSums(contingency_table)
  prob_count2<-rowSums(contingency_table)
  nmi<-2*I(contingency_table)/(H(prob_count1)+H(prob_count2))
  return(nmi)
}