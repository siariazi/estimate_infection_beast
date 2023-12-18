library(ape)
###################################
#main functions
###################################
setwd("C:/Siavash/Codes")

break_matrix <- function(mat){
  k = c()
  for(i in 1:dim(mat)[1]){
    if(mat[1,i]==0){
      k = c(k,i)
    }
  }
  # dropping k's rows and columns
  m1 = mat[-k,-k]
  # keeping k's rows and columns
  m2 = mat[k,k]
  output = list(m1,m2)
  return(output)
}

convert_newick <- function(mat){
  #num_zero = 0
  if(sqrt(length(mat)==1)){
    return(paste(":",toString(mat),sep = ""))
  }
  else if(sqrt(length(mat))==2){
    new_mat = mat - min(mat)
    dv = new_mat[which(!new_mat==0)]
    return(paste("(",":",toString(dv[1]),",",":",toString(dv[2]),")",":",toString(min(mat)),sep = ""))
  }
  
  else if(sqrt(length(mat))>2){
    branch_length = min(mat)
    newm = mat - branch_length
    out = break_matrix(newm)
    return(paste("(",convert_newick(out[[1]]),",",convert_newick(out[[2]]),")",":",
                 toString(branch_length),sep = ""))
  }
}

toNewick <- function(dis_matrix){
  out <- convert_newick(dis_matrix)
  return(paste("(",out,")",";",sep = ""))
}

###################################
#checking
###################################
m2 = as.matrix(read.table(file = "data2.txt",sep = ",",header = FALSE))
m4 = as.matrix(read.table(file = "data4.txt",sep = ",",header = FALSE))
m = as.matrix(read.table(file = "data.txt",sep = ",",header = FALSE))
m = as.matrix(read.table(file = "distance_matrix.csv",header = FALSE))
text = toNewick(m)
mytree <- read.tree(text = toNewick(m))
plot(mytree)

################
r1 = c(0.5,0.2)
r2 = c(0.2,0.5)
mat1 = rbind(r1,r2)
text1 = toNewick(mat1)

r1 = c(0.9,0.2,0.5)
r2 = c(0.2,0.9,0.2)
r3 = c(0.5,0.2,0.9)
mat1 = rbind(r1,r2,r3)
text1 = toNewick(mat1)
mytree <- read.tree(text = toNewick(mat1))
plot(mytree)


r1 = c(1.0, 0.2, 0.5, 0.5)
r2 = c(0.2, 1.0, 0.2, 0.2)
r3 = c(0.5, 0.2, 1.0, 0.9)
r4 = c(0.5, 0.2, 0.9, 1.0)
mat2 = rbind(r1,r2,r3,r4)
text2 = toNewick(mat2)
mytree <- read.tree(text = toNewick(mat2))
plot(mytree)


r1 = c(1.0, 0.2, 0.5, 0.5)
r2 = c(0.2, 1.0, 0.2, 0.2)
r3 = c(0.5, 0.2, 1.0, 0.9)
r4 = c(0.5, 0.2, 0.9, 1.0)
mat2 = rbind(r1,r2,r3,r4)
text2 = toNewick(mat2)
mytree <- read.tree(text = toNewick(mat2))
plot(mytree)

r1 = c(1.6, 0.2, 0.5, 0.5, 0.5)
r2 = c(0.2, 1.6, 0.2, 0.2, 0.2)
r3 = c(0.5, 0.2, 1.6, 0.9, 0.9)
r4 = c(0.5, 0.2, 0.9, 1.6, 1.0)
r5 = c(0.5, 0.2, 0.9, 1.0, 1.6)
mat2 = rbind(r1,r2,r3,r4,r5)
text2 = toNewick(mat2)
mytree <- read.tree(text = toNewick(mat2))
plot(mytree)

# this matrix is decomposable:
r1 = c(0.29000000000000004, 0.043, 0.199, 0.2, 0.043, 0.199)
r2 = c(0.043, 0.29900000000000004, 0.043, 0.043, 0.267, 0.043)
r3 = c(0.199, 0.043, 0.29900000000000004, 0.199, 0.043, 0.29000000000000004)
r4 = c(0.2, 0.043, 0.199, 0.29900000000000004, 0.043, 0.199)
r5 = c(0.043, 0.267, 0.043, 0.043, 0.29900000000000004, 0.043)
r6 = c(0.199, 0.043, 0.29000000000000004, 0.199, 0.043, 0.29900000000000004)
mat2 = rbind(r1,r2,r3,r4,r5,r6)
text2 = toNewick(mat2)
mytree <- read.tree(text = toNewick(mat2))
plot(mytree)

# this is not:
r1 = c(0.29000000000000004, 0.043, 0.199, 0.2, 0.043, 0.199, 0.199)
r2 = c(0.043, 0.32300000000000006, 0.043, 0.043, 0.267, 0.043, 0.043)
r3 = c(0.199, 0.043, 0.32300000000000006, 0.199, 0.043, 0.29000000000000004, 0.29000000000000004)
r4 = c(0.2, 0.043, 0.199, 0.29900000000000004, 0.043, 0.199, 0.199)
r5 = c(0.043, 0.267, 0.043, 0.043, 0.32300000000000006, 0.043, 0.043)
r6 = c(0.199, 0.043, 0.29000000000000004, 0.199, 0.043, 0.32300000000000006, 0.29900000000000004)
r7 = c(0.199, 0.043, 0.29000000000000004, 0.199, 0.043, 0.29900000000000004, 0.32300000000000006)
mat2 = rbind(r1,r2,r3,r4,r5,r6,r7)
text2 = toNewick(mat2)
mytree <- read.tree(text = toNewick(mat2))
plot(mytree)
