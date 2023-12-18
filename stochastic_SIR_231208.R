# a script based on stochastic_SIR_231205.R that goes through a loop to simulate mutliple trees
library(ggplot2)
library(ape)

testPars <- c(0.5,0.05,0.01,0.001,100,1,20) # beta, gamma, psi, sigma, kappa, i0, mT
beta <- testPars[1]
gamma <- testPars[2]
psi <- testPars[3]
sigma <- testPars[4]
kappa <- testPars[5]
i0 <- testPars[6]
mT <- testPars[7]

# a function that returns the indices of a vector given a condition 
find_indices <- function(lst, condition) {
  indices <- which(condition(lst))
  return(indices)
}

event <- function(e,deltat){
  treeMtrx <<- diag(nrow(treeMtrx)) * deltat * alive + treeMtrx
  # infection event
  if (e == 1) {
    ind <- sample(which(state == 1), size = 1)  # pick parent index
    
    # Update treeMtrx by adding row and column
    newRow <- treeMtrx[ind, ]
    treeMtrx <<- rbind(treeMtrx, newRow)
    newCol <- c(newRow, treeMtrx[ind, ind])
    treeMtrx <<- cbind(treeMtrx, newCol)
    
    # Update state and alive vectors
    state <<- c(state, 1)
    alive <<- c(alive, 1)
    
    # Update epiState
    newState <- c(epiState[nrow(epiState), ] + c(deltat, -1, 1, 0))
    epiState <<- rbind(epiState, newState)
  }
  # recovery
  else if (e == 2){
    ind = sample(which(state == 1), size = 1)   # pick parent index
    state[ind] <<- 0
    alive[ind] <<- 0
    # Update epiState
    newState <- c(epiState[nrow(epiState), ] + c(deltat, 0, -1, 1))
    epiState <<- rbind(epiState, newState)
  }
  # sampling
  else if (e == 3){
    ind = sample(which(state == 1), size = 1)
    state[ind] <<- -1
    alive[ind] <<- 0
    # Update epiState
    newState <- c(epiState[nrow(epiState), ] + c(deltat, 0, -1, 1))
    epiState <<- rbind(epiState, newState)
  }
  # waning immunity 
  else if (e == 4){
    # Update epiState
    newState <- c(epiState[nrow(epiState), ] + c(deltat, 1, 0, -1))
    epiState <<- rbind(epiState, newState)
  }
  #Update to present day *empty event*
  else if (e == 0){
    # Update epiState
    newState <- c(epiState[nrow(epiState), ] + c(deltat, 0, 0, 0))
    epiState <<-  rbind(epiState, newState)
  }
  else {
    cat("error in event id", "\n")
  }
}

gillespie <- function(){
  
  # setting the states
  t <- 0
  S <- epiState[1,2]
  I <- epiState[1,3]
  R <- epiState[1,4]
  # calculating rates
  rates <- c(beta/kappa * S * I, gamma * I, psi * I, sigma * R)
  
  # computing total rate
  totalRate <- sum(rates)
  
  # computing deltat
  deltat <- round(rexp(1, rate = totalRate), 3)
  
  # picking the event proportional to its weight (rate)
  choices <- 1:length(rates)
  e <- sample(choices, size = 1, prob = rates/totalRate)
  
  # Perform the Gillespie algorithm steps
  while (t + deltat < mT) {
    #perform event
    event(e,deltat)
    t <- t + deltat
    #pick new deltat
    S <- epiState[nrow(epiState),2]
    I <- epiState[nrow(epiState),3]
    R <- epiState[nrow(epiState),4]
    rates <- c(beta/kappa*S*I,gamma*I,psi*I,sigma*R)
    totalRate <- sum(rates)
    if (totalRate == 0) {
      deltat <- mT - t
      e <- 0
    }
    else {
      deltat <- round(rexp(1, rate = totalRate), 3)
      choices <- 1:length(rates)
      e <- sample(choices, size = 1, prob = rates/totalRate)
    }
  }
}



# Plotting the results
#ggplot(epiData, aes(x = time)) +
#  geom_line(aes(y = S, color = "Susceptible")) +
#  geom_line(aes(y = I, color = "Infectious")) +
#  geom_line(aes(y = R, color = "Recovered")) +
#  labs(x = "Time", y = "Population", color = "Status") +
#  theme_minimal()

sampledTree <- function(){
  
  # Find indices where state equals -1 (sampled)
  inds <- which(state == -1)
  # Create the sampTree sub-matrix from treeMtrx based on the indices
  sampTree <- treeMtrx[inds, inds]
  
  return(sampTree)
}

# functions to convert a distance matrix to Newick tree
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
  if(sqrt(length(mat)==1)){
    return(paste("xAz:",toString(mat),sep = ""))
  }
  else if(sqrt(length(mat))==2){
    new_mat = mat - min(mat)
    dv = new_mat[which(!new_mat==0)]
    return(paste("(xAz:",toString(dv[1]),",xAz:",toString(dv[2]),"):",toString(min(mat)),sep = ""))
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
  out <- paste("(",out,");",sep = "")
  out <- addLabel(out)
  return(out)
}
addLabel <- function(text){
  j <- 1
  textl <- strsplit(text, "")[[1]]
  labelList <- c()
  
  for (i in seq_along(textl)) {
    if (textl[i] == "A") {
      jChar <- as.character(j)
      textl <- append(textl, jChar, after = i)
      labelList <- c(labelList, paste0("A", j))
      j <- j + 1
    }
  }
  
  labelList <- c(labelList, "A0")
  treeTxtL <- paste(textl, collapse = "")
  return(treeTxtL)
}

########################################
i = 1
# Create an empty data frame with specified column names
trees <- data.frame(iteration = numeric(), size = numeric())

# removing the first row that was duplicated
minSize = 5

# How many I want to put per that size
iter = 2
setwd("/Users/siavashriazi/Desktop/SFU/iqtree-2.2.2.6-MacOSX") 
epiFileDir = '/Users/siavashriazi/Desktop/SFU/beast_alignments/'
while (i <= iter){
  treeMtrx <- matrix(0, nrow = i0, ncol = i0)
  state <- rep(1, i0)
  alive <- rep(1, i0)
  # initializing the states, I duplicate the first row and then I delete it 
  epiState <- matrix(rep(c(0, kappa - i0, i0, 0),2), nrow = 2, ncol = 4, byrow = TRUE)
  gillespie()
  smpTree <- sampledTree()
  if (dim(smpTree)[1] > minSize){
    epiState2 <- epiState[-1,]
    # converting epiState to dataframe
    epiData <- as.data.frame(epiState2)
    colnames(epiData) <- c("time","S","I","R")
    rownames(epiData) <- seq(0,dim(epiState2)[1]-1)
    # naming the SIR (epiData) file name based on the iteration and minSize
    epiFileName = paste(minSize,"_SIR_",i,".csv",sep = "")
    epiFilePath = paste(epiFileDir, epiFileName, sep = "")
    #write.csv(epiData,epiFilePath)
    # recording the tree size in its data frame
    trees[i,]$iteration <- i
    trees[i,]$size <- dim(smpTree)[1]
    #epiStateList[[i]] <- epiState2
    treeTxt <- toNewick(smpTree)
    # Replace ';' with 'xA0z;'
    treeTxt2 <- gsub(";", "xA0z;", treeTxt)
    RootedFileName <- paste(minSize,"_tree_",i,".nwk",sep = "")
    writeLines(treeTxt2, RootedFileName)
    i = i+1
  }
  else{}
}

# naming the trees file name based on the min size 
treesFileName = paste(minSize,"_trees.csv",sep = "")
row.names(trees) <- seq(1,iter)
treesFilePath = paste(epiFileDir,treesFileName,sep = "")
write.csv(trees,treesFilePath)

# after this point, the trees should be converted by iqtree and then 
# use python script to convert_fasta_231209.py to add date to the alignemnt
# then the reuslt of the python code should be oppened by BEAUti and then xml_modify.R