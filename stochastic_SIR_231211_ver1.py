#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 19:17:52 2023

@author: siavashriazi
"""

import random as rnd
import numpy as np
import matplotlib.pyplot as plt
import copy
import os
from ete3 import Tree
import re

testPars = [0.5,0.05,0.01,0.001,100,1,20] #beta,gamma,psi,sigma,kappa,i0,T
testPars = [0.5,0.05,0.01,0.001,300,1,50] #beta,gamma,psi,sigma,kappa,i0,T

# Here we simulate a tree from a stochastic SIR model to do so we follow a number of variables.
# treeMtrx is a matrix whoes elemnt i,j gives the shared ancestry between linegaes i and j
# state is a vector that indicates the state of each lineage. See states below.  Note that the alive vector is like the state vector but just denotes who is alive and who has either recovered or been sampled.  There should be a better way to make this alive vector but I don't know how.
# t is the time variable (forward in time)
# epiState is a matrix whoes rows give [t, S, I, R] for each event in the stocahstic process
# States: 1: infected, 0:recovered, -1:sampled
# Events: 1: transmission, 2: recovery, 3:sampling, 4:immunity loss, 0: non-event to itterate to present day
# The class sim tree has four functions:
# 1. init sets the inital parameters
# 2. event updates the tree and epidemic given that event e has occured after Deltat time
# 3. gillespie Simulates a sequence of events arising from an stochastic SIR model
# 4. sampledTree Returns the sampled treeMtrx and a vector of sampling times and birth times measured FORWARD in time

# a function that calculate the time to the origin for each lineage
def calculate_lineage_time(node, time_to_root):
    if node.is_leaf():
        lineage_times[node.name] = time_to_root
    else:
        for child in node.children:
            calculate_lineage_time(child, time_to_root + node.dist)

# function that returns indices of elements in a list given a condition             
def find_indices(lst, condition):
    return [i for i, elem in enumerate(lst) if condition(elem)]

# a function that rounds the number in newick text
def round_tree(newick_tree):
    # Define a regular expression pattern to find ": <number>," patterns
    pattern = r':\s*(-?\d+\.\d+)(?=[,)])'

    # Use a lambda function as the replacement to round the number to 3 decimal places
    rounded_tree = re.sub(pattern, lambda match: ': {:.3f}'.format(round(float(match.group(1)), 3)), newick_tree)

    # Remove spaces before and after colons
    rounded_tree = rounded_tree.replace(' :', ':').replace(': ', ':')
    return rounded_tree

class simTree:
    def __init__(self,pars):
        #setting parameters
        self.beta=pars[0]
        self.gamma=pars[1]
        self.psi=pars[2]
        self.sigma = pars[3]
        self.kappa = pars[4]
        self.i0 = pars[5]
        self.T=pars[6]
        
        #setting state
        self.treeMtrx = np.zeros((self.i0, self.i0))

        self.state = [1] * self.i0
        self.alive = [1] * self.i0
        self.epiState=np.array([[0,self.kappa-self.i0,self.i0,0,0]])
        #self.gillespie()#simulate tree
    def event(self,e,Deltat):
        #Add delta t to infected lineages
        self.treeMtrx=np.identity(len(self.treeMtrx))*Deltat*self.alive+self.treeMtrx
        #print(self.treeMtrx)
        if e==1: #infection
            ind=rnd.choice(find_indices(self.state, lambda x: x==1)) #pick parent
            #update tree matrix, state vector, alive vector
            self.treeMtrx=np.vstack((self.treeMtrx,self.treeMtrx[ind])) #add row to tree mtrx
            col=np.transpose(np.hstack((self.treeMtrx[ind],self.treeMtrx[ind,ind])))
            self.treeMtrx=np.vstack((np.transpose(self.treeMtrx),col))#adding column
            #print(self.treeMtrx)
            self.state=self.state+[1]
            self.alive=self.alive+[1]
            #update epiState
            self.epiState=np.vstack((self.epiState,self.epiState[-1]+[Deltat,-1,1,0,0]))
        elif e==2:#recovery
            ind=rnd.choice(find_indices(self.state, lambda x: x==1))# pick lineage to die
            self.state[ind]=0
            self.alive[ind]=0
            #update epiState
            self.epiState=np.vstack((self.epiState,self.epiState[-1]+[Deltat,0,-1,1,0]))
        elif e==3:#samplint
            ind=rnd.choice(find_indices(self.state, lambda x: x==1))# pick lineage for sampling
            self.state[ind]=-1
            self.alive[ind]=0
            #update epiState
            self.epiState=np.vstack((self.epiState,self.epiState[-1]+[Deltat,0,-1,1,1]))
        elif e==4:#waning
            #update epiState
            self.epiState=np.vstack((self.epiState,self.epiState[-1]+[Deltat,1,0,-1,0]))
        elif e==0: #Update to present day *empty event*
            self.epiState=np.vstack((self.epiState,self.epiState[-1]+[Deltat,0,0,0,0]))
        else:
            print("ERROR in event id")
    def gillespie(self):
        #initialize
        t=0
        S=self.epiState[-1,1]
        I=self.epiState[-1,2]
        R=self.epiState[-1,3]
        rates=[self.beta/self.kappa*S*I,self.gamma*I,self.psi*I,self.sigma*R]
        totalRate = sum(rates)
        Deltat=round(np.random.exponential(scale=1/totalRate),3)
        e=rnd.choices(np.linspace(1,len(rates),len(rates)), weights=rates)[0]
        while t+Deltat<self.T:
            #perform event
            self.event(e,Deltat)
            t+=Deltat
            #pick new deltat
            S=self.epiState[-1,1]
            I=self.epiState[-1,2]
            R=self.epiState[-1,3]
            rates=[self.beta/self.kappa*S*I,self.gamma*I,self.psi*I,self.sigma*R]
            totalRate = sum(rates)
            if totalRate==0:
                Deltat=self.T-t
                e=0
            else:
                Deltat=round(np.random.exponential(scale=1/totalRate),3)
                e=rnd.choices(np.linspace(1,len(rates),len(rates)), weights=rates)[0]
        #Last step
        self.event(0,self.T-t)
        self.sampledTree()
    def sampledTree(self):
        # Extracts the sampled tree
        # Extracts the observed sampling times recoded FORWARD in time(yVec)
        # Extracts the observed birth times recoded FORWARD in time (xVec)
        inds=find_indices(self.state, lambda x: x==-1)
        self.sampTree=self.treeMtrx[inds][:,inds]
        self.yVec=np.diagonal(self.sampTree)# sampling times are the diagonal
        # birth times are the (non-duplicated) off diagonals greater than 0
        temp2=np.reshape(np.triu(self.sampTree, k=1),len(self.sampTree)*len(self.sampTree))
        temp2=[x for x in temp2 if x > 0]
        self.xVec=np.array(list(dict.fromkeys(temp2)))
    # a recursive matrix to break the matrix 
    def convert_newick(self,mat):
        if np.shape(mat)[0] == 1:
            #return(":"+str(mat[0][0]))
            return "xAz:" + str(mat[0][0])
        elif np.shape(mat)[0] == 2:
            new_mat = mat - np.amin(mat)
            # dv collects non zero elements of the new mat 
            dv = new_mat[np.nonzero(new_mat)]
            #return("(:"+str(dv[0])+",:"+str(dv[1])+"):"+str(np.amin(mat)))
            return "(xAz:" + str(dv[0]) + ",xAz:" + str(dv[1]) + "):" + str(np.amin(mat))
        elif np.shape(mat)[0] > 2:
            branch_length =  np.amin(mat)
            # substracting min value of all elements
            newm = mat - branch_length
            out = self.break_matrix(newm)
            return "(" + self.convert_newick(out[0])  + "," + self.convert_newick(out[1]) + "):" + str(branch_length)

    # break matrix breaks the matrix to two matrices
    def break_matrix(self,mat):
        mat2 = copy.deepcopy(mat)
        k = []
        for i in range(np.shape(mat2)[0]):
            if mat2[0][i] == 0:
                k.append(i)
            #print(i)
        m1 = np.delete(mat2,k,1)
        m1 = np.delete(m1,k,0)
        m2 = mat[np.ix_(k,k)]
        output = [m1,m2]
        return output

    # toNweick outputs the final result
    def toNewick(self):
        out = self.convert_newick(self.sampTree)
        self.treeTxt = "("+out+");"
        self.treeTxt = round_tree(self.treeTxt)
        #self.treeTxt = "("+out+");"
    
    def add_label(self):
        j = 1
        textl = list(self.treeTxt)
        new_textl = []
        label_list = []

        for char in textl:
            new_textl.append(char)
            if char == 'A':
                new_textl.append(str(j))
                label_list.append("A" + str(j))
                j += 1
                
        label_list.append("A0")
        self.treeTxtL = ''.join(new_textl)


minSize = 40
maxSize = 110
i = 1
it = 1

while i < it+1:
    
    tree1=simTree(testPars)
    tree1.gillespie()
    
    if minSize <= np.shape(tree1.sampTree)[0] < maxSize:
        
        tree1.toNewick()
        tree1.add_label()
        newickTree = tree1.treeTxtL
        i+=1 
    else:
        pass
    
    
tree1=simTree(testPars)
tree1.gillespie()       
tree1.toNewick()
#round_tree(tree1.treeTxt)
np.shape(tree1.sampTree)
tree1.add_label()
newickTree = tree1.treeTxtL
newickTree
# plotting stochastic simulation
sPlot, = plt.plot(tree1.epiState[:,0],tree1.epiState[:,1], color='orange', label="Suscepible")
iPlot, = plt.plot(tree1.epiState[:,0],tree1.epiState[:,2], color='red', label="Infected")
rPlot, = plt.plot(tree1.epiState[:,0],tree1.epiState[:,3], color='purple', label="Recovered")
plt.legend(handles=[sPlot,iPlot, rPlot])
#plt.plot(tree1.epiState[:,0],tree1.epiState[:,1],tree1.epiState[:,0],tree1.epiState[:,2],tree1.epiState[:,0],tree1.epiState[:,3])
plt.xlabel('forward time (t)')
plt.ylabel('Conuts')
plt.show()

fileDirectory = '/Users/siavashriazi/Desktop/SFU/iqtree-2.2.2.6-MacOSX/'
writeFileDir = '/Users/siavashriazi/Desktop/SFU/beast_alignments/'
treeFileName = f"{minSize}_tree_{it}.nwk"
alignmentFileName = f"{minSize}_tree_{it}.phy"
alignmentWriteName = f"{minSize}_tree_{it}.fasta"
treeWriteName = f"{minSize}_tree_{it}.nwk"

newickTree = newickTree.replace('xA0z;', ';')
# Parse the Newick tree, it doesn't work because it has the root, but we need the root for iqtree
t = Tree(newickTree)
print(t)

# Calculate the time to root for each lineage
lineage_times = {}

calculate_lineage_time(t, 0)  # Start with the root node and initial time of 0

lineage_times['xA0z'] = 0
            
# Print the lineage times
for lineage, time in lineage_times.items():
    print(f"Lineage {lineage}: Time to Root = {time}")
    

tree2 =  newickTree.replace(';', 'xA0z;')   
# writing the tree to a text file
treeFilePath = os.path.join(fileDirectory, treeFileName)
file1 = open(treeFilePath, "w")
file1.write(tree2)
file1.close()

################################################################################
# after converting the tree to alignment by iqtree, the file needs to be trimmed
alignmentFilePath = fileDirectory + alignmentFileName
file2 = open(alignmentFilePath, "r")
content = file2.read()
file2.close()
# deleting all the content before A
content2 = content[content.index("x"):]
# removing whitepsaces
content2 = content2.strip()

# Split the string based on '\nxA0z'
parts = content2.split('\nxA0z')

# Take the part before '\nxA0z'
content2= parts[0]

# Define a character to insert between the sequence name and the sequence
insert_char = '|'
# Define the character to add after the sequence name
suffix_char = 'y'

# Iterate over lineage times and update the content2 string
for lineage, time in lineage_times.items():
    if lineage in content2:
        content2 = content2.replace(lineage, f'{lineage}{insert_char}{time:.3f}{suffix_char}', 1)
# Iterate over lineage times and update the content2 string


content2 = content2.replace('x','>')
content2 = content2.replace('z','')
content2 = content2.replace('y','\n')
alignmentWritePath = writeFileDir + alignmentWriteName
file3 = open(alignmentWritePath, "w")
file3.write(content2)
file3.close()


# after converting the alignment to xml file using beauti2 I write a clean version of newick tree in BEAST folder
tree3 = tree2.replace('x','')
tree3 = tree3.replace('z','')
# Find the last occurrence of ':'
lastColonIndex = tree3.rfind(':')
# Remove everything after the last ':'
tree4 = tree3[:lastColonIndex]
tree5 = tree4[1:]
treeWritePath = writeFileDir + treeWriteName
file4 = open(treeWritePath, "w")
file4.write(tree5)
file4.close()
