#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 17:10:14 2023

@author: siavashriazi
A script that takes the alignemnt results from iqtree and stochastic_SIR_122308.R and convert it to 
fasta files with dates, clean the newick tree and save all of them in beast_alignemnt 
"""

from ete3 import Tree

# Calculate the time to root for each lineage
lineage_times = {}

def calculate_lineage_time(node, time_to_root):
    if node.is_leaf():
        lineage_times[node.name] = time_to_root
    else:
        for child in node.children:
            calculate_lineage_time(child, time_to_root + node.dist)


fileDirectory = '/Users/siavashriazi/Desktop/SFU/iqtree-2.2.2.6-MacOSX/'
writeFileDir = '/Users/siavashriazi/Desktop/SFU/beast_alignments/'
minSize = 5

for it in range(1,4):
    treeFileName = f"{minSize}_tree_{it}.nwk"
    # Define the file path
    treeFilePath = fileDirectory + treeFileName

    # Read the Newick tree from the file
    with open(treeFilePath, 'r') as file:
        newickTree = file.readline().strip()
        
    # Parse the Newick tree, it doesn't work because it has the root, but we need the root for iqtree
    # Remove the trailing label "xA0z"
    newickTree = newickTree.replace('xA0z;', ';')

    tree = Tree(newickTree)
    #print(tree)

    # Calculate the time to root for each lineage
    lineage_times = {}

    calculate_lineage_time(tree, 0)  # Start with the root node and initial time of 0

    alignmentFileName = f"{minSize}_tree_{it}.phy"
    alignmentFilePath = fileDirectory + alignmentFileName
    # after converting the tree to alignment by iqtree, the file needs to be trimmed
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
            content2 = content2.replace(lineage, f'{lineage}{insert_char}{time:.6f}{suffix_char}', 1)


    content2 = content2.replace('x','>')
    content2 = content2.replace('z','')
    content2 = content2.replace('y','\n')
    fastaFileName = f"{minSize}_tree_{it}.fasta"
    fastaFilePath = writeFileDir + fastaFileName 
    file3 = open(fastaFilePath, "w")
    file3.write(content2)
    file3.close()


    # after converting the alignment to xml file using beauti2 I write a clean version of newick tree in BEAST folder
    tree2 = newickTree.replace('x','')
    tree3 = tree2.replace('z','')
    # Find the last occurrence of ':'
    lastColonIndex = tree3.rfind(':')
    # Remove everything after the last ':'
    tree4 = tree3[:lastColonIndex]
    tree5 = tree4[1:]
    newTreeFileName = f"{minSize}_tree_{it}.txt"
    newTreeFilePath = writeFileDir + newTreeFileName
    file4 = open(newTreeFilePath, "w")
    file4.write(tree5)
    file4.close()
    
# after this point the files should be open with BEAUti be named like minSize_iter.xml and then 
# be clean by xml_modify.R