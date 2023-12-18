# estimate_infection_beast
Files for generating an epidemic with a SIR model, making a tree(s), converting tree to fasta file and use it in BEAST to estimate the infection rate with phylodynamic BD-SIR. \
stochastic_SIR_231208.R: a R file to generate trees in a loop \
Iqtree code: \ \
cd /Users/siavashriazi/Desktop/SFU/iqtree-2.2.2.6-MacOSX \
bin/iqtree2 --alisim 15_tree_1 -m JC -t 15_tree_1.nwk \
\
Loop through values of i from 1 to 3 \
for i in {1..3}; do \
    # Formulate the filenames using the value of i \
    alisim_file="15_tree_$i" \
    nwk_file="15_tree_$i.nwk" \
    # Run the command with the modified filenames \
    bin/iqtree2 --alisim "$alisim_file" -m JC -t "$nwk_file" \
done \ 
\
convert_fasta_231209.py: After converting trees to fasta file with iqtree, this file adds the dates to fasta files, clean them and save the final fasta files and newick trees \
stochastic_SIR_231211_ver1.py: a python file that generates a single tree, and after converting the tree to a fasta file it adds the dates and trim the fasta file. \
xml_modify.R: a script that after loading the fasta file, modify the .xml file to remove the tree operators to keep the tree constant. \
base_xml.R: a script to go update a template .xml file with data and from a fasta file, this script doesn't remove tree operators. \
xml_all_modify: a combination of xml_modify.R and base_xml.R to update a template .xml file with data from a fasta file and remove the tree operators to keep the tree constant. \
plot_BDSIR.R: a script from BEAST2 website that analysis the result of a .log file from running phylodynamic BD-SIR pakcage. \
220612_convert_dism_newick.R: file from Yexuan to convert a distance matrix to a newick tree. 
