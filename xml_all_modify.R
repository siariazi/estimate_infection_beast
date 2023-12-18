# this script is a combination of base_xml.R that changes taxa and tree according to fasta file and newick.nwk
# and xml_modify.R that changes the operators to keep the tree constant 
library(xml2)
library(Biostrings)
library(dplyr)
setwd("/Users/siavashriazi/Desktop/SFU/beast_alignments")

minSize <- 10
iter <- 1
origin_value <- 20.0  # Update this with the actual value

file_name <- paste(minSize,"_tree_",iter,sep = "")
fasta_file_name <- paste(file_name,".fasta",sep = "") 
tree_file_name <- paste(file_name,".nwk",sep = "")
xml_template_file <- "template3.xml"
output_xml_file <- paste(minSize,"_",iter,".xml",sep="")


# Function to read a file as a single string
read_file_as_string <- function(file_path) {
  lines <- readLines(file_path)
  fasta_string <- paste(lines, collapse = "\n")
  return(fasta_string)
}

# Function to read Newick tree from file
read_newick_tree <- function(file_path) {
  tree_string <- readLines(file_path, warn = FALSE)
  # Assuming the file contains a single line with the Newick tree
  return(tree_string[1])
}

# Function to parse fasta string and return sequence IDs, taxon values, dates, and sequence values
parse_fasta <- function(fasta_string) {
  lines <- unlist(strsplit(fasta_string, "\n"))
  id_lines <- grep("^>", lines, value = TRUE)
  sequence_lines <- grep("^[^>]", lines, value = TRUE)
  
  ids <- sapply(id_lines, function(x) sub("^>([^|]+).*", "\\1", x))
  taxa_dates <- sapply(id_lines, function(x) {
    matches <- regmatches(x, regexec("^>([^|]+)\\|(.+)", x))
    return(matches[[1]][-1])  # Remove the full match and return the capture groups
  })
  values <- sapply(sequence_lines, function(x) gsub(" ", "", x))
  
  return(list(ids = ids, taxa = taxa_dates[1, ], dates = taxa_dates[2, ], values = values))
}

# Function to update XML with fasta data
update_xml_with_fasta <- function(fasta_data, newick_tree, xml_file, output_file) {
  # Load the XML file
  doc <- read_xml(xml_file)
  
  # Update the <data> id attribute to 'file_name'
  data_node <- xml_find_first(doc, "//data")
  
  # Remove all existing sequence nodes
  xml_find_all(doc, "//sequence") %>% xml_remove()
  
  # Add new sequence nodes based on fasta data
  for (i in seq_along(fasta_data$ids)) {
    xml_add_child(data_node, "sequence", 
                  id = paste0(fasta_data$taxa[i], "|", fasta_data$dates[i]), 
                  spec = "Sequence", 
                  taxon = fasta_data$taxa[i], 
                  totalcount = "4", 
                  value = fasta_data$values[i])
  }
  
  # Construct the value string for the trait node
  trait_value <- paste(paste(fasta_data$taxa, fasta_data$dates, sep="="), collapse=",")
  
  # Find the trait node and update its value
  trait_node <- xml_find_first(doc, "//trait[@id='dateTrait.t:template']")
  xml_set_attr(trait_node, "value", trait_value)
  
  # Update the Newick tree
  init_node <- xml_find_first(doc, "//init[@id='NewickTree.t:template']")
  xml_set_attr(init_node, "newick", newick_tree)
  
  # Update the origin value
  parameter_node <- xml_find_first(doc, "//parameter[@id='originEs.t:template']")
  xml_text(parameter_node) <- as.character(origin_value)
  
  # Save the updated XML
  return(doc)
}
# Read fasta file as string and parse data
fasta_string <- read_file_as_string(fasta_file_name)
fasta_data <- parse_fasta(fasta_string)

# Read Newick tree from file
newick_tree <- read_newick_tree(tree_file_name)
# Update XML template with fasta data
doc = update_xml_with_fasta(fasta_data, newick_tree, xml_template_file, output_xml_file)
# Find all sequence nodes
sequence_nodes <- xml_find_all(doc, ".//sequence")

# Loop through each sequence node and modify the taxon attribute
for (node in sequence_nodes) {
  taxon_attr <- xml_attr(node, "taxon")
  taxon_attr <- gsub("\\|.*", "", taxon_attr)  # Remove everything after '|'
  xml_set_attr(node, "taxon", taxon_attr)
}

# Find the trait node containing the values
trait_node <- xml_find_first(doc, ".//trait[@traitname='date']")

# Extract and modify the value attribute
value_attr <- xml_attr(trait_node, "value")
value_attr <- gsub("\\|.*?=", "=", value_attr, perl = TRUE)  # Remove everything after '|'
xml_set_attr(trait_node, "value", value_attr)

# Define a list of operator IDs that should have their weight set to 0
operators_to_modify <- c(
  "FrequenciesExchanger.s:template",
  "FrequenciesExchangerX.s:template",
  "BDSIR_serialUniformOperator.t:template",
  "BDSIR_serialnarrow.t:template",
  "BDSIR_serialwide.t:template",
  "BDSIR_serialWilsonBalding.t:template",
  "SIR_treeRoot_operatorEs.t:template",
  "SIR_tree_operatorEs.t:template",
  "strictClockUpDownOperator.c:template",
  "strictClockUpDownOperatorX.c:template"
)



# second part, setting operators to 0 
# Find and update the weight attribute for specified operators
all_operator_nodes <- xml_find_all(doc, ".//operator")

for (node in all_operator_nodes) {
  node_id <- xml_attr(node, "id")
  if (node_id %in% operators_to_modify || grepl("SIR_subtreeslide_operator", node_id)) {
    xml_set_attr(node, "weight", "0.0")
  } else if (node_id %in% c("SIR_treeRoot_operatorEs.t:template", "SIR_tree_operatorEs.t:template")) {
    affecting_operator <- xml_find_first(node, ".//affectingOperator")
    if (!is.null(affecting_operator)) {
      xml_set_attr(affecting_operator, "weight", "0.0")
    }
  }
}

# Save the modified XML
write_xml(doc,paste(file_name,"_2.xml",sep = "") )  # Replace with your desired output file path
