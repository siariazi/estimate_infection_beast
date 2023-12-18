library(xml2)
setwd("/Users/siavashriazi/Desktop/SFU/beast_alignments")
# the name of alignment2 should be ch
# Load the XML file
minSize = 10
iter = 1
file_name = paste(minSize,"_",iter,sep = "")
tree_name = paste(minSize,"_tree_",iter,sep = "")
xml_file <- paste(file_name,".xml",sep = "")  # Replace with your file path
doc <- read_xml(xml_file)

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
  paste("FrequenciesExchanger.s:",tree_name,sep=""),
  paste("FrequenciesExchangerX.s:",tree_name,sep=""),
  paste("BDSIR_serialUniformOperator.t:",tree_name,sep=""),
  paste("BDSIR_serialnarrow.t:",tree_name,sep=""),
  paste("BDSIR_serialwide.t:",tree_name,sep=""),
  paste("BDSIR_serialWilsonBalding.t:",tree_name,sep=""),
  paste("SIR_treeRoot_operatorEs.t:",tree_name,sep=""),
  paste("SIR_tree_operatorEs.t:",tree_name,sep=""),
  paste("strictClockUpDownOperator.c:",tree_name,sep = ""),
  paste("strictClockUpDownOperatorX.c:",tree_name,sep = "")
)

# Find and update the weight attribute for specified operators
all_operator_nodes <- xml_find_all(doc, ".//operator")

for (node in all_operator_nodes) {
  node_id <- xml_attr(node, "id")
  if (node_id %in% operators_to_modify || grepl("SIR_subtreeslide_operator", node_id)) {
    xml_set_attr(node, "weight", "0.0")
  } else if (node_id %in% c(paste("SIR_treeRoot_operatorEs.t:",tree_name,sep=""), paste("SIR_tree_operatorEs.t:",tree_name,sep=""))) {
    affecting_operator <- xml_find_first(node, ".//affectingOperator")
    if (!is.null(affecting_operator)) {
      xml_set_attr(affecting_operator, "weight", "0.0")
    }
  }
}

# Save the modified XML
write_xml(doc,paste(file_name,"_2.xml",sep = "") )  # Replace with your desired output file path
