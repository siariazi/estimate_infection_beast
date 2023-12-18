library(xml2)
library(Biostrings)
library(dplyr)
setwd("/Users/siavashriazi/Desktop/SFU/beast_alignments")

minSize <- 10
iter <- 1
file_name <- paste(minSize,"_tree_",iter,sep = "")
fasta_file_name <- paste(file_name,".fasta",sep = "") 
tree_file_name <- paste(file_name,".nwk",sep = "")
xml_template_file <- "template.xml"
output_xml_file <- paste(minSize,"_",iter,".xml",sep="")

# Function to read a file as a single string
read_file_as_string <- function(file_path) {
  # Read lines from file
  lines <- readLines(file_path)
  
  # Combine lines into a single string
  fasta_string <- paste(lines, collapse = "\n")
  
  return(fasta_string)
}

# Read fasta file as string
fasta_string <- read_file_as_string(fasta_file_name)


# Parse fasta IDs
fasta_ids <- parse_fasta(fasta_string)

# Function to read a file as a single string
read_file_as_string <- function(file_path) {
  # Read lines from file
  lines <- readLines(file_path)
  
  # Combine lines into a single string
  fasta_string <- paste(lines, collapse = "\n")
  
  return(fasta_string)
}

# Function to parse fasta string and return sequence IDs and values
parse_fasta <- function(fasta_string) {
  # Split the string into lines and filter out lines that start with '>'
  lines <- unlist(strsplit(fasta_string, "\n"))
  id_lines <- grep("^>", lines, value = TRUE)
  sequence_lines <- grep("^[^>]", lines, value = TRUE)
  
  # Extract IDs and sequences
  ids <- sapply(id_lines, function(x) sub("^>([^ ]+).*", "\\1", x))
  values <- sapply(sequence_lines, function(x) gsub(" ", "", x)) # Remove spaces
  
  return(list(ids = ids, values = values))
}

# Function to update XML with fasta data
update_xml_with_fasta <- function(fasta_data, xml_file, output_file) {
  # Load the XML file
  doc <- read_xml(xml_file)
  
  # Remove all existing sequence nodes
  xml_find_all(doc, "//sequence") %>% xml_remove()
  
  # Add new sequence nodes based on fasta data
  data_node <- xml_find_first(doc, "//data") # Assuming there's only one <data> node
  for (i in seq_along(fasta_data$ids)) {
    new_node <- xml_add_child(data_node, "sequence", id = fasta_data$ids[i], spec = "Sequence", taxon = fasta_data$ids[i], totalcount = "4", value = fasta_data$values[i])
  }
  
  # Update the root (or relevant parent) node's id to 'file_name'
  xml_set_attr(xml_root(doc), "id", "file_name")
  
  # Save the updated XML
  write_xml(doc, output_file)
}

# Read fasta file as string and parse data
fasta_string <- read_file_as_string(fasta_file_name)
fasta_data <- parse_fasta(fasta_string)

# Update XML template with fasta data
update_xml_with_fasta(fasta_data, xml_template_file, output_xml_file)
#############################################################
