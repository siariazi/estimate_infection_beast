#library(devtools)
#install_github("snailvet/OOPidemic", auth_token = "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")
directory <- "/Users/siavashriazi/SFU/beast_alignments"
setwd(directory)
library(OOPidemic)
library(Biostrings)

# set up a reference strain with a randomized genome
ref_strain <- ReferenceStrain$new(
  name = "ref_strain",
  g_len = 1000,
  mut_rate = 0.016,
  dna = TRUE
)

# set up a group 
group <- Group$new(
  id = 1,
  ref_strain = ref_strain,
  init_inf = 1,
  init_sus = 1000,
  inf_rate = 0.35,
  inc_shape = 0, # turns off exposure compartment
  find_infector_method = "serial",
  si_shape = 5, si_rate = 2,
  rec_shape = 8, rec_rate = 1 # mean is: shape/rate oddly when I increase shape epidemic goes up! 
)

# set up a lab to take samples 
lab <- Lab$new()

# run simluation
group$run_simulation(lab)

sample_host_ids <- vapply(lab$wg_sequences, function(wgs) wgs$host$id, numeric(1L))
inf_host_ids <- vapply(group$recovered_hosts(), function(r_host) r_host$id, numeric(1L))

all(sample_host_ids %in% inf_host_ids)

cd_plot <- plot_compartment_dynamics(group)

cd_plot

fasta_df <- lab$fasta_df()
head(fasta_df)

metadata_df <- lab$metadata_df()
head(metadata_df)

hostdata_df <- lab$hostdata_df()
head(hostdata_df)

# Concatenate the name and collection_time with a '|' separator
combined_names <- paste(fasta_df$name, metadata_df$collection_time, sep = "|")

# Convert the dataframe to a DNAStringSet or AAStringSet object
sequences <- DNAStringSet(fasta_df$genome)

# Set the names of the sequences
names(sequences) <- combined_names

tree_size = dim(metadata_df)[1] 
tree_size
# Write to a FASTA file
writeXStringSet(sequences, filepath = paste("tree_",tree_size,".fasta",sep = ""))

# writing a csv file
write.csv(metadata_df,file=paste("tree_",tree_size,"_meta.csv",sep = ""))
