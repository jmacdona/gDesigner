# (c) Michael Crone and James T. MacDonald, 2017
import gRNA_design.env
import gRNA_design.offinder_functions
import gRNA_design.melting_functions
import gRNA_design.pairfold_functions
import getopt, sys

# example command line: python grna_design.py -p Paths.txt -s Settings.txt -i r2o_634395117458942835 -o output.fa

paths = gRNA_design.env.Paths.Instance()
settings = gRNA_design.env.Settings.Instance()

paths_filename = ""
settings_filename = ""
input_r2o_dir = ""
output_filename = ""

letters = 'i:p:s:o:'
opts, params = getopt.getopt(sys.argv[1:], letters)
for o, p in opts:
    if o == "-p":
        paths_filename = p
    if o == "-s":
        settings_filename = p
    if o == "-i":
        input_r2o_dir = p
    if o == "-o":
        output_filename = p

if paths_filename == "" or settings_filename == "" or input_r2o_dir == "" or output_filename == "":
    sys.exit('Error: -p -s -i -o options are required to be set')

# initialise settings
paths.read_file(paths_filename)
settings.read_file(settings_filename)
settings.output_root_dir = input_r2o_dir    # output to same r2o directory

input_seq_filename = input_r2o_dir + "/" + "sequencesFinalFile.fa"
#output_hits_filename = "sdfasdfads"

# run offinder hits search in forward direction
output_fwd_hits_filename = gRNA_design.offinder_functions.run_offinder(input_seq_filename, "fwd")
output_fwd_filt_fasta_file = settings.output_root_dir + "/" + settings.crispr_type + "_" + "fwd" + "_" \
                           + "offinder_filtered.fa"
gRNA_design.offinder_functions.offinder_filter_sequences(input_seq_filename, output_fwd_hits_filename, output_fwd_filt_fasta_file, "fwd")

# run offinder hits search in reverse direction
output_rev_hits_filename = gRNA_design.offinder_functions.run_offinder(output_fwd_filt_fasta_file, "rev")
output_rev_filt_fasta_file = settings.output_root_dir + "/" + settings.crispr_type + "_" + "rev" + "_" \
                           + "offinder_filtered.fa"
gRNA_design.offinder_functions.offinder_filter_sequences(output_fwd_filt_fasta_file, output_rev_hits_filename, output_rev_filt_fasta_file, "rev")

# run MELTING
output_melting_filt_fasta_file = settings.output_root_dir + "/" + settings.crispr_type + "_" \
                           + "melting_filtered.fa"
gRNA_design.melting_functions.filter(output_rev_filt_fasta_file, output_melting_filt_fasta_file)

# pairfold network elimination
gRNA_design.pairfold_functions.filter_net(output_melting_filt_fasta_file, output_filename)






