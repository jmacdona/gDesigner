# gDesigner - a pipeline for the of design synthetic gRNAs for CRISPRi

gDesigner is a Python pipeline for designing synthetic gRNAs for CRISPRi applications. The pipeline takes an input list of DNA sequences in FASTA format, and after various screening steps, returns a list of orthogonal gRNAs. We have used sequences generated using [R2oDNA Designer](http://www.r2odna.com/), but a suitable library of FASTA sequences from any source could be used (e.g. from [NUPACK design](http://www.nupack.org/design) algorithm). Example R2oDNA Designer job specifications input files (see the example_r2o_jobSpecifications/ directory) generated FASTA sequences are included with this repository (see the directory ADD).

The associated publication for this software can be seen [here](http://) (gDesigner: computational design of synthetic gRNAs for Cpf-1-based transcriptional repression in mammalian cells (2018). Michael Crone, James T. MacDonald &ast;, Paul S. Freemont &ast;, and Velia Siciliano &ast;).

## Requirements

### Operating systems

We have only tested this software on UNIX-style operating systems  - Red Hat Enterprise Linux 7, and MacOS. The software will currently likely not run on Windows operating systems due to the difference in the way paths are specified and they way shell works.

### Python

The code is written for Python 2.7+ and requires the NumPy numerical, Biopython bioinformatics and NetworkX graph libraries.

### MELTING

MELTING is an open source GPL licenced program for calculating the free energies of nucleic acid duplexes, including DNA:RNA duplexes. The current stable version (v5.1.1) be obtained from [here](http://www.ebi.ac.uk/biomodels/tools/melting/).

### Cas-OFFinder

Cas-OFFinder is an open source BSD licenced software tool for searching genomes for potential off-target hits. It can be obtained from [here](http://www.rgenome.net/cas-offinder/portable) and the source code can be obtained from the [GitHub repository](https://github.com/snugel/cas-offinder).

### MultiRNAFold 2.0

The PairFold program from the MulitRNAFold software suite is used to calculate pairwise gRNA:gRNA free energies and can be obtained from [here](http://www.rnasoft.ca/download.html).

### Genomes to be scanned for off-target sequences.

Any genome sequences in FASTA format to be scanned for potential off-targets should be placed in the genomes/ directory. The E. coli BL21 and MG1655 genomes sequences are provided in the repository because they are reasonably small. Larger genomes, such as as the human genome, will have to downloaded and placed in this folder by the user as they are too large to be placed in this repository. Any FASTA files in the genomes/ will be automatically scanned.

## Installation

Once all the software dependencies are installed the Python pipeline needs to know where to find all the executables. The Paths.txt file in the main directory is a template example that can be edited for your particular installation locations and follows a two column tab-separated key/value format (be careful not to leave any extraneous extra whitespace characters in this file).

For example this Paths.txt file:
```
grna_design_install_dir	/home/jmacdona/dev/gDesigner/
offinder_exe	/home/jmacdona/dev/gDesigner/external/bin/cas-offinder
pairfold_exe	/home/jmacdona/dev/gDesigner/external/MultiRNAFold-2.0/pairfold
melting_exe	/home/jmacdona/dev/gDesigner/external/RNAMelting/executable/melting
```
tells the software the the gDesigner software is installed in the directory  /home/jmacdona/dev/gDesigner/, the offinder binary executable can be found at /home/jmacdona/dev/gDesigner/external/bin/cas-offinder, etc.

## Running the pipeline

The software is set up to take the output from R2oDNA Designer as the input into this pipeline but any FASTA format input sequences can be used. Here, we will assume you are using R2oDNA Designer generated sequences.
 
Firstly, unpack the emailed zip file from R2oDNA Designer into an appropriate directory, e.g. in the installation directory of gDesigner. You now have a directory called something like r2o_XXXXXXXXXXXXXXXXXX where the X's are some sequence of random numbers.

The program can then be run with the command line:

```
python grna_design.py -p Paths.txt -s Settings.txt -i r2o_XXXXXXXXXXXXXXXXXX -o grna_output.fa
```

where the Path.txt file is described above, the -i argument is location of the unpacked R2oDNA directory from the unpacked zip file, the -o argument is the output FASTA file for the final set of orthogonal sequences, and the Settings.txt contains the constraints the sequences need to satisfy in a tab-separated key/value format. Settings parameters are summarised in the following table:

| Parameter  | Value type | Description |
| ------------- | ------------- | ------------- |
| crispr_type | String (lbcpf1 or ascpf1 or cas9) | The CRISPR system being used |
| offinder_compute_type | Character (C or G or A) | The Cas-OFFinder processor type (depends on your Cas-OFFinder installation, we recommend using graphics card acceleration using CUDA if possible)|
| min_offinder_score_fwd | Float (0 to 1) | Off-target score in the forward direction |
| min_offinder_score_rev | Float (0 to 1) | Off-target score in the reverse strand direction (probably unnecessary for most use cases so can be set to 0) |
| offinder_max_mismatch | Integer | Cas-Offinder to return only hits with less than or equal to this number of mismatches|
| melting_kconc | Float | K+ concentration for use in MELTING RNA:DNA hybrid free energy calculations|
| melting_mgconc | Float | Mg2+ concentration for use in MELTING RNA:DNA hybrid free energy calculations |
| melting_naconc | Float | Na+ for use in MELTING RNA:DNA hybrid free energy calculations |
| melting_GFE_fwd_min | Float | MELTING RNA:DNA hybrid free energy minimum cutoff |
| melting_GFE_fwd_max | Float | MELTING RNA:DNA hybrid free energy maximum cutoff |
| melting_GFE_rev_min | Float | MELTING RNA:DNA hybrid free energy minimum cutoff in the reverse strand direction (probably unimportant for most use cases so can be set to a wide range) |
| melting_GFE_rev_max | Float | MELTING RNA:DNA hybrid free energy maximum cutoff in the reverse strand direction (probably unimportant for most use cases so can be set to a wide range) |
| max_submatch | Integer | Maximum length of substring match between different gRNA spacer sequences (used in maximum independent set search step) |
| pf_inter_cutoff | Float | Minimum Pairfold folding free energy between different gRNA molecules (used in maximum independent set search step) |
| sw_cutoff_for | Float | Smith-Waterman (ednafull matrix) maximum score cutoff (used in maximum independent set search step) |
| sw_cutoff_rev | Float | Smith-Waterman (ednafull matrix) maximum score cutoff in the reverse strand direction (used in maximum independent set search step) |

Some examples Settings files are provided in the distribution (e.g. see Settings_lbcpf1.txt).

## Examples

To filter synthetic gRNA for lbCpf1 use this command line (all example input files are provided in this repository):

```
python grna_design.py -p Paths.txt -s Settings_lbcpf1.txt -i r2o_example_lbcpf1 -o output_grna_lbcpf1.fa
```

Where r2o_example_lbcpf1 is an unzipped directory output from [R2oDNA Designer](http://www.r2odna.com/) and the final set of sequences is written to the file output_grna_lbcpf1.fa.

## License

This software is distributed under the GNU GPL license, version 3.

&copy; Michael Crone and James T. MacDonald
