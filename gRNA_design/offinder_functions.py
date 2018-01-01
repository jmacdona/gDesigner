# (c) Michael Crone and James T. MacDonald, 2017
import env
import subprocess
from Bio import SeqIO
import sys
import numpy as np
import csv


def run_offinder(input_seq_filename, direction):
    paths = env.Paths.Instance()
    settings = env.Settings.Instance()

    inputf = open(input_seq_filename, 'r')

    formatted_input_file = settings.output_root_dir + "/" + settings.crispr_type + "_" + direction + "_" \
                           + "offinder.seq"
    output_hits_filename = settings.output_root_dir + "/" + settings.crispr_type + "_" + direction + "_" \
                           + "offinder.hits"

    genomes_loc = paths.grna_design_install_dir + "/genomes/"

    var_length = settings.spacer_pos_end - settings.spacer_pos_start
    var_Ns = "N" * var_length
    pattern = ""
    if (settings.pam_five_prime == 1):
        pattern = settings.pam_seq + var_Ns
    else:
        pattern = var_Ns + settings.pam_seq

    pam_len = len(settings.pam_seq)
    pam_Ns = "N" * pam_len

    form_in_f = open(formatted_input_file, 'w+')
    form_in_f.write(genomes_loc + "\n")
    form_in_f.write(pattern + "\n")

    num_seqs = 0

    if direction == "fwd":
        for record in SeqIO.parse(inputf, "fasta"):
            num_seqs += 1
            if (settings.pam_five_prime == 1):
                settings.reformat_spacer_pos_start = pam_len
                settings.reformat_spacer_pos_end = settings.reformat_spacer_pos_start + var_length
                form_in_f.write(pam_Ns + str(record.seq[settings.spacer_pos_start:settings.spacer_pos_end])
                                + " " + str(settings.offinder_max_mismatch) + "\n")
            else:
                settings.reformat_spacer_pos_start = 0
                settings.reformat_spacer_pos_end = var_length
                form_in_f.write(str(record.seq[settings.spacer_pos_start:settings.spacer_pos_end]) + pam_Ns
                                + " " + str(settings.offinder_max_mismatch) + "\n")
    elif direction == "rev":
        for record in SeqIO.parse(inputf, "fasta"):
            num_seqs += 1
            if (settings.pam_five_prime == 1):
                settings.reformat_spacer_pos_start = pam_len
                settings.reformat_spacer_pos_end = settings.reformat_spacer_pos_start + var_length
                form_in_f.write(pam_Ns + str(
                    record.seq[settings.spacer_pos_start:settings.spacer_pos_end].reverse_complement()) + " "
                                + str(settings.offinder_max_mismatch) + "\n")
            else:
                settings.reformat_spacer_pos_start = 0
                settings.reformat_spacer_pos_end = var_length
                form_in_f.write(str(record.seq[settings.spacer_pos_start:settings.spacer_pos_end].reverse_complement())
                                + pam_Ns + " " + str(settings.offinder_max_mismatch) + "\n")
    else:
        sys.exit('Error: unrecognised direction: ' + direction + "\n")

    form_in_f.close()

    print "INFO: Generated cas-offinder input file with: " + str(num_seqs) + " sequences"

    cmd_args = []

    cmd_args.append(paths.offinder_exe)
    cmd_args.append(formatted_input_file)
    cmd_args.append(settings.offinder_compute_type)
    cmd_args.append(output_hits_filename)

    print("cas_offinder command is: " + str(cmd_args))
    result = subprocess.check_output(cmd_args)
    print result

    return output_hits_filename


# score a sequence
def offinder_scoring_fun(mismatches, weights):
    n_ = len(weights) - 1
    if len(mismatches) == 0:
        score = 100
        return score
    else:
        score = 100 * (1 - weights[mismatches]).prod()
        if len(mismatches) > 1:
            mean_pairwise = float(np.sum(mismatches[1:]) - np.sum(mismatches[:-1])) / (len(mismatches) - 1)
            # changed 19 to 22 to represent number of bases
            mpw_factor = ((float((n_ - mean_pairwise)) / n_) * 4 + 1)
            scl_factor = pow(len(mismatches), 2)
            score = score / (mpw_factor * scl_factor)
            score = max([score, 0])
    return score


def offinder_filter_sequences(in_fasta_file, offinder_hits_file, output_fasta_file, direction):

    paths = env.Paths.Instance()
    settings = env.Settings.Instance()

    aim_score = 0.8
    if direction == "fwd":
        aim_score = settings.min_offinder_score_fwd
    elif direction == "rev":
        aim_score = settings.min_offinder_score_rev
    else:
        sys.exit('Error: unrecognised direction: ' + direction + "\n")

    weights = settings.offinder_score_weights  # np.array([])

    sequences = []
    mismatches = []

    with open(offinder_hits_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for seq, acc, geno, mm, s, nummiss in reader:
            mmarray = []
            icount = 0
            for char in mm[4:27]:
                if char.islower():
                    mmarray.append(icount)
                icount += 1
                score = offinder_scoring_fun(mmarray, weights)  # scoring_fun(mmarray)
            sequences.append(seq)
            mismatches.append(score)

    outputfinal = open(output_fasta_file, 'w')

    num_seqs = 0

    for record in SeqIO.parse(in_fasta_file, "fasta"):
        countmm = 0
        score = 0
        addscore = 0
        totalscore = 0
        mcount = 0
        for row in sequences:
            fasta_var_seq = str(record.seq[settings.spacer_pos_start:settings.spacer_pos_end])
            if (direction == "rev"):
                fasta_var_seq = str(record.seq[settings.spacer_pos_start:settings.spacer_pos_end].reverse_complement())
            off_var_seq = row[settings.reformat_spacer_pos_start:settings.reformat_spacer_pos_end]
            #print fasta_var_seq + " " + off_var_seq
            if fasta_var_seq == off_var_seq:
                addscore = addscore + mismatches[mcount]
                countmm += 1
            mcount += 1
        totalscore = 100 / (100 + addscore)
        if totalscore > settings.min_offinder_score_fwd:
            #print(str(record.id) + '\t' + str(totalscore) + '\t' + str(countmm))
            outputfinal.write(">" + str(record.id) + "\n")
            outputfinal.write( str(record.seq) + "\n")
            num_seqs += 1

    print "INFO: Output: " + str(num_seqs) + " filtered sequences"




