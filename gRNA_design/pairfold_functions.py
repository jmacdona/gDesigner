# (c) Michael Crone and James T. MacDonald, 2017
import env
import subprocess
from Bio import SeqIO
#import sys
#import math
import networkx as nx
import networkx.algorithms as nx_alg
from Bio import pairwise2
#from Bio.SubsMat import MatrixInfo as matlist


def has_submatch(seq1, seq2, max_submatch):
    if len(seq1.seq) < max_submatch:
        return False
    for ii in range(0, len(seq1.seq)-max_submatch+1):
        sub_seq = seq1.seq[ii:ii+max_submatch]
        mi = seq2.seq.find(sub_seq)
        #print str(sub_seq) + " " + str(seq2.seq) + " " + str(ii) + "\n"
        if mi >= 0:
            return True
    return False


def get_align_score(seq1, seq2):
    # not EDNAFULL but close enough
    # match, mismatch, gap_open, gap_ext
    return pairwise2.align.localms(seq1, seq2, 5, -4, -10, -0.5, score_only=True)


def filter_net(input_fasta_file, output_fasta_file):
    paths = env.Paths.Instance()
    settings = env.Settings.Instance()


    pf_inter_cutoff = settings.pf_inter_cutoff
    sw_cutoff_for = settings.sw_cutoff_for    #35
    sw_cutoff_rev = settings.sw_cutoff_rev    #35
    max_submatch = settings.max_submatch     #10

    inputf = open(input_fasta_file, "r")
    outputf = open(output_fasta_file, "w")

    seq_dict = {}
    seq_list = []
    G = nx.Graph()
    G_inv = nx.Graph()

    print "reading sequences"
    for record in SeqIO.parse(inputf, "fasta"):
        seq_dict[record.id] = record
        seq_list.append(record)
        G.add_node(record.id)
        G_inv.add_node(record.id)

    print "iterating over sequences"
    for ii in range(0, len(seq_list)):
        seq_i = seq_list[ii]
        #print id_i
        #print seq_i
        var_seq_i_rec = seq_i[settings.spacer_pos_start:settings.spacer_pos_end]
        var_seq_i = str(var_seq_i_rec.seq)   #str(seq_i.seq[settings.spacer_pos_start:settings.spacer_pos_end])
        for jj in range(ii+1, len(seq_list)):
            seq_j = seq_list[jj]
            var_seq_j_rec = seq_j[settings.spacer_pos_start:settings.spacer_pos_end]
            var_seq_j = str(var_seq_j_rec.seq)  #str(seq_j.seq[settings.spacer_pos_start:settings.spacer_pos_end])
            submatched = has_submatch(var_seq_i_rec, var_seq_j_rec, max_submatch) or \
                           has_submatch(var_seq_i_rec, var_seq_j_rec.reverse_complement(), max_submatch)

            is_connected = submatched

            # only do Smith-Waterman calculation if not already connected
            if is_connected != True:
                sw_score_for = get_align_score(var_seq_i_rec.seq, var_seq_j_rec.seq)
                sw_score_rev = get_align_score(var_seq_i_rec.seq, var_seq_j_rec.reverse_complement().seq)
                #print str(sw_score_for) + " " + str(sw_score_rev)
                if sw_score_for > sw_cutoff_for or sw_score_rev > sw_cutoff_rev:
                    is_connected = True

            # only do expensive Pairfold calculation if not already connected
            if is_connected != True:
                pf_temp = 37
                pf_cmd = paths.pairfold_exe + " " + var_seq_i + " " + var_seq_j + " -t " + str(pf_temp) + " -m RNA | grep MFE | cut -f3"
                print pf_cmd
                mfe = subprocess.check_output(pf_cmd, shell=True, stderr=subprocess.STDOUT)
                mfe = float(mfe)
                #print(mfe)
                if mfe < pf_inter_cutoff:
                    is_connected = True

            if is_connected:
                G.add_edge(seq_i.id, seq_j.id)
            else:
                G_inv.add_edge(seq_i.id, seq_j.id)

    print "Read " + str(G.size()) + " edges"
    print "Read " + str(G_inv.size()) + " non-edges"

    cl = list(nx_alg.find_cliques(G_inv))

    best_index = -1
    best_size = -1

    #print cl
    for ii in range(len(cl)):
        if len(cl[ii]) > best_size:
            best_size = len(cl[ii])
            best_index = ii

    print "maximum clique is of size: " + str(best_size)
    print cl[best_index]

    for lab in (cl[best_index]):
        record = seq_dict[lab]
        # print record
        # output as FASTA
        SeqIO.write(record, outputf, "fasta")
