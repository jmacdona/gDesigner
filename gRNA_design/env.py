# (c) Michael Crone and James T. MacDonald, 2017
import utils
import sys
import numpy as np

@utils.Singleton
class Paths:
    def __init__(self):
        self.grna_design_install_dir = ""
        self.offinder_exe = ""
        self.pairfold_exe = ""
        self.melting_exe = ""
        print 'Paths singleton object created'

    def read_file(self, input_filename):
        print("Readings paths...")
        inp_f = open(input_filename,"r")
        for line in inp_f:
            cols = line.rstrip().split('\t')
            if len(cols) == 2:
                var = cols[0]
                val = cols[1]
                if (var == "grna_design_install_dir"):
                    self.grna_design_install_dir = val
                elif (var == "offinder_exe"):
                    self.offinder_exe = val
                elif (var == "pairfold_exe"):
                    self.pairfold_exe = val
                elif (var == "melting_exe"):
                    self.melting_exe = val
                else:
                    sys.exit('Error: unrecognised setting option: ' + var)
            else:
                print "INFO: ignoring line: " + line + "\n"


@utils.Singleton
class Settings:
    def __init__(self):
        self.output_root_dir = ""
        self.crispr_type = ""     # lbcpf1 or ascpf1 or cas9
        self.pam_seq = ""     # TTTN or
        self.pam_five_prime = 0   # 0 (PAM on 5') or 1 (PAM on 3')
        self.spacer_pos_start = 0   # start position of spacer sequence in input sequences
        self.spacer_pos_end = 20    # end position of spacer sequence in input sequences (Exclusive - follow python)
        self.reformat_spacer_pos_start = 0  # start position of spacer sequence in offinder formatted sequences
        self.reformat_spacer_pos_end = 20   # end position of spacer sequence in offinder formatted sequences
        self.offinder_max_mismatch = 5
        self.offinder_compute_type = "C"    #(C: using CPUs, G: using GPUs, A: using accelerators)
        self.min_offinder_score_fwd = 0.8
        self.min_offinder_score_rev = 0.8
        self.offinder_score_weights = np.array([])

        #Physiological Concentrations of electrolytes for MELTING calcs
        self.melting_kconc = 0.200      # 'K=0.200'
        self.melting_mgconc = 0.003     # 'Mg=0.003'
        self.melting_naconc = 0.008     # 'Na=0.008'

        # ['GibbsFwd'])>-130000:
        self.melting_GFE_fwd_min = -130000
        # ['GibbsFwd'])<-120000:
        self.melting_GFE_fwd_max = -120000
        # ['GibbsRev'])>-130000:
        self.melting_GFE_rev_min = -130000
        # ['GibbsRev'])<-120000:
        self.melting_GFE_rev_max = -120000

        self.max_submatch = 10

        self.pf_inter_cutoff = -12

        self.sw_cutoff_for = 35
        self.sw_cutoff_rev = 35

        # scoring based on FnCpf1 in Fonfara et al
        self.ascpf1_score_weights = np.array([
            0.63, 0.59, 0.78, 0.70, 0.67, 0.72, 0.61, 0.43, 0.31, 0.23,
            0.13, 0.32, 0.21, 0.00, 0.35, 0.21, 0.08, 0.25, 0.45, 0.38,
            0.34, 0.36, 0.00])
        # scoring based on FnCpf1 in Fonfara et al
        self.lbcpf1_score_weights = np.array([
            0.63, 0.59, 0.78, 0.70, 0.67, 0.72, 0.61, 0.43, 0.31, 0.23,
            0.13, 0.32, 0.21, 0.00, 0.35, 0.21, 0.08, 0.25, 0.45, 0.38,
            0.34, 0.36, 0.00])
        self.cas9_score_weights = np.array([
            0.000, 0.000, 0.014, 0.000, 0.000, 0.395, 0.317, 0.000, 0.389, 0.079,
            0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583])
        print 'Settings singleton object created'

    def read_file(self, input_filename):
        print("Reading settings...")
        inp_f = open(input_filename,"r")
        for line in inp_f:
            cols = line.rstrip().split('\t')
            if len(cols) == 2:
                var = cols[0]
                val = cols[1]
                if (var == "crispr_type"):
                    self.crispr_type = val
                    if (val != "lbcpf1" and val != "ascpf1" and val != "cas9" and val != "other"):
                        sys.exit("ERROR: unknown crispr_type: '"+ val + "'\n")
                    elif (val == "other"):
                        print("INFO: user defined other crispr_type")
                    elif (val == "lbcpf1"):
                        self.pam_seq = "TTTN"
                        self.pam_five_prime = 1
                        self.spacer_pos_start = 21
                        self.spacer_pos_end = 44
                        self.offinder_max_mismatch = 7
                        self.offinder_score_weights = self.lbcpf1_score_weights
                    elif (val == "ascpf1"):
                        self.pam_seq = "TTTN"
                        self.pam_five_prime = 1
                        self.spacer_pos_start = 20
                        self.spacer_pos_end = 43
                        self.offinder_max_mismatch = 7
                        self.offinder_score_weights = self.ascpf1_score_weights
                    elif (val == "cas9"):
                        self.pam_seq = "NRG"
                        self.pam_five_prime = 0
                        self.spacer_pos_start = 0
                        self.spacer_pos_end = 20
                        self.offinder_max_mismatch = 5
                        self.offinder_score_weights = self.cas9_score_weights
                elif (var == "offinder_compute_type"):
                    if val == "C" or val == "G" or val == "A":
                        self.offinder_compute_type = val
                    else:
                        sys.exit("ERROR: unknown offinder_compute_type: '"+ val + "'\n")
                elif (var == "min_offinder_score_fwd"):
                    self.min_offinder_score_fwd = float(val)
                elif (var == "min_offinder_score_rev"):
                    self.min_offinder_score_rev = float(val)
                elif (var == "offinder_max_mismatch"):
                    self.offinder_max_mismatch = int(val)
                elif (var == "melting_kconc"):
                    self.melting_kconc = float(val)
                elif (var == "melting_mgconc"):
                    self.melting_mgconc = float(val)
                elif (var == "melting_naconc"):
                    self.melting_naconc = float(val)
                elif (var == "melting_GFE_fwd_min"):
                    self.melting_GFE_fwd_min = float(val)
                elif (var == "melting_GFE_fwd_max"):
                    self.melting_GFE_fwd_max = float(val)
                elif (var == "melting_GFE_rev_min"):
                    self.melting_GFE_rev_min = float(val)
                elif (var == "melting_GFE_rev_max"):
                    self.melting_GFE_rev_max = float(val)
                elif (var == "max_submatch"):
                    self.max_submatch = int(val)
                elif (var == "pf_inter_cutoff"):
                    self.pf_inter_cutoff = float(val)
                elif (var == "sw_cutoff_for"):
                    self.sw_cutoff_for = float(val)
                elif (var == "sw_cutoff_rev"):
                    self.sw_cutoff_rev = float(val)
                else:
                    sys.exit('Error: unrecognised setting option: ' + var)

            else:
                print "INFO: ignoring line: " + line + "\n"



