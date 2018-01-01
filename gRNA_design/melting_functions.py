# (c) Michael Crone and James T. MacDonald, 2017
import env
import subprocess
from Bio import SeqIO
import sys
import math

#inputfile = open('sequences.txt',"rU")
#outputfile = open('output.txt',"w")

def dnameltoutput(sequence):
    paths = env.Paths.Instance()
    settings = env.Settings.Instance()

    kconc = 'K=' + str(settings.melting_kconc)       #0.200'
    mgconc = 'Mg=' + str(settings.melting_mgconc)    #0.003'
    naconc = 'Na=' + str(settings.melting_naconc)    #0.008'

    rnaseq=sequence.transcribe()
    dnaseq=rnaseq.back_transcribe()
    dnaseq=dnaseq.complement()
    meltinput = paths.melting_exe + ' -S ' + str(rnaseq) + ' -C ' + str(dnaseq) + ' -H rnadna -P 0.001 -E ' \
                + kconc + ' ' + mgconc + ' ' + naconc
    print(meltinput)
    output = subprocess.check_output(meltinput, shell=True)
    #print output
    enthpos = output.find('cal/mol ( ',)
    enthalpy=float(output[enthpos+10:enthpos+18].replace(',',''))
    entrpos = output.find('cal/mol-K ( ',)
    entropy=round(float(output[entrpos+12:entrpos+18].replace(',','')),1)
    tmpos = output.find('ture : ',)
    tm=round(float(output[tmpos+7:tmpos+8].replace(',','.')),1)
    deltag = round(enthalpy - 310.15*entropy,1) # 37 C = 310.15 K
    #Add Kd
    rt = 8.3144598*310.15
    kd = math.exp(deltag/(rt))*1
    return deltag, kd

def filter(input_fasta_file, output_fasta_file):
    paths = env.Paths.Instance()
    settings = env.Settings.Instance()

    inputf = open(input_fasta_file, "r")
    outputf = open(output_fasta_file, "w")

    num_seqs = 0
    for record in SeqIO.parse(inputf, "fasta"):
        sequence = record.seq[settings.spacer_pos_start:settings.spacer_pos_end]     #[21:44]
        deltag, kd = dnameltoutput(sequence)
        deltagrev, kdrev = dnameltoutput(sequence.reverse_complement())
        #print(record.id + '\t' + str(deltag) + '\t' + str(kd) + '\t' + str(deltagrev) + '\t'+ str(kdrev))

        if deltag > settings.melting_GFE_fwd_min:
            if deltag < settings.melting_GFE_fwd_max:
                if deltagrev > settings.melting_GFE_rev_min:
                    if deltagrev < settings.melting_GFE_rev_max:
                        #print()
                        num_seqs += 1
                        outputf.write(">" + str(record.id) + "\n")
                        outputf.write( str(record.seq) + "\n")

        #outputf.write(record.id + '\t' + str(deltag) + '\t' + str(kd) + '\t' + str(deltagrev) + '\t'+ str(kdrev) + "\n")
    inputf.close()
    outputf.close()

    print "INFO: Output: " + str(num_seqs) + " filtered sequences"



#dnarnamelt(inputfile,outputfile)
