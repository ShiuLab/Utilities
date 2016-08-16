"""
PURPOSE:
Use biopython's mapping function to map motifs to a DNA sequence. Improves on previous lab motif mapping scripts by setting a 
maximum precision desired when calculating the score distribution, the computational cost of doing this without a limit 
grows exponentially with the length of the motif. 


INPUT:
  -m         PWM file (ending in .pwm) or directory with .pwm files
  -fasta     Fasta file with ONE sequence only.
  -AT        ATbias in organism (Arabidopsis = 0.33)
  -GC        GCbias in organism (Arabidopsis = 0.17)
  -p         Threshold value for false positives (Default = 1e-5)
  -out       Name for output file

OUTPUT:
  -out    Results: Sequence - Motif - Motif_consensus - Hit_Seq - Hit_position - Hit_score -Threshold

AUTHOR: Christina Azodi

REVISIONS:   Submitted 8/10/2016

"""

import os,sys, math
from Bio import motifs, SeqIO
from Bio.Seq import Seq
import numpy as np
import time

def make_pssm(m, ATbias, GCbias, p):
  """Convert pwm into pssm (position-specific scoring matrices) to use for mapping 
  and calulate thresholds based on lenght and information content of motif"""
  
  #Add pseudocount so that the probablity never becomes zero which will result in pssm = -inf.
  pwm = m.counts.normalize(pseudocounts = 0.001)  
  consensus = pwm.consensus
  print("Motif consensus sequence: " + str(consensus))

  # Make position-specific scoring matrices & generage stats
  background = {'A': float(ATbias), 'C': float(GCbias), 'G': float(GCbias), 'T': float(ATbias)}
  pssm = pwm.log_odds(background)
  print("Maximum score obtainable from motif: %4.2f" % pssm.max)
  print("Minimum score obtainable from motif: %4.2f" % pssm.min) 
  print("Mean pssm score (a measure of information content of motif compared to background): %s" % str(pssm.mean(background)))

  # Determine the threshold for mapping 
  # Note: Since the space for a score distribution grows exponentially with motif length, 
        # this uses an approximation with a given precision to keep computation cost manageable
  distribution = pssm.distribution(background=background, precision=10**4)

  # Set the threshold to p
  threshold = distribution.threshold_fpr(p)
  print("Threshold for matching (fpr=0.00001): %3f" % threshold)

  # Alternative threshold method: balance FN with FP (100 false negatives for every 1 false positive)
  #threshold = distribution.threshold_balanced(100) 
  #print("Threshold for matching (fnr/fpr=100): %3f" % threshold)
  
  return(pssm, consensus, threshold)

def mapping(seq_file, pssm, threshold, consensus, motif_name, motif_len):
  """Map the pssm to the sequence in fasta format"""
  seq_name = seq_file.strip().split("/")[-1]
  for record in SeqIO.parse(seq_file, "fasta"):
    seq = Seq(str(record.seq), m.alphabet)
    local_array = np.array([])
    for position, score in pssm.search(seq, threshold=threshold):
      real_seq = record.seq[position:position+motif_len]
      out.write("%s\t%s\t%s\t%s\t%s\t%4.5f\t%4.5f\n" % (seq_name, motif_name, consensus, real_seq, str(position), score, threshold))

  #return(x)

if __name__=="__main__":
  start_time = time.time()

  # Set default values
  ATbias = 0.33
  GCbias = 0.17
  p = 0.00001
  include = "all"

  #Define arguments
  for i in range (1,len(sys.argv),2):
    if sys.argv[i] == '-m':             # pwm sequence file (TAMO)
      f = sys.argv[i+1]
    if sys.argv[i] == '-AT':            # ATbias in organism (Arabidopsis = 0.33)
      ATbias = sys.argv[i+1]
    if sys.argv[i] == '-GC':            # GCbias in organism (Arabidopsis = 0.17)
      GCbias = sys.argv[i+1]
    if sys.argv[i] == '-fasta':         # FASTA file of the sequence
      seq_file = sys.argv[i+1]
    if sys.argv[i] == '-p':             # threshold (default = 1e-5)
      p = sys.argv[i+1]
    if sys.argv[i] == '-out':           # Results output name
      out_file = sys.argv[i+1] 
    if sys.argv[i] == '-include':       # List of file names from pwm directory to include (default = all)
      include = sys.argv[i+1]   
  
  if len(sys.argv) <= 1:
    print(__doc__)
    exit()

  #threshold= math.pow(10,-float(sys.argv[2])) # P-val, 5 means -log(n,10)
  out = open(out_file, 'w')
  out.write("Sequence\tMotif\tMotif_consensus\tHit_Seq\tHit_position\tHit_score\tThreshold\n")

  # Specify if you only want to include certain pwm's in the directory
  if include == "all":
    include_if = os.listdir(f)
  else:
    include_if = open(include, 'r').read().splitlines()


  if f[-4:] == ".pwm":
    #Make PSSM for mapping
    motif_name = f.strip().split("/")[-1]
    m = motifs.read(open(f), "pfm")
    motif_len = len(m)
    pssm, consensus, threshold = make_pssm(m, ATbias, GCbias, p)
    hits = mapping(seq_file, pssm, threshold, consensus, motif_name, motif_len)

  else:
    for motif_file in os.listdir(f):
      if motif_file in include_if:
        motif_name = motif_file.strip().split("/")[-1]
        print(motif_name)
        path = f + motif_file
        m = motifs.read(open(path), "pfm")
        motif_len = len(m)
        pssm, consensus, threshold = make_pssm(m, ATbias, GCbias, p)
        hits = mapping(seq_file, pssm, threshold, consensus, motif_name, motif_len)
      else:
        pass
  end_time = time.time()
  print("Run time: " + str(float(end_time)-float(start_time)))

