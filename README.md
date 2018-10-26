# Utilities
Scripts useful for manipulation and basic analysis of common file types

## Scripts
  -FastaManager.py: get_sequences, getseq2, get_group_seq, get_stretch2, get_gc, fasta_to_oneline, fasta_to_phylip, fasta_to_stockholm, oneline_to_fasta, parse_desc, simplify, prefix, size_filter, compare_lists, cleanup, get_sizes, divide, divide1seq, index_names, index_pairs, change_names, rename, rename_all, del_redun_names, del_redun_seq, check_redun, get_sp, parse_ensembl_fasta, convert_header, delete, concat, locate, indiv, mask, count, get_longest, format
  
  -GFFUtil.py: Includes sort, merge, sort_merge, merge_depth, length, lengths, min_len, len_prcntl, compare_lens, prefix_id, loc_id, loc_id_keep, coords2gff, mask, overlap, and overlap+.

  -Map_PWM_biopy.py: Map PWM or dictionary of PWMs to a sequence in a fasta file.
  
## Running R on hpcc
  1. Load the following modules:
      
      module load GCC/7.3.0-2.30
      
      module load OpenMPI/3.1.1
      
      module load R/3.5.1-X11-20180604
      
  2. Run the R script:
  
         R --vanilla --slave --args <arg 1> <arg 2> > script.R
