# DEV NOTES ##
# General
# -These functions are derived from the overlap functions in 
# /mnt/home/lloydjo1/scripts/gff_manager.py
#
# -General imports are used in more than one function
#
# Version History
# - (6/10/2016) Initial version development
# - (6/13/2016) Add overlap by bisection for target features a single position

#######################
### GENERAL IMPORTS ###
#######################

import sys
import math

#################
### FUNCTIONS ###
#################

## 0. General Functions ##
# Functions use in more than one pipeline keep. Should include only very
# basic processes not already in python

def FilterList(in_list,indexes,keep=True):
    ''' Takes a in_list and filteres by a list of indexes

    Input:
         -in_list = the unmodified input list
         -indexes = the list of indexes to fliter the list by
         -keep = flag that control whether to keep (True) or exclude (False) the indexes
 
    Output:
         -out_list = the filtered list 
      
    '''

    if keep == True:
        out_list = [in_list[i] for i in range(len(in_list)) if i in indexes] 
        return out_list
    elif keep == False:
        out_list = [in_list[i] for i in range(len(in_list)) if i not in indexes]
        return out_list
    else:
        print "Warning: Unrecogized value " + str(keep) + " for flag keep. Returning the original list"
        return in_list

## 1. Indexed Overlapping Functions ##
# Functions for finding overlaps between two sets of features with length (i.e.
# two list of genes) using an indexed dictioanry

def FindSeqMax(lines,seq_idx,stop_idx,seq_max_dict,sep="\t"):
    ''' For a set features, for each sequence (chromosome,scaffold,contig, etc)
        find the largest index for each sequence
        
        Input:
             -lines = The lines for the feature file
             -seq_idx = The index of the sequence column in the feature file
             -stop_idx = The index of the stop coordinate column in the feature file 
             -max_dict = The dictionary which we will record the max index in

        Optional:
             -sep = Seperator character for file parsing
        
        Output:
             -None, dictionary is updated in the function

    '''
    # Get a list of sequences and stop-coords by sequence    
    SeqCoords = [[ln.split(sep)[seq_idx],int(ln.split(sep)[stop_idx])] for ln in lines ]
    Seqs = set([v[0] for v in SeqCoords])

    # Update dictionary
    for s in Seqs:
        current_max = max([v[1] for v in SeqCoords if v[0] == s])
        # If the sequence is not new, check if update is needed
        if s in seq_max_dict.keys():
            previous_max = seq_max_dict[s]
            if current_max > previous_max:
                seq_max_dict[s] = current_max
        # If the sequence is new, add it to the dict
        else:
            seq_max_dict[s] = current_max

def IndexSequences(seq_max_dict,span):
    ''' For a set of sequences and their lengths, build and indexed
    dictionary of positions within that dictionary 

    Inputs:
          -seq_max_dict = A Dictionary of sequences and their lengths (i.e. max vlaues)
          -span = the size of indexes in the dictionary

    Ouputs:
          -seq_index_dict = A dictionary of indexed sequences 
 
    Note:
          -IMPORTANT: The start of an index will overlap with the start of the previous index
           (i.e. [0_10000][10000_20000]). Therefore, the inculsion test should be start <= position < end
    '''
    # Define out output dictionary
    seq_index_dict = {}

    # For each sequence
    for seq in seq_max_dict.keys():
        seq_index_dict[seq] = {}
        max = seq_max_dict[seq]
        max_index = math.floor(max/span) # Will return the divisor, rounded down to nearest int
        max_index = int(max_index + 1)   # Add one to ensure last index goes beyond max
        
        # Define indexes
        for index_start in range(0,max_index*span,span):
            index_end = index_start + span
            seq_index_dict[seq][str(index_start) + "_" + str(index_end)] = set()
    return seq_index_dict 

def GetIndexes(start,stop,span):
    ''' Get the set indexes covere by a feature
        
        Input:
             -start = start postion of the feature
             -stop = stop position of the feature
             -span = width of indexes

        Ouput:
             index_list = list of indxes in start_stop format
    
    '''
    # Convert start/stop to int
    start = int(start)
    stop = int(stop)

    # Get first span before and beyond the feature
    ind_start = int(math.floor(start/span)*span)
    ind_stop = int((math.floor(stop/span)+1)*span)

    # Define index list as all indexes betw
    index_list = [str(i) + "_" + str(i+span) for i in range(ind_start,ind_stop,span)]

    return index_list

def AddToIndex(file_vars,seq_index_dict,span,sep="\t"):
    ''' Add the features from a file to an indexed postional dictionary

        Input:
             -file_vars = A list defining the following features of the feature file:
                    [lines, sequence_col#, start_col#, stop_col#] 
             -seq_index_dict = A dictionary of indexed sequences
             -span = the size of an indexed region in the overlap dict

        Optional:
             -sep = Seperator character for file parsing

        Output:
             -None, dictionary is upated in the script

        Subfunctions:
             -GetIndex(start,stop,span)
    '''

    # Get variables for our feature file
    [lines,seq_idx,start_idx,stop_idx] = file_vars

    # Make sure the indexes are ints
    seq_idx = int(seq_idx)
    start_idx = int(start_idx)
    stop_idx = int(stop_idx)
    
    # For each line
    for ln in lines:
        # Split the ln
        split_ln = ln.split(sep)

        # Get the values we need
        seq = split_ln[seq_idx]
        start = split_ln[start_idx]
        stop = split_ln[stop_idx]
 
        # Get everything else
        remainder = ",".join([split_ln[i] for i in range(len(split_ln)) if i not in [seq_idx,start_idx,stop_idx]])
        
        # Get indexes of the feature
        dict_indexes = GetIndexes(int(start),int(stop),span)

        # Add feature to each index it overlaps
        for idx in dict_indexes:
            try:
                seq_index_dict[seq][idx].add((start,stop,remainder))
            except KeyError:
                print "Warning: Index " + idx + " for feature at " + seq + ": " + start + " - " + stop + " is not in the indexed sequence dictionary"

def TestOverlap(start1,stop1,start2,stop2):
    ''' Test if two features overlap based on their position

        Input:
             -start1 = start of feature 1
             -stop1 = stop of feature 1
             -start2 = start of feature 2
             -stop2 = stop of feature 2
       
        Ouptut:
             -bool = True if overlap, false if not

    '''
        
    # Make sure we have integers
    start1 = int(start1)
    stop1 = int(stop1)
    start2 = int(start2)
    stop2 = int(stop2)

    # Determine if there is an overlap
    overlap = False
    if start2 >= start1 and start2 <= stop1:
        overlap = True
    elif start1 >= start2 and start1 <= stop2:
        overlap = True
    return overlap

def FindOverlapsByIndexing(ref,target,span,sep="\t"):
    ''' Find overlaps between features in a reference file
        and a target file using an indexed dictionary

        Input:
             -ref = A list defining the following features of the reference file:
                    [lines, sequence_col#, start_col#, stop_col#]
             -target = A list definining the following features of the target file:
                    [lines, sequence_col#,start_col#, stop_col#]
             -span = the size of an indexed region in the overlap dict
        
        Optional:
             -sep = Seperator character for file parsing

        Output:
             -overlap_dict = a dictionary of overlaps where:
                  > the dictionary keys define sequence, start position, and stop position of the reference features
                  > the dictionary values are a list where:
                       -The first value is the remainder of the reference feature
                       -Subsequent values are tuples containing the start, stop, and remainder of target features
                       -If there are no overlaps, there will be one tuple of (NA,NA,NA)
                       -The last value is the number of overlaps
                  > remainder is a comma seperate list of all other values for a feature (i.e.e not seq, start, stop)
             -sorted_keys = a list of keys of overlap_dict ordered by according to their appearance in the ref file

        Subfunctions:
             -FindSeqMax(lines,seq_idx,stop_idx,seq_max_dict)
             -IndexSequenecs(seq_max_dict,span)
             -AddToIndex(file_vars,seq_index_dict,span)
             -TestOverlap(start1,stop1,start2,stop2)

        Notes:
             -IMPORTANT: The function EXPLICITLY treat each line in the reference 
             and target file as a seperate features because multiple features
             may have the same coordiantes (i.e. reads)
             -IMPORTANT: The start and stop coordinates MUST be ordered beforehand.
             This script assumes that the start coordindate ALWAYS comes before 
             the stop coordinate on the defined sequenec
             -IMPORTANT: The lines for each file should be striped of the "\n" character before
             being passed through this function
    '''

    # Breakdown our variable list, cause we don't need all of them all the time
    ref_lines = ref[0]
    tar_lines = target[0]    
    [ref_seq,ref_start,ref_stop] = [int(i) for i in ref[1:]]
    [tar_seq,tar_start,tar_stop] = [int(i) for i in target[1:]]

    # Determine the size of each sequence
    seq_max_dict = {}
    FindSeqMax(ref_lines,ref_seq,ref_stop,seq_max_dict,sep)
    FindSeqMax(tar_lines,tar_seq,tar_stop,seq_max_dict,sep)

    # Max an indexed dictionary 
    seq_index_dict = IndexSequences(seq_max_dict,span)

    # Add the target featurs to the 
    AddToIndex(target,seq_index_dict,span,sep)

    # For each referrence feature...
    sorted_keys = []
    overlap_dict = {}
    for ln_idx in range(len(ref_lines)): # IMPORTANT, we use index here because there may be multiple featurs with the same location
        # Split ln
        split_ln = ref_lines[ln_idx].split(sep)
 
        # Get the values we need
        seq = split_ln[ref_seq]
        start = split_ln[ref_start]
        stop = split_ln[ref_stop]

        # Get everything else
        remainder = ",".join([split_ln[i] for i in range(len(split_ln)) if i not in [ref_seq,ref_start,ref_stop]])

        # Define entry in overlap dictionary and add to sorted_keys
        overlap_dict[(seq,start,stop,ln_idx)] = [remainder]
        sorted_keys.append((seq,start,stop,ln_idx))

        # Get indexes of the feature
        dict_indexes = GetIndexes(start,stop,span)

        # Get features in overlapping indecies
        comparison_features = []
        for index in dict_indexes:
            try:
                comparison_features.extend(list(seq_index_dict[seq][index]))
            except KeyError:
                print "Warning: Index " + index + " for feature at " + seq + ": " + start + " - " + stop + " is not in the indexed sequence dictionary"
                  
        comparison_features = list(set(comparison_features)) # De-dup the comparison list

        # Test for overlap
        count_overlaps = 0
        for comp in comparison_features:
            comp_start = comp[0]
            comp_stop = comp[1]
            overlap = TestOverlap(start,stop,comp_start,comp_stop)
            # If there is an overlap, add it to our dictionary
            if overlap == True:
                count_overlaps = count_overlaps + 1
                overlap_dict[(seq,start,stop,ln_idx)].append(comp)

        # If the there were not overlaps, create a ("NA","NA","NA") field
        if count_overlaps == 0:
            overlap_dict[(seq,start,stop,ln_idx)].append(("NA","NA","NA"))
            overlap_dict[(seq,start,stop,ln_idx)].append(count_overlaps)
        else:
            overlap_dict[(seq,start,stop,ln_idx)].append(count_overlaps)
    
    return overlap_dict, sorted_keys

## 2. Bisect Overlap Functions ##
# Functions for find the set of single point features (i.e. SNPs, methylated sites)
# with a feature with length (i.e. a gene)

def UpdateSeqDict(lines,seq_idx,seq_dict,sep="\t"):
    '''Add sequences in a feature file to the sequence dictionary

    Input:
          -lines = The lines from the feature file 
          -seq_idx = The index of sequence column in the feature file
          -seq_dict = The dictionary to be updated

    Optional:
         -sep = Seperator character for file parsing

    Ouput:
          -None, dictionary is updated in the function

    '''

    # Get a list of sequenecs
    Seqs = set([ln.split(sep)[seq_idx] for ln in lines])
    
    # Update dictionary
    for s in Seqs:
        if not s in seq_dict:
            seq_dict[s] = [[],[]]

def SortPairedLists(listA,listB):
    ''' Sort a pair of lists by the values in the first list

    Input:
         listA = first list
         listB = second list
    
    Ouput:
         sort_listA = first list sorted by first list
         sort_listB = second list sorted by first list

    '''
    sort_listA = [A for (A,B) in sorted(zip(listA,listB), key=lambda pair: pair[0])]
    sort_listB = [B for (A,B) in sorted(zip(listA,listB), key=lambda pair: pair[0])]
    return [sort_listA,sort_listB]


def AddtoSeq(file_vars,seq_dict,sep="\t"):
    ''' Add single postion features from a file to sequence dictionary

    Input:
         -file_vars = A list definfing the following features of the features file:
                      [lines, sequence_col#, position_col#]
         -seq_dict = A dictionary of sequence to be filled

    Optional:
         -sep = Seperator chracter for file parsing

    Output:
          -None, dictionary is updated in the script

    Subfunctions:
          -FilterList(in_list,indexes,keep=True)
          -SortPairedLists(list1,list2)
    '''

    # Get variables for out featulre file
    [lines,seq_idx,posit_idx] = file_vars
    
    # Make sure indexes are ints
    seq_idx = int(seq_idx)
    posit_idx = int(posit_idx)
    
    # For each sequence in the sequence dictionary
    for seq in seq_dict.keys():
        # Get indexes for the current seq
        
        # Get features on the current seq
        seq_indexes = [i for i in range(len(lines)) if lines[i].split("\t")[seq_idx] == seq]

        # Get position and remainder
        position_list = [int(lines[i].split("\t")[posit_idx]) for i in seq_indexes] #Need intergers for sorting and comparison
        remainder_list = [",".join(FilterList(lines[i].split(sep),[seq_idx,posit_idx],keep=False)) for i in seq_indexes]
        
        # Sort lists and add to dict
        [sort_position_list, sort_remainder_list] = SortPairedLists(position_list,remainder_list)
        seq_dict[seq][0] = sort_position_list
        seq_dict[seq][1] = sort_remainder_list

def BisectFeature(start,stop,lists):
    '''Find the position defineded by list overlapped by the region defined by start top stop

    Input:
         -start = the start loci of the region
         -stop = the end loci of the region
         -list = the position [0] and remainder [1] values of the features 

    Ouput: 
         -filt_list = the lists in list filtered for overlap in the region

    Note:
         -This function requires the bisect library
    '''
    import bisect
    # Find where the ergion start/stop fit in the sorted position list
    start_index = bisect.bisect_left(lists[0],start)
    stop_index = bisect.bisect_right(lists[0],stop)

    # Alter return based on the poistion of the stop index
    if stop_index > len(lists[0]) - 1: # If the stop index is beyond the list
        return [lists[0][start_index:],lists[1][start_index:]]
    else: 
        return [lists[0][start_index:stop_index],lists[1][start_index:stop_index]]        

def FindOverlapByBisection(ref,target,sep="\t"):
    '''Find overlaps between a set of features with length (in ref) and a 
    set of single point features (in target)

    Input:
          -ref = A list defining the following features of the reference file:
                    [lines, sequence_col#, start_col#, stop_col#]
          -target = A list definining the following features of the target file:
                    [lines, sequence_col#, position_col#]

    Optional:
         -sep = Seperator character for file parsing

    Output:
         -overlap_dict = a dictionary of overlaps where:
              > the dictionary keys define sequence, start position, and stop position of the reference features
              > the dictionary values are a list where:
                   -The first value is the remainder of the reference feature
                   -Subsequent values are tuples containing the position and remainder of target features
                   -If there are no overlaps, there will be one tuple of (NA,NA)
                   -The last value is the number of overlaps
              > remainder is a comma seperate list of all other values for a feature (i.e.e not seq, start, stop)
         -sorted_keys = a list of keys of overlap_dict ordered by according to their appearance in the ref file
    
    Subfunctions:
         -UpdateSeqDict(lines,seq_idx,seq_dictsep="\t")
         -AddtoSeq(file_vars,seq_dict,sep="\t")
         -BisectFeature(start,stop,lists)

    Notes:
         -See Notes for "FindOverlapsByIndexing"

    '''

    # Breakdown variable list, not all values will be needed all the time
    ref_lines = ref[0]
    tar_lines = target[0]
    [ref_seq,ref_start,ref_stop] = [int(i) for i in ref[1:]]
    [tar_seq,tar_posit] = [int(i) for i in target[1:]]
 
    # Make a dictionary of sequences
    seq_dict = {}
    UpdateSeqDict(ref_lines,ref_seq,seq_dict,sep)
    UpdateSeqDict(tar_lines,tar_seq,seq_dict,sep)
    
    # Add target features to the sequene dictionary
    AddtoSeq(target,seq_dict,sep)

    # For each referrence feature
    sorted_keys = []
    overlap_dict = {}
    for ln_idx in range(len(ref_lines)): # IMPORTANT, we use index here because there may be multiple features with the same location
         
        # Split ln
        split_ln = ref_lines[ln_idx].split(sep)

        # Get the values we need
        seq = split_ln[ref_seq]
        start = split_ln[ref_start]
        stop =  split_ln[ref_stop]
       
        # Get everyhing else
        remainder = ",".join([split_ln[i] for i in range(len(split_ln)) if i not in [ref_seq,ref_start,ref_stop]])

        # Define entry in overlap dictionary and add to sorted_keys
        overlap_dict[(seq,start,stop,ln_idx)] = [remainder]
        sorted_keys.append((seq,start,stop,ln_idx))

        # Get the seq of the feature
        seq_features = seq_dict[seq]

        # Get the overlaping position
        overlapping_features = BisectFeature(int(start),int(stop),seq_features)
        overlap_number = len(overlapping_features[0])

        # Update overlap_dict
        if overlap_number == 0:
            overlap_dict[(seq,start,stop,ln_idx)].append(("NA","NA"))
        else:
            overlap_dict[(seq,start,stop,ln_idx)].extend([(str(pair[0]),pair[1]) for pair in zip(overlapping_features[0],overlapping_features[1])])
        overlap_dict[(seq,start,stop,ln_idx)].append(overlap_number)
 
    return overlap_dict, sorted_keys

## 3. Generic Utility Functions ##
# Generic utility functions for reading and writing files related to overlap

def ReadCtlFileForOverlap(file):
    '''Reads a control file for running finding overlaps

       Input:
            -file = Path the control file

       Ouput:
            -argument_dict = a dictionary of input values

    '''

    # Read control file and parse values
    argument_dict = {}
    ctl_lines = [ln.strip() for ln in open(file,"r").readlines()]
    for ln in ctl_lines:
        [key,value] = ln.split("\t")
        argument_dict[key] = value
        
    # Parse important values
    try:
        ref_values = argument_dict["Ref"]
        ref_values = ref_values.split(",")
        argument_dict["Ref"] = ref_values
    except Exception as e:
        print "Error: Comma seperated list of values is required for the Ref argument"
        raise e
    try:
        tar_values = argument_dict["Tar"]
        tar_values = tar_values.split(",")
        argument_dict["Tar"] = tar_values
    except Exception as e:
        print "Error: Comma seperated list of values is required for the Target argument"
        raise e
    try:
        span = argument_dict["Span"]
        argument_dict["Span"] = int(span)
    except Exception as e:
        print "Error: An interger value is require for the Span argument"
        raise e        

    return argument_dict

def FilterLines(lines,keep_idxs,sep="\t"):
    ''' Takes a list of lines and filters to only keep certains columns

    Input:
          -lines = a list of text lines to be filtered
          -keep_idxs = a list of indexes of the columns to keep after filtering

    Optional:
          -sep = the character that seperates columns in lines

    Output:
          -filtered_lines = a list of lines filtered for columns defined by keep_idxs

    Subfunctions:
          - FilterList(in_list,indexes,keep=True)

    Notes:
          -IMPORTANT: The text lines in "lines" should be stripped of the "\n" character before being
          passed to this function. It assumes no "\n" is present and will not at one if the
          final field is dropped
          -IMPORTANT: Recall that python indexes from 0, so keep_idxs should be adjusted BEFORE passing
          to this function 
    '''

    split_lines = [ln.split(sep) for ln in lines]
    filtered_lines = [sep.join(FilterList(split_ln,keep_idxs)) for split_ln in split_lines]
    return filtered_lines

def WriteOverlapLines(overlap_dict,sorted_keys,outfile):
    ''' Write overlaps from a overlap dictionary to an output file

        Input:
              -overlap_dict = an overlap dictionary generated by "FindOverlapsByIndexing"
              -sorted_keys = a list of keys generated by "FindOverlapsByIndexing"
              -outfile = the name of the outfile

        Output:
              -Writes a file with the following format:
                   > Each line records: Sequence RefStart RefStop RefRemainder TarStart Tar StopRemainder
                   > The same reference sequence may appear on multiple lines
                   > The number of overlaps per referrence is not recorded here (see "WriteOverlapNumber")
    '''
    
    # Make a list of lines to write  
    outlines = []  
    for key in sorted_keys:
        # Define a baseline for each ref feature

        remainder = overlap_dict[key][0]
        base_ln = "\t".join(list(key[0:3])) + "\t" + remainder

        # Write a line for each overlap
        overlaps = overlap_dict[key][1:-1]
        for ov in overlaps:
            outlines.append(base_ln + "\t" + "\t".join(list(ov)) + "\n")

    #Write the lines to file
    output = open(outfile,"w")
    output.write("".join(outlines))
    output.close()
    
def WriteOverlapNumber(overlap_dict,sorted_keys,outfile):
    ''' Write number of overlap for each reference features to an output file
        
        Inputs:
              -overlap_dict = an overlap dictionary generated by "FindOverlapsByIndexing"
              -sorted_keys = a list of keys generated by "FindOverlapsByIndexing"
              -outfile = the name of the outfile

        Output:
              -Writes a file with the following format:
                   > Each line records: Sequence RefStart RefStop RefRemainder #ofOverlaps
 
    '''

    # Make a list of lines to write
    outlines = []
    for key in sorted_keys:
        # Define a baseline for each ref feature
        remainder = overlap_dict[key][0]
        number = overlap_dict[key][-1]
        outlines.append("\t".join(list(key[0:3])) + "\t" + remainder + "\t" + str(number) + "\n")
 
    #Write the lines to file
    output = open(outfile,"w")
    output.write("".join(outlines))
    output.close()
