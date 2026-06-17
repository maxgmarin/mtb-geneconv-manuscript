#!/usr/bin/env python3

# Authors: Max Marin (mgmarin@g.harvard.edu),
# Purpose: Script for updating a BAM file for better visualization in IGV

# For each SAM entry it does the following:
#        a) Change the MQ to 60,
#        b) Change the SOFT clipping to HARD clipping,
#        c) Identify if read/genome is reverse complimented (+ or - strand)
#        d) Update the QuerySeq of the alignment, by subsetting the aligned genome/read based on HardClip lengths and +/- strand.


__doc__ = """Script for updating SAM files output by minimap2's self-alignment function (for visualization of sequence mismatches between homologous regions)"""

import sys
import argparse
import pysam
from Bio import SeqIO



def convert_SoftToHardClip(cigar_Tup):
    """ Function for updating all Soft-Clipping SAM CIGAR entries to Hard-Clipping """
        
    if cigar_Tup[0] != 4:
        return cigar_Tup
    else:
        cigar_List = list(cigar_Tup)
        cigar_List[0] = 5
                
        return tuple(cigar_List)


def updatePySam_Aln_MM2SelfAlign_ForIGVviz(i_Aln, i_SeqRecords):
    '''
    Function that updates a SAM entry (PySam alignment object). It does the following:
        a) Change the MQ to 60,
        b) Change the SOFT clipping to HARD clipping,
        c) Identify if read/genome is reverse complimented (+ or - strand)
        d) Update the QuerySeq of the alignment, by subsetting the aligned genome/read based on HardClip lengths and +/- strand.
    '''
    ## Step 1: Update the CIGAR string (Soft to Hard Clipping)
    
    # Create a new CIGAR where all SoftClipping becomes HardClipping
    i_Aln_UpdatedCigar_NoSoftClip = [ convert_SoftToHardClip(x) for x in i_Aln.cigar ]
    
    # Update the ALN Cigar
    i_Aln.cigar = i_Aln_UpdatedCigar_NoSoftClip
    
    ## Step 2: Check what the 
    if (i_Aln.cigar[0][0] == 4) or (i_Aln.cigar[0][0] == 5):
        FirstClip_Len = i_Aln.cigar[0][1]
    else:
        FirstClip_Len = 0
        
        
    if (i_Aln.cigar[-1][0] == 4) or (i_Aln.cigar[-1][0] == 5):
        LastClip_Len = i_Aln.cigar[-1][1]
    else:
        LastClip_Len = 0
        

    
    ## Step 3: Change MQ to 60
    i_Aln.mapq = 60 
    
    ## Step 4: Check if aligned sequence is reverse complimented (or not)
    i_Aln_IsRevComp = i_Aln.is_reverse

    
    
    ## Step 5: Slice the aligned genome (The query seq) based on the clipping lengths
    
    QuerySeq_Name = i_Aln.qname
    
    Query_SeqRecord = i_SeqRecords[QuerySeq_Name]
    Query_SeqRecord_RevComp = i_SeqRecords[QuerySeq_Name].reverse_complement()
    
    if not i_Aln_IsRevComp: # Clip normally if on the + strand
        if LastClip_Len != 0:
            New_QuerySeq = str(Query_SeqRecord.seq[FirstClip_Len : -LastClip_Len])  
        else:
            New_QuerySeq = str(Query_SeqRecord.seq[FirstClip_Len :]) 
            
    else: # Invert clipping if on the - strand
        if LastClip_Len != 0:
            New_QuerySeq = str( Query_SeqRecord_RevComp.seq[FirstClip_Len : -LastClip_Len])
        else:
            New_QuerySeq = str(Query_SeqRecord_RevComp.seq[FirstClip_Len :]) 

    
    i_Aln.query_sequence = New_QuerySeq
        
    return i_Aln
    




def main():

    # 1) Parse user arguments
    parser = argparse.ArgumentParser(
        description="")

    parser.add_argument('--input_sam', type=str, 
                        help="Input SAM/BAM file from Minimap2's self-alignment function")


    parser.add_argument('--aligned_sequence_fasta', type=str, 
                        help="FASTA file containing the genome(s)/sequence(s) used as input to minimap2's self-alignment function")

    parser.add_argument('--output_sam', type=str,
                        help="Updated SAM file where the query sequence has been filled in")

    args = parser.parse_args()

    # 2) Define input and output paths
    input_SAMorBAM = args.input_sam
    output_Updated_SAM_PATH = args.output_sam

    aligned_Seq_FA_PATH = args.aligned_sequence_fasta


    # 3) Parse the aligned sequences (from FASTA)
    input_SeqRecordsDict = {}


    print("List of all sequence IDs in FASTA:")
    for index, record in enumerate(SeqIO.parse(aligned_Seq_FA_PATH, "fasta")):
        
        record_ID = record.id
        print(record_ID)
        
        input_SeqRecordsDict[record_ID] = record
        

    # 4) Process the input BAM and output to the new SAM file 
    input_MM2_SA_PySam = pysam.AlignmentFile(input_SAMorBAM, "r")

    In_Header = input_MM2_SA_PySam.header

    with pysam.AlignmentFile(output_Updated_SAM_PATH, "w", header = In_Header) as outf:
        for  in_Aln in input_MM2_SA_PySam.fetch():
        
            out_Aln = updatePySam_Aln_MM2SelfAlign_ForIGVviz(in_Aln, input_SeqRecordsDict)

            outf.write(out_Aln)


if __name__ == "__main__":
    sys.exit(main())

