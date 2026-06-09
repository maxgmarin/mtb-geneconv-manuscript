# _____.py

import pandas as pd
import numpy as np

from dna_features_viewer import GraphicFeature, GraphicRecord, CircularGraphicRecord, BiopythonTranslator

from Bio import SeqIO

import bioframe as bf

import matplotlib.pyplot as plt

import matplotlib.ticker as mticker

import subprocess

from typing import Tuple, Dict, Any, Optional

from tqdm import tqdm



from .kmer import read_kmers_from_file, build_kmers, hash_kmers, hash_kmers_ToUnqNP, hash_kmers_ToSet, jaccard_containment_FromSets

from .general import Rv_dist

from .genomeviz_utils import generate_SNP_GraphicFeatures





def insertSNPs_IntoRef(in_SNPs_DF, In_Seq):
    Out_Seq_List = list(In_Seq)

    for i, row in in_SNPs_DF.iterrows():
        pos_0based = row["Pos_1based"] - 1
        Parent_Nt = row["Parent_Call"]
        Child_Nt = row["Child_Call"]
 
        Out_Seq_List[pos_0based] = Child_Nt
        #print(f"Changed the {pos_0based} to {Child_Nt}")

    Out_Seq_Mutated = "".join(Out_Seq_List)
    return Out_Seq_Mutated


def insert_ParentSNPs_IntoRef(in_SNPs_DF, In_Seq):
    Out_Seq_List = list(In_Seq)

    for i, row in in_SNPs_DF.iterrows():
        pos_0based = row["Pos_1based"] - 1
        Parent_Nt = row["Parent_Call"]
        Child_Nt = row["Child_Call"]
 
        Out_Seq_List[pos_0based] = Parent_Nt

    Out_Seq_Mutated = "".join(Out_Seq_List)
    return Out_Seq_Mutated



def extract_EventInfo(event_row):
    
    # A) Pull out event info for event of interest

    EventInfo = {}
    
    EventInfo["RE_EventID"] = event_row["EventID"]
    EventInfo["RE_Ovrlap_Genes"] = event_row["Overlap_Genes"]
    EventInfo["RE_Start"] = event_row["start_0based"]
    EventInfo["RE_End"] = event_row["end_1based"]

    EventInfo["RE_SNP_Count"] = event_row["snp_count"]

    EventInfo["RE_Parent_Node"] = event_row["Parent_Node"]
    EventInfo["RE_Child_Node"] = event_row["Child_Node"]

    EventInfo["RE_Lineage"] = event_row["Lineage"]
    EventInfo["taxa_List"] = event_row["taxa_List"]

    RE_Start, RE_End = event_row["start_0based"], event_row["end_1based"]
    EventInfo["RE_Coords"] = f"NC_000962.3:{RE_Start}-{RE_End}"

    RE_Child_Node = EventInfo["RE_Child_Node"]
    if RE_Child_Node[:4] == "Node":
        EventInfo["RE_TermNode"] = False
    else:
        EventInfo["RE_TermNode"] = True

    return EventInfo
    


def Compare_RE_Vs_Paralog_Hashes(i_Recomb_Hashes_NonUNQ, i_p_Q_Hashes):

    listOf_KmerMatchStatus = []

    for i, i_hash in enumerate(i_Recomb_Hashes_NonUNQ):
        
        # Check if k-mer hash is present in the paralog of interest
        i_k_Pres = int(i_hash in i_p_Q_Hashes)

        listOf_KmerMatchStatus.append(i_k_Pres)

    return np.array(listOf_KmerMatchStatus)





def compare_EventSeq_kmers_To_Paralogs(i_EventInfo_Row,
                                       i_AllEvent_SNPs,
                                       i_HomologyRegions_DF,
                                       i_RefSeq,
                                       i_RefSeq_K11Hashes = "None"):
                                         #i_Ref_GraphicRecord):
    
    #print(i_EventInfo_Row)
    # 1) Pull out the event specific info (Event_152)
    i_EventInfo = extract_EventInfo(i_EventInfo_Row)

    RE_Start = i_EventInfo["RE_Start"]
    RE_End = i_EventInfo["RE_End"]
    RE_Coords = i_EventInfo["RE_Coords"]
    RE_EventID =  i_EventInfo["RE_EventID"]
                                     
    #print(RE_EventID)                            
    i_Pad_Size = 600 
    
    i_RE_Start_Pad = RE_Start - i_Pad_Size
    i_RE_End_Pad = RE_End + i_Pad_Size
    k_size = 11

                                     
    # 2) Query SNPs associated with the event (From Base-Reconstruction by Gubbins)

    Event_SNPs = i_AllEvent_SNPs.query(f" (EventID == '{RE_EventID}') ")    

    # 3) Infer the EVENT seq and the ANCESTRAL seq (SNPs only) in WHOLE ref genome

    i_Rv_Seq_WiEventSNPs = insertSNPs_IntoRef(Event_SNPs, i_RefSeq)
    i_Rv_Seq_WiAncestralSNPs = insert_ParentSNPs_IntoRef(Event_SNPs, i_RefSeq)


                                     
                                     
    # 4) Get all paralogs (Homologous sequences) to events of interest.

    # E) Parse Homologous regions (Query homologous regions)
    sub_DF_Overlap_Hm_Regions = bf.select(i_HomologyRegions_DF, RE_Coords, cols = ("Target_Name", "Target_Start", "Target_End")) 
    
    Num_Overlapping_Hm_Regions = sub_DF_Overlap_Hm_Regions.shape[0]


    # 5) Create Dict of sequence info for ancestral seq + paralogous seqs



    # A) Get k-mer info from REseq-PADDED
    
    i_RecombSeq_FromSNPsInRv_PADDED = i_Rv_Seq_WiEventSNPs[i_RE_Start_Pad : i_RE_End_Pad]
    
    RecombSeq_PADDED_Kmers = build_kmers(str(i_RecombSeq_FromSNPsInRv_PADDED), 11)
    RecombSeq_PADDED_Hashes = np.array( hash_kmers(RecombSeq_PADDED_Kmers) )
                                     
    #print(RecombSeq_PADDED_Hashes.shape[0])

    #print("Length of PADDED event region w/ EVENT SNPs ", len(i_RecombSeq_FromSNPsInRv_PADDED))

                                     
    # Infer ANCESTRAL-seq w/ padding
    i_AncestralSeq_FromSNPsInRv_PADDED        = i_Rv_Seq_WiAncestralSNPs[i_RE_Start_Pad : i_RE_End_Pad]
    
    i_AncestralSeq_FromSNPsInRv_PADDED_Kmers  = build_kmers( str(i_AncestralSeq_FromSNPsInRv_PADDED),
                                                           k_size)
    
    i_AncestralSeq_FromSNPsInRv_PADDED_Hashes = hash_kmers_ToUnqNP(i_AncestralSeq_FromSNPsInRv_PADDED_Kmers  )  

    #print("Length of PADDED event region w/ ANC SNPs ", len(i_AncestralSeq_FromSNPsInRv_PADDED))
                   
                                     
    i_dictOf_Seqs_ToCompare = {}
    i_dictOf_Seqs_ToCompare['Ancestral_Seq'] = {}
    i_dictOf_Seqs_ToCompare['Ancestral_Seq']["Seq"] = i_AncestralSeq_FromSNPsInRv_PADDED
    
    i_dictOf_Seqs_ToCompare['Ancestral_Seq']["Kmers"] = i_AncestralSeq_FromSNPsInRv_PADDED_Kmers
    
    i_dictOf_Seqs_ToCompare['Ancestral_Seq']["Hashes"] = i_AncestralSeq_FromSNPsInRv_PADDED_Hashes 

    # C) Get k-mer info for each paralog
    
    for i, row in sub_DF_Overlap_Hm_Regions.iterrows():
    
        p_QueryGeneNames = row["QueryOverlap_Genes"]
        p_Q_Chr = row["Query_Name"]
        p_Q_Start = row["Query_Start"]
        p_Q_End = row["Query_End"]
        
    
        # If query is from H37Rv, subset seq from H37Rv genome
        if p_Q_Chr == "NC_000962.3":
            p_Q_Seq = str( i_RefSeq[p_Q_Start: p_Q_End] )
            
            #p_Name = f"{p_QueryGeneNames}-{p_Q_Chr}-{p_Q_Start}-{p_Q_End}"
            p_Name = f"{p_QueryGeneNames}-H37Rv-{p_Q_Start}-{p_Q_End}"
            
            p_Kmers = build_kmers(str(p_Q_Seq), k_size)
            
            p_Hashes = hash_kmers_ToUnqNP(p_Kmers)
                                  
            i_dictOf_Seqs_ToCompare[p_Name] = {}
            
            i_dictOf_Seqs_ToCompare[p_Name]["Seq"] = p_Q_Seq
            i_dictOf_Seqs_ToCompare[p_Name]["Kmers"] = p_Kmers
            i_dictOf_Seqs_ToCompare[p_Name]["Hashes"] = p_Hashes   
    
            #print(i, p_QueryGeneNames, p_Q_Start, p_Q_End, len(p_Q_Seq))  


    # 6) Compare ALL sequences against the recombinant sequence
    
    KmerComp_Results_V2 = {}
    
    for p_Name, SeqInfoDict in i_dictOf_Seqs_ToCompare.items():
    
        KmerComp_Results_V2[p_Name] = {}
        
        #print(p_Name, p_Q_Seq)
        #print(SeqInfoDict)
        i_p_Q_Hashes = SeqInfoDict["Hashes"]
    
        RESeq_Vs_Paralog_NP = Compare_RE_Vs_Paralog_Hashes(RecombSeq_PADDED_Hashes,
                                                           i_p_Q_Hashes)
        
        KmerComp_Results_V2[p_Name]["Kmatch_PADDED"] = RESeq_Vs_Paralog_NP 

        
    # Subset Anc-seq k-mer match for JUST the recomb region
    i_Pad_Shift = 10
    
    # Subset for JUST the recomb region
    
    AncSeq_Kmatch_NP =  KmerComp_Results_V2["Ancestral_Seq"]["Kmatch_PADDED"]
    
    AncSeq_Kmatch_RERegion = AncSeq_Kmatch_NP[i_Pad_Size - i_Pad_Shift : - (i_Pad_Size - i_Pad_Shift*2)]
    
    KmerComp_Results_V2["Ancestral_Seq"]["Kmatch_RERegion"] = AncSeq_Kmatch_RERegion
    
    
    RERegion_DiffToAnc_Mask = np.invert(AncSeq_Kmatch_RERegion.astype(bool))
    RERegion_DiffToAnc_Mask_Int = RERegion_DiffToAnc_Mask.astype(int)
    
    RERegion_SameToAnc_Mask = AncSeq_Kmatch_RERegion.astype(bool)
    
    KmerComp_Results_V2["Ancestral_Seq"]["RecombVsAncestral_Diff_Mask"] = RERegion_DiffToAnc_Mask
    KmerComp_Results_V2["Ancestral_Seq"]["RecombVsAncestral_Diff_Int"] = RERegion_DiffToAnc_Mask_Int
    KmerComp_Results_V2["Ancestral_Seq"]["RecombVsAncestral_Same_Mask"] = RERegion_SameToAnc_Mask
    
    KmerComp_Results_V2["Ancestral_Seq"]["Kmatch_Score"] = 0
    KmerComp_Results_V2["Ancestral_Seq"]["Kmatch_Num"] = 0
    KmerComp_Results_V2["Ancestral_Seq"]["Kmatch_NumEvaluated"] = 0

                                     
    # Now process info for the paralogs (All sequence that are not the ancestral seq)
    for p_Name, SeqInfoDict in i_dictOf_Seqs_ToCompare.items():
        
        if p_Name != 'Ancestral_Seq':
            p_Kmatch_NP = KmerComp_Results_V2[p_Name]["Kmatch_PADDED"]
            
            p_Kmatch_RERegion = p_Kmatch_NP[i_Pad_Size - i_Pad_Shift : - (i_Pad_Size - i_Pad_Shift*2)]
            
            p_Kmatch_RERegion_AncDiffPos = p_Kmatch_RERegion[RERegion_DiffToAnc_Mask] # Note the subsetted "RERegion" specific NP array should be shifted over by 5 idxs. (?)
    
            p_Kmatch_RERegion_MaskedByAncDiff = np.where(RERegion_SameToAnc_Mask, 0, p_Kmatch_RERegion)
    
            KmerComp_Results_V2[p_Name]["Kmatch_RERegion"] = AncSeq_Kmatch_RERegion
            KmerComp_Results_V2[p_Name]["Kmatch_RERegion_DiffToAncPos"] = p_Kmatch_RERegion_AncDiffPos
            KmerComp_Results_V2[p_Name]["Kmatch_RERegion_MaskedDiffKmers"] = p_Kmatch_RERegion_MaskedByAncDiff

            # print(p_Kmatch_RERegion_AncDiffPos)
            if p_Kmatch_RERegion_AncDiffPos.shape[0] > 0:
                match_Score_AncDiffPos = p_Kmatch_RERegion_AncDiffPos.sum() / p_Kmatch_RERegion_AncDiffPos.shape[0]
            
            else: 
                match_Score_AncDiffPos = 0
                
            #print(f"{p_Name} has a k-mer match of {match_Score_AncDiffPos} to the recomb event: {RE_EventID}")
            
            KmerComp_Results_V2[p_Name]["Kmatch_Score"] = match_Score_AncDiffPos
            KmerComp_Results_V2[p_Name]["Kmatch_Num"] = p_Kmatch_RERegion_AncDiffPos.sum()
            KmerComp_Results_V2[p_Name]["Kmatch_NumEvaluated"] = p_Kmatch_RERegion_AncDiffPos.shape[0]



    # Compare CHANGE k-mers to entire H37Rv genome

    if isinstance(i_RefSeq_K11Hashes, np.ndarray):
        KmerComp_Results_V2["Entire_Rv_Seq"] = {}

        RESeq_Vs_AllRvSeq_NP = Compare_RE_Vs_Paralog_Hashes(RecombSeq_PADDED_Hashes,
                                                            i_RefSeq_K11Hashes)
        
        KmerComp_Results_V2["Entire_Rv_Seq"]["Kmatch_PADDED"] = RESeq_Vs_AllRvSeq_NP 

        p_Kmatch_NP = RESeq_Vs_AllRvSeq_NP
            
        p_Kmatch_RERegion = p_Kmatch_NP[i_Pad_Size - i_Pad_Shift : - (i_Pad_Size - i_Pad_Shift*2)]
        
        p_Kmatch_RERegion_AncDiffPos = p_Kmatch_RERegion[RERegion_DiffToAnc_Mask] # Note the subsetted "RERegion" specific NP array should be shifted over by 5 idxs. (?)

        p_Kmatch_RERegion_MaskedByAncDiff = np.where(RERegion_SameToAnc_Mask, 0, p_Kmatch_RERegion)

        KmerComp_Results_V2["Entire_Rv_Seq"]["Kmatch_RERegion"] = AncSeq_Kmatch_RERegion
        KmerComp_Results_V2["Entire_Rv_Seq"]["Kmatch_RERegion_DiffToAncPos"] = p_Kmatch_RERegion_AncDiffPos
        KmerComp_Results_V2["Entire_Rv_Seq"]["Kmatch_RERegion_MaskedDiffKmers"] = p_Kmatch_RERegion_MaskedByAncDiff

        if p_Kmatch_RERegion_AncDiffPos.shape[0] > 0:
            match_Score_AncDiffPos = p_Kmatch_RERegion_AncDiffPos.sum() / p_Kmatch_RERegion_AncDiffPos.shape[0]
        else: 
            match_Score_AncDiffPos = 0
        
        KmerComp_Results_V2["Entire_Rv_Seq"]["Kmatch_Score"] = match_Score_AncDiffPos
        KmerComp_Results_V2["Entire_Rv_Seq"]["Kmatch_Num"] = p_Kmatch_RERegion_AncDiffPos.sum()
        KmerComp_Results_V2["Entire_Rv_Seq"]["Kmatch_NumEvaluated"] = p_Kmatch_RERegion_AncDiffPos.shape[0]
                                     

    
    return KmerComp_Results_V2, i_dictOf_Seqs_ToCompare





def compare_EventSeq_kmers_To_Paralogs_V2(i_EventInfo_Row,
                                         i_AllEvent_SNPs,
                                         i_HomologyRegions_DF,
                                         i_RefSeq,
                                         i_RefSeq_K11Hashes = "None"):
                                         #i_Ref_GraphicRecord):
    
    #print(i_EventInfo_Row)
    # 1) Pull out the event specific info (Event_152)
    i_EventInfo = extract_EventInfo(i_EventInfo_Row)

    RE_Start = i_EventInfo["RE_Start"]
    RE_End = i_EventInfo["RE_End"]
    RE_Coords = i_EventInfo["RE_Coords"]
    RE_EventID =  i_EventInfo["RE_EventID"]
                                     
    #print(RE_EventID)                            
    i_Pad_Size = 600 
    
    i_RE_Start_Pad = RE_Start - i_Pad_Size
    i_RE_End_Pad = RE_End + i_Pad_Size
    k_size = 11

                                     
    # 2) Query SNPs associated with the event (From Base-Reconstruction by Gubbins)

    Event_SNPs = i_AllEvent_SNPs.query(f" (EventID == '{RE_EventID}') ")    

    # 3) Infer the EVENT seq and the ANCESTRAL seq (SNPs only) in WHOLE ref genome

    i_Rv_Seq_WiEventSNPs = insertSNPs_IntoRef(Event_SNPs, i_RefSeq)
    i_Rv_Seq_WiAncestralSNPs = insert_ParentSNPs_IntoRef(Event_SNPs, i_RefSeq)


                                     
                                     
    # 4) Get all paralogs (Homologous sequences) to events of interest.

    # E) Parse Homologous regions (Query homologous regions)
    sub_DF_Overlap_Hm_Regions = bf.select(i_HomologyRegions_DF, RE_Coords, cols = ("Target_Name", "Target_Start", "Target_End")) 
    
    Num_Overlapping_Hm_Regions = sub_DF_Overlap_Hm_Regions.shape[0]


    # 5) Create Dict of sequence info for ancestral seq + paralogous seqs



    # A) Get k-mer info from REseq-PADDED
    
    i_RecombSeq_FromSNPsInRv_PADDED = i_Rv_Seq_WiEventSNPs[i_RE_Start_Pad : i_RE_End_Pad]
    
    RecombSeq_PADDED_Kmers = build_kmers(str(i_RecombSeq_FromSNPsInRv_PADDED), 11)
    RecombSeq_PADDED_Hashes = np.array( hash_kmers(RecombSeq_PADDED_Kmers) )
                                     
    #print(RecombSeq_PADDED_Hashes.shape[0])

    #print("Length of PADDED event region w/ EVENT SNPs ", len(i_RecombSeq_FromSNPsInRv_PADDED))

              
    # Infer ANCESTRAL-seq w/ padding
    i_AncestralSeq_FromSNPsInRv_PADDED        = i_Rv_Seq_WiAncestralSNPs[i_RE_Start_Pad : i_RE_End_Pad]
    
    i_AncestralSeq_FromSNPsInRv_PADDED_Kmers  = build_kmers( str(i_AncestralSeq_FromSNPsInRv_PADDED),
                                                           k_size)
    
    i_AncestralSeq_FromSNPsInRv_PADDED_Hashes = hash_kmers_ToUnqNP(i_AncestralSeq_FromSNPsInRv_PADDED_Kmers  )  

    #print("Length of PADDED event region w/ ANC SNPs ", len(i_AncestralSeq_FromSNPsInRv_PADDED))
                   
                                     
    i_dictOf_Seqs_ToCompare = {}
    i_dictOf_Seqs_ToCompare['Ancestral_Seq'] = {}
    i_dictOf_Seqs_ToCompare['Ancestral_Seq']["Seq"] = i_AncestralSeq_FromSNPsInRv_PADDED
    
    i_dictOf_Seqs_ToCompare['Ancestral_Seq']["Kmers"] = i_AncestralSeq_FromSNPsInRv_PADDED_Kmers
    
    i_dictOf_Seqs_ToCompare['Ancestral_Seq']["Hashes"] = i_AncestralSeq_FromSNPsInRv_PADDED_Hashes 

    # C) Get k-mer info for each paralog
    
    for i, row in sub_DF_Overlap_Hm_Regions.iterrows():
    
        p_QueryGeneNames = row["QueryOverlap_Genes"]
        p_Q_Chr = row["Query_Name"]
        p_Q_Start = row["Query_Start"]
        p_Q_End = row["Query_End"]
        
    
        # If query is from H37Rv, subset seq from H37Rv genome
        if p_Q_Chr == "NC_000962.3":
            p_Q_Seq = str( i_RefSeq[p_Q_Start: p_Q_End] )
            
            #p_Name = f"{p_QueryGeneNames}-{p_Q_Chr}-{p_Q_Start}-{p_Q_End}"
            p_Name = f"{p_QueryGeneNames}-H37Rv-{p_Q_Start}-{p_Q_End}"
            
            p_Kmers = build_kmers(str(p_Q_Seq), k_size)
            
            p_Hashes = hash_kmers_ToUnqNP(p_Kmers)
                                  
            i_dictOf_Seqs_ToCompare[p_Name] = {}
            
            i_dictOf_Seqs_ToCompare[p_Name]["Seq"] = p_Q_Seq
            i_dictOf_Seqs_ToCompare[p_Name]["Kmers"] = p_Kmers
            i_dictOf_Seqs_ToCompare[p_Name]["Hashes"] = p_Hashes   
    
            #print(i, p_QueryGeneNames, p_Q_Start, p_Q_End, len(p_Q_Seq))  


    # 6) Compare ALL sequences against the recombinant sequence
    
    KmerComp_Results_V2 = {}
    
    for p_Name, SeqInfoDict in i_dictOf_Seqs_ToCompare.items():
    
        KmerComp_Results_V2[p_Name] = {}
        
        #print(p_Name, p_Q_Seq)
        #print(SeqInfoDict)
        i_p_Q_Hashes = SeqInfoDict["Hashes"]
    
        RESeq_Vs_Paralog_NP = Compare_RE_Vs_Paralog_Hashes(RecombSeq_PADDED_Hashes,
                                                           i_p_Q_Hashes)
        
        KmerComp_Results_V2[p_Name]["Kmatch_PADDED"] = RESeq_Vs_Paralog_NP 

        
    # Subset Anc-seq k-mer match for JUST the recomb region
    i_Pad_Shift = 10
    
    # Subset for JUST the recomb region
    
    AncSeq_Kmatch_NP =  KmerComp_Results_V2["Ancestral_Seq"]["Kmatch_PADDED"]
    
    AncSeq_Kmatch_RERegion = AncSeq_Kmatch_NP[i_Pad_Size - i_Pad_Shift : - (i_Pad_Size - i_Pad_Shift*2)]
    
    KmerComp_Results_V2["Ancestral_Seq"]["Kmatch_RERegion"] = AncSeq_Kmatch_RERegion
    
    
    RERegion_DiffToAnc_Mask = np.invert(AncSeq_Kmatch_RERegion.astype(bool))
    RERegion_DiffToAnc_Mask_Int = RERegion_DiffToAnc_Mask.astype(int)
    
    RERegion_SameToAnc_Mask = AncSeq_Kmatch_RERegion.astype(bool)
    
    KmerComp_Results_V2["Ancestral_Seq"]["RecombVsAncestral_Diff_Mask"] = RERegion_DiffToAnc_Mask
    KmerComp_Results_V2["Ancestral_Seq"]["RecombVsAncestral_Diff_Int"] = RERegion_DiffToAnc_Mask_Int
    KmerComp_Results_V2["Ancestral_Seq"]["RecombVsAncestral_Same_Mask"] = RERegion_SameToAnc_Mask
    
    KmerComp_Results_V2["Ancestral_Seq"]["Kmatch_Score"] = 0
    KmerComp_Results_V2["Ancestral_Seq"]["Kmatch_Num"] = 0
    KmerComp_Results_V2["Ancestral_Seq"]["Kmatch_NumEvaluated"] = 0

                                     
    # Now process info for the paralogs (All sequence that are not the ancestral seq)
    for p_Name, SeqInfoDict in i_dictOf_Seqs_ToCompare.items():
        
        if p_Name != 'Ancestral_Seq':
            p_Kmatch_NP = KmerComp_Results_V2[p_Name]["Kmatch_PADDED"]
            
            p_Kmatch_RERegion = p_Kmatch_NP[i_Pad_Size - i_Pad_Shift : - (i_Pad_Size - i_Pad_Shift*2)]
            
            p_Kmatch_RERegion_AncDiffPos = p_Kmatch_RERegion[RERegion_DiffToAnc_Mask] # Note the subsetted "RERegion" specific NP array should be shifted over by 5 idxs. (?)
    
            p_Kmatch_RERegion_MaskedByAncDiff = np.where(RERegion_SameToAnc_Mask, 0, p_Kmatch_RERegion)
    
            KmerComp_Results_V2[p_Name]["Kmatch_RERegion"] = AncSeq_Kmatch_RERegion
            KmerComp_Results_V2[p_Name]["Kmatch_RERegion_DiffToAncPos"] = p_Kmatch_RERegion_AncDiffPos
            KmerComp_Results_V2[p_Name]["Kmatch_RERegion_MaskedDiffKmers"] = p_Kmatch_RERegion_MaskedByAncDiff

            # print(p_Kmatch_RERegion_AncDiffPos)
            if p_Kmatch_RERegion_AncDiffPos.shape[0] > 0:
                match_Score_AncDiffPos = p_Kmatch_RERegion_AncDiffPos.sum() / p_Kmatch_RERegion_AncDiffPos.shape[0]
            
            else: 
                match_Score_AncDiffPos = 0
                
            #print(f"{p_Name} has a k-mer match of {match_Score_AncDiffPos} to the recomb event: {RE_EventID}")
            
            KmerComp_Results_V2[p_Name]["Kmatch_Score"] = match_Score_AncDiffPos
            KmerComp_Results_V2[p_Name]["Kmatch_Num"] = p_Kmatch_RERegion_AncDiffPos.sum()
            KmerComp_Results_V2[p_Name]["Kmatch_NumEvaluated"] = p_Kmatch_RERegion_AncDiffPos.shape[0]




    RecombVsAncSeq_Dict = {}


    i_RecombSeq_Hashes = hash_kmers_ToSet(RecombSeq_PADDED_Kmers  )  
    i_AncSeq_Hashes = hash_kmers_ToSet(i_AncestralSeq_FromSNPsInRv_PADDED_Kmers  )  

    Hashes_UnqToRecombVsAncestral = i_RecombSeq_Hashes - i_AncSeq_Hashes
    Kmers_UnqToRecombVsAncestral = set(RecombSeq_PADDED_Kmers) - set(i_AncestralSeq_FromSNPsInRv_PADDED_Kmers)

    # print(len(RecombSeq_PADDED_Kmers),
    #     len(i_AncestralSeq_FromSNPsInRv_PADDED_Kmers),
    #     len(Hashes_UnqToRecombVsAncestral),
    #     len(Kmers_UnqToRecombVsAncestral))  

    RecombVsAncSeq_Dict["Hashes_UnqToRecombVsAncestral"] = Hashes_UnqToRecombVsAncestral
    RecombVsAncSeq_Dict["Kmers_UnqToRecombVsAncestral"] = Kmers_UnqToRecombVsAncestral

    KmerComp_Results_V2["Ancestral_Seq"]["RecombVsAncSeq_Dict"] = RecombVsAncSeq_Dict
    # Compare CHANGE k-mers to entire H37Rv genome

    if not isinstance(i_RefSeq_K11Hashes, str):
        KmerComp_Results_V2["Entire_Rv_Seq"] = {}


        RESeq_Vs_AllRvSeq_NP = Compare_RE_Vs_Paralog_Hashes(RecombSeq_PADDED_Hashes, #np.unique(np.array(list(Hashes_UnqToRecombVsAncestral))),
                                                            np.array(list(i_RefSeq_K11Hashes))) #i_RefSeq_K11Hashes)
        
        KmerComp_Results_V2["Entire_Rv_Seq"]["Kmatch_PADDED"] = RESeq_Vs_AllRvSeq_NP 

        p_Kmatch_NP = RESeq_Vs_AllRvSeq_NP
            
        p_Kmatch_RERegion = p_Kmatch_NP[i_Pad_Size - i_Pad_Shift : - (i_Pad_Size - i_Pad_Shift*2)]
        
        p_Kmatch_RERegion_AncDiffPos = p_Kmatch_RERegion[RERegion_DiffToAnc_Mask] # Note the subsetted "RERegion" specific NP array should be shifted over by 5 idxs. (?)

        p_Kmatch_RERegion_MaskedByAncDiff = np.where(RERegion_SameToAnc_Mask, 0, p_Kmatch_RERegion)

        KmerComp_Results_V2["Entire_Rv_Seq"]["Kmatch_RERegion"] = AncSeq_Kmatch_RERegion
        KmerComp_Results_V2["Entire_Rv_Seq"]["Kmatch_RERegion_DiffToAncPos"] = p_Kmatch_RERegion_AncDiffPos
        KmerComp_Results_V2["Entire_Rv_Seq"]["Kmatch_RERegion_MaskedDiffKmers"] = p_Kmatch_RERegion_MaskedByAncDiff

        if p_Kmatch_RERegion_AncDiffPos.shape[0] > 0:
            match_Score_AncDiffPos = p_Kmatch_RERegion_AncDiffPos.sum() / p_Kmatch_RERegion_AncDiffPos.shape[0]
        else: 
            match_Score_AncDiffPos = 0
        
        KmerComp_Results_V2["Entire_Rv_Seq"]["Kmatch_Score"] = match_Score_AncDiffPos
        KmerComp_Results_V2["Entire_Rv_Seq"]["Kmatch_Num"] = p_Kmatch_RERegion_AncDiffPos.sum()
        KmerComp_Results_V2["Entire_Rv_Seq"]["Kmatch_NumEvaluated"] = p_Kmatch_RERegion_AncDiffPos.shape[0]

        KmerComp_Results_V2["Entire_Rv_Seq"]["RecombVsAncSeq_Dict"] = RecombVsAncSeq_Dict
        
        if len(Hashes_UnqToRecombVsAncestral) > 0:
            KmerComp_Results_V2["Entire_Rv_Seq"]["JC_RecombKmers_ToAllRvSeq"] = jaccard_containment_FromSets(Hashes_UnqToRecombVsAncestral,
                                                                                                             i_RefSeq_K11Hashes )
        else:
            KmerComp_Results_V2["Entire_Rv_Seq"]["JC_RecombKmers_ToAllRvSeq"] = np.nan

    
    return KmerComp_Results_V2, i_dictOf_Seqs_ToCompare




























# def check_snp_overlap(i_Event_SNPs, i_Event_HmRegion_Combo_Var_DF):
#     # Calculate the total number of SNPs associated with the event
#     RE_SNP_Count = len(i_Event_SNPs)
    
#     # Merging the two dataframes on specific columns to find overlapping SNPs
#     df3 = pd.merge(i_Event_SNPs, i_Event_HmRegion_Combo_Var_DF, how='inner', right_on='Query_End', left_on='Pos_1based')
    
#     # Filtering rows where the base overlap conditions are met
#     df3_BaseOvrlap = df3[ ((df3["Parent_Call"] == df3["Ref"]) & (df3["Child_Call"] == df3["Alt"])) |
#                            ((df3["Parent_Call"] == df3["Alt"]) & (df3["Child_Call"] == df3["Ref"])) ]
    
#     # Counting the number of unique positions with overlaps
#     NumOvrlap = df3_BaseOvrlap["Pos_1based"].nunique()
    
#     # Checking for a specific condition to issue a warning
#     if NumOvrlap > RE_SNP_Count:
#         print("WARNING, # SNPs match greater than possible Event SNPs")
    
#     # Return the number of overlapping SNPs
#     return NumOvrlap


def check_snp_overlap_V2(i_Event_SNPs, i_Event_HmRegion_Combo_Var_DF):
    # Calculate the total number of SNPs associated with the event
    RE_SNP_Count = len(i_Event_SNPs)
    
    # Merging the two dataframes on specific columns to find overlapping SNPs
    df3 = pd.merge(i_Event_SNPs, i_Event_HmRegion_Combo_Var_DF, how='inner', right_on='Target_End', left_on='Pos_1based')
    
    # Filtering rows where the base overlap conditions are met
    df3_BaseOvrlap = df3[ ((df3["Parent_Call"] == df3["Ref"]) & (df3["Child_Call"] == df3["Alt"])) |
                           ((df3["Parent_Call"] == df3["Alt"]) & (df3["Child_Call"] == df3["Ref"])) ]
    
    # Counting the number of unique positions with overlaps
    NumOvrlap = df3_BaseOvrlap["Pos_1based"].nunique()
    
    # Checking for a specific condition to issue a warning
    if NumOvrlap > RE_SNP_Count:
        print("WARNING, # SNPs match greater than possible Event SNPs")
    
    # Return the number of overlapping SNPs
    return NumOvrlap







############## Recomb Event K-mer Viz Functions ##############



def genViz_RecombEvent_Vs_Paralogs_V3(Genome_Graphic_Record, 
                                       i_KmerComp_Results_Dict,
                                       i_Event_SNPs,
                                       i_Target_Coords,
                                       i_Viz_Coords,
                                       i_Recomb_Coords,
                                       ShowMatch_ToAllRvKmers = True):
                                           
    i_RE_Start_Pad, i_RE_End_Pad =  i_Target_Coords
    
    Viz_Start, Viz_End = i_Viz_Coords
                                           
    i_RE_Start, i_RE_End = i_Recomb_Coords


    All_SeqIDs = i_KmerComp_Results_Dict.keys()
    
    All_SeqIDs = [seqID for seqID in All_SeqIDs if seqID != "Recomb_Seq"]


    if not ShowMatch_ToAllRvKmers:
        All_SeqIDs = [seqID for seqID in All_SeqIDs if seqID != "Entire_Rv_Seq"]
                                           
                                           
    # A) Setup subplots
    #NumPlots = len(All_SeqIDs) + 1
    NumPlots = len(All_SeqIDs)

    fig, axs = plt.subplots(NumPlots, 1, figsize=(15, 2 * NumPlots), sharex=True)

    
    # B) Pull out coord info
    RE_Start_PadShifted = i_RE_Start_Pad + 5
    RE_End_PadShifted = i_RE_End_Pad - 5
    
    #Start_Coord = i_RE_Start_Pad
    #End_Coord = i_RE_End_Pad
    
    # C) Create graphic records for viz

    ### C.1) Make and plot cropped gene annotations
    Genome_Graphic_Record_cropped = Genome_Graphic_Record.crop((Viz_Start, Viz_End + 1))
    Genome_Graphic_Record_cropped.plot(strand_in_label_threshold=5, ax = axs[0])

    #Mtb_H37rv_Graphic_Record_cropped_WiSNPs = AddSNPs_ToGraphicRecord(Genome_Graphic_Record_cropped, i_Event_SNPs)   
    #Mtb_H37rv_Graphic_Record_cropped_WiSNPs.plot(strand_in_label_threshold=5, ax = axs[0], plot_sequence = False)

    ### Create graphic record containing features objects for each SNP in the event ###
    i_Event_SNP_GFeats = generate_SNP_GraphicFeatures(i_Event_SNPs)

    i_Event_SNPs_GraphicRecord = GraphicRecord(sequence_length= 4411532, features = i_Event_SNP_GFeats)

    i_Event_SNPs_GraphicRecord_cropped = i_Event_SNPs_GraphicRecord.crop((Viz_Start, Viz_End + 1))

    sf = mticker.ScalarFormatter(useOffset=False)
    sf.set_scientific(False)
    axs[0].xaxis.set_major_formatter(sf)
    
    i_Event_SNPs_GraphicRecord_cropped.plot(strand_in_label_threshold=5, ax = axs[0])



    #Genome_Graphic_Record_cropped.plot(strand_in_label_threshold=5, ax = axs[0], plot_sequence = False)

                                           
    # D) Highlight regions of the recomb event that have changed k-mers
    
    pos_Coords_RERegion = np.arange(i_RE_Start - 5 , i_RE_End + 15 ) #+ 1)
    
    i_RERegion_DiffToAnc_Mask_Int = i_KmerComp_Results_Dict["Ancestral_Seq"]["RecombVsAncestral_Diff_Int"]
    
    axs[0].fill_between(pos_Coords_RERegion,
                        np.minimum(i_RERegion_DiffToAnc_Mask_Int, 0.5) * 100,
                        -1.5,
                        #y2 = -10.5
                        alpha = 0.15, color="red", linewidth = 0, label = "Recombinant region (Changed K-mers)")
                                           
    axs[0].set_title(f'Recombinant sequence w/ mutations highlighted')

    
    # E) Create a k-mer match plot for each paralog
    
    pos_Coords_1 = np.arange(RE_Start_PadShifted, RE_End_PadShifted) #+ 1)
                               
    #for i, p_Name in enumerate(All_SeqIDs):
    for i, p_Name in enumerate(All_SeqIDs[1:]):

        i_Kmatch_NP = i_KmerComp_Results_Dict[p_Name]["Kmatch_PADDED"]
    
        axs[i + 1].fill_between(pos_Coords_1, 
                                i_Kmatch_NP,
                                alpha=1, color="#FDDB91", label = 'Match to reference locus')
    
        if p_Name != "Ancestral_Seq":
            i_Kmatch_MaskedDiffKmers = i_KmerComp_Results_Dict[p_Name]["Kmatch_RERegion_MaskedDiffKmers"]
        
            axs[i + 1].fill_between(pos_Coords_RERegion, 
                                    i_Kmatch_MaskedDiffKmers,
                                    alpha = 0.95, color="red", label = 'Match to recombinant "k-mer Only"') 

        
        Kmatch_Score =  i_KmerComp_Results_Dict[p_Name]["Kmatch_Score"]
        Kmers_NumMatch = i_KmerComp_Results_Dict[p_Name]["Kmatch_Num"]
        Kmers_NumEvaluated = i_KmerComp_Results_Dict[p_Name]["Kmatch_NumEvaluated"] 

        Kmatch_Score_Rounded = str(round(Kmatch_Score, 4))
    

        axs[i + 1].set_title(f'{p_Name} --- k-mer match to putative recomb seq = {Kmatch_Score_Rounded} ({Kmers_NumMatch}/{Kmers_NumEvaluated})')

        
        axs[i + 1].set(ylabel='k-mer match')
        
        axs[i + 1].legend(loc="upper right", bbox_to_anchor=(1.3, 1.05))
        axs[i + 1].set_ylim(0, 1)
        
    axs[-1].set(xlabel='H37rv Ref Position (bp)')


    Anc_Kmatch_NP = i_KmerComp_Results_Dict["Ancestral_Seq"]["Kmatch_PADDED"]

    #print(pos_Coords_1.shape)
    #print(Anc_Kmatch_NP.shape)

    axs[0].fill_between(pos_Coords_1, 
                        Anc_Kmatch_NP * 100,
                        -0.5,
                        alpha=0.15, color="#FDDB91", label = 'Match to reference (non-recombinant)locus')


    axs[0].legend(loc = "upper right", 
                  bbox_to_anchor=(1.3, 1.05))
    
    #axs[0].set_ylim(-0.5, 2)

    RE_Event_ID = i_Event_SNPs["EventID"].values[0]
    RE_NumSNPs = i_Event_SNPs.shape[0]
                                           
    fig.suptitle(f'{RE_Event_ID} - {RE_NumSNPs} SNPs - k-mer analysis relative to paralogs', fontsize=20)

    return fig, axs
                                           





def process_and_visualize_event(event_id,
                                dict_of_kmer_analysis,
                                gubbins_snps_df,
                                gre_df,
                                graphic_record,
                                ShowMatch_ToAllRvKmers = False,
                                zoom_out_dist=100,
                                pad_size=600,
                                pad_shift=10):
           
    # Retrieve Kmer Comparison Results for the given event
    kmer_comp_results = dict_of_kmer_analysis[event_id]

    # Query the SNP data for the event
    event_snps = gubbins_snps_df.query(f"EventID == '{event_id}'")

    # Extract event information from GRE dataframe
    event_info = extract_EventInfo(gre_df.query(f"EventID == '{event_id}'").iloc[0])

    # Calculate visualization parameters
    re_start, re_end = event_info['RE_Start'], event_info['RE_End']
    viz_start = re_start - zoom_out_dist
    viz_end = re_end + zoom_out_dist

    # Calculate padding for visualization boundaries
    re_start_pad = re_start - pad_size
    re_end_pad = re_end + pad_size

    # Generate the visualization
    fig, axs = genViz_RecombEvent_Vs_Paralogs_V3(graphic_record, 
                                                 kmer_comp_results,
                                                 event_snps,
                                                 (re_start_pad, re_end_pad),
                                                 (viz_start, viz_end),
                                                 (re_start, re_end),
                                                 ShowMatch_ToAllRvKmers)
    return fig, axs










































def getEventCoords_RvToAsm(i_EventInfo, i_LiftOvr_RvToAsm_PAF):

    EventCoords_TabDelim = "'" +  "\\t".join(["NC_000962.3",
                                              str(i_EventInfo["RE_Start"]),
                                              str(i_EventInfo["RE_End"]) ]) + "'"

    LiftOver_CMD = f'echo -e {EventCoords_TabDelim} | paftools.js liftover {i_LiftOvr_RvToAsm_PAF} -'      
    #print("LiftOver_CMD:", LiftOver_CMD)
    
    LiftOver_Out = subprocess.run( [LiftOver_CMD], shell=True, stdout=subprocess.PIPE).stdout
    LiftOver_List = LiftOver_Out.decode().split("\t")
    #print(LiftOver_List)
    #Asm_Event_Chr = LiftOver_List[0]
    Asm_Event_Start = int(LiftOver_List[1])
    Asm_Event_End = int(LiftOver_List[2])

    return Asm_Event_Start, Asm_Event_End


def getEventCoords_RvToAsm_WithPad(i_EventInfo, i_LiftOvr_RvToAsm_PAF, PadLen = 600):

    EventCoords_TabDelim = "'" +  "\\t".join(["NC_000962.3",
                                              str(i_EventInfo["RE_Start"] - PadLen),
                                              str(i_EventInfo["RE_End"] + PadLen) ]) + "'"

    LiftOver_CMD = f'echo -e {EventCoords_TabDelim} | paftools.js liftover {i_LiftOvr_RvToAsm_PAF} -'      
    #print("LiftOver_CMD:", LiftOver_CMD)
    
    LiftOver_Out = subprocess.run( [LiftOver_CMD], shell=True, stdout=subprocess.PIPE).stdout
    LiftOver_List = LiftOver_Out.decode().split("\t")
    #print(LiftOver_List)
    #Asm_Event_Chr = LiftOver_List[0]
    Asm_Event_Start = int(LiftOver_List[1])
    Asm_Event_End = int(LiftOver_List[2])

    return Asm_Event_Start, Asm_Event_End


def get_RecombSeq_FromTargetAsm(i_TarAsm_FA, i_LiftOvr_RvToAsm_PAF, i_EventInfo):

    # Read in the entire "target assembly" genome sequence (1 contig)
    i_TarAsm_Seq = SeqIO.read(i_TarAsm_FA, "fasta").seq

    # Use MM2 alignment + paftools.js liftover to get the corresponding coordinates in the genome of interest      
    i_TarAsm_Start, i_TarAsm_End = getEventCoords_RvToAsm(i_EventInfo, i_LiftOvr_RvToAsm_PAF)

    #print("Inferred coords in the Target Asm", i_TarAsm_Start, i_TarAsm_End)
    # Use liftover coordinates to subset for recomb sequence in 
    i_RecombSeq_FromTarAsm = str( i_TarAsm_Seq[i_TarAsm_Start: i_TarAsm_End])
    
    return i_RecombSeq_FromTarAsm


def compare_EventSeq_kmers_To_Paralogs_ProvideRecombSeq(i_EventInfo_Row,
                                                         i_AllEvent_SNPs,
                                                         i_SampleID_ToPaths_Dict,
                                                         i_HomologyRegions_DF,
                                                         i_RefSeq,):
                                         #i_Ref_GraphicRecord):
    
    #print(i_EventInfo_Row)
    # 1) Pull out the event specific info (Event_152)
    i_EventInfo = extract_EventInfo(i_EventInfo_Row)

    RE_Start = i_EventInfo["RE_Start"]
    RE_End = i_EventInfo["RE_End"]
    RE_Coords = i_EventInfo["RE_Coords"]
    RE_EventID =  i_EventInfo["RE_EventID"]
                                     
    #print(RE_EventID)                            
    i_Pad_Size = 600 
    
    i_RE_Start_Pad = RE_Start - i_Pad_Size
    i_RE_End_Pad = RE_End + i_Pad_Size
    k_size = 11

                      

    RE_TaxaList = i_EventInfo["taxa_List"]
    RE_FirstTaxa_Sample = RE_TaxaList.split(",")[0].split("'")[1]
    

    # Read in the entire "target assembly" genome sequence (1 contig)
    i_TarAsm_FASTA = i_SampleID_ToPaths_Dict[RE_FirstTaxa_Sample]['LR_Asm_FA']
    i_TarAsm_Seq = SeqIO.read(i_TarAsm_FASTA, "fasta").seq
    
    i_RvToAsm_PAF = i_SampleID_ToPaths_Dict[RE_FirstTaxa_Sample]['H37RvToAsm_ForLiftOff_PAF']

    i_EventSeq_Asm1Extracted = get_RecombSeq_FromTargetAsm(i_TarAsm_FASTA, i_RvToAsm_PAF, i_EventInfo)


    # Use MM2 alignment + paftools.js liftover to get the corresponding coordinates in the genome of interest      
    i_TarAsm_Start, i_TarAsm_End = getEventCoords_RvToAsm_WithPad(i_EventInfo, i_RvToAsm_PAF, PadLen = i_Pad_Size)
    
    i_RecombSeq_PADDED_FromTarAsm = str( i_TarAsm_Seq[i_TarAsm_Start: i_TarAsm_End])

    # 2) Query SNPs associated with the event (From Base-Reconstruction by Gubbins)
    
    Event_SNPs = i_AllEvent_SNPs.query(f" (EventID == '{RE_EventID}') ")    

    # 3) Infer the EVENT seq and the ANCESTRAL seq (SNPs only) in WHOLE ref genome

    i_Rv_Seq_WiEventSNPs = insertSNPs_IntoRef(Event_SNPs, i_RefSeq)
    i_Rv_Seq_WiAncestralSNPs = insert_ParentSNPs_IntoRef(Event_SNPs, i_RefSeq)


                                     
                                     
    # 4) Get all paralogs (Homologous sequences) to events of interest.

    # E) Parse Homologous regions (Query homologous regions)
    sub_DF_Overlap_Hm_Regions = bf.select(i_HomologyRegions_DF, RE_Coords, cols = ("Target_Name", "Target_Start", "Target_End")) 
    
    Num_Overlapping_Hm_Regions = sub_DF_Overlap_Hm_Regions.shape[0]


    # 5) Create Dict of sequence info for ancestral seq + paralogous seqs



    # A) Get k-mer info from REseq-PADDED
    
    #i_RecombSeq_FromSNPsInRv_PADDED = i_Rv_Seq_WiEventSNPs[i_RE_Start_Pad : i_RE_End_Pad]


    RecombSeq_PADDED_Kmers = build_kmers(str(i_RecombSeq_PADDED_FromTarAsm), 11)
    RecombSeq_PADDED_Hashes = np.array( hash_kmers(RecombSeq_PADDED_Kmers) )
                                     
    #print(RecombSeq_PADDED_Hashes.shape[0])

    #print("Length of PADDED event region w/ EVENT SNPs ", len(i_RecombSeq_FromSNPsInRv_PADDED))

              
    # Infer ANCESTRAL-seq w/ padding
    i_AncestralSeq_FromSNPsInRv_PADDED        = i_Rv_Seq_WiAncestralSNPs[i_RE_Start_Pad : i_RE_End_Pad]
    
    i_AncestralSeq_FromSNPsInRv_PADDED_Kmers  = build_kmers( str(i_AncestralSeq_FromSNPsInRv_PADDED),
                                                           k_size)
    
    i_AncestralSeq_FromSNPsInRv_PADDED_Hashes = hash_kmers_ToUnqNP(i_AncestralSeq_FromSNPsInRv_PADDED_Kmers  )  

    #print("Length of PADDED event region w/ ANC SNPs ", len(i_AncestralSeq_FromSNPsInRv_PADDED))
                   
                                     
    i_dictOf_Seqs_ToCompare = {}
    i_dictOf_Seqs_ToCompare['Ancestral_Seq'] = {}
    i_dictOf_Seqs_ToCompare['Ancestral_Seq']["Seq"] = i_AncestralSeq_FromSNPsInRv_PADDED
    
    i_dictOf_Seqs_ToCompare['Ancestral_Seq']["Kmers"] = i_AncestralSeq_FromSNPsInRv_PADDED_Kmers
    
    i_dictOf_Seqs_ToCompare['Ancestral_Seq']["Hashes"] = i_AncestralSeq_FromSNPsInRv_PADDED_Hashes 

    # C) Get k-mer info for each paralog
    
    for i, row in sub_DF_Overlap_Hm_Regions.iterrows():
    
        p_QueryGeneNames = row["QueryOverlap_Genes"]
        p_Q_Chr = row["Query_Name"]
        p_Q_Start = row["Query_Start"]
        p_Q_End = row["Query_End"]
        
    
        # If query is from H37Rv, subset seq from H37Rv genome
        if p_Q_Chr == "NC_000962.3":
            p_Q_Seq = str( i_RefSeq[p_Q_Start: p_Q_End] )
            
            #p_Name = f"{p_QueryGeneNames}-{p_Q_Chr}-{p_Q_Start}-{p_Q_End}"
            p_Name = f"{p_QueryGeneNames}-H37Rv-{p_Q_Start}-{p_Q_End}"
            
            p_Kmers = build_kmers(str(p_Q_Seq), k_size)
            
            p_Hashes = hash_kmers_ToUnqNP(p_Kmers)
                                  
            i_dictOf_Seqs_ToCompare[p_Name] = {}
            
            i_dictOf_Seqs_ToCompare[p_Name]["Seq"] = p_Q_Seq
            i_dictOf_Seqs_ToCompare[p_Name]["Kmers"] = p_Kmers
            i_dictOf_Seqs_ToCompare[p_Name]["Hashes"] = p_Hashes   
    
            #print(i, p_QueryGeneNames, p_Q_Start, p_Q_End, len(p_Q_Seq))  


    # 6) Compare ALL sequences against the recombinant sequence
    
    KmerComp_Results_V2 = {}
    
    for p_Name, SeqInfoDict in i_dictOf_Seqs_ToCompare.items():
    
        KmerComp_Results_V2[p_Name] = {}
        
        #print(p_Name, p_Q_Seq)
        #print(SeqInfoDict)
        i_p_Q_Hashes = SeqInfoDict["Hashes"]
    
        RESeq_Vs_Paralog_NP = Compare_RE_Vs_Paralog_Hashes(RecombSeq_PADDED_Hashes,
                                                           i_p_Q_Hashes)
        
        KmerComp_Results_V2[p_Name]["Kmatch_PADDED"] = RESeq_Vs_Paralog_NP 

        
    # Subset Anc-seq k-mer match for JUST the recomb region
    i_Pad_Shift = 10
    
    # Subset for JUST the recomb region
    
    AncSeq_Kmatch_NP =  KmerComp_Results_V2["Ancestral_Seq"]["Kmatch_PADDED"]
    
    AncSeq_Kmatch_RERegion = AncSeq_Kmatch_NP[i_Pad_Size - i_Pad_Shift : - (i_Pad_Size - i_Pad_Shift*2)]
    
    KmerComp_Results_V2["Ancestral_Seq"]["Kmatch_RERegion"] = AncSeq_Kmatch_RERegion
    
    
    RERegion_DiffToAnc_Mask = np.invert(AncSeq_Kmatch_RERegion.astype(bool))
    RERegion_DiffToAnc_Mask_Int = RERegion_DiffToAnc_Mask.astype(int)
    
    RERegion_SameToAnc_Mask = AncSeq_Kmatch_RERegion.astype(bool)
    
    KmerComp_Results_V2["Ancestral_Seq"]["RecombVsAncestral_Diff_Mask"] = RERegion_DiffToAnc_Mask
    KmerComp_Results_V2["Ancestral_Seq"]["RecombVsAncestral_Diff_Int"] = RERegion_DiffToAnc_Mask_Int
    KmerComp_Results_V2["Ancestral_Seq"]["RecombVsAncestral_Same_Mask"] = RERegion_SameToAnc_Mask
    
    KmerComp_Results_V2["Ancestral_Seq"]["Kmatch_Score"] = 0
    KmerComp_Results_V2["Ancestral_Seq"]["Kmatch_Num"] = 0
    KmerComp_Results_V2["Ancestral_Seq"]["Kmatch_NumEvaluated"] = 0

                                     
    # Now process info for the paralogs (All sequence that are not the ancestral seq)
    for p_Name, SeqInfoDict in i_dictOf_Seqs_ToCompare.items():
        
        if p_Name != 'Ancestral_Seq':
            p_Kmatch_NP = KmerComp_Results_V2[p_Name]["Kmatch_PADDED"]
            
            p_Kmatch_RERegion = p_Kmatch_NP[i_Pad_Size - i_Pad_Shift : - (i_Pad_Size - i_Pad_Shift*2)]
            
            p_Kmatch_RERegion_AncDiffPos = p_Kmatch_RERegion[RERegion_DiffToAnc_Mask] # Note the subsetted "RERegion" specific NP array should be shifted over by 5 idxs. (?)
    
            p_Kmatch_RERegion_MaskedByAncDiff = np.where(RERegion_SameToAnc_Mask, 0, p_Kmatch_RERegion)
    
            KmerComp_Results_V2[p_Name]["Kmatch_RERegion"] = AncSeq_Kmatch_RERegion
            KmerComp_Results_V2[p_Name]["Kmatch_RERegion_DiffToAncPos"] = p_Kmatch_RERegion_AncDiffPos
            KmerComp_Results_V2[p_Name]["Kmatch_RERegion_MaskedDiffKmers"] = p_Kmatch_RERegion_MaskedByAncDiff

            # print(p_Kmatch_RERegion_AncDiffPos)
            if p_Kmatch_RERegion_AncDiffPos.shape[0] > 0:
                match_Score_AncDiffPos = p_Kmatch_RERegion_AncDiffPos.sum() / p_Kmatch_RERegion_AncDiffPos.shape[0]
            
            else: 
                match_Score_AncDiffPos = 0
                
            #print(f"{p_Name} has a k-mer match of {match_Score_AncDiffPos} to the recomb event: {RE_EventID}")
            
            KmerComp_Results_V2[p_Name]["Kmatch_Score"] = match_Score_AncDiffPos
            KmerComp_Results_V2[p_Name]["Kmatch_Num"] = p_Kmatch_RERegion_AncDiffPos.sum()
            KmerComp_Results_V2[p_Name]["Kmatch_NumEvaluated"] = p_Kmatch_RERegion_AncDiffPos.shape[0]




    RecombVsAncSeq_Dict = {}


    i_RecombSeq_Hashes = hash_kmers_ToSet(RecombSeq_PADDED_Kmers  )  
    i_AncSeq_Hashes = hash_kmers_ToSet(i_AncestralSeq_FromSNPsInRv_PADDED_Kmers  )  

    Hashes_UnqToRecombVsAncestral = i_RecombSeq_Hashes - i_AncSeq_Hashes
    Kmers_UnqToRecombVsAncestral = set(RecombSeq_PADDED_Kmers) - set(i_AncestralSeq_FromSNPsInRv_PADDED_Kmers)

    # print(len(RecombSeq_PADDED_Kmers),
    #     len(i_AncestralSeq_FromSNPsInRv_PADDED_Kmers),
    #     len(Hashes_UnqToRecombVsAncestral),
    #     len(Kmers_UnqToRecombVsAncestral))  

    RecombVsAncSeq_Dict["Hashes_UnqToRecombVsAncestral"] = Hashes_UnqToRecombVsAncestral
    RecombVsAncSeq_Dict["Kmers_UnqToRecombVsAncestral"] = Kmers_UnqToRecombVsAncestral


    
    return KmerComp_Results_V2, i_dictOf_Seqs_ToCompare






















def EventMapping_Part1_compare_all_events_to_HmMap_Paralogs(
    i_GRE_DF: pd.DataFrame,
    i_HmPair_DF: pd.DataFrame,
    i_Gubbins_SNPs_EventOnly_Anno_DF: pd.DataFrame,
    i_HmMapAln_Var_DF: pd.DataFrame,
    i_H37Rv_Seq,
    *,
    genome_chr: str = "NC_000962.3",
    max_hm_regions: int = 15,
    round_ndigits: int = 4,
    verbose: bool = True,
) -> Tuple[pd.DataFrame, Dict[str, Dict[str, Any]]]:
    """
    Map each recombination event to overlapping homologous regions and compute:
      - SNP overlap fraction (P_SNPMatch)
      - k-mer match fraction (P_KmerMatch)
      - genomic distance between event and homolog (Dist_Event_To_Hm)

    Returns
    -------
    RE_MatchToHmRegions_DF : pd.DataFrame
        Concatenated per-event matches.
    dictOf_FullKmerAnalysis_PerEvent : dict
        { EventID : i_KmerCompResults } from compare_EventSeq_kmers_To_Paralogs_V2
    """
    # ---------- validation ----------
    def _require_cols(df: pd.DataFrame, need: set, df_name: str):
        missing = need - set(df.columns)
        if missing:
            raise ValueError(f"{df_name} is missing required columns: {sorted(missing)}")

    _require_cols(
        i_GRE_DF,
        {
            "EventID", "Overlap_Genes", "start_0based", "end_1based",
            "snp_count", "Parent_Node", "Child_Node", "Lineage"
        },
        "i_GRE_DF",
    )
    _require_cols(
        i_HmPair_DF,
        {
            "Target_Name", "Target_Start", "Target_End",
            "QueryOverlap_Genes", "Query_Start", "Query_End"
        },
        "i_HmPair_DF",
    )
    _require_cols(
        i_Gubbins_SNPs_EventOnly_Anno_DF,
        {"EventID"},
        "i_Gubbins_SNPs_EventOnly_Anno_DF",
    )
    _require_cols(
        i_HmMapAln_Var_DF,
        {"Target_Name", "Target_Start", "Target_End", "Query_Name", "Query_Start", "Query_End"},
        "i_HmMapAln_Var_DF",
    )

    if not isinstance(max_hm_regions, int) or max_hm_regions <= 0:
        raise ValueError("max_hm_regions must be a positive integer.")

    # ---------- core ----------
    listOf_Event_To_HmRegion_DFs = []
    dictOf_FullKmerAnalysis_PerEvent = {}

    for _, event_row in tqdm(i_GRE_DF.iterrows(), total=len(i_GRE_DF)):
        RE_EventID = event_row["EventID"]
        RE_Ovrlap_Genes = event_row["Overlap_Genes"]
        RE_Start = int(event_row["start_0based"])
        RE_End = int(event_row["end_1based"])
        Total_RE_SNP_Count = int(event_row["snp_count"])

        RE_Coords = f"{genome_chr}:{RE_Start}-{RE_End}"

        # 1) homologous regions overlapping the event
        sub_DF_Overlap_Hm_Regions = bf.select(
            i_HmPair_DF, RE_Coords, cols=("Target_Name", "Target_Start", "Target_End")
        )

        if sub_DF_Overlap_Hm_Regions.empty:
            if verbose:
                print(f"No overlapping paralogous regions for {RE_EventID}, Genes = {RE_Ovrlap_Genes}")
            continue

        # 2) event SNPs
        Event_SNPs = i_Gubbins_SNPs_EventOnly_Anno_DF.query("EventID == @RE_EventID")
        if Event_SNPs.empty and verbose:
            print(f"Warning: no SNPs found for EventID={RE_EventID}")

        # 3) k-mer comparison against paralogs
        i_KmerCompResults, _ = compare_EventSeq_kmers_To_Paralogs_V2(
            event_row, i_Gubbins_SNPs_EventOnly_Anno_DF, i_HmPair_DF, i_H37Rv_Seq
        )
        dictOf_FullKmerAnalysis_PerEvent[RE_EventID] = i_KmerCompResults

        # 4) per homolog calculations
        rows = []
        for _, HmRegionRow in sub_DF_Overlap_Hm_Regions.head(max_hm_regions).iterrows():
            Hm_Tar_Genes = HmRegionRow["QueryOverlap_Genes"]
            Hm_Tar_Start = int(HmRegionRow["Query_Start"])
            Hm_Tar_End = int(HmRegionRow["Query_End"])

            Hm_Target_Coords = f"{genome_chr}:{Hm_Tar_Start}-{Hm_Tar_End}"
            Hm_Tar_ID = f"{Hm_Tar_Genes}-H37Rv-{Hm_Tar_Start}-{Hm_Tar_End}"

            km = i_KmerCompResults.get(Hm_Tar_ID)
            if km is None:
                if verbose:
                    print(f"Warning: no k-mer result for {RE_EventID} vs {Hm_Tar_ID}; skipping.")
                continue

            Kmatch_Score = km.get("Kmatch_Score", np.nan)
            NumKmers_Match = km.get("Kmatch_Num", np.nan)
            NumKmers_Total = km.get("Kmatch_NumEvaluated", np.nan)

            # SNP overlap
            All_HmVar_Ovrlap_RE_DF = bf.select(
                i_HmMapAln_Var_DF, RE_Coords, cols=("Target_Name", "Target_Start", "Target_End")
            )
            Event_HmRegion_Combo_Var_DF = bf.select(
                All_HmVar_Ovrlap_RE_DF, Hm_Target_Coords, cols=("Query_Name", "Query_Start", "Query_End")
            )

            Num_SNPs_Ovrlap = check_snp_overlap_V2(Event_SNPs, Event_HmRegion_Combo_Var_DF)

            P_SNP_Match = Num_SNPs_Ovrlap / Total_RE_SNP_Count if Total_RE_SNP_Count > 0 else np.nan

            middle_paralog = (Hm_Tar_Start + Hm_Tar_End) / 2.0
            middle_event = (RE_Start + RE_End) / 2.0
            Dist_To_Paralog = Rv_dist(middle_paralog, middle_event)

            rows.append(
                (
                    RE_EventID, RE_Start, RE_End,
                    RE_Ovrlap_Genes, Hm_Tar_ID,
                    Hm_Tar_Genes, Hm_Tar_Start, Hm_Tar_End,
                    Num_SNPs_Ovrlap, Total_RE_SNP_Count, P_SNP_Match,
                    NumKmers_Match, NumKmers_Total, Kmatch_Score, Dist_To_Paralog
                )
            )

        if not rows:
            continue

        RE_To_HmRegion_Match_DF = pd.DataFrame(rows, columns=[
            "EventID","start_0based","end_1based","Overlap_Genes","HomologTargetID",
            "Hm_Overlap_Genes","Target_Start","Target_End",
            "NumMatch_SNPs","Total_SNPs","P_SNPMatch",
            "NumMatch_Kmers","Total_Kmers","P_KmerMatch","Dist_Event_To_Hm"
        ]).sort_values(["P_SNPMatch", "P_KmerMatch"], ascending=False)

        listOf_Event_To_HmRegion_DFs.append(RE_To_HmRegion_Match_DF)

    # ---------- concat + post ----------
    if not listOf_Event_To_HmRegion_DFs:
        empty_cols = [
            "EventID","start_0based","end_1based","Overlap_Genes","HomologTargetID",
            "Hm_Overlap_Genes","Target_Start","Target_End",
            "NumMatch_SNPs","Total_SNPs","P_SNPMatch",
            "NumMatch_Kmers","Total_Kmers","P_KmerMatch","Dist_Event_To_Hm"
        ]
        return pd.DataFrame(columns=empty_cols), dictOf_FullKmerAnalysis_PerEvent

    RE_MatchToHmRegions_DF = pd.concat(listOf_Event_To_HmRegion_DFs, ignore_index=True)

    for col in ("P_SNPMatch", "P_KmerMatch"):
        if col in RE_MatchToHmRegions_DF.columns:
            RE_MatchToHmRegions_DF[col] = RE_MatchToHmRegions_DF[col].round(round_ndigits)

    for col in ("Overlap_Genes", "Hm_Overlap_Genes"):
        if col in RE_MatchToHmRegions_DF.columns:
            RE_MatchToHmRegions_DF[col] = RE_MatchToHmRegions_DF[col].fillna("None")

    return RE_MatchToHmRegions_DF, dictOf_FullKmerAnalysis_PerEvent






def EventMapping_Part2_AnnotateEvents_ByBestMatch(
    i_GRE_DF: pd.DataFrame,
    i_RE_MatchToHmRegions_DF: pd.DataFrame,
    i_Gubbins_SNPs_EventOnly_Anno_DF: pd.DataFrame,
    i_Mtb_HM_Var_SNPs_DF: pd.DataFrame,
    *,
    genome_chr: str = "NC_000962.3",
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Annotate GRE events with paralog match summaries.

    Now builds HM_Var_SNPs_TrimUnq_DF internally from i_Mtb_HM_Var_SNPs_DF via:
        TarCol = ['Query_Name','Query_Start','Query_End','Ref','Alt','SNP']
        HM_Var_SNPs_TrimUnq_DF = i_Mtb_HM_Var_SNPs_DF[TarCol].drop_duplicates()
    """
    # --------- validation ----------
    def _require_cols(df: pd.DataFrame, need: set, name: str):
        missing = need - set(df.columns)
        if missing:
            raise ValueError(f"{name} is missing required columns: {sorted(missing)}")

    _require_cols(i_GRE_DF, {"EventID", "start_0based", "end_1based", "Child_Node"}, "i_GRE_DF")
    _require_cols(
        i_RE_MatchToHmRegions_DF,
        {
            "EventID","start_0based","end_1based","Overlap_Genes","HomologTargetID",
            "Hm_Overlap_Genes","Target_Start","Target_End",
            "NumMatch_SNPs","Total_SNPs","P_SNPMatch",
            "NumMatch_Kmers","Total_Kmers","P_KmerMatch","Dist_Event_To_Hm"
        },
        "i_RE_MatchToHmRegions_DF",
    )
    _require_cols(
        i_Gubbins_SNPs_EventOnly_Anno_DF,
        {"EventID", "Child_Call", "Pos_1based"},
        "i_Gubbins_SNPs_EventOnly_Anno_DF",
    )
    _require_cols(
        i_Mtb_HM_Var_SNPs_DF,
        {"Query_Name", "Query_Start", "Query_End", "Ref", "Alt", "SNP"},
        "i_Mtb_HM_Var_SNPs_DF",
    )

    # Build trimmed + unique HM SNPs DF
    TarCol = ['Query_Name', 'Query_Start', 'Query_End', 'Ref', 'Alt', 'SNP']
    HM_Var_SNPs_TrimUnq_DF = i_Mtb_HM_Var_SNPs_DF[TarCol].drop_duplicates()

    out_rows = []

    for _, event_row in tqdm(i_GRE_DF.iterrows(), total=len(i_GRE_DF)):
        RE_EventID = event_row["EventID"]
        RE_Start = int(event_row["start_0based"])
        RE_End = int(event_row["end_1based"])
        RE_Child_Node = str(event_row["Child_Node"])
        RE_Coords = f"{genome_chr}:{RE_Start}-{RE_End}"

        row = dict(event_row)
        row["IsTermNode"] = not RE_Child_Node.startswith("Node")

        # Event SNPs
        Event_SNPs = i_Gubbins_SNPs_EventOnly_Anno_DF.query("EventID == @RE_EventID").copy()
        if Event_SNPs.empty:
            if verbose:
                print(f"Warning: no SNPs for EventID={RE_EventID}")
            row["Freq_SNP_FoundInAnyPR"] = 0.0
        else:
            Event_SNPs["Pos_0based"] = Event_SNPs["Pos_1based"] - 1
            Gub_cols = ("Child_Call", "Pos_0based", "Pos_1based")
            Paf_cols = ("Alt", "Query_Start", "Query_End")
            Event_SNPs_CtByPRSNP_DF = bf.count_overlaps(
                Event_SNPs, HM_Var_SNPs_TrimUnq_DF, cols1=Gub_cols, cols2=Paf_cols
            ).rename(columns={"count": "NumPRSNP"})
            freq = (Event_SNPs_CtByPRSNP_DF.query("NumPRSNP > 0").shape[0] / Event_SNPs.shape[0]) if Event_SNPs.shape[0] > 0 else 0.0
            row["Freq_SNP_FoundInAnyPR"] = freq

        # All homolog matches for this event
        RE_To_HmRegion_Match_DF = i_RE_MatchToHmRegions_DF.query("EventID == @RE_EventID")
        row["NumHmTargets"] = int(RE_To_HmRegion_Match_DF.shape[0])

        if RE_To_HmRegion_Match_DF.shape[0] > 0:
            Max_KmerJC_ToHm = float(RE_To_HmRegion_Match_DF["P_KmerMatch"].max())
            Max_SNPmatch_ToHm = float(RE_To_HmRegion_Match_DF["P_SNPMatch"].max())
        else:
            Max_KmerJC_ToHm = 0.0
            Max_SNPmatch_ToHm = 0.0

        # K-mer top matches
        if RE_To_HmRegion_Match_DF.shape[0] > 0:
            RE_To_HmRegion_MaxKmerMatch_Only_DF = RE_To_HmRegion_Match_DF.query("P_KmerMatch == @Max_KmerJC_ToHm")
            row["NumHm_AtMaxKmerMatch"] = int(RE_To_HmRegion_MaxKmerMatch_Only_DF.shape[0])
            row["Top_KmerMatch_HomologTarIDs"]  = "-".join(RE_To_HmRegion_MaxKmerMatch_Only_DF["HomologTargetID"].astype(str).tolist())
            row["Top_KmerMatch_HomologGeneIDs"] = "-".join(RE_To_HmRegion_MaxKmerMatch_Only_DF["Hm_Overlap_Genes"].fillna("None").astype(str).tolist())
            row["Max_KmerMatch_ToHm"] = Max_KmerJC_ToHm
            row["DistToHm_TopKmatch"] = float(RE_To_HmRegion_MaxKmerMatch_Only_DF["Dist_Event_To_Hm"].values[0]) if RE_To_HmRegion_MaxKmerMatch_Only_DF.shape[0] == 1 else float("nan")
            row["TotalKmers_Eval"] = int(RE_To_HmRegion_Match_DF["Total_Kmers"].iloc[0])
        else:
            row["NumHm_AtMaxKmerMatch"] = 0
            row["Top_KmerMatch_HomologTarIDs"] = ""
            row["Top_KmerMatch_HomologGeneIDs"] = ""
            row["Max_KmerMatch_ToHm"] = float("nan")
            row["DistToHm_TopKmatch"] = float("nan")
            row["TotalKmers_Eval"] = 0

        # SNP top matches
        if RE_To_HmRegion_Match_DF.shape[0] > 0:
            RE_To_HmRegion_MaxSNPMatch_Only_DF = RE_To_HmRegion_Match_DF.query("P_SNPMatch == @Max_SNPmatch_ToHm")
            row["NumHm_AtMaxSNPMatch"] = int(RE_To_HmRegion_MaxSNPMatch_Only_DF.shape[0])
            row["Top_SNPMatch_HomologTarIDs"]  = "-".join(RE_To_HmRegion_MaxSNPMatch_Only_DF["HomologTargetID"].astype(str).tolist())
            row["Top_SNPMatch_HomologGeneIDs"] = "-".join(RE_To_HmRegion_MaxSNPMatch_Only_DF["Hm_Overlap_Genes"].fillna("None").astype(str).tolist())
            row["Max_SNPMatch_ToHm"] = Max_SNPmatch_ToHm
            row["DistToHm_TopSNPmatch"] = float(RE_To_HmRegion_MaxSNPMatch_Only_DF["Dist_Event_To_Hm"].values[0]) if RE_To_HmRegion_MaxSNPMatch_Only_DF.shape[0] == 1 else float("nan")
        else:
            row["NumHm_AtMaxSNPMatch"] = 0
            row["Top_SNPMatch_HomologTarIDs"] = ""
            row["Top_SNPMatch_HomologGeneIDs"] = ""
            row["Max_SNPMatch_ToHm"] = float("nan")
            row["DistToHm_TopSNPmatch"] = float("nan")

        out_rows.append(row)

    GRE_AnnoByMatch_DF = pd.DataFrame(out_rows)

    GRE_AnnoByMatch_DF["Max_KmerMatch_ToHm"] = GRE_AnnoByMatch_DF["Max_KmerMatch_ToHm"].fillna(0)

    return GRE_AnnoByMatch_DF




def RunEventMapping_to_HmMapParalogs_V1(
    i_GRE_DF: pd.DataFrame,
    i_HmPair_DF: pd.DataFrame,
    i_Gubbins_SNPs_EventOnly_Anno_DF: pd.DataFrame,
    i_Mtb_HM_Var_SNPs_DF: pd.DataFrame,
    i_H37Rv_Seq,
    *,
    genome_chr: str = "NC_000962.3",
    max_hm_regions: int = 15,
    round_ndigits: int = 4,
    verbose: bool = True,
) -> Tuple[pd.DataFrame, Dict[str, Dict[str, Any]], pd.DataFrame]:
    """
    Wrapper to:
      1) Match recombination events to homologous regions (SNP + kmer overlap).
      2) Annotate recombination events with paralog comparison info.

    Returns
    -------
    i_GRE_AnnoByMatch_DF : pd.DataFrame
        GRE_DF annotated with paralog comparison information.
    i_GRE_MatchToHmRegions_DF : pd.DataFrame
        Matches between GRE events and homologous regions.
    i_FullKmerAnalysis_PerEvent : dict
        Kmer comparison results for each event.
    """
    # Step 1: run match_events_to_homologs
    i_GRE_MatchToHmRegions_DF, i_FullKmerAnalysis_PerEvent = EventMapping_Part1_compare_all_events_to_HmMap_Paralogs(
        i_GRE_DF,
        i_HmPair_DF,
        i_Gubbins_SNPs_EventOnly_Anno_DF,
        i_Mtb_HM_Var_SNPs_DF,
        i_H37Rv_Seq,
        genome_chr=genome_chr,
        max_hm_regions=max_hm_regions,
        round_ndigits=round_ndigits,
        verbose=verbose,
    )

    # Step 2: run annotate_events_with_paralog_matches
    i_GRE_AnnoByMatch_DF = EventMapping_Part2_AnnotateEvents_ByBestMatch(
        i_GRE_DF,
        i_GRE_MatchToHmRegions_DF,
        i_Gubbins_SNPs_EventOnly_Anno_DF,
        i_Mtb_HM_Var_SNPs_DF,
        genome_chr=genome_chr,
        verbose=verbose,
    )

    return i_GRE_AnnoByMatch_DF, i_GRE_MatchToHmRegions_DF, i_FullKmerAnalysis_PerEvent










def getEventKmerMatchInfo_TopParalogMatches(i_Events_DF_Trim_DF,
                                            i_Events_MatchToEachParalog_DF):
    """
    For each mapped GCE event, extract the homologous-region matches
    that have the maximum k-mer match score.

    Parameters
    ----------
    i_Events_DF_Trim_DF : pd.DataFrame
        Must contain columns ["EventID", "Max_KmerMatch_ToHm"].
    i_Events_MatchToHm_DF : pd.DataFrame
        Must contain columns including ["EventID", "P_KmerMatch"].

    Returns
    -------
    pd.DataFrame
        mGCEToParalogs_TopKmerMatches_Info_DF
    """

    listOfDFs = []

    for _, row in i_Events_DF_Trim_DF.iterrows():
        event_id = row["EventID"]
        max_match = row["Max_KmerMatch_ToHm"]

        # subset rows matching this event ID *and* max kmer match score
        subset_df = i_Events_MatchToEachParalog_DF.query(
            f"P_KmerMatch == {max_match} and EventID == '{event_id}'"
        ).drop(["NumMatch_SNPs", "P_SNPMatch"], axis=1, errors="ignore")

        listOfDFs.append(subset_df)

    # concatenate all
    i_EventToTopParalogMatches_df = pd.concat(listOfDFs, ignore_index=True)

    # assign ref genome name
    i_EventToTopParalogMatches_df["seqname"] = "NC_000962.3"

    return i_EventToTopParalogMatches_df




def distribute_mGCE_donor_acceptor_counts_UnqParalogRegions(
    in_HM_MergedRegions_Anno_DF: pd.DataFrame,
    in_Mapped_GCE_DF: pd.DataFrame,
    in_Putative_GCE_DF: pd.DataFrame,
    i_HM_Var_SNPs_DF: pd.DataFrame,
    i_mGCEToParalogs_TopKmerMatches_Info_DF: pd.DataFrame,
):
    """
    Distribute donor/acceptor counts across homologous regions and compute
    per-region event/SNP summaries.

    Parameters
    ----------
    in_HM_MergedRegions_Anno_DF : pd.DataFrame
        HM merged regions annotation with at least:
        ["HmRegionID", "Chr", "Start", "End"].
    in_Mapped_GCE_DF : pd.DataFrame
        Mapped GCE events, with at least:
        ["EventID", "seqname", "start_0based", "end_1based"].
    in_Putative_GCE_DF : pd.DataFrame
        All putative GC events with at least:
        ["seqname", "start_0based", "end_1based"].
    i_HM_Var_SNPs_DF : pd.DataFrame
        HM-region SNPs with at least:
        ["Query_Name", "Query_Start", "Query_End", "Ref", "Alt", "SNP"].
    i_mGCEToParalogs_TopKmerMatches_Info_DF : pd.DataFrame
        Top k-mer GCE→paralog matches with at least:
        ["EventID", "HomologTargetID"].

    Returns
    -------
    HmRegions_DonorAndEventCt_DF : pd.DataFrame
        Copy of in_HM_MergedRegions_Anno_DF with additional summary columns.
    Overall_Donor_Count : float
        Total donor-normalized count across all regions.
    """

    # Build trimmed + unique HM SNPs DF
#    UnqSNPs_TarCol = ['Query_Name', 'Query_Start', 'Query_End', 'Ref', 'Alt', 'SNP']
    UnqSNPs_TarCol = ['Target_Name', 'Target_Start', 'Target_End', 'Ref', 'Alt', 'SNP']

    i_HM_Var_SNPs_TrimUnq_DF = i_HM_Var_SNPs_DF[UnqSNPs_TarCol].drop_duplicates()

    unique_HmRegionID = in_HM_MergedRegions_Anno_DF["HmRegionID"].unique()

    HmRegion_Donor_NormCt_dict    = {key: 0.0 for key in unique_HmRegionID}
    HmRegion_Acceptor_NormCt_dict = {key: 0.0 for key in unique_HmRegionID}

    Overall_Donor_Count = 0.0

    # Loop over mapped GCEs
    for _, row in in_Mapped_GCE_DF.iterrows():
        event_ID = row["EventID"]

        # Define event range
        event_Start = int(row["start_0based"])
        event_End   = int(row["end_1based"])
        event_Range = f"NC_000962.3:{event_Start}-{event_End}"

        # Step 1: Get the list of all BEST Donor matches for each event
        i_EventToDonor_TopMatches_DF = (
            i_mGCEToParalogs_TopKmerMatches_Info_DF
            .query("EventID == @event_ID")
        )

        Num_PossibleDonors = i_EventToDonor_TopMatches_DF.shape[0]
        if Num_PossibleDonors == 0:
            # No donor match for this event, skip donor/acceptor distribution
            continue

        weight = 1.0 / Num_PossibleDonors

        # Step 2: Add putative DONOR counts to each homologous region
        for _, match_row in i_EventToDonor_TopMatches_DF.iterrows():
            HmTargetID = match_row["HomologTargetID"]

            PutDonor_Start = HmTargetID.split("-")[-2]
            PutDonor_End   = HmTargetID.split("-")[-1]
            PutDonor_Range = f"NC_000962.3:{PutDonor_Start}-{PutDonor_End}"

            i_Donor_HmRegion = bf.select(
                in_HM_MergedRegions_Anno_DF,
                PutDonor_Range,
                cols=["Chr", "Start", "End"],
            )

            if i_Donor_HmRegion.shape[0] == 0:
                print("WARNING: No donor HM region for:", PutDonor_Range)
                continue

            if i_Donor_HmRegion.shape[0] > 1:
                print("BIGGGG ERRROR! Multiple donor HM regions for:", PutDonor_Range)

            i_Donor_HmRegionID = i_Donor_HmRegion["HmRegionID"].iloc[0]
            HmRegion_Donor_NormCt_dict[i_Donor_HmRegionID] += weight
            Overall_Donor_Count += weight

        # Step 3: Add putative ACCEPTOR counts to each homologous region
        # (acceptor is the event location itself; original code repeats per donor)
        for _, _ in i_EventToDonor_TopMatches_DF.iterrows():
            i_Acceptor_HmRegion = bf.select(
                in_HM_MergedRegions_Anno_DF,
                event_Range,
                cols=["Chr", "Start", "End"],
            )

            if i_Acceptor_HmRegion.shape[0] == 0:
                print("WARNING: No acceptor HM region for:", event_Range)
                continue

            if i_Acceptor_HmRegion.shape[0] > 1:
                print("BIGGGG ERRROR for ACCEPTOR COUNTING!", event_Range)

            i_Acceptor_HmRegionID = i_Acceptor_HmRegion["HmRegionID"].iloc[0]
            HmRegion_Acceptor_NormCt_dict[i_Acceptor_HmRegionID] += weight

    print("Overall_Donor_Count:", Overall_Donor_Count)

    # Start from a copy of the HM region annotations
    HmRegions_DonorAndEventCt_DF = in_HM_MergedRegions_Anno_DF.copy()

    # Add normalized donor counts
    HmRegions_DonorAndEventCt_DF["Norm_Donor_Count"] = (
        HmRegions_DonorAndEventCt_DF["HmRegionID"].map(HmRegion_Donor_NormCt_dict)
    )
    # If you want acceptor counts in the table as well, uncomment:
    # HmRegions_DonorAndEventCt_DF["Norm_Acceptor_Count"] = (
    #     HmRegions_DonorAndEventCt_DF["HmRegionID"].map(HmRegion_Acceptor_NormCt_dict)
    # )

    RE_CoordCols       = ("seqname", "start_0based", "end_1based")
    HmRegion_CoordCols = ("Chr", "Start", "End")

    # Count overlapping MAPPED GC Events per region
    HmRegions_DonorAndEventCt_DF = bf.count_overlaps(
        HmRegions_DonorAndEventCt_DF,
        in_Mapped_GCE_DF,
        cols1=HmRegion_CoordCols,
        cols2=RE_CoordCols,
    ).rename(columns={"count": "N_Events_Mapped"})

    # Count overlapping PUTATIVE GC Events per region
    HmRegions_DonorAndEventCt_DF = bf.count_overlaps(
        HmRegions_DonorAndEventCt_DF,
        in_Putative_GCE_DF,
        cols1=HmRegion_CoordCols,
        cols2=RE_CoordCols,
    ).rename(columns={"count": "N_Events_Putative"})

    # Skew and normalization metrics
    HmRegions_DonorAndEventCt_DF["EventSkew"] = (
        HmRegions_DonorAndEventCt_DF["N_Events_Mapped"]
        - HmRegions_DonorAndEventCt_DF["Norm_Donor_Count"]
    )

    HmRegions_DonorAndEventCt_DF["DonAccCt"] = (
        HmRegions_DonorAndEventCt_DF["Norm_Donor_Count"]
        + HmRegions_DonorAndEventCt_DF["N_Events_Mapped"]
    )

    # Avoid divide-by-zero with replace({0: pd.NA})
    HmRegions_DonorAndEventCt_DF["EventSkew_Norm"] = (
        HmRegions_DonorAndEventCt_DF["EventSkew"]
        / HmRegions_DonorAndEventCt_DF["DonAccCt"].replace({0: pd.NA})
    )

    HmRegions_DonorAndEventCt_DF["EventSkew_Norm"] = HmRegions_DonorAndEventCt_DF["EventSkew_Norm"].fillna(0)

    
    HmRegions_DonorAndEventCt_DF["GCE_FractionMapped"] = (
        HmRegions_DonorAndEventCt_DF["N_Events_Mapped"]
        / HmRegions_DonorAndEventCt_DF["N_Events_Putative"].replace({0: pd.NA})
    )

    HmRegions_DonorAndEventCt_DF["GCE_FractionMapped"] = HmRegions_DonorAndEventCt_DF["GCE_FractionMapped"].fillna(0)
    
    # Add counts of ALL PR-SNPs per HmRegion
    HmRegions_DonorAndEventCt_DF = bf.count_overlaps(
        HmRegions_DonorAndEventCt_DF,
        i_HM_Var_SNPs_DF,
        cols1=("Chr", "Start", "End"),
        cols2=("Target_Name", "Target_Start", "Target_End"),
    ).rename(columns={"count": "Num_PR_SNPs_All"})

    # Add counts of UNIQUE PR-SNPs per HmRegion
    HmRegions_DonorAndEventCt_DF = bf.count_overlaps(
        HmRegions_DonorAndEventCt_DF,
        i_HM_Var_SNPs_TrimUnq_DF,
        cols1=("Chr", "Start", "End"),
        cols2=("Target_Name", "Target_Start", "Target_End"),
    ).rename(columns={"count": "Num_PR_SNPs_Unq"})

    return HmRegions_DonorAndEventCt_DF, Overall_Donor_Count




def compute_HmPair_EventToDonor_counts(
    in_HmPair_DF: pd.DataFrame,
    in_Mapped_GCE_DF: pd.DataFrame,
    in_mGCEToParalogs_TopKmerMatches_Info_DF: pd.DataFrame,
):
    """
    Compute donor assignment counts for each homologous Query–Target region pair.

    Parameters
    ----------
    in_HmPair_DF : pd.DataFrame
        Homology pair alignments. Must contain:
        ["Query_Start", "Query_End", "Target_Start", "Target_End",
         "QueryCoords", "TargetCoords", "QueryToTarget_ID"]
    in_Mapped_GCE_DF : pd.DataFrame
        Mapped GCE events with:
        ["EventID", "seqname", "start_0based", "end_1based"]
    in_mGCEToParalogs_TopKmerMatches_Info_DF : pd.DataFrame
        Top k-mer paralog matches per event. Must contain:
        ["EventID", "HomologTargetID"]

    Returns
    -------
    HmPair_EventToDonor_Ct_Dict : dict
        { QueryToTarget_ID → donor-normalized count }
    Overall_Donor_Count : float
        Total distributed weight across all assignments
    """

    HmPair_EventToDonor_Ct_Dict = {
        key: 0.0 for key in in_HmPair_DF["QueryToTarget_ID"].unique()
    }

    Overall_Donor_Count = 0.0

    # Iterate over all true mapped GC events
    for _, row in in_Mapped_GCE_DF.iterrows():

        event_ID = row["EventID"]
        event_Chr = "NC_000962.3"
        event_Start = int(row["start_0based"])
        event_End   = int(row["end_1based"])
        event_Range = f"{event_Chr}:{event_Start}-{event_End}"

        # Get all top donor matches for this event
        i_EventToDonor_TopMatches_DF = (
            in_mGCEToParalogs_TopKmerMatches_Info_DF
            .query("EventID == @event_ID")
        )

        Num_PossibleDonors = i_EventToDonor_TopMatches_DF.shape[0]
        if Num_PossibleDonors == 0:
            continue

        weight = 1.0 / Num_PossibleDonors

        # Which hm-pair rows overlap the "Event"?
        i_HmPair_MatchEvent = bf.select(
            in_HmPair_DF,
            event_Range,
            cols=["Target_Name", "Target_Start", "Target_End"]
        )

        # Loop through K-mer donor matches
        for _, match_row in i_EventToDonor_TopMatches_DF.iterrows():

            HmTargetID = match_row["HomologTargetID"]
            PutDonor_Start = HmTargetID.split("-")[-2]
            PutDonor_End   = HmTargetID.split("-")[-1]
            PutDonor_Start = int(PutDonor_Start)
            PutDonor_End   = int(PutDonor_End)

            # Match donor coordinates to the hm-pair rows
            i_HmPair_MatchEventAndDonor = i_HmPair_MatchEvent.query(
                "Query_Start == @PutDonor_Start and Query_End == @PutDonor_End"
            )

            if i_HmPair_MatchEventAndDonor.shape[0] > 1:
                print("\n❗ Multiple Hm-Pair alignments found!")
                print("Event:", event_ID, event_Range)
                print("Donor:", f"{PutDonor_Start}-{PutDonor_End}")
                print(i_HmPair_MatchEventAndDonor["QueryToTarget_ID"].values)

            elif i_HmPair_MatchEventAndDonor.shape[0] == 0:
                print("\n❗ No matching Hm-Pair alignment found!")
                print("Event:", event_ID, event_Range)
                print("Donor:", f"{PutDonor_Start}-{PutDonor_End}")
                continue  # nothing to increment

            # Retrieve ID of this query–target pair
            pair_id = i_HmPair_MatchEventAndDonor["QueryToTarget_ID"].values[0]

            # Increment dictionary
            HmPair_EventToDonor_Ct_Dict[pair_id] += weight
            Overall_Donor_Count += weight

    return HmPair_EventToDonor_Ct_Dict, Overall_Donor_Count







def summarize_HmPair_MappedEvent_counts(
    in_HmPair_DF: pd.DataFrame,
    in_Mapped_GCE_DF: pd.DataFrame,
    in_Putative_GCE_DF: pd.DataFrame,
    in_mGCEToParalogs_TopKmerMatches_Info_DF: pd.DataFrame,
    bin_start: int = 0,
    bin_end: int = 2300000,
    bin_size: int = 50000,
) -> pd.DataFrame:
    """
    Add event count, distance binning, gene labels, normalization, and
    putative GCE overlap counts to the HmPair table.

    Parameters
    ----------
    in_HmPair_DF : pd.DataFrame
        HM pair alignments with at least:
        ["QueryToTarget_ID", "Dist_Middles",
         "QueryOverlap_Genes", "TargetOverlap_Genes",
         "Query_Name", "Query_Start", "Query_End"].
    in_Putative_GCE_DF : pd.DataFrame
        Putative GCE events with coords columns:
        ["seqname", "start_0based", "end_1based"].
    in_HmPair_EventToDonor_Ct_Dict : dict
        Mapping {QueryToTarget_ID -> donor-normalized event count}.
    bin_start, bin_end, bin_size : int
        Parameters for distance binning of 'Dist_Middles'.

    Returns
    -------
    HmPair_EventCt_DF : pd.DataFrame
        Copy of in_HmPair_DF with additional columns:
        - "EventsMapped"
        - "Dist_Bin"
        - "Dist_Bin_Mid"
        - "DonorGenes"
        - "AcceptorGenes"
        - "Normalized_EventsMapped"
        - "N_Ovrlap_pGCE"
    """


    # 1) Compute donor counts per HmPair
    i_HmPair_EventToDonor_Ct_Dict, i_Overall_Donor_Count = compute_HmPair_EventToDonor_counts(
        in_HmPair_DF = in_HmPair_DF,
        in_Mapped_GCE_DF = in_Mapped_GCE_DF,
        in_mGCEToParalogs_TopKmerMatches_Info_DF = in_mGCEToParalogs_TopKmerMatches_Info_DF,
    )

    # (Optional) log this if useful
    print("Overall_Donor_Counts inferred across all paralog-pairs:", i_Overall_Donor_Count)
    
    HmPair_EventCt_DF = in_HmPair_DF.copy()

    # Map donor-normalized counts
    HmPair_EventCt_DF["EventsMapped"] = (
        HmPair_EventCt_DF["QueryToTarget_ID"].map(i_HmPair_EventToDonor_Ct_Dict)
    )

    # Distance bins
    bins = range(bin_start, bin_end + bin_size, bin_size)
    HmPair_EventCt_DF["Dist_Bin"] = pd.cut(HmPair_EventCt_DF["Dist_Middles"], bins=bins)
    HmPair_EventCt_DF["Dist_Bin_Mid"] = HmPair_EventCt_DF["Dist_Bin"].apply(
        lambda x: x.mid if pd.notnull(x) else pd.NA
    )

    # Donor / acceptor gene labels
    HmPair_EventCt_DF["DonorGenes"] = HmPair_EventCt_DF["QueryOverlap_Genes"]
    HmPair_EventCt_DF["AcceptorGenes"] = HmPair_EventCt_DF["TargetOverlap_Genes"]

    # Normalize event counts (0–1 scale) and cast to float
    HmPair_EventCt_DF["EventsMapped"] = HmPair_EventCt_DF["EventsMapped"].astype(float)
    max_events = HmPair_EventCt_DF["EventsMapped"].max()
    HmPair_EventCt_DF["Normalized_EventsMapped"] = (
        HmPair_EventCt_DF["EventsMapped"] / max_events if pd.notnull(max_events) else pd.NA
    )

    # Add count of overlapping PUTATIVE GC Events
    RE_CoordCols = ("seqname", "start_0based", "end_1based")
    Hm_Query_Cols = ["Query_Name", "Query_Start", "Query_End"]
    Hm_Target_Cols = ["Target_Name", "Target_Start", "Target_End"]

    HmPair_EventCt_DF = bf.count_overlaps(
        HmPair_EventCt_DF,
        in_Putative_GCE_DF,
        cols1=Hm_Target_Cols,
        cols2=RE_CoordCols,
    ).rename(columns={"count": "N_Ovrlap_pGCE_ToTarget"})

    # Fill NA gene fields with "_"
    HmPair_EventCt_DF["AcceptorGenes"] = HmPair_EventCt_DF["AcceptorGenes"].fillna("_")
    HmPair_EventCt_DF["DonorGenes"] = HmPair_EventCt_DF["DonorGenes"].fillna("_")

    return HmPair_EventCt_DF, i_HmPair_EventToDonor_Ct_Dict




