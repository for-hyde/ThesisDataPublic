import pandas as pd
from pysmf.core.bamobject import BamObject
from pysmf.core.sparsematrix import SparseMatrix
from pysmf.utils.regions import expand_region, find_context_sites
    

def call_SMF(bam: BamObject, region: tuple, context: str="GC", span: int=250):
    """
    Call single-molecule methylation features in a specified genomic region.
    Args:
        bam (BamObject): An instance of BamObject for BAM file operations.
        region (tuple): A tuple in the format of (chromosome, start, end) specifying the genomic region.
        context (str): The methylation context to analyze (default is "GC").
        span (int): The number of base pairs to extend on either side of the region (default is 250).
    Returns:
        pd.DataFrame: A DataFrame where rows represent reads and columns represent genomic sites with methylation status.
                      Methylation status is encoded as 1 (methylated), -1 (unmethylated/converted), and 0 (not covered).
    """
    #Handle if region has directional info (Add for single site)
    if len(region) == 3:
        region = (*region, "+") #Default to forward strand

    #Expand region by span
    #Make a single call of function for this. 
    new_region = expand_region(region, span)

    # Find Reference contexts
    sites = find_context_sites(bam.ref_fa,new_region, context = context)

    #Tracking reads and methylation info
    read_ids = []
    methylation_info = {site:[] for site in sites}
    cl = len(context)   #len of context sequence
    
    #Define region to parse
    c = new_region[0]
    s, e = sites[0], sites[-1]+cl 
    i = 0

    for read in bam.fetch(c,s,e):
        #Check if read has already been processed.
        read_id = read.query_name
        if read_id in read_ids:
            continue
        read_ids.append(read_id)

        #Initialize 0 for all sites.
        for key in methylation_info.keys():
            methylation_info[key].append(0)

        #Get sequences
        seq = read.get_forward_sequence()
        
        for index, reference in read.get_aligned_pairs():
            #Check that reference is context site and aligns to sequence.
            if reference not in sites or not index:
                continue
            #Call reference with each identified site.
            #Sequence different from reference - Converted
            #This site can be improved to ensure conversions are canonical.
            re = int(index+cl)
            if seq[index:re] != context:
                methylation_info[reference][i] = 1
            else:
                methylation_info[reference][i] = -1

        #Processing of mates, only occurs if mate is mapped and neither of the pair have been added yet. 
        if read.mate_is_mapped:
            mate = bam.bam.mate(read)
            #Ensure reads are not duplicated. 
            m_read_id = mate.query_name
            if m_read_id in read_ids:
                continue
            read_ids.append(m_read_id)

            mate_seq = mate.get_forward_sequence()
            for index, reference in mate_seq.get_aligned_pairs():
                #Parse for sequences
                if reference not in sites or not index:
                    continue
                re = int(index+cl)
                #Sequence different from reference - Converted
                if mate_seq[index:re] != context and methylation_info[reference][i] != -1:
                    methylation_info[reference][i] = 1
                else:
                    methylation_info[reference][i] = -1

        i += 1  #Increment read index
    return pd.DataFrame(methylation_info)