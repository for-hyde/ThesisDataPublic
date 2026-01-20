def expand_region(region, span):
    """
    Support function for expanding the region around a given site. Called in the call_SMF and call_SMF_average functions.
    Args:
        region (tuple): A tuple in the format of (chromosome, start, end) specifying the genomic region. 
        span (int): A integer describing the distance from the region to expand. 
    Returns:
        region (tuple): Same as initial region, but with expanded limits.
    """
    c,s_orig,e_orig,dir = region
    s = s_orig - span
    e = e_orig + span
    return (c,s,e,dir)

def find_context_sites(ref, region, context = "GC"):
    """
    Support function for finding all instances of a given context in a region. Called in the call_SMF and call_SMF_average functions.
    Args:
        ref (pysam.FastaFile): A pysam FastaFile object for reference sequence operations.
        region (tuple): A tuple in the format of (chromosome, start, end) specifying the genomic region. 
        context (str): The methylation context to analyze (default is "GC").
    Returns:
        sites (list): A list of all positions in the region that match the given context.
    """
    c,s,e,_ = region
    seq = ref.fetch(c,s,e)
    sites = []
    cl = len(context)
    for i in range(len(seq)-cl+1):
        if seq[i:i+cl] == context:
            sites.append(s+i)
    return sites

