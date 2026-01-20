class BamObject:
    """"
    A class to handle BAM file operations using pysam. Initialization 
    requires a BAM file and a reference FASTA file.
    """
    def __init__(self, bamfile, ref_fa):
        import pysam

        self.bam = pysam.AlignmentFile(bamfile, "rb")
        self.ref_fa = pysam.FastaFile(ref_fa)

    def fetch(self, chrom, start, end):
        """Fetch reads from BAM file in the specified region."""
        return self.bam.fetch(chrom, start, end)

    def count_reads(self, chrom, start, end):
        """Fetch the count of reads in the specified region."""
        return self.bam.count(chrom, start, end)
    
    def pileup(self, chrom, start, end, truncate=True):
        """Generate a pileup of reads in the specified region.
        If truncate is True, only include reads that fully overlap the region.
        """
        return self.bam.pileup(chrom, start, end, truncate=truncate)

    def fetch_reference_sequence(self, chrom, start, end):
        """Fetch reference sequence from the FASTA file in the specified region."""
        return self.ref_fa.fetch(chrom, start, end)

    def close(self):
        self.bam.close()
        self.ref_fa.close()

class BamObjectFactory:
    @staticmethod
    def create(bamfile, ref_fa):
        return BamObject(bamfile, ref_fa)
