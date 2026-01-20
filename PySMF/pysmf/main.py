import pysmf
import os
import pandas as pd

bam_path = "/Volumes/sd17l002/p/SingleMolculeFootprinting/data/Outputs/Biscuit/ERR4165149_DE_.sorted.bam"
ref_path = "/Volumes/sd17l002/p/SingleMolculeFootprinting/data/ReferenceGenome/GRCm39.genome.fa"
motif_files_path = "/Volumes/sd17l002/p/SingleMolculeFootprinting/data/RawData/Promoter_data/"
motif_files = [motif_files_path+file for file in os.listdir(motif_files_path)]

bam = pysmf.BamObject(bam_path, ref_path)

df = pd.read_csv(motif_files[0], sep=" ")

regions=[(r["seqnames.x"], r["start.x"], r["end.x"], r["strand.x"])for idx, r in df.iterrows()]
max = regions[985]

smf = pysmf.call_SMF(bam, max, context="GC", span=250)