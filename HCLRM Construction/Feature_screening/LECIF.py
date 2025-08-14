#!/usr/bin/env python3
import numpy as np
import pandas as pd
import pyBigWig

# ——— 1. read original gene pairs ———
df = pd.read_csv(
    "/data/yangchangbo/ZSQ/HCLRs/xgboost_species_TF_cor_phastCons100_30_HCLR.txt",
    sep="\t", dtype=str
)

# ——— 2. GTF:gene→(chr,start,end)  ———
gtf = "/data/yangchangbo/ZSQ/HCLRs/gencode.v48.annotation.gtf"
gene_coords = {}
with open(gtf) as f:
    for ln in f:
        if ln.startswith('#'): continue
        cols = ln.split('\t')
        if cols[2] != "gene": continue
        attrs = cols[8]
        # extract gene_name
        import re
        m = re.search(r'gene_name "([^"]+)"', attrs)
        if not m: continue
        name = m.group(1)
        chrom = cols[0]
        start = int(cols[3]) - 1  # convert to 0-based
        end = int(cols[4])
        gene_coords[name] = (chrom, start, end)

# ——— 3.  bigWig to get Mean_LECIF for genes ———
bw = pyBigWig.open("/data/yangchangbo/ZSQ/HCLRs/hg38.LECIFv1.1.bw")

def get_lecif(gene):
    # 1. check if gene has coordinates
    coord = gene_coords.get(gene)
    if coord is None:
        # gene not found, return nan
        return np.nan

    chrom, start, end = coord
    # 2. check if interval is valid
    if not (isinstance(chrom, str) and end > start and start >= 0):
        return np.nan

    # 3. try to get values from bigWig
    try:
        vals = bw.values(chrom, start, end, numpy=True)
    except RuntimeError:
        # any RuntimeError (e.g. chromosome name mismatch) is considered missing
        return np.nan

    # 4. if the whole interval is nan or the array is empty, return nan
    if vals is None or len(vals) == 0 or np.all(np.isnan(vals)):
        return np.nan

    # 5. return the mean value
    return np.nanmean(vals)

# ——— 4. calculate Mean_LECIF for target and source ———
df['target_LECIF'] = df['target'].apply(get_lecif)
df['source_LECIF'] = df['source'].apply(get_lecif)

# ——— 5. write new file ———
out = "/data/yangchangbo/ZSQ/HCLRs/xgboost_species_TF_cor_phastCons100_30_LECIF_HCLR.txt"
df.to_csv(out, sep='\t', index=False)
print("Wrote:", out)