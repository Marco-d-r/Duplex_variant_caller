import pybedtools
import json

def get_density_mask(vcf_path, window, max_vars, genome):
    vcf = pybedtools.BedTool(vcf_path)
    clusters = vcf.cluster(d=window)
    df = clusters.to_dataframe()
    
    if df.empty:
        return pybedtools.BedTool("", from_string=True)
    
    cluster_col = df.columns[-1]
    bad_ids = set(df[df[cluster_col].value_counts().reindex(df[cluster_col]).values > max_vars][cluster_col])
    
    bad_regions = [f"{f.chrom}\t{f.start}\t{f.end}" for f in clusters if int(f[-1]) in bad_ids]
    if not bad_regions:
        return pybedtools.BedTool("", from_string=True)
        
    return pybedtools.BedTool("\n".join(bad_regions), from_string=True).slop(b=window, genome=genome).merge()

# Inputs from Snakemake
wes = pybedtools.BedTool(snakemake.input.wes)
blacklist = pybedtools.BedTool(snakemake.input.blacklist)
vcf = snakemake.input.vcf
d_cfg = snakemake.config["density"]

# 1. Generate Masks (Ensure these functions return sorted BedTools)
m20 = get_density_mask(vcf, d_cfg["window_small"], d_cfg["max_vars_small"], d_cfg["genome_build"])
m250 = get_density_mask(vcf, d_cfg["window_large"], d_cfg["max_vars_large"], d_cfg["genome_build"])

# 2. Combine all exclusions
# Explicitly sort after concatenating to ensure merge() won't fail
all_masks = blacklist.cat(m20).cat(m250).sort().merge()

# 3. Calculate Final Target
# Subtract can also occasionally result in unsorted output depending on the BED version
clean_targets = wes.subtract(all_masks).sort()
clean_targets.saveas(snakemake.output.bed)

# 4. Record Stats (FIXED: Add .sort() before total_coverage)
# total_coverage calls merge(), so we must sort the intersection first
stats = {
    "original_wes_size": wes.total_coverage(),
    "total_masked_bases": all_masks.intersect(wes).sort().total_coverage(),
    "density_mask_20bp_size": m20.intersect(wes).sort().total_coverage(),
    "density_mask_250bp_size": m250.intersect(wes).sort().total_coverage(),
    "final_target_size": clean_targets.total_coverage()
}

with open(snakemake.output.stats, 'w') as f:
    json.dump(stats, f)