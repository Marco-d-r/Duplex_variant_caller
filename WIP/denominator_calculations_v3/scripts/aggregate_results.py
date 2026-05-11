import json
import pandas as pd

# 1. Load the Masking Stats (the denominator "Space" stats)
with open(snakemake.input.mask_stats) as f:
    report = json.load(f)

# 2. Initialize counters for the BAM "Evidence" stats
callable_sum = 0
raw_obs_sum = 0
marginal_costs = {}

# 3. Aggregate all chromosome JSONs
for j_file in snakemake.input.chrom_jsons:
    with open(j_file) as f:
        c_data = json.load(f)
        
        # Sum up the main counts
        callable_sum += c_data.get("callable_bases_sum", 0)
        raw_obs_sum += c_data.get("total_observations", 0)
        
        # Sum up all marginal gain keys
        for key, val in c_data.items():
            if "marginal_gain" in key:
                marginal_costs[key] = marginal_costs.get(key, 0) + val

# 4. Finalize the report object
report.update({
    "sample": snakemake.wildcards.sample,
    "total_raw_observations": raw_obs_sum,
    "final_callable_bases": callable_sum,
    **marginal_costs
})

# 5. Save as CSV
df = pd.DataFrame([report])
df.to_csv(snakemake.output.csv, index=False)