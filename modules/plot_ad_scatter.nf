#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process PLOT_AD_SCATTER {
    tag "plot_ad_scatter"
    publishDir "${params.outdir}/genotype_qc", mode: 'copy'

    input:
    path ad_gt_tsv

    output:
    path "ad_scatter.png"
    path "ad_scatter.filtered.tsv"

    script:
    def max_points = params.ad_plot_max_points ?: 0
    """
    python - <<'PY'
    import pandas as pd
    import matplotlib.pyplot as plt

    # Read AD and GT data
    df = pd.read_csv("${ad_gt_tsv}", sep="\\t", header=None, names=["CHROM","POS","REF","ALT","SAMPLE","GT","AD"])
    
    # Split AD into REF and ALT counts
    ad = df["AD"].astype(str).str.split(",", expand=True)
    df["REF_COUNT"] = pd.to_numeric(ad[0], errors="coerce")
    df["ALT_COUNT"] = pd.to_numeric(ad[1], errors="coerce")

    # Remove missing data
    df = df.dropna(subset=["REF_COUNT","ALT_COUNT","GT"])

    # Simplify GT values and filter to common genotypes
    df["GT_SIMPLE"] = df["GT"].astype(str).str.replace("|","/", regex=False)
    df = df[df["GT_SIMPLE"].isin(["0/0","0/1","1/0","1/1","./."])]

    # Normalize heterozygous GT representation
    df.loc[df["GT_SIMPLE"].isin(["1/0"]), "GT_SIMPLE"] = "0/1"

    # Downsample if too many points
    max_points = int(${max_points})
    if max_points and len(df) > max_points: df = df.sample(n=max_points, random_state=1)
    
    # Save filtered data
    df.to_csv("ad_scatter.filtered.tsv", sep="\\t", index=False)
    
    # setting colors for each genotype
    colors = {"0/0":"#1f77b4", "0/1":"#ff7f0e", "1/1":"#2ca02c", "./.":"#7f7f7f"}

    # Plotting
    plt.figure(figsize=(7,7))
    for gt, sub in df.groupby("GT_SIMPLE"):
        plt.scatter(sub["REF_COUNT"], sub["ALT_COUNT"], s=15, alpha=0.5,
                    c=colors.get(gt, "#000000"), label=gt, linewidths=0)
    plt.xlabel("Reference allele read count (AD[0])")
    plt.ylabel("Alternate allele read count (AD[1])")
    plt.title("Allele depth scatter colored by genotype")
    plt.legend(title="GT", markerscale=2)
    plt.grid(True, alpha=0.2)
    plt.tight_layout()
    plt.savefig("ad_scatter.png", dpi=200)
    PY
    """
}
