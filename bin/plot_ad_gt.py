#!/usr/bin/env python3
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

def main():
    parser = argparse.ArgumentParser(description="Plot Ref vs Alt AD counts")
    parser.add_argument("--inputs", nargs='+', required=True, help="TSV file")
    parser.add_argument("--output", required=True, help="Output PNG filename")
    args = parser.parse_args()

    # Merge all input files
    df_list = []
    for f in args.inputs:
        try:
            # Create a temporary DataFrame to check if the file is empty
            temp_df = pd.read_csv(f, sep='\t')
            if not temp_df.empty:
                df_list.append(temp_df)
        except Exception as e:
            print(f"Warning: Could not read {f}: {e}", file=sys.stderr)

    if not df_list:
        print("No data found to plot.", file=sys.stderr)
        sys.exit(1)

    df = pd.concat(df_list, ignore_index=True)

    # Clean AD column: remove rows with missing AD
    df = df[df['AD'] != '.']
    
    # Split AD into Ref and Alt counts
    # Multi-allelic (e.g., 10,5,0) cases only use the first two (Ref and First Alt)
    # Split by comma and convert to numeric
    ad_split = df['AD'].str.split(',', expand=True)
    
    # If there are not at least two columns after split, exit with error
    if ad_split.shape[1] < 2:
        print("Error: AD column does not contain Ref,Alt format.", file=sys.stderr)
        sys.exit(1)

    df['Ref_Count'] = pd.to_numeric(ad_split[0], errors='coerce')
    df['Alt_Count'] = pd.to_numeric(ad_split[1], errors='coerce')
    
    # Remove rows with NaN in Ref_Count or Alt_Count
    df = df.dropna(subset=['Ref_Count', 'Alt_Count'])

    # Convert GT to categorical labels
    def classify_gt(gt_str):
        if pd.isna(gt_str):
            return "Unknown"
        
        # Phased and unphased genotypes are treated the same
        gt_norm = gt_str.replace('|', '/')
        
        if gt_norm == '0/0':
            return 'HomRef'
        elif gt_norm in ['0/1', '1/0']:
            return 'Het'
        elif gt_norm == '1/1':
            return 'HomAlt'
        elif './.' in gt_norm:
            return 'Missing'
        else:
            # Handle multi-allelic or complex genotypes
            return 'Other'

    # Convert GT to categorical labels
    df['GT_Label'] = df['GT'].apply(classify_gt)


    # Plot setting
    plt.figure(figsize=(6, 5))
    sns.set_style("whitegrid")

    # Scatter plot
    sns.scatterplot(
        data=df, 
        x='Ref_Count', 
        y='Alt_Count', 
        hue='GT_Label', 
        alpha=0.5,
        s=15,
        palette='viridis' # Color palette
    )

    plt.title('Allele Count: Ref vs Alt by Genotype')
    plt.xlabel('Reference Allele Count')
    plt.ylabel('Alternate Allele Count')
    plt.legend(title='Genotype', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # X=Y line
    max_val = max(df['Ref_Count'].max(), df['Alt_Count'].max())
    plt.plot([0, max_val], [0, max_val], ls="--", c=".3", label="1:1 Ratio")

    plt.tight_layout()
    plt.savefig(args.output, dpi=500)
    print(f"Plot saved to {args.output}")

if __name__ == "__main__":
    main()