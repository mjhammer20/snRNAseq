# Imports
import argparse
import scanpy as sc
import pandas as pd
from anndata import AnnData
import gc
import json

# Function to Extract Nearest Genes from GWAS Summary Stats
def extract_nearest_genes(gwas_summary_stats: str) -> pd.DataFrame:
    """
    Extract nearest genes from GWAS summary stats file

    Args:
        gwas_summary_stats: Path to GWAS summary stats file

    Returns:
        DataFrame with nearest genes information
    """
    gwas_df = pd.read_csv(gwas_summary_stats, sep='\t')
    # Extract relevant columns (e.g., SNP, nearest gene)
    nearest_genes = gwas_df['Nearest Gene'].unique().tolist()

    # Append GBA to the list of nearest genes
    nearest_genes.append('GBA')
    
    return nearest_genes

def extract_marker_genes(mmc_extended_results_fp: str, gene_mapping_df: pd.DataFrame) -> list[str]:
    """
    Extract marker genes from marker genes file

    Args:
        mmc_extended_results_fp: Path to marker genes JSON file
        gene_mapping_df: DataFrame containing gene mapping information

    Returns:
        List of marker genes mapped to unique gene names
    """
    # Load the JSON file containing the extended results
    with open(mmc_extended_results_fp, "r") as f:
        extended_results = json.load(f)
    
    # Extract marker genes from the extended results
    marker_genes = set().union(*extended_results["marker_genes"].values())

    # Map marker genes to unique gene names using the gene mapping DataFrame
    mapped_marker_genes = list(set(gene_mapping_df[gene_mapping_df["gene_id"].isin(marker_genes)]['gene_name']))

    return mapped_marker_genes

# Function to run feature selection on the AnnData object for each region
def select_features(
        adata: AnnData, 
        batch_key: str,  
        n_top_genes: int, 
        nearest_genes: list[str],
        markers: list[str]
        ) -> tuple[AnnData, pd.DataFrame, pd.DataFrame]:
    """
    Do feature selection and add PCA

    Args:
        adata: AnnData object to process
        batch_key: Key in adata.obs to use for batch information in HVG selection
        n_top_genes: Number of HVG genes to keep [3000]
        nearest_genes: List of nearest genes to retain in the dataset regardless of whether they are selected as HVGs
        markers: List of marker genes to retain in the dataset regardless of whether they are selected as HVGs
    Returns:
        adata: Processed AnnData object with HVGs and PCA
        full_features: DataFrame with metadata for all genes
        reduced_features: DataFrame with metadata for HVG, nearest, and marker genes
    """

    # Identify highly variable genes using the Pearson residuals method
    hvgs_full = sc.experimental.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        batch_key=batch_key,
        flavor="pearson_residuals",
        check_values=True,
        subset=False,
        inplace=False,
    )

    # Sort genes by how often they selected as hvg within each batch and break ties with median rank of residual variance across batches
    hvgs_full.sort_values(
        ["highly_variable_nbatches", "highly_variable_rank"],
        ascending=[False, True],
        na_position="last",
        inplace=True,
    )

    # Save the full set of features before filtering to HVGs, nearest genes, and markers
    full_features = adata.var.copy()

    # Filter to HVGs, nearest genes, and markers
    hvgs_full = hvgs_full.iloc[: n_top_genes].index.to_list()
    genes_to_keep = set(hvgs_full + nearest_genes + markers)
    adata = adata[:, adata.var.index.isin(genes_to_keep)].copy()

    # Save the filtered features
    reduced_features = adata.var.copy()

    return (adata, full_features, reduced_features)


def main(args: argparse.Namespace):

    # Extract arguments
    adata_input_prefix = args.adata_input_prefix
    batch_key = args.batch_key
    regions = args.regions.split()
    gwas_summary_stats = args.gwas_summary_stats
    mmc_results_prefix = args.mmc_results_prefix
    gene_mapping_file = args.gene_mapping_file
    n_top_genes = args.n_top_genes
    adata_output_prefix = args.adata_output_prefix
    output_all_genes_prefix = args.output_all_genes_prefix
    output_hvg_genes_prefix = args.output_hvg_genes_prefix

    # Set CPUs to use for parallel computing
    sc.settings.n_jobs = -1

    # Extract nearest genes from GWAS summary stats
    nearest_genes = extract_nearest_genes(gwas_summary_stats)
    print(f"Extracted {len(nearest_genes)} unique nearest genes from GWAS summary stats")

    # Load the gene mapping file
    gene_mapping_df = pd.read_csv(gene_mapping_file, sep='\t')

    # Iterate through regions and process each one
    for region in regions:
        print(f"Processing region: {region}")

        # Extract marker genes for the current region
        mmc_extended_results_fp = f"{mmc_results_prefix}_{region}.extended_results.json"
        markers = extract_marker_genes(mmc_extended_results_fp = mmc_extended_results_fp, gene_mapping_df = gene_mapping_df)
        print(f"Extracted {len(markers)} unique marker genes for region {region}")

        # Load adata
        adata_fp = f"{adata_input_prefix}_{region}.h5ad"
        print(f"Loading AnnData object for region {region} from {adata_fp}...")
        adata = sc.read_h5ad(adata_fp)

        # Select Features
        print(f"Selecting features for region {region}...")
        adata, full_features, reduced_features = select_features(adata = adata, batch_key = batch_key, n_top_genes = n_top_genes, nearest_genes = nearest_genes, markers = markers)

        # Save the feature selected adata
        adata_output_fp = f"{adata_output_prefix}_{region}.h5ad"
        print(f"Saving feature selected AnnData object for region {region} to {adata_output_fp}...")
        adata.write_h5ad(filename=adata_output_fp, compression="gzip")

        # Save the full and reduced feature metadata to CSV files
        output_all_genes_fp = f"{output_all_genes_prefix}_{region}.csv"
        print(f"Saving full feature metadata for region {region} to {output_all_genes_fp}...")
        full_features.to_csv(output_all_genes_fp, index=True)

        output_hvg_genes_fp = f"{output_hvg_genes_prefix}_{region}.csv"
        print(f"Saving reduced feature metadata for region {region} to {output_hvg_genes_fp}...")
        reduced_features.to_csv(output_hvg_genes_fp, index=True)

        # Clean up memory
        del adata
        gc.collect()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Normalize and identify features (HVG)"
    )
    parser.add_argument(
        "--adata-input-prefix",
        type=str,
        required=True,
        help="AnnData file prefix for a dataset",
    )
    parser.add_argument(
        "--regions",
        type=str,
        required=True,
        help="Space separated list of regions to process",
    )
    parser.add_argument(
        "--gwas-summary-stats",
        type=str,
        required=True,
        help="Path to GWAS summary stats file",
    )
    parser.add_argument(
        "--mmc-results-prefix",
        type=str,
        required=True,
        help="MMC results file prefix for marker gene extraction",
    )
    parser.add_argument(
        "--gene-mapping-file",
        type=str,
        required=True,
        help="Path to gene mapping file for mapping marker genes to unique gene names",
    )
    parser.add_argument(
        "--batch-key",
        type=str,
        required=True,
        help="Key in AnnData object for batch information",
    )
    parser.add_argument(
        "--n-top-genes",
        type=int,
        required=True,
        help="Number of HVG genes to keep [3000]",
    )
    parser.add_argument(
        "--adata-output-prefix",
        type=str,
        required=True,
        help="Output file prefix to save AnnData object to",
    )
    parser.add_argument(
        "--output-all-genes-prefix",
        type=str,
        required=True,
        help="Output file prefix to save feature metadata (full genes)",
    )
    parser.add_argument(
        "--output-hvg-genes-prefix",
        type=str,
        required=True,
        help="Output file prefix to save HVG metadata (highly variable genes)",
    )

    args = parser.parse_args()
    main(args)
