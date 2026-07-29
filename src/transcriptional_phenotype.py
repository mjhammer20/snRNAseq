# Imports
import scanpy as sc
import argparse
import pandas as pd
from pathlib import Path
import gc


# Function to read the MMC results from the CSV file and return a pandas DataFrame
def read_csv_results(csv_results_path: str | Path) -> pd.DataFrame:
    """
    Read the results file and return a pandas DataFrame. The first 4 lines of the file are header information and should be skipped.

    Args:
        csv_results_path: Path to the CSV file containing the MMC results. The first 4 lines of the file are header information and should be skipped.
    
    Returns:
        DataFrame containing the MMC results
    """
    # Read the CSV file, skipping the first 4 lines of header information
    print(f"Reading MMC results from {csv_results_path}...")
    
    results = pd.read_csv(csv_results_path, header=4)
    results["cell_id"] = results["cell_id"].astype(str)
    results.set_index("cell_id", inplace=True)
    results.index.name = None
    
    return results

# Function to transform the MMC results and assign cell types based on the correlation coefficients and bootstrapping probabilities.
def transform_mmc_results(mmc_results: pd.DataFrame):
    """
    Transform the MMC results to assign cell types based on the correlation coefficients and bootstrapping probabilities. If either the correlation coefficient or bootstrapping probability is less than 0.5, the cell type will be assigned as "Unknown".

    Args:
        mmc_results: DataFrame containing the MMC results with columns for class and subclass names, correlation coefficients, and bootstrapping probabilities.
    Returns:
        DataFrame with columns for cell type, phenotype, correlation coefficient, bootstrapping probability, class name, subclass name, and supertype name.
    """

    # Extract cell type mappings for 3 highest levels of the taxonomy (supertype, type, subtype)
    mmc_results["cell_supertype"] = mmc_results.iloc[:,1]
    mmc_results["prob_supertype"] = mmc_results.iloc[:,2]
    mmc_results["rho_supertype"] = mmc_results.iloc[:,4]
    mmc_results["cell_type"] = mmc_results.iloc[:,6]
    mmc_results["prob_type"] = mmc_results.iloc[:,7]
    mmc_results["rho_type"] = mmc_results.iloc[:,9]
    mmc_results["cell_subtype"] = mmc_results.iloc[:,11]
    mmc_results["prob_subtype"] = mmc_results.iloc[:,12]
    mmc_results["rho_subtype"] = mmc_results.iloc[:,14]

    # Change the phenotype to unknown if the correlation or bootstrap probability < 0.5
    mmc_results.loc[mmc_results["rho_supertype"] < 0.5, "cell_supertype"] = "Unknown"
    mmc_results.loc[mmc_results["prob_supertype"] < 0.5, "cell_supertype"] = "Unknown"
    mmc_results.loc[mmc_results["rho_type"] < 0.5, "cell_type"] = "Unknown"
    mmc_results.loc[mmc_results["prob_type"] < 0.5, "cell_type"] = "Unknown"
    mmc_results.loc[mmc_results["rho_subtype"] < 0.5, "cell_subtype"] = "Unknown"
    mmc_results.loc[mmc_results["prob_subtype"] < 0.5, "cell_subtype"] = "Unknown"

    # Keep only the relevant columns for merging with the adata obs
    mmc_results = mmc_results[["cell_supertype", "cell_type", "cell_subtype", "rho_supertype", "prob_supertype", "rho_type", "prob_type", "rho_subtype", "prob_subtype"]]

    return mmc_results


def main(args: argparse.Namespace):

    # Extract arguments
    adata_input_prefix = args.adata_input_prefix
    mmc_results_prefix = args.mmc_results_prefix
    adata_output_prefix = args.adata_output_prefix
    output_cell_types_file_prefix = args.output_cell_types_file_prefix
    regions = args.regions.split()

    # Iterate over regions and process each one
    for region in regions:
        print(f"Processing region: {region}")
        
        # Load MMC results
        mmc_results_fp = f"{mmc_results_prefix}_{region}.results.csv"
        print(f"Loading MMC results for region {region} from {mmc_results_fp}...")
        mmc_results = read_csv_results(mmc_results_fp)

        # Transform the MMC results to assign cell types based on the correlation coefficients and bootstrapping probabilities
        print("Transforming MMC results...")
        transformed_mmc_results = transform_mmc_results(mmc_results)
        
        # Save the results to parquet file
        output_cell_types_fp = f"{output_cell_types_file_prefix}_{region}.parquet"
        print(f"Saving cell type assignments for region {region} to {output_cell_types_fp}...")
        transformed_mmc_results.to_parquet(output_cell_types_fp, compression="gzip")

        # Load the adata
        adata_input_fp = f"{adata_input_prefix}_{region}.h5ad"
        print(f"Loading AnnData object for region {region} from {adata_input_fp}...")
        adata = sc.read_h5ad(adata_input_fp)

        # Merge the results with the adata obs
        print("Merging MMC results with AnnData object...")
        adata.obs = adata.obs.merge(transformed_mmc_results, left_index=True, right_index=True)
        
        # Save the adata
        adata_output_fp = f"{adata_output_prefix}_{region}.h5ad"
        print(f"Saving AnnData object with assigned cell types to {adata_output_fp}...")
        adata.write_h5ad(filename=adata_output_fp, compression="gzip")

        # Clean up memory
        del adata
        del mmc_results
        del transformed_mmc_results
        gc.collect()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Assign cell type to high-fidelity mappings/transcriptional phenotype (MMC)")
    parser.add_argument(
        "--adata-input-prefix",
        type=str,
        required=True,
        help="Prefix for AnnData object files",
    )
    parser.add_argument(
        "--mmc-results-prefix",
        type=str,
        required=True,
        help="Path to MMC results CSV",
    )
    parser.add_argument(
        "--regions",
        type=str,
        required=True,
        help="Space-separated list of regions to process.",
    )
    parser.add_argument(
        "--adata-output-prefix",
        type=str,
        required=True,
        help="Output file to save AnnData object to",
    )
    parser.add_argument(
        "--output-cell-types-file-prefix",
        type=str,
        required=True,
        help="Output file to write cell types to",
    )

    args = parser.parse_args()
    main(args)
