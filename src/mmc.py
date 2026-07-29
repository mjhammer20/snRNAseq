# Imports
import pathlib
from cell_type_mapper.cli.map_to_on_the_fly_markers import OnTheFlyMapper
from cell_type_mapper.cli.from_specified_markers import FromSpecifiedMarkersRunner
import cell_type_mapper.test_utils.cache_wrapper as cache_module
import argparse


def run_mmc(
          adata_input_fp: str, 
          mmc_precomputed_stats: str, 
          mmc_marker_genes: str, 
          mmc_gene_mapper_db: str, 
          extended_results_json: str, 
          csv_results: str, 
          log_file: str, 
          n_processors: int, 
          max_gb: int, 
          rng_seed: int, 
          chunk_size: int, 
          n_runners_up: int):
    '''
    Run Allen Institute's Cell Type Mapper (MMC) on a given AnnData object with specified parameters and output paths.
    
    args:
        adata_input_fp (str): File path to input AnnData object (h5ad file)
        mmc_precomputed_stats (str): File path to MMC precomputed stats HD5F file
        mmc_marker_genes (str): File path to MMC marker genes JSON file (optional; if not provided, MMC will compute markers on the fly)
        mmc_gene_mapper_db (str): File path to MMC gene mapper SQLite database
        extended_results_json (str): File path to output extended results JSON file
        csv_results (str): File path to output CSV results file
        log_file (str): File path to output log file
        n_processors (int): The number of independent worker processes to spin up
        max_gb (int): Maximum GB of memory to use (will be converted to bytes and passed to MMC)
        rng_seed (int): Random seed for reproducibility
        chunk_size (int): Number of cells to process in a chunk (tune based on memory constraints)
        n_runners_up (int): Number of runner-up cell type predictions to include in output

    '''

    # Config definition for MMC
    config = {
        "precomputed_stats": {
            "path": mmc_precomputed_stats
        },
        "type_assignment": {
            "bootstrap_factor": 0.5,
            "bootstrap_iteration": 100,
            "normalization": "raw",
            "rng_seed": rng_seed,
            "chunk_size": chunk_size,
            "n_runners_up": n_runners_up
        },
        "query_path": adata_input_fp,
        "extended_result_path": extended_results_json,
        "csv_result_path": csv_results,
        "log_path": log_file,
        "cloud_safe": True,
        "verbose_csv": True,
        "max_gb": max_gb,
        "gene_mapping": {
            "db_path": mmc_gene_mapper_db
        }
    }

    # If marker genes are provided, use FromSpecifiedMarkersRunner; otherwise, use OnTheFlyMapper to compute markers and map cell types
    if mmc_marker_genes:
        print("Running FromSpecifiedMarkersRunner")
        config["type_assignment"]["n_processors"] = n_processors
        config["query_markers"] = {
            "serialized_lookup": mmc_marker_genes
        }
        runner = FromSpecifiedMarkersRunner(args=[], input_data=config)
    else:
        print("Running OnTheFlyMapper")
        config["reference_markers"] = {
            "log2_fold_min_th": 0.5
        }
        config["query_markers"] = {
            "n_per_utility": 15,
            "genes_at_a_time": 1
        }
        runner = OnTheFlyMapper(args=[], input_data=config)

    runner.run()

def main(args: argparse.Namespace):
    
    # Extract arguments
    mmc_marker_genes = args.mmc_marker_genes
    mmc_precomputed_stats = args.mmc_precomputed_stats
    mmc_gene_mapper_db = args.mmc_gene_mapper_db
    adata_input_prefix = args.adata_input_prefix
    regions = args.regions.split()
    n_processors = args.n_processors
    max_gb = args.max_gb
    rng_seed = args.rng_seed
    chunk_size = args.chunk_size
    n_runners_up = args.n_runners_up
    output_prefix = args.output_prefix

    # Iterate over regions and run MMC for each region's stratified dataset
    for region in regions:

        # Define adata input file path for this region
        adata_input_fp = f"{adata_input_prefix}_{region}.h5ad"
        print(f"Processing region: {region} with input file: {adata_input_fp}")

        # Define output file paths for this region
        extended_results_json = f"{output_prefix}_{region}.extended_results.json"
        print(f"EXTENDED_RESULTS: {extended_results_json}")
        csv_results = f"{output_prefix}_{region}.results.csv"
        print(f"CSV_RESULTS: {csv_results}")
        log_file = f"{output_prefix}_{region}.log.txt"
        print(f"LOG_FILE: {log_file}")

        # Run MMC for this region
        run_mmc(
            adata_input_fp=adata_input_fp,
            mmc_precomputed_stats=mmc_precomputed_stats,
            mmc_marker_genes=mmc_marker_genes,
            mmc_gene_mapper_db=mmc_gene_mapper_db,
            extended_results_json=extended_results_json,
            csv_results=csv_results,
            log_file=log_file,
            n_processors=n_processors,
            max_gb=max_gb,
            rng_seed=rng_seed,
            chunk_size=chunk_size,
            n_runners_up=n_runners_up
        )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Map cell types with Allen Institute's Cell Type Mapper")
    parser.add_argument(
        "--adata-input-prefix",
        type=str,
        required=True,
        help="Prefix for AnnData object for stratified datasets",
    )
    parser.add_argument(
        "--regions",
        type=str,
        required=True,
        help="Path to text file with list of regions to process (one region per line)",
    )
    parser.add_argument(
        "--mmc-precomputed-stats",
        type=str,
        required=True,
        help="Path to MMC precomputed stats HD5F file",
    )
    parser.add_argument(
        "--mmc-marker-genes",
        type=str,
        required=False,
        help="Path to MMC marker genes JSON file",
    )
    parser.add_argument(
        "--mmc-gene-mapper-db",
        type=str,
        required=True,
        help="Path to MMC gene mapper SQLite database",
    )
    parser.add_argument(
        "--output-prefix",
        type=str,
        required=True,
        help="Prefix used to name output files",    
    )
    parser.add_argument(
        "--n-processors",
        type=int,
        required=True,
        help="The number of independent worker processes to spin up",
    )
    parser.add_argument(
        "--max-gb",
        type=float,
        required=True,
        help="Maximum GB of memory to use (will be converted to bytes and passed to MMC)",
    )
    parser.add_argument(
        "--rng-seed",
        type=int,
        required=True,
        help="Random seed for reproducibility",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        required=True,
        help="Number of cells to process in a chunk (tune based on memory constraints)",
    )
    parser.add_argument(
        "--n-runners-up",
        type=int,
        required=True,
        help="Number of runner-up cell type predictions to include in output",
    )

    args = parser.parse_args()
    main(args)

