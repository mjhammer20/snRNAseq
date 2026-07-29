# Imports
import os
import argparse
import anndata as ad
import scvi
import gc
import multiprocessing as mp
import torch
from concurrent.futures import ProcessPoolExecutor, as_completed

# Function to integrate data with scVI
def integrate_with_scvi(
    adata: ad.AnnData,
    batch_key: str,
    latent_key: str,
) -> tuple[ad.AnnData, scvi.model.SCVI]:
    """
    Fit scVI model to AnnData object

    Args:
        adata: AnnData object to integrate with scVI
        batch_key: Key in adata.obs for batch information
        latent_key: Key to save the scVI latent representation to in adata.obsm
    
    Returns:
        Tuple of (integrated AnnData object, scVI model)
    """

    # Fixed parameters
    n_latent = 30
    n_layers = 2
    train_size = 0.85
    scvi_epochs = 300
    accelerator = "gpu"
    dispersion = "gene-cell" 
    gene_likelihood = "zinb"
    early_stopping = True
    early_stopping_patience = 20

    # Set parameters based on numerical instabilities
    threshold_cells = 3.05e6 # No. of cells in Sep 2025 PMDBS sc cohort (Lee, Hardy, Hafler, Jakobsson, Scherzer)
    if adata.n_obs > threshold_cells:
        plan_kwargs = {"lr": 1e-4}
        gradient_clip_val = 5.0
        print(f"AnnData object contains {adata.n_obs} which is > {threshold_cells}")
        print(f"--- Using learning rate: {plan_kwargs}")
        print(f"--- Using gradient clipping: {gradient_clip_val}")
    else:
        # Defaults
        plan_kwargs = {"lr": 1e-3} 
        gradient_clip_val = None
        print(f"AnnData object contains {adata.n_obs} which is < {threshold_cells}")
        print(f"--- Using default learning rate: {plan_kwargs}")
        print(f"--- Using default gradient clipping: {gradient_clip_val}")

    # Adjust batch_size based on dataset size
    # Use smaller batch size or drop_last=True to avoid single-sample batches
    # Calculate optimal batch size to avoid last batch with only 1 sample
    n_samples = adata.n_obs
    batch_size = 128  # default
    for candidate_size in [256, 128, 64, 32, 16]:
        remainder = n_samples % candidate_size
        # Accept if divides evenly OR last batch has >= 2 samples
        if remainder == 0 or remainder >= 2:
            batch_size = candidate_size
            break

    # Setup AnnData for scVI
    noise = ["doublet_score", "pct_counts_mt", "pct_counts_rb"]
    categorical_covariate_keys = None
    scvi.model.SCVI.setup_anndata(
        adata,
        layer=None,
        batch_key=batch_key,
        continuous_covariate_keys=noise,
        categorical_covariate_keys=categorical_covariate_keys,
    )

    # Initialize scVI model
    model = scvi.model.SCVI(
        adata,
        n_layers=n_layers,
        n_latent=n_latent,
        dispersion=dispersion,
        gene_likelihood=gene_likelihood,
    )

    # Train the model
    print(f"Training scVI model with batch size {batch_size}, learning rate {plan_kwargs['lr']}, and gradient clipping {gradient_clip_val}")
    model.train(
        train_size=train_size,
        max_epochs=scvi_epochs,
        early_stopping=early_stopping,
        early_stopping_patience=early_stopping_patience,
        accelerator=accelerator,
        batch_size=batch_size,
        gradient_clip_val=gradient_clip_val,
        plan_kwargs=plan_kwargs,
    )

    # Save the latent representation to adata.obsm
    adata.obsm[latent_key] = model.get_latent_representation()

    return (adata, model)

# Function for processing each region
def process_region(
        region: str, 
        adata_input_prefix: str, 
        adata_output_prefix: str,
        output_scvi_dir: str, 
        batch_key: str, 
        latent_key: str,
        gpu_id: int
):
    """
    Train scVI model and save integrated AnnData for a single region

    Args:
        region: Region to process
        adata_input_prefix: Input file prefix for AnnData object
        adata_output_prefix: Output file prefix for integrated AnnData object
        output_scvi_dir: Output directory to save scVI model
        batch_key: Key in adata.obs for batch information
        latent_key: Key to save the scVI latent representation to in adata.obsm
        gpu_id: GPU ID to use for this process

    Returns: 
        None
    """
    # Set CUDA device for this worker
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)

    print(f"GPU {gpu_id} - Processing region: {region}")

    # Load data
    adata_input = f"{adata_input_prefix}_{region}.h5ad"
    print(f"GPU {gpu_id} - Loading AnnData for region {region} from {adata_input}")
    adata = ad.read_h5ad(adata_input)

    # Process data
    print(f"GPU {gpu_id} - Integrating data for region {region} with scVI")
    adata, model = integrate_with_scvi(adata=adata, batch_key=batch_key, latent_key=latent_key)
    
    # Save the scVI model
    region_scvi_dir = f"{output_scvi_dir}/{region}"
    os.makedirs(region_scvi_dir, exist_ok=True)
    print(f"GPU {gpu_id} - Saving scVI model for region {region} to {region_scvi_dir}")
    model.save(region_scvi_dir, overwrite=True)

    # Save the integrated AnnData object
    adata_output = f"{adata_output_prefix}_{region}.h5ad"
    print(f"GPU {gpu_id} - Saving integrated AnnData for region {region} to {adata_output}")
    adata.write_h5ad(filename=adata_output, compression="gzip")

    # Clean up memory
    del adata
    del model
    gc.collect()

    print(f"GPU {gpu_id} - Region {region} complete!")


def main(args: argparse.Namespace):
    
    # Extract arguments    
    latent_key = args.latent_key
    batch_key = args.batch_key
    adata_input_prefix = args.adata_input_prefix
    adata_output_prefix = args.adata_output_prefix
    output_scvi_dir = args.output_scvi_dir
    regions = args.regions.split()

    # Create output directory
    os.makedirs(output_scvi_dir, exist_ok=True)

    # Number of GPUs
    n_gpus = 4

    # Prevent torch from using too many threads which can cause issues in multiprocessing
    torch.set_num_threads(1)
        
    # Set scVI settings before loading data
    scvi.settings.dl_num_workers = 0  # Disable multiprocessing in data loader
    
    # Parallelize across regions using ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=n_gpus) as executor:
        futures = {}
        for idx, region in enumerate(regions):
            gpu_id = idx % n_gpus  # Round-robin GPU assignment
            future = executor.submit(
                process_region,
                region=region, adata_input_prefix=adata_input_prefix, adata_output_prefix=adata_output_prefix,
                output_scvi_dir=output_scvi_dir, batch_key=batch_key, latent_key=latent_key, gpu_id=gpu_id
            )
            futures[future] = region
        
        # Wait for tasks as they complete (not in order)
        for future in as_completed(futures):
            region = futures[future]
            try:
                future.result()
                print(f"✓ Region {region} completed")
            except Exception as e:
                print(f"✗ Region {region} failed: {e}")
                import traceback
                traceback.print_exc()

    print("All regions processed!")
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run scVI integration")
    parser.add_argument(
        "--latent-key",
        type=str,
        required=True,
        help="Latent key to save the scVI latent to",
    )
    parser.add_argument(
        "--batch-key",
        type=str,
        required=True,
        help="Key in AnnData object for batch information",
    )
    parser.add_argument(
        "--adata-input-prefix",
        type=str,
        required=True,
        help="AnnData object prefix for a dataset",
    )
    parser.add_argument(
        "--adata-output-prefix",
        type=str,
        required=True,
        help="Output file prefix to save AnnData object to",
    )
    parser.add_argument(
        "--regions",
        type=str,
        required=True,
        help="Regions to process, space separated",
    )
    parser.add_argument(
        "--output-scvi-dir",
        type=str,
        required=True,
        help="Output folder to save scVI model",
    )

    args = parser.parse_args()
    main(args)
