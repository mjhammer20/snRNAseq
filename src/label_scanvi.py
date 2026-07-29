# Imports
import os
import argparse
import scvi
import gc
import torch
import anndata as ad
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed

# Function to label cells with scANVI
def label_with_scanvi(
    adata: ad.AnnData,
    model: scvi.model.SCVI,
    num_workers: int,
    latent_key: str,
    cell_type_label_key: str,
    predictions_key: str, 
    predictions_score_key: str
) -> tuple[ad.AnnData, scvi.model.SCANVI]:
    """
    Fit scANVI model to AnnData object

    Args:
        adata: AnnData object to label with scANVI
        model: Pre-trained scVI model to initialize scANVI from
        num_workers: Number of workers for data loading during training
        latent_key: Key to save the scANVI latent representation to in adata.obsm
        predictions_key: Key to save the scANVI predictions to in adata.obs
        predictions_score_key: Key to save the scANVI predictions scores to in adata.obs

    Returns:
        Tuple of (labeled AnnData object, trained scANVI model)
    """
    # Fixed parameters
    scanvi_epochs = 300
    accelerator = "gpu"
    early_stopping = True
    early_stopping_patience = 20
    train_size = 0.85

    # Adjust batch_size based on dataset size
    n_samples = adata.n_obs
    batch_size = 128
    for candidate_size in [256, 128, 64, 32, 16]:
        remainder = n_samples % candidate_size
        if remainder == 0 or remainder >= 2:
            batch_size = candidate_size
            break

    # Set parameters based on numerical instabilities
    threshold_cells = 3.05e6
    
    if adata.n_obs > threshold_cells:
        plan_kwargs = {"lr": 5e-5}
        gradient_clip_val = 200.0
    else:
        plan_kwargs = {"lr": 5e-4}
        gradient_clip_val = 200.0

    # Initialize scANVI model
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        model,
        adata=adata,  # Train on ALL cells, including rare types
        labels_key=cell_type_label_key,
        unlabeled_category="Unknown",
    )

    # Train the model
    try:
        print("Training scANVI model")
        scanvi_model.train(
            accelerator=accelerator,
            max_epochs=scanvi_epochs,
            early_stopping=early_stopping,
            early_stopping_patience=early_stopping_patience,
            train_size=train_size,
            batch_size=batch_size,
            datasplitter_kwargs={"num_workers": num_workers},
            gradient_clip_val=gradient_clip_val,
            plan_kwargs=plan_kwargs,
        )

    except ValueError as e:
        if "nan" in str(e).lower():
            print(f"[WARNING] NaN detected during training: {e}")
            del scanvi_model
            gc.collect()
            torch.cuda.empty_cache()

            scanvi_model = scvi.model.SCANVI.from_scvi_model(
                model,
                adata=adata,
                labels_key=cell_type_label_key,
                unlabeled_category="Unknown",
            )
            print("Re-training with even more conservative learning rate...")
            plan_kwargs = {"lr": 1e-4}
            gradient_clip_val = 500.0
            scanvi_model.train(
                accelerator=accelerator,
                max_epochs=scanvi_epochs,
                early_stopping=early_stopping,
                early_stopping_patience=early_stopping_patience,
                train_size=train_size,
                batch_size=batch_size,
                datasplitter_kwargs={"num_workers": num_workers},
                gradient_clip_val=gradient_clip_val,
                plan_kwargs=plan_kwargs,
            )
        else:
            raise e

    # Generate predictions and latent representation
    print("Generating scANVI predictions")
    adata.obsm[latent_key] = scanvi_model.get_latent_representation(adata)
    adata.obs[predictions_key] = scanvi_model.predict(adata)
    adata.obs[predictions_score_key] = scanvi_model.predict(adata, soft=True).max(axis=1).values
    return (adata, scanvi_model)


# Function to process a single region with scANVI
def process_region(
        region: str,
        adata_input_prefix: str,
        scvi_outputs_dir_prefix: str,
        adata_output_prefix: str,
        scanvi_outputs_dir_prefix: str,
        output_cell_types_file_prefix: str,
        latent_key: str,
        cell_type_label_key: str,
        predictions_key: str,
        predictions_score_key: str,
        num_workers: int,
        gpu_id: int
):
    """
    Train scANVI model and save labeled AnnData for a single region

    Args:
        region: Region to process
        adata_input_prefix: Input file prefix for AnnData object
        scvi_outputs_dir_prefix: Directory containing scVI models
        adata_output_prefix: Output file prefix for labeled AnnData object
        scanvi_outputs_dir_prefix: Output directory to save scANVI model
        output_cell_types_file_prefix: Output file prefix for cell types
        latent_key: Key to save the scANVI latent representation
        predictions_key: Key to save the scANVI predictions
        num_workers: Number of workers for data loading
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

    region_scvi_dir = f"{scvi_outputs_dir_prefix}/{region}"
    print(f"GPU {gpu_id} - Loading scVI model for region {region} from {region_scvi_dir}")
    model = scvi.model.SCVI.load(
        dir_path=region_scvi_dir,
        adata=adata,
    )
    
    # Predict cell types with scANVI
    print(f"GPU {gpu_id} - Predicting cell types for region {region} with scANVI")
    adata, scanvi_model = label_with_scanvi(adata=adata, model=model, num_workers=num_workers, latent_key=latent_key, cell_type_label_key=cell_type_label_key, predictions_key=predictions_key, predictions_score_key=predictions_score_key)

    # Save scanvi model
    region_scanvi_dir = f"{scanvi_outputs_dir_prefix}/{region}"
    os.makedirs(region_scanvi_dir, exist_ok=True)
    print(f"GPU {gpu_id} - Saving scANVI model for region {region} to {region_scanvi_dir}") 
    scanvi_model.save(region_scanvi_dir, overwrite=True)

    # Save labeled AnnData object
    adata_output = f"{adata_output_prefix}_{region}.h5ad"
    adata.write_h5ad(filename=adata_output, compression="gzip")

    # Save cell types to parquet
    output_cell_types_file = f"{output_cell_types_file_prefix}_{region}.parquet"
    adata.obs[[predictions_key, predictions_score_key]].to_parquet(output_cell_types_file, compression="gzip")

    # Clean up memory
    del adata
    del model
    del scanvi_model
    gc.collect()

    print(f"GPU {gpu_id} - Region {region} complete!")


def main(args: argparse.Namespace):
    # Extract arguments
    regions = args.regions.split()
    latent_key = args.latent_key_scanvi
    cell_type_label_key = args.cell_type_label_key
    predictions_key = args.predictions_key
    predictions_score_key = f"{predictions_key}_score"
    adata_input_prefix = args.adata_input_prefix
    scvi_outputs_dir_prefix = args.scvi_outputs_dir_prefix
    adata_output_prefix = args.adata_output_prefix
    scanvi_outputs_dir_prefix = args.scanvi_outputs_dir_prefix
    output_cell_types_file_prefix = args.output_cell_types_file_prefix

    # Create output directory
    os.makedirs(scanvi_outputs_dir_prefix, exist_ok=True)

    # Number of GPUs
    n_gpus = 4

    # Prevent torch from using too many threads which can cause issues in multiprocessing
    torch.set_num_threads(1)
        
    # Set number of workers for scANVI training
    num_workers = 0  # Disable multiprocessing in data loader to avoid nested multiprocessing
    scvi.settings.dl_num_workers = num_workers  
    
    # Parallelize across regions using ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=n_gpus) as executor:
        futures = {}
        for idx, region in enumerate(regions):
            gpu_id = idx % n_gpus  # Round-robin GPU assignment
            future = executor.submit(
                process_region,
                region=region, adata_input_prefix=adata_input_prefix, scvi_outputs_dir_prefix=scvi_outputs_dir_prefix,
                adata_output_prefix=adata_output_prefix, scanvi_outputs_dir_prefix=scanvi_outputs_dir_prefix, output_cell_types_file_prefix=output_cell_types_file_prefix,
                latent_key=latent_key, cell_type_label_key=cell_type_label_key, predictions_key=predictions_key, predictions_score_key=predictions_score_key, num_workers=num_workers, gpu_id=gpu_id
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
    parser = argparse.ArgumentParser(description="Leverage cell-type from MMC to assign the rest of the cells with scANVI")
    parser.add_argument(
        "--latent-key-scanvi",
        type=str,
        required=True,
        help="Latent key to save the scANVI latent to",
    )
    parser.add_argument(
        "--cell-type-label-key",
        type=str,
        required=True,
        help="Key in adata.obs to use for cell type labels for scANVI training",
    )
    parser.add_argument(
        "--predictions-key",
        type=str,
        required=True,
        help="scANVI cell type predictions column name in AnnData object",
    )
    parser.add_argument(
        "--adata-input-prefix",
        type=str,
        required=True,
        help="AnnData object for a dataset",
    )
    parser.add_argument(
        "--scvi-outputs-dir-prefix",
        type=str,
        required=True,
        help="Saved scVI outputs folder",
    )
    parser.add_argument(
        "--regions",
        type=str,
        required=True,
        help="Regions to process, space separated",
    )
    parser.add_argument(
        "--adata-output-prefix",
        type=str,
        required=True,
        help="Output file to save AnnData object to",
    )
    parser.add_argument(
        "--scanvi-outputs-dir-prefix",
        type=str,
        required=True,
        help="Output folder to save scANVI model",
    )
    parser.add_argument(
        "--output-cell-types-file-prefix",
        type=str,
        required=True,
        help="Output file to write cell types to",
    )

    args = parser.parse_args()
    main(args)
