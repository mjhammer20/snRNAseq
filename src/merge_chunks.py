# Imports 
import scanpy as sc
from anndata import AnnData
from anndata import concat as ad_concat
import argparse
import gc


# Merge filtered chunks back into single AnnData object
def merge_chunks(chunk_files: list, output_file: str):
    """Merge filtered chunks back into single AnnData object."""
    print(f"Loading {len(chunk_files)} chunks...")

    adata_merged = None
    for i, f in enumerate(chunk_files):
        
        print(f"  Loading chunk {i+1}: {f}")
        if adata_merged is None:
            adata_merged = sc.read_h5ad(f)
        else:
            adata = sc.read_h5ad(f)
            adata_merged = ad_concat([adata_merged, adata], merge="same", uns_merge="same", index_unique="_")

            # Clean up memory
            del adata
            gc.collect()
    
    print(f"Merged to {adata_merged.n_obs} cells")
    adata_merged.write_h5ad(output_file, compression="gzip")
    print(f"Saved: {output_file}")


# Main
if __name__=='__main__':
    output_dir = "/mnt/output/output"
    chunk_files = [
        f"{output_dir}/Putamen/merged_cleaned_filtered_22.h5ad",
        f"{output_dir}/Putamen/merged_cleaned_filtered_32.h5ad"
    ]
    output_file = f"{output_dir}/merged_cleaned_filtered_Putamen.h5ad"
    merge_chunks(chunk_files, output_file)