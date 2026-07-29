import scanpy as sc
from anndata import AnnData
from pathlib import Path
import gc

def split_adata(adata_input: str, output_dir: str, n_chunks: int = 2):
    """
    Split AnnData object into chunks and save each chunk separately.
    
    Args:
        adata_input: Path to input .h5ad file
        output_dir: Directory to save chunk files
        n_chunks: Number of chunks to split into (default: 2)
    """
    print(f"Loading {adata_input}...")
    adata_path = f'{output_dir}/{adata_input}'
    adata = sc.read_h5ad(adata_path, backed='r')
    
    n_obs = adata.n_obs
    chunk_size = n_obs // n_chunks
    
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    for i in range(n_chunks):
        start_idx = i * chunk_size
        # Last chunk gets remaining cells
        end_idx = n_obs if i == n_chunks - 1 else (i + 1) * chunk_size
        
        print(f"Chunk {i+1}/{n_chunks}: cells {start_idx}-{end_idx}")
        
        # Extract chunk
        adata_chunk = adata[start_idx:end_idx].to_memory()
        
        # Save chunk
        chunk_path = f"{output_dir}/merged_cleaned_unfiltered_1_chunk_{i+1}.h5ad"
        adata_chunk.write_h5ad(chunk_path, compression="gzip")
        print(f"  Saved: {chunk_path}")

        # Clean up memory
        del adata_chunk
        gc.collect()

    print("Done splitting!")

if __name__ == "__main__":
    
    adata_input = "merged_cleaned_unfiltered_1.h5ad"
    output_dir = "/mnt/output/output"
    n_chunks = 2

    split_adata(adata_input, output_dir, n_chunks)