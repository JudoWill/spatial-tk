#!/usr/bin/env python3
"""
Data I/O operations for Xenium spatial clustering tool.

This module handles loading Xenium spatial datasets from .zarr format,
concatenating multiple samples, and saving processed results.
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import anndata as ad
import pandas as pd
import spatialdata as sd


def load_sample_metadata(csv_path: str) -> pd.DataFrame:
    """
    Load sample metadata from CSV file.
    
    The CSV must contain at minimum:
    - sample: Sample name/identifier
    - path: Path to .zarr file
    
    Additional columns (e.g., status, location) are preserved as metadata.
    
    Args:
        csv_path: Path to CSV file with sample information
        
    Returns:
        DataFrame with sample metadata
        
    Raises:
        ValueError: If required columns are missing
    """
    logging.info(f"Loading sample metadata from {csv_path}")
    
    df = pd.read_csv(csv_path)
    
    # Validate required columns
    required_cols = ["sample", "path"]
    missing_cols = set(required_cols) - set(df.columns)
    if missing_cols:
        raise ValueError(f"CSV missing required columns: {missing_cols}")
    
    # Get optional metadata columns
    metadata_cols = [col for col in df.columns if col not in required_cols]
    
    logging.info(f"Loaded {len(df)} samples")
    if metadata_cols:
        logging.info(f"  Metadata columns: {', '.join(metadata_cols)}")
    
    return df


def load_spatial_datasets(sample_df: pd.DataFrame) -> List[Tuple[str, sd.SpatialData]]:
    """
    Load Xenium spatial datasets from .zarr files.
    
    Args:
        sample_df: DataFrame with 'sample' and 'path' columns
        
    Returns:
        List of tuples (sample_name, SpatialData object)
    """
    logging.info("Loading spatial datasets")
    
    spatial_data_list = []
    
    for idx, row in sample_df.iterrows():
        sample_name = row["sample"]
        zarr_path = row["path"]
        
        logging.info(f"  Loading {sample_name} from {zarr_path}")
        
        try:
            sdata = sd.read_zarr(zarr_path)
            spatial_data_list.append((sample_name, sdata))
            
            # Log basic info about the dataset
            if sdata.table is not None:
                n_cells = sdata.table.n_obs
                n_genes = sdata.table.n_vars
                logging.info(f"    {n_cells} cells × {n_genes} genes")
            else:
                logging.warning(f"    No expression table found in {sample_name}")
                
        except Exception as e:
            logging.error(f"  Failed to load {sample_name}: {e}")
            raise
    
    logging.info(f"Successfully loaded {len(spatial_data_list)} spatial datasets")
    return spatial_data_list


def concatenate_spatial_data(
    spatial_data_list: List[Tuple[str, sd.SpatialData]],
    sample_df: pd.DataFrame
) -> sd.SpatialData:
    """
    Concatenate multiple SpatialData objects into one, preserving metadata.
    
    Args:
        spatial_data_list: List of (sample_name, SpatialData) tuples
        sample_df: DataFrame with sample metadata
        
    Returns:
        Concatenated SpatialData object with metadata in .table.obs
    """
    logging.info("Concatenating spatial datasets")
    
    if len(spatial_data_list) == 0:
        raise ValueError("No spatial datasets to concatenate")
    
    # Get table accessor (handle both .table and .tables API)
    def get_table(sdata):
        if hasattr(sdata, 'tables') and len(sdata.tables) > 0:
            return list(sdata.tables.values())[0]
        elif hasattr(sdata, 'table'):
            return sdata.table
        return None
    
    if len(spatial_data_list) == 1:
        sample_name, sdata = spatial_data_list[0]
        logging.info("Single sample - no concatenation needed")
        
        # Add metadata to the table
        table = get_table(sdata)
        if table is not None:
            metadata_cols = [col for col in sample_df.columns if col not in ["sample", "path"]]
            sample_metadata = sample_df[sample_df["sample"] == sample_name].iloc[0]
            
            table.obs["sample"] = sample_name
            for col in metadata_cols:
                table.obs[col] = sample_metadata[col]
        
        return sdata
    
    # Extract SpatialData objects and their names
    sdata_dict = {name: sdata for name, sdata in spatial_data_list}
    
    # Concatenate using spatialdata's concatenate function
    # Pass as dict to handle duplicate label names
    try:
        concatenated_sdata = sd.concatenate(
            sdata_dict,
            region_key="region",
            instance_key="instance_id",
            concatenate_tables=True
        )
        
        # Get the concatenated table
        table = get_table(concatenated_sdata)
        
        if table is not None:
            # Add sample names - extract from region key
            # spatialdata.concatenate adds element prefixes (e.g., "cell_circles-Drexel-Pos")
            # so we need to extract just the sample name part
            if "region" in table.obs.columns:
                # Try to extract sample name from region
                # Region format: "element_name-sample_name"
                def extract_sample_name(region_str):
                    # Split on '-' and look for matching sample names
                    for sample_name in sample_df["sample"].values:
                        if str(region_str).endswith(str(sample_name)):
                            return sample_name
                    # Fallback: return the region as is
                    return region_str
                
                table.obs["sample"] = table.obs["region"].apply(extract_sample_name)
            
            # Add additional metadata from CSV
            metadata_cols = [col for col in sample_df.columns if col not in ["sample", "path"]]
            if metadata_cols:
                logging.info(f"  Adding metadata columns: {', '.join(metadata_cols)}")
                
                # Create a mapping from sample name to metadata
                metadata_dict = {}
                for col in metadata_cols:
                    metadata_dict[col] = sample_df.set_index("sample")[col].to_dict()
                
                # Add metadata to obs
                for col in metadata_cols:
                    table.obs[col] = table.obs["sample"].map(metadata_dict[col])
            
            total_cells = table.n_obs
            total_genes = table.n_vars
            logging.info(f"Concatenation complete: {total_cells} total cells × {total_genes} genes")
        else:
            logging.warning("No table found after concatenation")
        
        return concatenated_sdata
        
    except Exception as e:
        logging.error(f"Failed to concatenate spatial datasets: {e}")
        raise


def save_spatial_data(sdata: sd.SpatialData, output_path: Path, overwrite: bool = False) -> None:
    """
    Save SpatialData object to .zarr format.
    
    Args:
        sdata: SpatialData object to save
        output_path: Path where .zarr will be saved
        overwrite: Whether to overwrite an existing store (required for inplace operations)
    """
    import shutil
    import tempfile
    
    logging.info(f"Saving spatial data to {output_path}")
    
    try:
        if overwrite and output_path.exists():
            # Workaround for SpatialData limitation: cannot overwrite a store that's currently in use
            # Save to temporary location first, then replace the original
            # See: https://github.com/scverse/spatialdata/discussions/520
            temp_dir = Path(tempfile.mkdtemp(prefix="spatialdata_tmp_", dir=output_path.parent))
            temp_path = temp_dir / output_path.name
            
            try:
                # Save to temporary location
                sdata.write(temp_path)
                
                # Remove original store
                shutil.rmtree(output_path)
                
                # Move temporary store to original location
                shutil.move(str(temp_path), str(output_path))
                
                # Clean up temporary directory
                temp_dir.rmdir()
                
                logging.info(f"Successfully saved spatial data (overwrite)")
            except Exception as e:
                # Clean up temporary directory on error
                if temp_dir.exists():
                    shutil.rmtree(temp_dir)
                raise
        else:
            # Normal save operation
            sdata.write(output_path, overwrite=overwrite)
            logging.info(f"Successfully saved spatial data")
    except Exception as e:
        logging.error(f"Failed to save spatial data: {e}")
        raise


def load_existing_spatial_data(zarr_path: Path) -> sd.SpatialData:
    """
    Load an existing processed SpatialData object.
    
    Args:
        zarr_path: Path to .zarr file
        
    Returns:
        SpatialData object
    """
    logging.info(f"Loading existing spatial data from {zarr_path}")
    
    try:
        sdata = sd.read_zarr(zarr_path)
        
        # Get table (handle both .table and .tables API)
        table = None
        if hasattr(sdata, 'tables') and len(sdata.tables) > 0:
            table = list(sdata.tables.values())[0]
        elif hasattr(sdata, 'table'):
            table = sdata.table
        
        if table is not None:
            logging.info(f"Loaded: {table.n_obs} cells × {table.n_vars} genes")
        
        return sdata
    except Exception as e:
        logging.error(f"Failed to load spatial data: {e}")
        raise

