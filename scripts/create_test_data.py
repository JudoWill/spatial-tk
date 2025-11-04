#!/usr/bin/env python3
"""
Script to create subsampled test data from full Xenium datasets.

This generates small (~1000 cells) test datasets for functional testing.
"""

import argparse
import logging
import sys
from pathlib import Path

import spatialdata as sd
import scanpy as sc


def setup_logging():
    """Configure logging."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def subsample_zarr(input_path: Path, output_path: Path, n_cells: int = 1000):
    """
    Subsample a .zarr file to a smaller number of cells.
    
    Args:
        input_path: Path to input .zarr file
        output_path: Path to output subsampled .zarr file
        n_cells: Number of cells to keep
    """
    logging.info(f"Loading {input_path}")
    sdata = sd.read_zarr(input_path)
    
    # Get table
    table = None
    if hasattr(sdata, 'tables') and len(sdata.tables) > 0:
        table = list(sdata.tables.values())[0]
    elif hasattr(sdata, 'table'):
        table = sdata.table
    
    if table is None:
        logging.error("No table found in spatial data")
        return False
    
    logging.info(f"  Original: {table.n_obs} cells × {table.n_vars} genes")
    
    # Subsample cells
    if table.n_obs > n_cells:
        sc.pp.subsample(table, n_obs=n_cells)
        logging.info(f"  Subsampled to: {table.n_obs} cells")
    else:
        logging.info(f"  Dataset already has ≤ {n_cells} cells, keeping all")
    
    # Update table in sdata
    if hasattr(sdata, 'tables') and len(sdata.tables) > 0:
        table_name = list(sdata.tables.keys())[0]
        sdata.tables[table_name] = table
    else:
        sdata.table = table
    
    # Save subsampled data
    logging.info(f"  Saving to {output_path}")
    sdata.write(output_path)
    
    return True


def create_test_csv(csv_path: Path, test_dir: Path):
    """
    Create a test CSV file pointing to subsampled data.
    
    Args:
        csv_path: Original CSV file path
        test_dir: Directory containing subsampled .zarr files
    """
    import pandas as pd
    
    logging.info(f"Creating test CSV from {csv_path}")
    
    # Read original CSV
    df = pd.read_csv(csv_path)
    
    # Update paths to point to subsampled files
    new_rows = []
    for _, row in df.iterrows():
        sample_name = row['sample']
        subsampled_name = f"subsampled_{sample_name}.zarr"
        subsampled_path = test_dir / subsampled_name
        
        if subsampled_path.exists():
            row['path'] = str(subsampled_path)
            new_rows.append(row)
        else:
            logging.warning(f"  Subsampled file not found for {sample_name}, skipping")
    
    if not new_rows:
        logging.error("No valid subsampled files found")
        return None
    
    # Create new DataFrame
    test_df = pd.DataFrame(new_rows)
    
    # Save test CSV
    test_csv_path = test_dir / "test_samples.csv"
    test_df.to_csv(test_csv_path, index=False)
    logging.info(f"  Saved test CSV to {test_csv_path}")
    
    return test_csv_path


def main():
    """Main execution."""
    parser = argparse.ArgumentParser(
        description="Create subsampled test data from Xenium datasets"
    )
    parser.add_argument(
        '--input-csv',
        required=True,
        help='Path to original samples CSV file'
    )
    parser.add_argument(
        '--output-dir',
        required=True,
        help='Directory to save subsampled test data'
    )
    parser.add_argument(
        '--n-cells',
        type=int,
        default=1000,
        help='Number of cells to keep per sample (default: 1000)'
    )
    
    args = parser.parse_args()
    
    setup_logging()
    
    logging.info("="*60)
    logging.info("Creating Test Data")
    logging.info("="*60)
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Read CSV
    import pandas as pd
    csv_path = Path(args.input_csv)
    
    if not csv_path.exists():
        logging.error(f"Input CSV not found: {csv_path}")
        sys.exit(1)
    
    df = pd.read_csv(csv_path)
    
    # Subsample each dataset
    success_count = 0
    for _, row in df.iterrows():
        sample_name = row['sample']
        input_zarr = Path(row['path'])
        
        if not input_zarr.exists():
            logging.warning(f"Input file not found: {input_zarr}, skipping")
            continue
        
        output_zarr = output_dir / f"subsampled_{sample_name}.zarr"
        
        logging.info(f"\nProcessing {sample_name}")
        if subsample_zarr(input_zarr, output_zarr, args.n_cells):
            success_count += 1
    
    # Create test CSV
    if success_count > 0:
        test_csv = create_test_csv(csv_path, output_dir)
        
        logging.info("\n" + "="*60)
        logging.info(f"Test data creation complete")
        logging.info(f"  Subsampled {success_count} datasets")
        logging.info(f"  Output directory: {output_dir}")
        if test_csv:
            logging.info(f"  Test CSV: {test_csv}")
        logging.info("="*60)
    else:
        logging.error("No datasets were successfully subsampled")
        sys.exit(1)


if __name__ == "__main__":
    main()

