#!/usr/bin/env python3
"""
Main CLI entry point for xenium_process.

This module provides the command-line interface for the Xenium spatial
transcriptomics processing pipeline.
"""

import argparse
import sys
import warnings

# Suppress warnings
warnings.filterwarnings('ignore', category=FutureWarning)

from xenium_process.commands import concat, normalize, cluster, annotate, differential
from xenium_process.utils.helpers import setup_logging


def create_parser() -> argparse.ArgumentParser:
    """
    Create the main argument parser with all subcommands.
    
    Returns:
        Configured ArgumentParser
    """
    parser = argparse.ArgumentParser(
        prog='xenium_process',
        description='Xenium Spatial Transcriptomics Processing Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Concatenate multiple samples
  xenium_process concat --input samples.csv --output merged.zarr
  
  # Normalize (in place to save space)
  xenium_process normalize --input merged.zarr --inplace --save-plots
  
  # Cluster with multiple resolutions
  xenium_process cluster --input merged.zarr --inplace --leiden-resolution 0.2,0.5,1.0
  
  # Annotate with markers
  xenium_process annotate --input merged.zarr --inplace --markers markers.csv
  
  # Differential analysis between groups
  xenium_process differential --input merged.zarr --output-dir results/ \\
      --groupby status --compare-groups HIV,NEG
  
  # Full pipeline (separate files)
  xenium_process concat --input samples.csv --output step1_concat.zarr
  xenium_process normalize --input step1_concat.zarr --output step2_normalized.zarr
  xenium_process cluster --input step2_normalized.zarr --output step3_clustered.zarr
  xenium_process annotate --input step3_clustered.zarr --output step4_annotated.zarr
  xenium_process differential --input step4_annotated.zarr --output-dir results/
        """
    )
    
    # Create subparsers
    subparsers = parser.add_subparsers(
        dest='command',
        help='Available commands',
        required=True
    )
    
    # Add concat subcommand
    concat_parser = subparsers.add_parser(
        'concat',
        help='Concatenate multiple Xenium .zarr files',
        description='Join multiple Xenium spatial datasets into a single .zarr file'
    )
    concat.add_arguments(concat_parser)
    concat_parser.set_defaults(func=concat.main)
    
    # Add normalize subcommand
    normalize_parser = subparsers.add_parser(
        'normalize',
        help='Normalize and preprocess data',
        description='Perform QC, filtering, normalization, and feature selection'
    )
    normalize.add_arguments(normalize_parser)
    normalize_parser.set_defaults(func=normalize.main)
    
    # Add cluster subcommand
    cluster_parser = subparsers.add_parser(
        'cluster',
        help='Perform clustering analysis',
        description='Run PCA, compute neighbors, UMAP, and Leiden clustering'
    )
    cluster.add_arguments(cluster_parser)
    cluster_parser.set_defaults(func=cluster.main)
    
    # Add annotate subcommand
    annotate_parser = subparsers.add_parser(
        'annotate',
        help='Annotate cell types',
        description='Perform cell type annotation using marker genes and/or ULM scoring'
    )
    annotate.add_arguments(annotate_parser)
    annotate_parser.set_defaults(func=annotate.main)
    
    # Add differential subcommand
    differential_parser = subparsers.add_parser(
        'differential',
        help='Differential expression analysis',
        description='Perform differential analysis between groups or find cluster markers'
    )
    differential.add_arguments(differential_parser)
    differential_parser.set_defaults(func=differential.main)
    
    return parser


def main():
    """Main entry point for the CLI."""
    # Setup logging
    setup_logging()
    
    # Parse arguments
    parser = create_parser()
    args = parser.parse_args()
    
    # Execute the appropriate command
    try:
        args.func(args)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

