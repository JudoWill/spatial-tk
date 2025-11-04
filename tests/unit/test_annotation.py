"""
Unit tests for annotation module.
"""

import pytest
import pandas as pd
from xenium_process.core import annotation


def test_load_marker_genes(test_markers_csv):
    """Test loading marker genes from CSV."""
    markers = annotation.load_marker_genes(str(test_markers_csv))
    
    # Check that markers were loaded
    assert isinstance(markers, dict)
    assert len(markers) > 0
    
    # Check structure
    for cell_type, genes in markers.items():
        assert isinstance(cell_type, str)
        assert isinstance(genes, list)
        assert len(genes) > 0


def test_annotate_with_markers(mock_adata_with_clusters, mock_markers):
    """Test cell type annotation with markers."""
    adata = mock_adata_with_clusters
    
    # Run annotation with tmin=1 for testing with small marker set
    cluster_key = 'leiden_res0p5'
    annotation_key = 'cell_type'
    
    adata = annotation.annotate_with_markers(
        adata,
        mock_markers,
        cluster_key=cluster_key,
        annotation_key=annotation_key,
        tmin=1
    )
    
    # Check that annotation was added
    assert annotation_key in adata.obs.columns
    
    # Check that all cells have an annotation
    assert adata.obs[annotation_key].notna().all()


def test_run_differential_expression(mock_adata_with_clusters):
    """Test differential expression analysis."""
    adata = mock_adata_with_clusters
    
    # Add normalized data
    import scanpy as sc
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    
    cluster_key = 'leiden_res0p5'
    adata = annotation.run_differential_expression(adata, cluster_key)
    
    # Check that differential expression was computed
    rank_key = f"rank_genes_{cluster_key}"
    assert rank_key in adata.uns


def test_run_differential_expression_resume(mock_adata_with_clusters):
    """Test that resume skips already computed DE."""
    adata = mock_adata_with_clusters
    
    # Add normalized data
    import scanpy as sc
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    
    cluster_key = 'leiden_res0p5'
    
    # Run once
    adata = annotation.run_differential_expression(adata, cluster_key)
    
    # Run again with resume
    adata = annotation.run_differential_expression(adata, cluster_key, resume=True)
    
    # Should still have results
    rank_key = f"rank_genes_{cluster_key}"
    assert rank_key in adata.uns

