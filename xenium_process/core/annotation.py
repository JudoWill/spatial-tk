#!/usr/bin/env python3
"""
Cell type annotation and enrichment scoring functions.

This module handles marker-based cell type annotation using decoupler,
differential expression analysis, and MLM score calculation for
multiple pathway/TF resources.
"""

import logging
from typing import Dict, List, Optional

import anndata as ad
import decoupler as dc
import pandas as pd
import scanpy as sc


def load_marker_genes(marker_path: str) -> Dict[str, List[str]]:
    """
    Load marker genes from CSV file.
    
    Args:
        marker_path: Path to CSV file with columns: cell_type, gene
        
    Returns:
        Dictionary mapping cell type to list of marker genes
    """
    logging.info(f"Loading marker genes from {marker_path}")
    
    df = pd.read_csv(marker_path)
    
    if not all(col in df.columns for col in ["cell_type", "gene"]):
        raise ValueError("Marker CSV must have 'cell_type' and 'gene' columns")
    
    # Group by cell type
    markers = df.groupby("cell_type")["gene"].apply(list).to_dict()
    
    total_markers = sum(len(genes) for genes in markers.values())
    logging.info(f"Loaded {len(markers)} cell types with {total_markers} total marker genes")
    
    return markers


def get_panglao_markers(
    organism: str = "human",
    min_sensitivity: float = 0.5,
    canonical_only: bool = True
) -> pd.DataFrame:
    """
    Get PanglaoDB markers with filtering.
    
    Best practices per: 
    https://decoupler.readthedocs.io/en/latest/notebooks/scell/rna_sc.html#panglaodb
    
    Args:
        organism: Organism name ('human' or 'mouse')
        min_sensitivity: Minimum sensitivity threshold (0-1)
        canonical_only: If True, only use canonical markers
        
    Returns:
        DataFrame with columns: source (cell_type), target (gene)
    """
    logging.info(f"Loading PanglaoDB markers (organism={organism}, min_sensitivity={min_sensitivity})")
    
    markers = dc.op.resource("PanglaoDB", organism=organism)
    
    # Apply filters
    filters = markers[organism].astype(bool)
    
    if canonical_only:
        filters &= markers["canonical_marker"].astype(bool)
    
    filters &= (markers[f"{organism}_sensitivity"].astype(float) > min_sensitivity)
    
    markers = markers[filters]
    
    # Remove duplicates
    markers = markers[~markers.duplicated(["cell_type", "genesymbol"])]
    
    # Rename columns to decoupler format
    markers = markers.rename(columns={"cell_type": "source", "genesymbol": "target"})
    
    logging.info(f"  Filtered to {len(markers)} PanglaoDB markers")
    
    return markers[["source", "target"]]


def annotate_with_markers(
    adata: ad.AnnData,
    markers: Dict[str, List[str]],
    cluster_key: str = "leiden",
    annotation_key: str = "cell_type",
    resume: bool = False,
    tmin: int = 2
) -> ad.AnnData:
    """
    Annotate clusters with cell types based on marker gene expression using
    decoupler's multivariate linear model (MLM) approach.
    
    This method uses enrichment analysis to test if marker gene collections
    are enriched in cells, similar to the approach in the Scverse tutorial.
    
    Args:
        adata: AnnData object
        markers: Dictionary mapping cell type to list of marker genes
        cluster_key: Key in adata.obs containing cluster assignments
        annotation_key: Key name for storing cell type annotations
        resume: If True, skip if annotation already exists
        tmin: Minimum number of targets per source (default: 2)
        
    Returns:
        AnnData object with cell type annotations added
    """
    if resume and annotation_key in adata.obs.columns:
        logging.info(f"Cell type annotation already exists (resuming)")
        return adata
    
    logging.info(f"Annotating cell types using decoupler MLM (cluster_key={cluster_key})")
    
    # Check which marker genes are present
    all_marker_genes = set()
    for genes in markers.values():
        all_marker_genes.update(genes)
    
    missing_genes = all_marker_genes - set(adata.var_names)
    if missing_genes:
        logging.info(f"Note: {len(missing_genes)} marker genes not found in dataset")
    
    # Convert markers dictionary to DataFrame format expected by decoupler
    # Format: columns 'source' (cell_type) and 'target' (gene)
    marker_rows = []
    for cell_type, genes in markers.items():
        for gene in genes:
            marker_rows.append({"source": cell_type, "target": gene})
    
    marker_df = pd.DataFrame(marker_rows)
    
    # Add weight column (all weights = 1)
    marker_df["weight"] = 1
    
    logging.info(f"Running MLM with {len(marker_df)} marker gene entries across {len(markers)} cell types")
    
    # Run multivariate linear model
    # This calculates enrichment scores for each cell type in each cell
    dc.mt.mlm(adata, net=marker_df, verbose=False, tmin=tmin)
    
    # Extract the MLM scores from adata.obsm
    # This creates a new AnnData-like object with cells x cell_types
    acts = dc.pp.get_obsm(adata, "score_mlm")
    
    # For each cluster, find the cell type with highest enrichment score
    # Use decoupler's rankby_group to get top scoring cell type per cluster
    enr = dc.tl.rankby_group(acts, groupby=cluster_key)
    
    # Get the top cell type (highest stat) for each cluster
    # Filter to positive stats only (stat > 0)
    annotation_dict = (
        enr[enr["stat"] > 0]
        .groupby("group", observed=True)
        .head(1)
        .set_index("group")["name"]
        .to_dict()
    )
    
    # Handle clusters that may not have positive enrichment scores
    all_clusters = adata.obs[cluster_key].unique()
    for cluster in all_clusters:
        if cluster not in annotation_dict:
            annotation_dict[cluster] = "Unknown"
    
    # Map cluster annotations to cells
    adata.obs[annotation_key] = adata.obs[cluster_key].map(annotation_dict)
    adata.obs[annotation_key] = adata.obs[annotation_key].astype("category")
    
    # Log annotation summary
    annotation_counts = adata.obs[annotation_key].value_counts()
    logging.info("Cell type annotation summary:")
    for cell_type, count in annotation_counts.items():
        logging.info(f"  {cell_type}: {count} cells")
    
    return adata


def calculate_mlm_scores(
    adata: ad.AnnData,
    use_panglao: bool = True,
    panglao_min_sensitivity: float = 0.5,
    tmin: int = 5,
    resume: bool = False
) -> ad.AnnData:
    """
    Pre-calculate MLM scores for multiple decoupler resources.
    
    Resources include:
    - hallmark: Hallmark gene sets
    - collectri: Transcription factor regulons
    - dorothea: TF activity inference
    - progeny: Pathway activity
    - PanglaoDB: Cell type markers (optional, filtered)
    
    Scores are stored in adata.obsm[f'score_mlm_{resource}']
    
    Args:
        adata: AnnData object with normalized data
        use_panglao: If True, include PanglaoDB markers
        panglao_min_sensitivity: Minimum sensitivity for PanglaoDB markers
        tmin: Minimum number of targets per source (default: 5)
        resume: If True, skip resources that already have scores
        
    Returns:
        AnnData object with MLM scores added to obsm
    """
    logging.info("Calculating MLM scores for pathway/TF resources")
    
    # Define resources using decoupler omnipath API (dc.op.*)
    resources = [
        ('hallmark', dc.op.hallmark(organism='human')),
        ('collectri', dc.op.collectri(organism='human')),
        ('dorothea', dc.op.dorothea(organism='human')),
        ('progeny', dc.op.progeny(organism='human'))
    ]
    
    if use_panglao:
        panglao = get_panglao_markers(
            organism="human",
            min_sensitivity=panglao_min_sensitivity,
            canonical_only=True
        )
        resources.append(('PanglaoDB', panglao))
    
    # Calculate MLM for each resource
    for name, resource in resources:
        obsm_key = f'score_mlm_{name}'
        
        if resume and obsm_key in adata.obsm:
            logging.info(f"  MLM scores for {name} already calculated (resuming)")
            continue
        
        logging.info(f"  Calculating MLM for {name}")
        
        try:
            # Run MLM
            dc.mt.mlm(data=adata, net=resource, tmin=tmin)
            
            # Store in named obsm key
            adata.obsm[obsm_key] = adata.obsm['score_mlm'].copy()
            
            # Get shape info
            n_sources = adata.obsm[obsm_key].shape[1] if len(adata.obsm[obsm_key].shape) > 1 else 1
            logging.info(f"    Calculated scores for {n_sources} sources")
            
        except Exception as e:
            logging.warning(f"  Failed to calculate MLM for {name}: {e}")
    
    logging.info("MLM score calculation complete")
    return adata


def run_differential_expression(
    adata: ad.AnnData,
    cluster_key: str,
    method: str = "wilcoxon",
    resume: bool = False
) -> ad.AnnData:
    """
    Run differential expression analysis to find marker genes for each cluster.
    
    Args:
        adata: AnnData object
        cluster_key: Key in adata.obs containing cluster assignments
        method: Statistical test to use (default: wilcoxon)
        resume: If True, skip if differential expression already computed
        
    Returns:
        AnnData object with differential expression results added
    """
    rank_key = f"rank_genes_{cluster_key}"
    
    if resume and "rank_genes_groups" in adata.uns and adata.uns.get("rank_genes_groups_key") == rank_key:
        logging.info(f"Differential expression already computed for {cluster_key} (resuming)")
        return adata
    
    logging.info(f"Running differential expression analysis for {cluster_key}")
    
    # Run rank_genes_groups
    sc.tl.rank_genes_groups(
        adata,
        groupby=cluster_key,
        method=method,
        use_raw=False,  # Use normalized data in .X
        key_added=rank_key,
        layer=None
    )
    
    # Store which key was used
    adata.uns["rank_genes_groups_key"] = rank_key
    
    n_clusters = adata.obs[cluster_key].nunique()
    logging.info(f"  Differential expression completed for {n_clusters} clusters")
    
    return adata


def save_differential_expression_results(
    adata: ad.AnnData,
    cluster_key: str,
    output_dir,
    n_genes: int = 100
) -> None:
    """
    Save differential expression results to CSV files.
    
    Args:
        adata: AnnData object with differential expression results
        cluster_key: Key in adata.obs containing cluster assignments
        output_dir: Directory to save output files
        n_genes: Number of top genes to save per cluster
    """
    rank_key = f"rank_genes_{cluster_key}"
    
    if rank_key not in adata.uns:
        logging.warning(f"  No differential expression results found for {cluster_key}")
        return
    
    logging.info(f"  Saving differential expression results for {cluster_key}")
    
    # Get the differential expression results as a DataFrame
    result = sc.get.rank_genes_groups_df(adata, group=None, key=rank_key)
    
    # Save all results
    de_dir = output_dir / "differential_expression"
    de_dir.mkdir(exist_ok=True)
    
    res_str = cluster_key.replace("leiden_res", "")
    
    # Save complete results
    all_results_path = de_dir / f"deg_all_clusters_res{res_str}.csv"
    result.to_csv(all_results_path, index=False)
    logging.info(f"    Saved all DE genes to {all_results_path}")
    
    # Save top N genes per cluster
    top_results_path = de_dir / f"deg_top{n_genes}_per_cluster_res{res_str}.csv"
    top_result = result.groupby('group').head(n_genes)
    top_result.to_csv(top_results_path, index=False)
    logging.info(f"    Saved top {n_genes} DE genes per cluster to {top_results_path}")

