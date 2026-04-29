#!/usr/bin/env python3
"""
Generate ROI subset zarrs from a single SpatialData zarr input.
"""

from __future__ import annotations

import argparse
import logging
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import spatialdata as sd
from spatialdata import bounding_box_query, transform

from spatial_tk.core.data_io import save_spatial_data
from spatial_tk.utils.helpers import get_table, set_table


def setup_logging(verbose: bool) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def iter_windows(
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
    frac: float,
    step_frac: float,
):
    width = (x_max - x_min) * frac
    height = (y_max - y_min) * frac
    if width <= 0 or height <= 0:
        return

    step_x = max(width * step_frac, 1e-6)
    step_y = max(height * step_frac, 1e-6)

    x_starts = np.arange(x_min, x_max - width + step_x, step_x)
    y_starts = np.arange(y_min, y_max - height + step_y, step_y)
    for xs in x_starts:
        for ys in y_starts:
            yield xs, xs + width, ys, ys + height


def iou(a: dict, b: dict) -> float:
    x0 = max(float(a["x_min"]), float(b["x_min"]))
    y0 = max(float(a["y_min"]), float(b["y_min"]))
    x1 = min(float(a["x_max"]), float(b["x_max"]))
    y1 = min(float(a["y_max"]), float(b["y_max"]))
    if x1 <= x0 or y1 <= y0:
        return 0.0
    inter = (x1 - x0) * (y1 - y0)
    area_a = (float(a["x_max"]) - float(a["x_min"])) * (float(a["y_max"]) - float(a["y_min"]))
    area_b = (float(b["x_max"]) - float(b["x_min"])) * (float(b["y_max"]) - float(b["y_min"]))
    denom = area_a + area_b - inter
    return float(inter / denom) if denom > 0 else 0.0


def get_xy_for_scanning(
    sdata: sd.SpatialData,
    spatial_key: str,
    target_coordinate_system: str,
) -> tuple[np.ndarray, np.ndarray]:
    adata = get_table(sdata)
    if adata is None:
        raise ValueError("No table found in SpatialData")

    spatial_attrs = adata.uns.get("spatialdata_attrs", {})
    annotated_region = spatial_attrs.get("region")
    if isinstance(annotated_region, (list, tuple)):
        annotated_region = annotated_region[0] if annotated_region else None

    if annotated_region in getattr(sdata, "shapes", {}):
        shape_element = transform(
            sdata.shapes[annotated_region],
            to_coordinate_system=target_coordinate_system,
        )
        centroids = shape_element.geometry.centroid
        return centroids.x.to_numpy(), centroids.y.to_numpy()

    if spatial_key not in adata.obsm:
        raise KeyError(f"Spatial key '{spatial_key}' missing in table.obsm")

    coords = np.asarray(adata.obsm[spatial_key])
    if coords.ndim != 2 or coords.shape[1] < 2:
        raise ValueError(f"Unexpected coordinate shape: {coords.shape}")
    return coords[:, 0], coords[:, 1]


def build_candidates(
    sample_name: str,
    source_path: Path,
    x: np.ndarray,
    y: np.ndarray,
    min_cells: int,
    max_cells: int,
    target_cells: int,
    window_fracs: list[float],
    step_frac: float,
    max_candidates: int,
) -> list[dict]:
    x_min, x_max = float(np.min(x)), float(np.max(x))
    y_min, y_max = float(np.min(y)), float(np.max(y))
    records: list[dict] = []
    idx = 0
    for frac in window_fracs:
        for wx0, wx1, wy0, wy1 in iter_windows(x_min, x_max, y_min, y_max, frac, step_frac):
            mask = (x >= wx0) & (x <= wx1) & (y >= wy0) & (y <= wy1)
            n_cells = int(mask.sum())
            if n_cells < min_cells or n_cells > max_cells:
                continue
            idx += 1
            records.append(
                {
                    "sample": sample_name,
                    "source_path": str(source_path),
                    "candidate_id": f"{sample_name}_cand_{idx:04d}",
                    "x_min": float(wx0),
                    "x_max": float(wx1),
                    "y_min": float(wy0),
                    "y_max": float(wy1),
                    "x_center": float((wx0 + wx1) / 2.0),
                    "y_center": float((wy0 + wy1) / 2.0),
                    "window_frac": float(frac),
                    "n_cells": int(n_cells),
                    "score": float(-abs(n_cells - target_cells) + frac * 100.0),
                }
            )

    records = sorted(records, key=lambda r: (r["score"], r["n_cells"]), reverse=True)
    return records[:max_candidates]


def export_roi(
    sdata: sd.SpatialData,
    row: dict,
    output_path: Path,
    coordinate_system: str,
    spatial_key: str,
    min_cells: int,
    max_cells: int,
    overwrite: bool,
    sample_name: str,
    status: str,
    location: str,
) -> tuple[bool, str, int]:
    roi = bounding_box_query(
        sdata,
        axes=("x", "y"),
        min_coordinate=[float(row["x_min"]), float(row["y_min"])],
        max_coordinate=[float(row["x_max"]), float(row["y_max"])],
        target_coordinate_system=coordinate_system,
        filter_table=True,
    )
    table = get_table(roi)
    if table is None:
        return False, "No table in ROI result", -1
    if spatial_key not in table.obsm:
        return False, f"Missing obsm['{spatial_key}']", int(table.n_obs)

    if "sample" not in table.obs:
        table.obs["sample"] = sample_name
    if status and "status" not in table.obs:
        table.obs["status"] = status
    if location and "location" not in table.obs:
        table.obs["location"] = location
    set_table(roi, table)

    n_cells = int(table.n_obs)
    if n_cells < min_cells or n_cells > max_cells:
        return False, f"ROI cell count out of range: {n_cells}", n_cells

    if output_path.exists():
        if overwrite:
            shutil.rmtree(output_path)
        else:
            return False, f"Output exists: {output_path}", n_cells

    try:
        save_spatial_data(roi, output_path, overwrite=overwrite)
        return True, "ok", n_cells
    except Exception as exc:
        return False, f"Failed to save ROI: {exc}", n_cells


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate ROI subset zarrs from one input zarr.")
    parser.add_argument("--input-zarr", type=Path, required=True, help="Input SpatialData .zarr path.")
    parser.add_argument("--output-dir", type=Path, required=True, help="Output root directory.")
    parser.add_argument("--sample-name", help="Sample name for outputs; defaults to input stem.")
    parser.add_argument("--status", default="", help="Optional status metadata to add (e.g. HIV/NEG).")
    parser.add_argument("--location", default="", help="Optional location metadata to add.")
    parser.add_argument("--n-rois", type=int, default=5, help="Number of ROI zarrs to export.")
    parser.add_argument("--min-cells", type=int, default=1000, help="Minimum cells per ROI.")
    parser.add_argument("--max-cells", type=int, default=5000, help="Maximum cells per ROI.")
    parser.add_argument("--target-cells", type=int, default=2500, help="Target cells for ranking.")
    parser.add_argument(
        "--window-fracs",
        default="0.15,0.2,0.25,0.3,0.35,0.4",
        help="Comma-separated bbox fractions of sample span.",
    )
    parser.add_argument("--step-frac", type=float, default=0.35, help="Sliding step fraction.")
    parser.add_argument("--max-candidates", type=int, default=300, help="Max candidate windows to keep.")
    parser.add_argument("--max-attempts", type=int, default=120, help="Max windows attempted for export.")
    parser.add_argument(
        "--target-coordinate-system",
        default="global",
        help="Coordinate system for bbox generation/query.",
    )
    parser.add_argument("--spatial-key", default="spatial", help="Coordinate key in table.obsm.")
    parser.add_argument("--min-center-distance", type=float, default=250.0, help="Min distance between ROI centers.")
    parser.add_argument("--max-iou", type=float, default=0.4, help="Max overlap between accepted ROIs.")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing output zarrs.")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logs.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    setup_logging(args.verbose)

    if not args.input_zarr.exists():
        raise FileNotFoundError(f"Input zarr not found: {args.input_zarr}")

    sample_name = args.sample_name or args.input_zarr.stem.replace(".zarr", "")
    window_fracs = [float(v.strip()) for v in args.window_fracs.split(",") if v.strip()]

    out_root = args.output_dir
    roi_dir = out_root / "rois"
    out_root.mkdir(parents=True, exist_ok=True)
    roi_dir.mkdir(parents=True, exist_ok=True)

    logging.info("Loading source zarr: %s", args.input_zarr)
    sdata = sd.read_zarr(args.input_zarr)
    if hasattr(sdata, "points") and sdata.points:
        logging.info("Dropping points in memory for faster bounding_box_query")
        sdata.points = {}

    x, y = get_xy_for_scanning(
        sdata=sdata,
        spatial_key=args.spatial_key,
        target_coordinate_system=args.target_coordinate_system,
    )
    candidates = build_candidates(
        sample_name=sample_name,
        source_path=args.input_zarr,
        x=x,
        y=y,
        min_cells=args.min_cells,
        max_cells=args.max_cells,
        target_cells=args.target_cells,
        window_fracs=window_fracs,
        step_frac=args.step_frac,
        max_candidates=args.max_candidates,
    )
    if not candidates:
        raise RuntimeError("No candidate ROI windows found. Try wider window-fracs or different min/max cells.")

    attempts: list[dict] = []
    accepted: list[dict] = []
    tries = 0
    for row in candidates:
        if len(accepted) >= args.n_rois or tries >= args.max_attempts:
            break
        tries += 1

        if any(
            np.hypot(float(row["x_center"]) - float(r["x_center"]), float(row["y_center"]) - float(r["y_center"]))
            < args.min_center_distance
            or iou(row, r) > args.max_iou
            for r in accepted
        ):
            continue

        roi_id = f"{sample_name}_roi_{len(accepted) + 1:02d}"
        out_path = roi_dir / f"{roi_id}.zarr"
        ok, msg, n_cells = export_roi(
            sdata=sdata,
            row=row,
            output_path=out_path,
            coordinate_system=args.target_coordinate_system,
            spatial_key=args.spatial_key,
            min_cells=args.min_cells,
            max_cells=args.max_cells,
            overwrite=args.overwrite,
            sample_name=sample_name,
            status=args.status,
            location=args.location,
        )
        record = dict(row)
        record.update(
            {
                "roi_id": roi_id,
                "roi_path": str(out_path),
                "export_ok": bool(ok),
                "export_message": msg,
                "export_n_cells": int(n_cells),
                "status": args.status,
                "location": args.location,
            }
        )
        attempts.append(record)

        if ok:
            accepted.append(record)
            logging.info("Exported %s (%d cells)", roi_id, n_cells)
        else:
            logging.warning("Skipped %s: %s", roi_id, msg)

    attempts_df = pd.DataFrame(attempts)
    attempts_csv = out_root / f"{sample_name}_roi_attempts.csv"
    attempts_df.to_csv(attempts_csv, index=False)

    manifest_cols = ["sample", "path", "status", "location", "roi_id", "n_cells"]
    if accepted:
        manifest = pd.DataFrame(
            [
                {
                    "sample": row["roi_id"],
                    "path": row["roi_path"],
                    "status": row["status"],
                    "location": row["location"],
                    "roi_id": row["roi_id"],
                    "n_cells": row["export_n_cells"],
                }
                for row in accepted
            ]
        )
    else:
        manifest = pd.DataFrame(columns=manifest_cols)

    manifest_csv = out_root / f"{sample_name}_roi_samples.csv"
    manifest.to_csv(manifest_csv, index=False)

    logging.info("Wrote attempts: %s", attempts_csv)
    logging.info("Wrote manifest: %s (%d rows)", manifest_csv, len(manifest))
    if len(accepted) < args.n_rois:
        logging.warning("Exported %d/%d ROIs. Consider relaxing constraints.", len(accepted), args.n_rois)


if __name__ == "__main__":
    main()
