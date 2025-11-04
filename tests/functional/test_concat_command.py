"""
Functional tests for concat command.
"""

import pytest
import subprocess
import sys
from pathlib import Path


def test_concat_command(test_samples_csv, tmp_zarr_cleanup):
    """Test concat command with test data."""
    if not test_samples_csv.exists():
        pytest.skip("Test samples CSV not found")
    
    output_path = tmp_zarr_cleanup / "concat_output.zarr"
    
    result = subprocess.run([
        sys.executable, '-m', 'xenium_process.cli',
        'concat',
        '--input', str(test_samples_csv),
        '--output', str(output_path)
    ], capture_output=True, text=True)
    
    # Check that command succeeded
    assert result.returncode == 0, f"Command failed: {result.stderr}"
    
    # Check that output file was created
    assert output_path.exists()
    
    # Verify the output is valid
    import spatialdata as sd
    sdata = sd.read_zarr(output_path)
    assert sdata is not None
    
    from xenium_process.utils.helpers import get_table
    table = get_table(sdata)
    assert table is not None
    assert table.n_obs > 0


def test_concat_with_downsampling(test_samples_csv, tmp_zarr_cleanup):
    """Test concat command with downsampling."""
    if not test_samples_csv.exists():
        pytest.skip("Test samples CSV not found")
    
    output_path = tmp_zarr_cleanup / "concat_downsampled.zarr"
    
    result = subprocess.run([
        sys.executable, '-m', 'xenium_process.cli',
        'concat',
        '--input', str(test_samples_csv),
        '--output', str(output_path),
        '--downsample', '0.5'
    ], capture_output=True, text=True)
    
    assert result.returncode == 0
    assert output_path.exists()

