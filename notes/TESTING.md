# Testing Guide

## Managing Disk Space During Tests

The functional tests create large `.zarr` files that can fill up `/tmp`. Here are the strategies in place to manage disk space:

### Automatic Solutions (Already Configured)

1. **Custom Temp Directory**: Tests use `.pytest_tmp/` in the project root instead of `/tmp`
   - Configure via `make test` (uses `--basetemp=.pytest_tmp`)
   - Or set `PYTEST_BASETEMP` environment variable

2. **Automatic Cleanup**: 
   - Fixtures automatically clean up after each test
   - `tmp_zarr_cleanup` fixture aggressively removes `.zarr` files
   - Garbage collection runs after each test

3. **No Retention**: pytest is configured to not retain any temp files

### Manual Solutions

#### Option 1: Use Custom Temp Directory (Recommended)

Run tests with a custom location on a partition with more space:

```bash
# In Makefile (already configured)
make test

# Or manually
pytest --basetemp=/path/to/larger/partition/.pytest_tmp

# Or set environment variable
export PYTEST_BASETEMP=/path/to/larger/partition/.pytest_tmp
pytest
```

#### Option 2: Clean Between Test Runs

```bash
# Clean test temp files
make clean-test

# Or manually
rm -rf .pytest_tmp/ .pytest_cache/
```

#### Option 3: Run Tests in Smaller Batches

```bash
# Run only unit tests (no large files)
make test-unit

# Run functional tests one at a time with cleanup
pytest tests/functional/test_concat_command.py && make clean-test
pytest tests/functional/test_normalize_command.py && make clean-test
pytest tests/functional/test_full_pipeline.py && make clean-test
```

#### Option 4: Use Smaller Test Data

Reduce the size of test datasets by creating even smaller subsampled data:

```bash
python scripts/create_test_data.py \
  --input-csv example.csv \
  --output-dir tests/test_data \
  --n-cells 100  # Even smaller (default is 500)
```

#### Option 5: Skip Functional Tests

```bash
# Run only unit tests (fast, no large files)
pytest tests/unit/ -v

# Or with make
make test-unit
```

#### Option 6: Monitor Disk Usage

```bash
# Check disk usage before tests
df -h /tmp
df -h .pytest_tmp

# Monitor during test run (in another terminal)
watch -n 1 'du -sh .pytest_tmp 2>/dev/null || echo "No temp files yet"'
```

### Pytest Configuration

The `pyproject.toml` includes these settings:

```toml
[tool.pytest.ini_options]
tmp_path_retention_count = 0
tmp_path_retention_policy = "none"
```

This ensures pytest doesn't retain temporary directories after test runs.

### Best Practices

1. **Always use `tmp_zarr_cleanup` fixture** for functional tests that create `.zarr` files
2. **Run `make clean-test`** before large test runs
3. **Monitor disk space** if running full test suite repeatedly
4. **Use test markers** to selectively run tests:

```bash
# Mark slow tests
@pytest.mark.slow
def test_full_pipeline():
    ...

# Skip slow tests
pytest -m "not slow"
```

### Troubleshooting

**Problem**: Still running out of space

**Solutions**:
1. Check if `.pytest_tmp/` is actually being used:
   ```bash
   ls -lah .pytest_tmp/
   ```

2. Verify basetemp is set:
   ```bash
   pytest --basetemp=.pytest_tmp --collect-only
   ```

3. Manually specify a location with more space:
   ```bash
   pytest --basetemp=$HOME/.pytest_tmp
   ```

4. Clean up aggressively during test run by adding to `conftest.py`:
   ```python
   @pytest.fixture(autouse=True, scope="function")
   def cleanup_zarr_after_each_test(request):
       yield
       import shutil
       import gc
       # Find and remove all .zarr in temp directories
       for zarr in Path("/tmp").glob("**/*.zarr"):
           shutil.rmtree(zarr, ignore_errors=True)
       gc.collect()
   ```

**Problem**: Tests are slow due to cleanup

**Solution**: The cleanup is necessary to prevent disk issues. Consider:
- Running tests in parallel with `pytest-xdist`: `pytest -n auto`
- Using ramdisk for temp (if available)
- Running unit and functional tests separately

