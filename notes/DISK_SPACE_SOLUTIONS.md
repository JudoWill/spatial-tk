# Disk Space Solutions for Testing

## Quick Reference

### Default behavior (already configured):
```bash
make test          # Uses .pytest_tmp/ instead of /tmp
make clean-test    # Clean up test files
```

### Key Changes Made:

1. **✅ Makefile updated**: All test targets use `--basetemp=.pytest_tmp`
2. **✅ conftest.py updated**: Automatic cleanup fixtures added
3. **✅ pyproject.toml updated**: Pytest configured to not retain temp files
4. **✅ .gitignore updated**: Excludes .pytest_tmp/
5. **✅ Functional tests updated**: Use `tmp_zarr_cleanup` fixture

### To use a different location:

```bash
# Temporary (single run)
pytest --basetemp=/data/tmp

# Permanent (in your shell rc)
export PYTEST_BASETEMP=$HOME/tmp/pytest
```

### Emergency cleanup:

```bash
# Clean pytest temp files
make clean-test

# Clean everything
make clean

# Manual cleanup
rm -rf .pytest_tmp/ .pytest_cache/
```

## File Locations

- **Configuration**: `pyproject.toml` (lines 84-87)
- **Cleanup fixtures**: `tests/conftest.py` (lines 15-125)
- **Makefile targets**: `makefile` (lines 43-63, 60)
- **Full guide**: `TESTING.md`
