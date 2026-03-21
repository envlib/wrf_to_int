# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

wrf_to_int converts WRF (Weather Research and Forecasting) output files to WPS intermediate files, specifically for metgrid.exe input when creating subdomains in different coordinate systems.

## Development Commands

```bash
# Install/sync dependencies (uses UV)
uv sync

# Run all tests
uv run pytest

# Run a single test file
uv run pytest wrf_to_int/tests/test_foo.py

# Run a single test
uv run pytest wrf_to_int/tests/test_foo.py::test_name

# Run tests with coverage
uv run pytest --cov --cov-report=xml

# Lint and typecheck (all three must pass in CI)
uv run ruff check .
uv run black --check --diff .
uv run mypy --install-types --non-interactive wrf_to_int

# Auto-format
uv run black .
uv run ruff check --fix .

# Build docs locally
uv run mkdocs serve
```

## Code Style

- **Line length:** 120 characters (black and ruff)
- **Black:** skip-string-normalization enabled (use single quotes)
- **Ruff:** relative imports are banned; use absolute imports from `wrf_to_int`
- **Ruff:** unused imports (F401) are not auto-removed
- **Target:** Python 3.10+ (tested on 3.10, 3.11, 3.12)

## Architecture

- `wrf_to_int/` — main package; version defined in `__init__.py`
- `wrf_to_int/tests/` — pytest tests
- Build backend: hatchling; version sourced from `wrf_to_int/__init__.py`
- Tests are excluded from sdist builds
- Docs use mkdocs-material with mkdocstrings (Google-style docstrings)
