# Gemini Context: `seq-smith`

## Project Overview

This project is a Python library for sequence alignment, with the core logic implemented in Rust for performance. It provides several alignment strategies, including global, local, "local_global", and overlap alignment. The library is designed for bioinformatics applications.

The project consists of three main parts:
1.  A Python package that serves as the main user interface.
2.  A Rust extension, built with `pyo3` and `maturin`, which contains the primary, high-performance alignment algorithms.
3.  A C implementation of the alignment algorithms, which appears to be a reference or legacy version. It is compiled into a shared library (`libalign.so`) and can be accessed via `ctypes` in `c_align.py`.

The main, tested, and recommended interface is the Rust-based extension.

## Building and Running

### Rust Extension (Primary)

The project uses `maturin` to build the Rust extension.

- **Build and install:**
  ```bash
  pip install .
  ```

- **Run tests:**
  ```bash
  pytest
  ```

### C Implementation (Reference)

The C implementation can be built using the provided `Makefile`.

- **Build the shared library:**
  ```bash
  make
  ```
  This will create `libalign.so` in the root directory.

- **Clean build artifacts:**
  ```bash
  make clean
  ```

## Development Conventions

- **Linting:** The project uses `ruff` for Python linting and `pre-commit` to enforce code style.
  - To run linters on all files:
    ```bash
    pre-commit run --all-files
    ```
- **Testing:** Tests are written using `pytest` and are located in the `tests/` directory. The tests primarily target the Rust implementation.
- **Dependencies:** Python dependencies are managed in `pyproject.toml`. Development dependencies are listed in `requirements-dev.txt`.
- **CI:** The GitHub Actions workflow in `.github/workflows/lint.yaml` runs `pre-commit` on every push.
