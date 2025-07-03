This is a Python and Nextflow repository for predicting whether TCRs (T-cell receptors) bind to p:MHCs (peptide:MHC complexes) using AlphaFold3. Please follow these guidelines when contributing:

## Code Standards

### Required Before Each Commit
- Run `black .` from the cloned project root (using the `tcrtrifold-experiments` conda environment) before committing any changes to ensure proper code formatting
- This will run black on all Python files to maintain consistent style

### Development Flow
- Test: `./scripts/test/<script_name>.sh`

## Repository Structure
- `data/`: Parquet files, AF3 inference outputs, and other data associated with triads (TCRs in complex with p:MHCs) and p:MHC complexes alone, organized on a per-dataset basis.
- `envs/`: Conda environments for different purposes in this repo. `env_runner.yaml` is the primary package code containing the package described by 'pyproject.toml' and its dependencies, often used by nextflow pipelines kicked off using the `nf-core` environment.
- `notebooks/`: Jupyter notebooks for analyzing results and making figures
- `results/`: Persistent storage for outputs of jupyter notebooks
- `scripts/`: Primary entrypoint for interacting with this repo. Contains slurm scripts (which should be run as bash scripts in the runner environment) that kick off the nextflow pipelines.
- `src/tcrtrifold/`: The python package code used by the pipelines in this repo, organized using the "src" layout.
- `workflows/`: Nextflow pipelines and their associated Python scripts for running AlphaFold3, formatting data, and extracting features.

## Key Guidelines
1. Follow Python best practices and idiomatic patterns
2. Use Numpy-style docstrings for Python methods and classes
3. Maintain existing code structure and organization
4. If issues are ambiguous, ask clarifying questions before attempting to solve them
5. While other datasets are present in this repo, `data/test` is the dataset intended to be used by the runner environment. Pipelines should be written in a way that works on test data but is generic enough to be run on other datasets.
