
## Extracting raw cognate triads from IEDB

The raw extracted data from IEDB is already available in this repo in `data/iedb_<class>/raw`

Alternatively, if you're interested in re-running our data extraction, clone our fork of the IEDB_IMMREP data repo and run the data extraction script.
```bash
git clone https://github.com/ljwoods2/IEDB_IMMREP.git
cd IEDB_IMMREP
git checkout new-categories
pip install .
chmod +x run.sh
./run.sh
```

## Re-running pipelines

```
conda env create -f nf-core.yaml
conda env create -f tcrtrifold-experiments.yaml

sbatch scripts/<dataset>/<pipeline_step>.sh
```


## Gemini location

Cloned repo and full inference data is available in 

```
/tgen_labs/altin/alphafold3/workspace/tcrtrifold-experiments
```