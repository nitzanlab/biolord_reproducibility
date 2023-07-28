## chemCPA reproducibility
1. Clone the chemCPA project from `https://github.com/theislab/chemCPA` (commit `a4a4ded0c3b949c64ff1ea51033be1b7c301c36b`; see `chemCPA` modifications)
2. `cd chemCPA` + `python setup.py install`
3. Download the sci-Plex 3 data to from [figshare](https://figshare.com/ndownloader/files/39324305) to `./data/perturbations/sciplex3`
4. Update  `chemCPA/paths.py` with the correct `PROJECT_DIR`.
5. Follow the instruction at `chemCPA` to create `rdkit2D_embedding_lincs_trapnell.parquet` and place it under `/chemCPA/embeddings/rdkit/data/embeddings/`. 
6. In `config_sciplex_rdkit_nonpretrained.yml` or (`config_sciplex_rdkit_pretrained.yml` for the pretrained version): <br>
   1. change `training.save_dir` to the desired directory. <br>
   2. change `project_root_dir` to the chemCPA project location. <br>
   3. For the pretrained version, change `pretrained_model_path`. The pretrained model we used was kindly provided by the authors of chemCPA. 
7. cd `/path/to/chemCPA/experiments/baseline_comparison` and run chemCPA training:
   1. non-pretrained: <br>
      ```{bash}
      python /path/to/biolord_reproducibility/scripts/comparison_script/manual_seml_sweep_nonpretrained.py
      ```
   2. pretrained: <br>
      ```{bash}
      python /path/to/biolord_reproducibility/scripts/comparison_script/manual_seml_sweep_pretrained.py
      ```
8. In `get_results.py`, set paths:
   ```{python}
   CONFIG_PATH = "/path/to/config.yml"
   CHECKPOINT_HASH = ""
   CHEMCPA_PATH = "/path/to/chemCPA/notebooks/"
   ```
9. to evaluate and save results: 
   ```{bash}
   python path/to/biolord_reproducibility/scripts/comparison_script/get_results.py -- type x
   ```
   with `x=pretrained/nonpretrained`.
10. Performance `.csv` files will be saved at `

### `chemCPA` modifications

1. Fix paths in `chemCPA/chemCPA/paths.py`
2. Add to `chemCPA/notebooks/utils.py`:
   1. Line 8: <br>
   ```
   import sys
   sys.path.append('/path/to/chemCPA')]
   from chemCPA.data import SubDataset, canonicalize_smiles, drug_names_to_once_canon_smiles
    ```
   2. line 466: `dataset.use_drugs_idx` as an argument to `compute_prediction`.
5. In `chemCPA/chemCPA/experiments_run.py`, change line 425 to `file_name = f"checkpoint.pt"` (your desired checkpoint name).
