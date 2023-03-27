## biolord reproducibility 

1. Run [1_perturbations_adamson_preprocessing.ipynb](https://github.com/nitzanlab/biolord_reproducibility/tree/main/notebooks/perturbations/adamson/1_perturbations_adamson_preprocessing.ipynb) or download the Adamson datasets (`pertrubseq_adamson` and `pertrubseq_adamson_single`)from [figshare](https://figshare.com/projects/biolord_datasets/160085) to `./data/perturbations/adamson/`

2. cd to the current folder:
    ```{bash}
   cd /biolord_reproducibility/scripts/biolord/adamson
   ```
3. run sweep:
    ```{bash}
       python base_experiment_adamson.py
   ```
   If you want to run over optimal parameters add the flag `--optimal_params True`