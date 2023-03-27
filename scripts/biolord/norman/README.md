## biolord reproducibility 

1. Download the data from [here](https://dataverse.harvard.edu/api/access/datafile/6894431). Move the uncompressed `norman2019.tar.gz` folder to `./data/perturbations/norman`. It should contain the subdirectory `data_pyg`. Move `essential_norman.pkl` and `go_essential_norman.csv` to `./norman`, So `./data/perturbations/norman` should contain `essential_norman.pkl, go_essential_norman.csv and norman2019/data_pyg/`. and 
2. Run [1_perturbations_norman_preprocessing.ipynb](https://github.com/nitzanlab/biolord_reproducibility/tree/main/notebooks/perturbations/norman/1_perturbations_norman_preprocessing.ipynb) or download the Norman datasets (`pertrubseq_norman` and `pertrubseq_norman_single`) from [figshare](https://figshare.com/projects/biolord_datasets/160085) to `./data/perturbations/norman/`

3. cd to the current folder:
    ```{bash}
   cd /biolord_reproducibility/scripts/biolord/norman
   ```
4. run sweep:
    ```{bash}
       python base_experiment_norman.py
   ```
   If you want to run over optimal parameters add the flag `--optimal_params True`