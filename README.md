# biolord's reproducibility repository

This code contains the code to reproduce the results shown in ["_Biological representation disentanglement of single-cell data_"](https://www.biorxiv.org/content/10.1101/2023.03.05.531195) (2023).

## Data availability
1. perturbation data (sci-Plex 3, [Srivatsan et al. 2020](https://www.science.org/doi/10.1126/science.aax6234)):
   1. Anndata file used for the initial [pre-processing notebook](https://github.com/nitzanlab/biolord_reproducibility/tree/main/notebooks/perturbations/sciplex3/1_perturbations_sciplex3_preprocessing.ipynb): [`sciplex_complete_middle_subset.h5ad`](https://f003.backblazeb2.com/file/chemCPA-datasets/sciplex_complete_middle_subset.h5ad); provided by  [Hetzel et al.](https://openreview.net/forum?id=vRrFVHxFiXJ).
   2. Processed anndata file used for the [analysis](https://github.com/nitzanlab/biolord_reproducibility/tree/main/notebooks/perturbations/sciplex3/2_perturbations_sciplex3_evaluation.ipynb): [`sciplex3_biolord.h5ad`](https://figshare.com/ndownloader/files/39324305).
2. perturbation data (Perturb-seq (1-gene), [Adamson et al. 2016](https://doi.org/10.1016/j.cell.2016.11.048)):
   1. Data used for the initial [pre-processing notebook](https://github.com/nitzanlab/biolord_reproducibility/tree/main/notebooks/perturbations/norman/1_perturbations_norman_preprocessing.ipynb): Downloaded from [here](https://dataverse.harvard.edu/api/access/datafile/6154417) using GEARS within the notebook.
   2. Processed anndata file used for training a biolord model [analysis](https://github.com/nitzanlab/biolord_reproducibility/tree/scripts/biolord/adamson/base_experiment_adamson.py): [`adamson_biolord.h5ad`](https://figshare.com/articles/dataset/perturbseq_adamson/22344214) and [`adamson_single_biolord.h5ad`](https://figshare.com/articles/dataset/perturbseq_adamson_single/22344445).
3. perturbation data (Perturb-seq (2-gene), [Norman et al. 2019](https://doi.org/10.1126/science.aax4438)):
   1. Data used for the initial [pre-processing notebook](https://github.com/nitzanlab/biolord_reproducibility/tree/main/notebooks/perturbations/norman/1_perturbations_norman_preprocessing.ipynb): Download the data from [here](https://dataverse.harvard.edu/api/access/datafile/6894431). Move the uncompressed `norman2019.tar.gz` folder to `./data/perturbations/norman`. It should contain the subdirectory `data_pyg`. Move `essential_norman.pkl` and `go_essential_norman.csv` to `./norman`, So `./data/perturbations/norman` should contain `essential_norman.pkl, go_essential_norman.csv and norman2019/data_pyg/`.
   2. Processed anndata file used for training a biolord model [analysis](https://github.com/nitzanlab/biolord_reproducibility/tree/scripts/biolord/norman/base_experiment_norman.py): [`norman_biolord.h5ad`](https://figshare.com/articles/dataset/perturbseq_nornan/22344253) and [`norman_single_biolord.h5ad`](https://figshare.com/articles/dataset/pertrubseq_norman_single/22344427).
4. spatio-temporal infection dataset ([Afriat et al.](https://www.nature.com/articles/s41586-022-05406-5)):
   1. Data used for the initial [pre-processing notebook](https://github.com/nitzanlab/biolord_reproducibility/tree/main/notebooks/spatio-temporal-infection/1_spatio-temporal-infection_preprocessing.ipynb): [10.5281/zenodo.7081862](https://zenodo.org/record/7081863#.Y5jhCuxBxBw) (provided by [Afriat et al.](https://www.nature.com/articles/s41586-022-05406-5))
   2. Processed anndata file used for the [infection analysis](https://github.com/nitzanlab/biolord_reproducibility/tree/main/notebooks/spatio-temporal-infection/2_spatio-temporal-infection_state.ipynb): [`adata_infected.h5ad`](https://figshare.com/ndownloader/files/39375713). 
   3. Processed anndata file used for the [abortive classification analysis](https://github.com/nitzanlab/biolord_reproducibility/tree/main/notebooks/spatio-temporal-infection/3_spatio-temporal-infection_abortive.ipynb): [`adata_abortive.h5ad`](https://figshare.com/ndownloader/files/39375752).

      

