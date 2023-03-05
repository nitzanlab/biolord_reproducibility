# biolord's reproducibility repository

This code contains the code to reproduce the results shown in "_Biological representation disentanglement of single-cell data_".

## Data availability
1. perturbation data (sci-Plex 3, [Srivatsan et al. 2020](https://www.science.org/doi/10.1126/science.aax6234)):
   1. Anndata file used for the initial [pre-processing notebook](https://github.com/nitzanlab/biolord_reproducibility/tree/main/notebooks/perturbation/1_perturbation_preprocessing.ipynb): [`sciplex_complete_middle_subset.h5ad`](https://f003.backblazeb2.com/file/chemCPA-datasets/sciplex_complete_middle_subset.h5ad); provided by  [Hetzel et al.](https://openreview.net/forum?id=vRrFVHxFiXJ).
   2. Processed anndata file used for the [analysis](https://github.com/nitzanlab/biolord_reproducibility/tree/main/notebooks/perturbation/2_perturbation_evaluation.ipynb): [`sciplex3_biolord.h5ad`](https://figshare.com/ndownloader/files/39324305).
2. spatio-temporal infection dataset ([Afriat et al.](https://www.nature.com/articles/s41586-022-05406-5)):
   1. Data used for the initial [pre-processing notebook](https://github.com/nitzanlab/biolord_reproducibility/tree/main/notebooks/spatio-temporal-infection/1_spatio-temporal-infection_preprocessing.ipynb): [10.5281/zenodo.7081862](https://zenodo.org/record/7081863#.Y5jhCuxBxBw) (provided by [Afriat et al.](https://www.nature.com/articles/s41586-022-05406-5))
   2. Processed anndata file used for the [infection analysis](https://github.com/nitzanlab/biolord_reproducibility/tree/main/notebooks/spatio-temporal-infection/2_spatio-temporal-infection_state.ipynb): [`adata_infected.h5ad`](https://figshare.com/ndownloader/files/39375713). 
   3. Processed anndata file used for the [abortive classification analysis](https://github.com/nitzanlab/biolord_reproducibility/tree/main/notebooks/spatio-temporal-infection/3_spatio-temporal-infection_abortive.ipynb): [`adata_abortive.h5ad`](https://figshare.com/ndownloader/files/39375752).

      

