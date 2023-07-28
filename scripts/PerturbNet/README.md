## PerturbNet reproducibility
1. Clone the PerturbNet project from `https://github.com/welch-lab/PerturbNet` (commit `d719212aab1eb7cc56d7413e611cdb67987f1aec`; In [pytorch_scvi/scvi_generate_z.p](https://github.com/welch-lab/PerturbNet/blob/main/pytorch_scvi/scvi_generate_z.py) change `model.model` to `model.module`).
2. Download the sci-Plex 3 data to from [figshare](https://figshare.com/ndownloader/files/39324305) to `./data/perturbations/sciplex3`.
3. Download the PerturbNet required files from [here](https://www.dropbox.com/sh/3nk6qk1653h2y1v/AACplqpkgt3LH9_JDrsk-Hg4a?dl=0) to `./data/perturbations/sciplex3`.
4. Run `1_perturbations_sciplex3_perturbnet-preprocess.ipynb` to create the reference `adata_train.h5ad`.
5. Update paths in  `perturbnet_scvi.py` and train an scVI model.
6. Update paths in  `perturbnet_cINN.py` and train the ciNN model.
7. Run `2_perturbations_sciplex3_perturbnet-analysis.ipynb`
