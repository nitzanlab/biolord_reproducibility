# import the modules
import sys

import os

import torch

import scanpy as sc
import numpy as np
import pandas as pd
from scipy import sparse
import scvi

sys.path.insert(0, "PATH_TO_PerturbNet") 
from pytorch_scvi.distributions import *
from pytorch_scvi.scvi_generate_z import *

from perturbnet.drug_perturb.util import *
from perturbnet.drug_perturb.cinn.modules.flow import *
from perturbnet.drug_perturb.chemicalvae.chemicalVAE import *

if __name__ == "__main__":

	# (1) load data
	## directories
	path_data = ""
	path_adata = ""
	path_chemvae_model = ""
	path_scvi_model_cinn = ""

	path_cinn_model_save = ""
	path_sciplex_onehot = ""
	path_chem_onehot = ""
	path_std_param = ""

	adata_train = sc.read(path_adata)
	input_ltpm_label = adata_train.obs.copy()

	## onehot
	data_sciplex_onehot = np.load(path_sciplex_onehot)
	data_chem_onehot = np.load(path_chem_onehot)

	data_trt = pd.read_csv(os.path.join(path_data, 'emb_named_chemvae_canonize.csv'))
	data_trt['Indices'] = list(range(data_trt.shape[0]))

	cell_embdata = input_ltpm_label.loc[:, ['treatment']].merge(data_trt, how = 'left', on = 'treatment')
	indices_onehot = list(cell_embdata['Indices'])

	data_sciplexKept_onehot = data_sciplex_onehot[indices_onehot]

	# cell type and dose covariates
	dose_pd = pd.get_dummies(list(adata_train.obs['dose'].astype(int).astype(str)))
	dose_onehot_data = dose_pd.values.astype('float64')

	cell_type_pd = pd.get_dummies(list(adata_train.obs['cell_type'].astype(str)))
	cell_onehot_data = cell_type_pd.values.astype('float64')

	dose_cell_onehot = np.concatenate((dose_onehot_data, cell_onehot_data), axis=1)

	# (2) load models
	## generation scvi

	scvi.data.setup_anndata(adata_train, layer="counts")
	scvi_model_cinn = scvi.model.SCVI.load(path_scvi_model_cinn, adata_train, use_cuda=False)
	scvi_model_de = scvi_predictive_z(scvi_model_cinn)

	device = 'cuda' if torch.cuda.is_available() else 'cpu'

	## ChemicalVAE
	model_chemvae = ChemicalVAE(n_char=data_chem_onehot.shape[2], max_len=data_chem_onehot.shape[1]).to(device)
	model_chemvae.load_state_dict(torch.load(path_chemvae_model, map_location=device))
	model_chemvae.eval()

	## standardization model
    std_model = Standardize(data_all = data_chem_onehot, model = model_chemvae, device = device)

	# cell type and dose covariates
	dose_pd = pd.get_dummies(list(input_ltpm_label['dose'].astype(int).astype(str)))
	dose_onehot_data = dose_pd.values.astype('float64')

	cell_type_pd = pd.get_dummies(list(input_ltpm_label['cell_type'].astype(str)))
	cell_onehot_data = cell_type_pd.values.astype('float64')

	## perturbnet

	## PCA
	if sparse.issparse(adata_train.X):
		usedata = adata_train.X.A
	else:
		usedata = adata_train.X

	if sparse.issparse(adata_train.layers['counts']):
		usedata_count = adata_train.layers['counts'].A
	else:
		usedata_count = adata_train.layers['counts']

	flow_model = ConditionalFlatCouplingFlow(conditioning_dim=203, # extra 7 columns from cell type and dose
											 # condition dimensions
											 embedding_dim=10,
											 conditioning_depth=2,
											 n_flows=20,
											 in_channels=10,
											 hidden_dim=1024,
											 hidden_depth=2,
											 activation="none",
											 conditioner_use_bn=True)

	model_c = Net2NetFlow_scVIChemStdStatesFlow(configured_flow=flow_model,
												first_stage_data=usedata_count,
												cond_stage_data = data_sciplexKept_onehot,
												model_con = model_chemvae,
												scvi_model = scvi_model_cinn,
												std_model = std_model,
												cell_type_onehot = cell_onehot_data,
												dose_onehot = dose_onehot_data)

	model_c.to(device=device)
	model_c.train(n_epochs=50, batch_size=128, lr=4.5e-6)
	#### save the model
	model_c.save(path_cinn_model_save)
