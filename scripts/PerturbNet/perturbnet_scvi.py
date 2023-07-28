import scanpy as sc
import scvi


if __name__ == "__main__":

	adata_path = "" # change to your setting
	model_path = "" # change to your setting

	adata = sc.read(adata_path)

	scvi.data.setup_anndata(adata, layer = "counts")

	model = scvi.model.SCVI(adata, n_latent=10)
	model.train(n_epochs=700)

	model.save(model_path)
