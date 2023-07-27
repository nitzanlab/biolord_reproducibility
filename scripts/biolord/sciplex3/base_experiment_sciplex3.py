import sys
sys.path.append("../../../")
from paths import DATA_DIR
import biolord
import argparse
import numpy as np
import scanpy as sc
import wandb
import os
from utils.utils_sciplex3 import compute_prediction, compute_baseline, create_df
from sciplex3_config import varying_arg

def train_dataset():
    """
    Read arguments if this script is called from a terminal.
    """
    parser = argparse.ArgumentParser(description=".")

    parser.add_argument("--checkpoint_freq", type=int, default=20) #100
    parser.add_argument("--biolord_epochs", type=int, default=500) #200
    parser.add_argument("--project_name", type=str, default="try") #200
    parser.add_argument("--params_seed", type=int)
    parser.add_argument("--batch_size", type=int, default=512) #200
    parser.add_argument("--model_name", type=str, default="model_name") #200

    args = parser.parse_args()

    np.random.seed(args.params_seed)

    varying_arg["adata"] = sc.read(os.path.join(DATA_DIR, "perturbations/sciplex3/sciplex3_biolord.h5ad"))
    varying_arg["ordered_attributes_keys"] = ["rdkit2d_dose"]
    varying_arg["categorical_attributes_keys"] = ["cell_type"]
    varying_arg["retrieval_attribute_key"] = None
    varying_arg["split_key"] = "split_ood"
    varying_arg["layer"] = None
    varying_arg["gene_likelihood"] = "normal"

    varying_arg["dataset_name"] = "sciplex3"
    
    module_params = {
        "decoder_width": varying_arg["decoder_width"],
        "decoder_depth": varying_arg["decoder_depth"],
        "attribute_nn_width":  varying_arg["attribute_nn_width"],
        "attribute_nn_depth": varying_arg["attribute_nn_depth"],
        "use_batch_norm": varying_arg["use_batch_norm"],
        "use_layer_norm": varying_arg["use_layer_norm"],
        "attribute_dropout_rate":  varying_arg["attribute_dropout_rate"],
        "unknown_attribute_noise_param": varying_arg["unknown_attribute_noise_param"],
        "seed": varying_arg["seed"],
        "n_latent_attribute_ordered": varying_arg["n_latent_attribute_ordered"],
        "n_latent_attribute_categorical": varying_arg["n_latent_attribute_categorical"],
        "loss_ae": varying_arg["loss_ae"],
        "reconstruction_penalty": varying_arg["reconstruction_penalty"],
        "unknown_attribute_penalty":varying_arg["unknown_attribute_penalty"],
    }
    
    trainer_params = {
        "n_epochs_warmup": 0,
        "latent_lr": varying_arg["latent_lr"],
        "latent_wd": varying_arg["latent_wd"],
        "decoder_lr": varying_arg["decoder_lr"],
        "decoder_wd": varying_arg["decoder_wd"],
        "attribute_nn_lr": varying_arg["attribute_nn_lr"],
        "attribute_nn_wd": varying_arg["attribute_nn_wd"],
        "step_size_lr": varying_arg["step_size_lr"],
        "cosine_scheduler": varying_arg["cosine_scheduler"],
        "scheduler_final_lr": varying_arg["scheduler_final_lr"]
    }

    wandb.init(project=args.project_name, entity="biolord", config = varying_arg)
    adata = varying_arg["adata"]


    biolord.Biolord.setup_anndata(
        adata,
        ordered_attributes_keys=varying_arg["ordered_attributes_keys"],
        categorical_attributes_keys=varying_arg["categorical_attributes_keys"],
        retrieval_attribute_key=varying_arg["retrieval_attribute_key"],
    )

    print(varying_arg["split_key"])

    model = biolord.Biolord(
        adata=adata,
        n_latent=varying_arg["n_latent"],
        model_name=args.model_name,
        module_params=module_params,
        train_classifiers=False,
        split_key=varying_arg["split_key"],
    )
    
    model.train(
        max_epochs=args.biolord_epochs,
        batch_size=args.batch_size,
        plan_kwargs=trainer_params,
        early_stopping=True,
        early_stopping_patience=20,
        check_val_every_n_epoch=args.checkpoint_freq,
        num_workers=1,
    )

    model.module.eval()

    idx_test_control = np.where(
        (adata.obs["split_ood"] == "test") & (adata.obs["control"] == 1)
    )[0]

    adata_test_control = adata[idx_test_control].copy()

    idx_ood = np.where(
        (adata.obs["split_ood"] == "ood")
    )[0]

    adata_ood = adata[idx_ood].copy()

    dataset_control = model.get_dataset(adata_test_control)
    dataset_ood = model.get_dataset(adata_ood)

    #evaluate predictions and baseline
    res_all = {}

    res_all["baseline"], _ = compute_baseline(
        model=model,
        adata=adata_ood,
        dataset=dataset_ood,
        dataset_control=dataset_control,
        use_DEGs=False,
        verbose=False,
    )

    res_all["biolord"], _ = compute_prediction(
        model=model,
        adata=adata_ood,
        dataset=dataset_ood,
        dataset_control=dataset_control,
        use_DEGs=False,
        verbose=False,
    )

    print("res_all[“biolord”]", res_all["biolord"])

    df_all = create_df(res_all)

    print("df_all", df_all)

    mean_all = df_all.groupby(by=["dose", "type"]).mean().reset_index()

    mean_all_numpy = mean_all.to_numpy()
    print("mean_all", mean_all)

    wandb.log({
        "mean_all_biolord/0.001": mean_all_numpy[1][2],
        "mean_all_biolord/0.01": mean_all_numpy[3][2],
        "mean_all_biolord/0.1": mean_all_numpy[5][2],
        "mean_all_biolord/1.0": mean_all_numpy[7][2],
    })

    wandb.log({
        "mean_all_baseline/0.001": mean_all_numpy[0][2],
        "mean_all_baseline/0.01": mean_all_numpy[2][2],
        "mean_all_baseline/0.1": mean_all_numpy[4][2],
        "mean_all_baseline/1.0": mean_all_numpy[6][2],
    })

if __name__ == "__main__":
    train_dataset()


