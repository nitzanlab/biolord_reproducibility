import os
import sys
import scanpy as sc
import anndata
import numpy as np
import pandas as pd
import argparse

import biolord
import pickle
import wandb


from utils_perturbations import (
    compute_metrics_cpu,
    bool2idx,
    repeat_n
)


sys.path.append("../../../")
from paths import DATA_DIR

adata = sc.read(
    DATA_DIR + "perturbations/norman/norman_biolord.h5ad",
)

adata_single = sc.read(
    DATA_DIR + "perturbations/norman/norman_single_biolord.h5ad",
)

with open(DATA_DIR + "perturbations/norman/no_perturb_subgroup_analysis.pkl", "rb") as f:
    no_perturb_subgroup = pickle.load(f)


def train_dataset():
    """
    Read arguments if this script is called from a terminal.
    """
    parser = argparse.ArgumentParser(description=".")
    parser.add_argument("--checkpoint_freq", type=int, default=20)
    parser.add_argument("--biolord_epochs", type=int, default=500)
    parser.add_argument("--project_name", type=str, default="norman") 
    parser.add_argument("--params_seed", type=int)
    parser.add_argument("--batch_size", type=int, default=32) 
    parser.add_argument("--model_name", type=str, default="norman")
    parser.add_argument("--optimal_params", default=False)

    args = parser.parse_args()
    
    np.random.seed(args.params_seed)
    if args.optimal_params:
        print("Using optimal params for train.")
        from norman_config_optimal import varying_arg
    else:
        print("Using random params for train.")
        from norman_config import varying_arg
    
    wandb.init(project=args.project_name, entity="biolord", config = varying_arg)
    wandb.log({'params_seed': args.params_seed})

    pert2neighbor = np.asarray([val for val in adata.uns["pert2neighbor"].values()])
    keep_idx = pert2neighbor.sum(0) > 0

    name_map = dict(adata.obs[["condition", "condition_name"]].drop_duplicates().values)
    ctrl = np.asarray(adata[adata.obs["condition"].isin(["ctrl"])].X.mean(0)).flatten()


    df_perts_expression = pd.DataFrame(adata.X.A, index=adata.obs_names, columns=adata.var_names)
    df_perts_expression["condition"] = adata.obs["condition"]
    df_perts_expression = df_perts_expression.groupby(["condition"]).mean()
    df_perts_expression = df_perts_expression.reset_index()

    single_perts_condition = []
    single_pert_val = []
    double_perts = []
    for pert in adata.obs["condition"].cat.categories:
        if len(pert.split("+")) == 1:
            continue
        elif "ctrl" in pert:
            single_perts_condition.append(pert)
            p1, p2 = pert.split("+")
            if p2 == "ctrl":
                single_pert_val.append(p1)
            else:
                single_pert_val.append(p2)
        else:
            double_perts.append(pert)
    single_perts_condition.append("ctrl")
    single_pert_val.append("ctrl")

    df_singleperts_expression = pd.DataFrame(df_perts_expression.set_index("condition").loc[single_perts_condition].values, index=single_pert_val)
    df_singleperts_emb = np.asarray([adata.uns["pert2neighbor"][p1][keep_idx] for p1 in df_singleperts_expression.index])
    df_singleperts_pca = sc.pp.pca(df_singleperts_emb)
    df_singleperts_condition = pd.Index(single_perts_condition)
    df_single_pert_val = pd.Index(single_pert_val)

    df_doubleperts_expression = df_perts_expression.set_index("condition").loc[double_perts].values
    df_doubleperts_condition = pd.Index(double_perts)

    np.random.seed(42)

    module_params = {
        "attribute_nn_width":  varying_arg["attribute_nn_width"],
        "attribute_nn_depth": varying_arg["attribute_nn_depth"],
        "use_batch_norm": varying_arg["use_batch_norm"],
        "use_layer_norm": varying_arg["use_layer_norm"],
        "attribute_dropout_rate":  varying_arg["attribute_dropout_rate"],
        "unknown_attribute_noise_param": varying_arg["unknown_attribute_noise_param"],
        "seed": varying_arg["seed"],
        "n_latent_attribute_ordered": varying_arg["n_latent_attribute_ordered"],
        "n_latent_attribute_categorical": varying_arg["n_latent_attribute_categorical"],
        "reconstruction_penalty": varying_arg["reconstruction_penalty"],
        "unknown_attribute_penalty": varying_arg["unknown_attribute_penalty"],
        "decoder_width": varying_arg["decoder_width"],
        "decoder_depth": varying_arg["decoder_depth"],
        "decoder_activation": varying_arg["decoder_activation"],
        "attribute_nn_activation": varying_arg["attribute_nn_activation"],
        "unknown_attributes": varying_arg["unknown_attributes"],
    }


    trainer_params = {
        "n_epochs_warmup": 0,
        "latent_lr": varying_arg["latent_lr"],
        "latent_wd": varying_arg["latent_wd"],
        "attribute_nn_lr": varying_arg["attribute_nn_lr"],
        "attribute_nn_wd": varying_arg["attribute_nn_wd"],
        "step_size_lr": varying_arg["step_size_lr"],
        "cosine_scheduler": varying_arg["cosine_scheduler"],
        "scheduler_final_lr": varying_arg["scheduler_final_lr"],
        "decoder_lr": varying_arg["decoder_lr"],
        "decoder_wd": varying_arg["decoder_wd"]
    }


    test_metrics_biolord_delta = {}
    test_metrics_biolord_delta_normalized = {}
    
    ordered_attributes_key = varying_arg["ordered_attributes_key"]

    biolord.Biolord.setup_anndata(
        adata_single,
        ordered_attributes_keys=[ordered_attributes_key],
        categorical_attributes_keys=None,
        retrieval_attribute_key=None,
    )

    for split_seed in range(1,6):
        test_metrics_biolord_delta[split_seed] = {}
        test_metrics_biolord_delta_normalized[split_seed] = {}

        train_idx = df_singleperts_condition.isin(adata[adata.obs[f"split{split_seed}"] == "train"].obs["condition"].cat.categories)
        train_condition_perts = df_singleperts_condition[train_idx]
        train_condition_perts_double = df_doubleperts_condition[df_doubleperts_condition.isin(adata[adata.obs[f"split{split_seed}"] == "train"].obs["condition"].cat.categories)]
        train_perts = df_single_pert_val[train_idx]

        model = biolord.Biolord(
            adata=adata_single,
            n_latent=varying_arg["n_latent"],
            model_name="norman",
            module_params=module_params,
            train_classifiers=False,
            split_key=f"split{split_seed}"
        )

        model.train(
            max_epochs=int(varying_arg["max_epochs"]),
            batch_size=32,
            plan_kwargs=trainer_params,
            early_stopping=True,
            early_stopping_patience=int(varying_arg["early_stopping_patience"]),
            check_val_every_n_epoch=5,
            num_workers=1,
            enable_checkpointing=False
        )
        adata_control = adata_single[adata_single.obs["condition"] == "ctrl"].copy()
        dataset_control = model.get_dataset(adata_control)

        dataset_reference = model.get_dataset(adata_single)

        n_obs = adata_control.shape[0]
        for ood_set in ["combo_seen0", "combo_seen1", "combo_seen2", "unseen_single"]:
            predictions_dict_delta = {}
            predictions_dict_mean = {}
            perts = adata[adata.obs[f"subgroup{split_seed}"] == ood_set].obs["condition"].cat.categories

            for i, pert in enumerate(perts):
                bool_de = adata.var_names.isin(
                            np.array(adata.uns["top_non_zero_de_20"][name_map[pert]])
                        )
                idx_de = bool2idx(bool_de)
                if pert in train_condition_perts:
                    idx_ref =  bool2idx(adata_single.obs["condition"] == pert)[0]
                    expression_pert = dataset_reference["X"][[idx_ref], :].mean(0).cpu().numpy()
                    test_preds_delta = expression_pert
                elif pert in train_condition_perts_double:
                    expression_pert = df_doubleperts_expression[df_doubleperts_condition.isin([pert]), :]
                    test_preds_delta = df_doubleperts_expression[df_doubleperts_condition.isin([pert]), :]
                elif "ctrl" in pert:
                    idx_ref =  bool2idx(adata_single.obs["condition"] == pert)[0]
                    expression_pert = dataset_reference["X"][[idx_ref], :].mean(0).cpu().numpy()

                    dataset_pred = dataset_control.copy()
                    dataset_pred[ordered_attributes_key] = repeat_n(dataset_reference[ordered_attributes_key][idx_ref, :], n_obs)
                    test_preds, _ = model.module.get_expression(dataset_pred)

                    test_preds_delta = test_preds.cpu().numpy()

                else:
                    expression_pert = df_doubleperts_expression[df_doubleperts_condition.isin([pert]), :]
                    test_preds_add = []
                    for p in pert.split("+"):
                        if p in train_perts:
                            test_predsp = df_singleperts_expression.values[df_single_pert_val.isin([p]), :]
                            test_preds_add.append(test_predsp[0, :])
                        else:
                            idx_ref =  bool2idx(adata_single.obs["perts_name"].isin([p]))[0]

                            dataset_pred = dataset_control.copy()
                            dataset_pred[ordered_attributes_key] = repeat_n(dataset_reference[ordered_attributes_key][idx_ref, :], n_obs)
                            test_preds, _ = model.module.get_expression(dataset_pred)
                            test_preds_add.append(test_preds.cpu().numpy())

                    test_preds_delta = test_preds_add[0] + test_preds_add[1] - ctrl

                predictions_dict_delta[name_map[pert]] = [
                        expression_pert.flatten(),
                        test_preds_delta.flatten(),
                        idx_de
                ]



            test_metrics_biolord_delta[split_seed][ood_set], _ = compute_metrics_cpu(predictions_dict_delta, ctrl=ctrl)

            test_metrics_biolord_delta_normalized[split_seed][ood_set] = {key_: val_ / no_perturb_subgroup[split_seed][ood_set][key_] for key_, val_ in test_metrics_biolord_delta[split_seed][ood_set].items()}
            
            print(test_metrics_biolord_delta_normalized[split_seed][ood_set])
            wandb.log({"%s/seed_%d_mse_de"%(ood_set,split_seed): test_metrics_biolord_delta[split_seed][ood_set]['mse_de'],
                       })

            wandb.log({"normalized_%s/seed_%d_mse_de"%(ood_set,split_seed): test_metrics_biolord_delta_normalized[split_seed][ood_set]['mse_de'],
                       })
            
    res_biolord_delta_normalized = {}
    for key in test_metrics_biolord_delta_normalized:
        res_biolord_delta_normalized[f"mse_de_seed{key}"] = pd.DataFrame(test_metrics_biolord_delta_normalized[key]).T["mse_de"]

    if args.optimal_params:
        pd.DataFrame(res_biolord_delta_normalized).to_csv(DATA_DIR + "perturbations/norman/biolord_normalized_mse_de_seeds.csv")
            
    for ood_set in ["combo_seen0", "combo_seen1", "combo_seen2", "unseen_single"]:
        print(ood_set)
        test_metrics_biolord_normalized_list = []
        for split_seed in range(1,6):
            print(f"test_metrics_biolord_delta_normalized[{split_seed}][{ood_set}]",test_metrics_biolord_delta_normalized[split_seed][ood_set])
            test_metrics_biolord_normalized_list.append(test_metrics_biolord_delta_normalized[split_seed][ood_set]['mse_de'])

        print(np.mean(test_metrics_biolord_normalized_list))
        wandb.log({
                "normalized_%s/mean_mse_de"%ood_set: np.mean(test_metrics_biolord_normalized_list),
                })

if __name__ == "__main__":
    train_dataset()
