import numpy as np
import pandas as pd
import torch
from torchmetrics import R2Score
from tqdm import tqdm


def bool2idx(x):
    """
    Returns the indices of the True-valued entries in a boolean array `x`
    """
    return np.where(x)[0]

def repeat_n(x, n):
    """
    Returns an n-times repeated version of the Tensor x,
    repetition dimension is axis 0
    """
    # copy tensor to device BEFORE replicating it n times
    device = "cuda" if torch.cuda.is_available() else "cpu"
    return x.to(device).view(1, -1).repeat(n, 1)

def compute_r2(y_true, y_pred):
    """
    Computes the r2 score for `y_true` and `y_pred`,
    returns `-1` when `y_pred` contains nan values
    """
    y_pred = torch.clamp(y_pred, -3e12, 3e12)
    metric = R2Score().to(y_true.device)
    metric.update(y_pred, y_true)  # same as sklearn.r2_score(y_true, y_pred)
    return metric.compute().item()

def compute_prediction(
    model,
    adata,
    dataset,
    cell_lines=None,
    dataset_control=None,
    use_DEGs=True,
    verbose=True
):
    pert_categories_index = pd.Index(adata.obs["cov_drug_dose_name"].values, dtype="category")
    allowed_cell_lines = []

    cl_dict = {
        torch.Tensor([0.]): "A549",
        torch.Tensor([1.]): "K562",
        torch.Tensor([2.]): "MCF7",
    }

    if cell_lines is None:
        cell_lines = ["A549", "K562", "MCF7"]

    print(cell_lines)
    layer = "X" if "X" in dataset else "layers"
    predictions_dict = {}
    drug_r2 = {}
    for cell_drug_dose_comb, category_count in tqdm(
        zip(*np.unique(pert_categories_index.values, return_counts=True))
    ):
        # estimate metrics only for reasonably-sized drug/cell-type combos
        if category_count <= 5:
            continue
        # doesn"t make sense to evaluate DMSO (=control) as a perturbation
        if (
            "dmso" in cell_drug_dose_comb.lower()
            or "control" in cell_drug_dose_comb.lower()
        ):
            continue

        # adata.var_names is the list of gene names
        # adata.uns["all_DEGs"] is a dict, containing a list of all differentiably-expressed
        # genes for every cell_drug_dose combination.
        bool_de = adata.var_names.isin(
            np.array(adata.uns["all_DEGs"][cell_drug_dose_comb])
        )
        idx_de = bool2idx(bool_de)

        # need at least two genes to be able to calc r2 score
        if len(idx_de) < 2:
            continue

        bool_category = pert_categories_index.get_loc(cell_drug_dose_comb)
        idx_all = bool2idx(bool_category)
        idx = idx_all[0]
        y_true = dataset[layer][idx_all, :].to(model.device)
        
                    
        dataset_comb = {}
        if dataset_control is None:
            n_obs = y_true.size(0).to(model.device)
            for key, val in dataset.items():
                dataset_comb[key] = val[idx_all].to(model.device)
        else:
            n_obs = dataset_control[layer].size(0)
            dataset_comb[layer] = dataset_control[layer].to(model.device)
            dataset_comb["ind_x"] = dataset_control["ind_x"].to(model.device)
            for key in dataset_control:
                if key not in [layer, "ind_x"]:
                    dataset_comb[key] = repeat_n(dataset[key][idx, :], n_obs)

        stop = False
        for tensor, cl in cl_dict.items():
            if (tensor == dataset["cell_type"][idx]).all():
                if cl not in cell_lines:
                    stop = True
        if stop:
            continue
            
        pred, _ = model.module.get_expression(dataset_comb)

        y_pred = pred.mean(0)
        y_true = y_true.mean(0)
        if use_DEGs:
            r2_m_de = compute_r2(y_true[idx_de].cuda(), y_pred[idx_de].cuda())
            print(f"{cell_drug_dose_comb}: {r2_m_de:.2f}") if verbose else None
            drug_r2[cell_drug_dose_comb] = r2_m_de
        else:
            r2_m = compute_r2(y_true.cuda(), y_pred.cuda())
            print(f"{cell_drug_dose_comb}: {r2_m:.2f}") if verbose else None
            drug_r2[cell_drug_dose_comb] = r2_m

        predictions_dict[cell_drug_dose_comb] = [y_true, y_pred, idx_de]
    return drug_r2, predictions_dict


def compute_baseline(
    model,
    adata,
    dataset,
    cell_lines=None,
    dataset_control=None,
    use_DEGs=True,
    verbose=True,
):
    pert_categories_index = pd.Index(adata.obs["cov_drug_dose_name"].values, dtype="category")
    allowed_cell_lines = []

    cl_dict = {
        torch.Tensor([0.]): "A549",
        torch.Tensor([1.]): "K562", 
        torch.Tensor([2.]): "MCF7",
    }
    
    cl_dict_op = {
        "A549":torch.Tensor([0.]),
        "K562": torch.Tensor([1.]),
        "MCF7": torch.Tensor([2.]),
    }

    if cell_lines is None:
        cell_lines = ["A549", "K562", "MCF7"]

    print(cell_lines)

    layer = "X" if "X" in dataset else "layers"
    predictions_dict = {}
    drug_r2 = {}
    for cell_drug_dose_comb, category_count in tqdm(
        zip(*np.unique(pert_categories_index.values, return_counts=True))
    ):
        # estimate metrics only for reasonably-sized drug/cell-type combos
        if category_count <= 5:
            continue

        # doesn"t make sense to evaluate DMSO (=control) as a perturbation
        if (
            "dmso" in cell_drug_dose_comb.lower()
            or "control" in cell_drug_dose_comb.lower()
        ):
            continue

        # adata.var_names is the list of gene names
        # adata.uns["all_DEGs"] is a dict, containing a list of all differentiably-expressed
        # genes for every cell_drug_dose combination.
        bool_de = adata.var_names.isin(
            np.array(adata.uns["all_DEGs"][cell_drug_dose_comb])
        )
        idx_de = bool2idx(bool_de)

        # need at least two genes to be able to calc r2 score
        if len(idx_de) < 2:
            continue

        bool_category = pert_categories_index.get_loc(cell_drug_dose_comb)
        idx_all = bool2idx(bool_category)
        idx = idx_all[0]
        y_true = dataset[layer][idx_all, :].to(model.device)
        
        cov_name = cell_drug_dose_comb.split("_")[0]
        cond = bool2idx(dataset_control["cell_type"] == cl_dict_op[cov_name])
        y_pred = dataset_control[layer][cond, :].to(model.device)

        stop = False
        for tensor, cl in cl_dict.items():
            if (tensor == dataset["cell_type"][idx]).all():
                if cl not in cell_lines:
                    stop = True
        if stop:
            continue
            
        y_pred = y_pred.mean(0)
        y_true = y_true.mean(0)
        if use_DEGs:
            r2_m_de = compute_r2(y_true[idx_de].cuda(), y_pred[idx_de].cuda())
            print(f"{cell_drug_dose_comb}: {r2_m_de:.2f}") if verbose else None
            drug_r2[cell_drug_dose_comb] = r2_m_de
        else:
            r2_m = compute_r2(y_true.cuda(), y_pred.cuda())
            print(f"{cell_drug_dose_comb}: {r2_m:.2f}") if verbose else None
            drug_r2[cell_drug_dose_comb] = r2_m

        predictions_dict[cell_drug_dose_comb] = [y_true, y_pred, idx_de]
    
    return drug_r2, predictions_dict


def create_df(res):
    dfs_ = []
    for key_, res_ in res.items():
        df_ = pd.DataFrame.from_dict(res_, orient="index", columns=["r2_de"])
        df_["type"] = key_
        dfs_.append(df_)
        
    df = pd.concat(dfs_)

    df["r2_de"] = df["r2_de"].apply(lambda x: max(x,0))
    df["cell_line"] = pd.Series(df.index.values).apply(lambda x: x.split("_")[0]).values
    df["drug"] = pd.Series(df.index.values).apply(lambda x: x.split("_")[1]).values
    df["dose"] = pd.Series(df.index.values).apply(lambda x: x.split("_")[2]).values
    df["dose"] = df["dose"].astype(float)

    df["combination"] = df.index.values
    df = df.reset_index()
    return df