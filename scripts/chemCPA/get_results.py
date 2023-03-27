CONFIG_PATH = "/path/to/config.yml"
CHECKPOINT_HASH = ""
CHEMCPA_PATH = "/path/to/chemCPA/notebooks/"

import sys
sys.path.append(CHEMCPA_PATH)
sys.path.append("../../")
from paths import DATA_DIR


# chemCPA utils file
from utils import load_config, load_dataset, load_smiles, load_model, compute_drug_embeddings, compute_pred, compute_pred_ctrl
from seml.config import generate_configs, read_config
import pandas as pd

from chemCPA.data import load_dataset_splits
import argparse
import wandb



def create_df(drug_r2_pretrained, type):

    df_pretrained = pd.DataFrame.from_dict(drug_r2_pretrained, orient="index", columns=["r2_de"])
    df_pretrained["type"] = type

    df = df_pretrained

    df["r2_de"] = df["r2_de"].apply(lambda x: max(x,0))
    df["cell_line"] = pd.Series(df.index.values).apply(lambda x: x.split("_")[0]).values
    df["drug"] = pd.Series(df.index.values).apply(lambda x: x.split("_")[1]).values
    df["dose"] = pd.Series(df.index.values).apply(lambda x: x.split("_")[2]).values
    df["dose"] = df["dose"].astype(float)

    df["combination"] = df.index.values

    df = df.reset_index()
    return df


#create args

parser = argparse.ArgumentParser()
parser.add_argument("--project_name", type=str, default="")
parser.add_argument("--file", type=str, default="")
parser.add_argument("--type", type=str, default="pretrained")
args = parser.parse_args()

print("args",args)

wandb.init(project="chemCPA" + args.project_name, entity="biolord", config = args)
wandb.config = args

seml_collection = "multi_task"
model_hash_pretrained = CHECKPOINT_HASH

dosages = [1e1, 1e2, 1e3, 1e4]
cell_lines = ["A549", "K562", "MCF7"]
use_DEGs = True

seml_config, slurm_config, experiment_config = read_config(CONFIG_PATH)

config = generate_configs(experiment_config)[0]

dataset, key_dict = load_dataset(config)

config["dataset"]["n_vars"] = dataset.n_vars
config["config_hash"] = model_hash_pretrained

canon_smiles_unique_sorted, smiles_to_pathway_map, smiles_to_drug_map = load_smiles(config, dataset, key_dict, True)
model, embedding_pretrained = load_model(config, canon_smiles_unique_sorted)


datasets = load_dataset_splits(
                    dataset_path=str(DATA_DIR) + "/perturbation/sciplex3_biolord.h5ad",
                    perturbation_key="condition",
                    dose_key="dose",
                    covariate_keys="cell_type",
                    smiles_key="SMILES",
                    degs_key="all_DEGs",
                    pert_category="cov_drug_dose_name",
                    split_key="split_ood",
                    return_dataset=False,
                    use_drugs_idx=True
                )


drug_r2_res, _ = compute_pred(
    model,
    datasets["ood"],
    genes_control=datasets["test_control"].genes,
    dosages=dosages,
    cell_lines=cell_lines,
    use_DEGs=False,
    verbose=False
)

df_all = create_df(drug_r2_res, args.type)
# save res
df_all.to_csv(str(DATA_DIR) + f"/perturbation/chemCPA_{type}_all.csv")


# calc mean
df_all_mean = df_all.groupby("dose", as_index=False)["r2_de"].mean()
print("df_all_mean", df_all_mean)

# calc median
df_all_med = df_all.groupby("dose", as_index=False)["r2_de"].median()

print("df_all_med", df_all_med)

wandb.log({
    "df_all_mean": df_all_mean["r2_de"].values[-1],
    "df_all_med": df_all_med["r2_de"].values[-1],
})

