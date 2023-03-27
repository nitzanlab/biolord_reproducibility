import numpy as np
import pandas as pd
import torch

from scipy.stats import pearsonr
from sklearn.metrics import mean_squared_error as mse


def bool2idx(x):
    """
    Returns the indices of the True-valued entries in a boolean array `x`
    """
    return np.where(x)[0]

def repeat_n(x, n):
    """combo_seen2
    Returns an n-times repeated version of the Tensor x,
    repetition dimension is axis 0
    """
    # copy tensor to device BEFORE replicating it n times
    device = "cuda" if torch.cuda.is_available() else "cpu"
    return x.to(device).view(1, -1).repeat(n, 1)


def compute_metrics(predictions_dict, ctrl=0):
    """
    Given results from a model run and the ground truth, compute metrics
    """
    metrics = {}
    metrics_pert = {}

    metric2fct = {
           "mse": mse,
           "pearson": pearsonr
    }
    
    for m in metric2fct.keys():
        metrics[m] = []
        metrics[m + "_de"] = []

    for pert in predictions_dict:

        metrics_pert[pert] = {}
    
        for m, fct in metric2fct.items():
            if m == "pearson":
                val = fct(
                    predictions_dict[pert][1].cpu().numpy()-ctrl.cpu().numpy(), 
                    predictions_dict[pert][0].cpu().numpy()-ctrl.cpu().numpy())[0]
                if np.isnan(val):
                    val = 0
            else:
                val = fct(predictions_dict[pert][1].cpu().numpy(), predictions_dict[pert][0].cpu().numpy())

            metrics_pert[pert][m] = val
            metrics[m].append(metrics_pert[pert][m])

       
        if pert != "ctrl":
            
            for m, fct in metric2fct.items():
                if m == "pearson":
                    val = fct(
                        predictions_dict[pert][1].cpu().numpy()[predictions_dict[pert][-1]]-ctrl.cpu().numpy()[predictions_dict[pert][-1]], 
                        predictions_dict[pert][0].cpu().numpy()[predictions_dict[pert][-1]]-ctrl.cpu().numpy()[predictions_dict[pert][-1]])[0]
                    if np.isnan(val):
                        val = 0
                else:
                    val = fct(
                        predictions_dict[pert][1].cpu().numpy()[predictions_dict[pert][-1]], 
                        predictions_dict[pert][0].cpu().numpy()[predictions_dict[pert][-1]]
                    )
                    
                metrics_pert[pert][m + "_de"] = val
                metrics[m + "_de"].append(metrics_pert[pert][m + "_de"])

        else:
            for m, fct in metric2fct.items():
                metrics_pert[pert][m + "_de"] = 0
    
    for m in metric2fct.keys():
        
        metrics[m] = np.mean(metrics[m])
        metrics[m + "_de"] = np.mean(metrics[m + "_de"])
    
    return metrics, metrics_pert

def compute_metrics(predictions_dict, ctrl=0):
    """
    Given results from a model run and the ground truth, compute metrics
    """
    metrics = {}
    metrics_pert = {}

    metric2fct = {
           "mse": mse,
           "pearson": pearsonr
    }
    
    for m in metric2fct.keys():
        metrics[m] = []
        metrics[m + "_de"] = []

    for pert in predictions_dict:

        metrics_pert[pert] = {}
            
        for m, fct in metric2fct.items():
            if m == "pearson":
                val = fct(predictions_dict[pert][1]-ctrl, predictions_dict[pert][0]-ctrl)[0]
                if np.isnan(val):
                    val = 0
            else:
                val = fct(predictions_dict[pert][1], predictions_dict[pert][0])

            metrics_pert[pert][m] = val
            metrics[m].append(metrics_pert[pert][m])

       
        if pert != "ctrl":
            
            for m, fct in metric2fct.items():
                if m == "pearson":
                    val = fct(
                        predictions_dict[pert][1][predictions_dict[pert][-1]]-ctrl[predictions_dict[pert][-1]], 
                        predictions_dict[pert][0][predictions_dict[pert][-1]]-ctrl[predictions_dict[pert][-1]])[0]
                    if np.isnan(val):
                        val = 0
                else:
                    val = fct(
                        predictions_dict[pert][1][predictions_dict[pert][-1]], 
                        predictions_dict[pert][0][predictions_dict[pert][-1]]
                    )
                    
                metrics_pert[pert][m + "_de"] = val
                metrics[m + "_de"].append(metrics_pert[pert][m + "_de"])

        else:
            for m, fct in metric2fct.items():
                metrics_pert[pert][m + "_de"] = 0
    
    for m in metric2fct.keys():
        
        metrics[m] = np.mean(metrics[m])
        metrics[m + "_de"] = np.mean(metrics[m + "_de"])
    
    return metrics, metrics_pert
