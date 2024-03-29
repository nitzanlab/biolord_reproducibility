import numpy as np

varying_arg = {
    "seed": 42,
    "unknown_attribute_noise_param": 0.2, 
    "use_batch_norm": False,
    "use_layer_norm": False, 
    "step_size_lr": 45, 
    "attribute_dropout_rate": 0.0, 
    "cosine_scheduler":True,
    "scheduler_final_lr":1e-5,
    "n_latent":32, 
    "n_latent_attribute_ordered": 32,
    "reconstruction_penalty": 10000.0,
    "attribute_nn_width": 64,
    "attribute_nn_depth" :2, 
    "attribute_nn_lr": 0.001, 
    "attribute_nn_wd": 4e-8,
    "latent_lr": 0.01,
    "latent_wd": 0.00001,
    "decoder_width": 32,
    "decoder_depth": 2,  
    "decoder_activation": True,
    "attribute_nn_activation": True,
    "unknown_attributes": False,
    "decoder_lr": 0.01,
    "decoder_wd": 0.01,
    "max_epochs":200,
    "early_stopping_patience": 200,
    "ordered_attributes_key": "perturbation_neighbors1",
    "n_latent_attribute_categorical": 16,
}

