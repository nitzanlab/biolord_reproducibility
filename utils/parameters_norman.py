module_params = {
    "decoder_width": 32,
    "decoder_depth": 2,
    "decoder_activation": True,
    "attribute_nn_width": 64,
    "attribute_nn_depth": 2,
    "attribute_nn_activation": True,
    "use_batch_norm": False,
    "use_layer_norm": False,
    "unknown_attribute_noise_param": 0.2,
    "seed": 42,
    "n_latent_attribute_ordered": 32,
    "n_latent_attribute_categorical": 16,
    "reconstruction_penalty": 10000.0,
    "unknown_attribute_penalty": 10000.0,
    "attribute_dropout_rate": 0.0,
    "unknown_attributes": False
}

trainer_params = {
    "n_epochs_warmup": 0,
    "latent_lr": 0.1,
    "latent_wd": 0.00001,
    "decoder_lr": 0.01,
    "decoder_wd": 0.01,
    "attribute_nn_lr": 0.001,
    "attribute_nn_wd": 4e-8,
    "step_size_lr": 45,
    "cosine_scheduler": True,
    "scheduler_final_lr": 1e-5,
}