import numpy as np

varying_arg = {
    "seed": np.random.choice([0, 1, 2, 3, 4]), "batch_size": np.random.choice([256, 512]),
    "unknown_attribute_noise_param": np.random.choice([0.1, 0.5, 1, 2, 5, 10, 20]),
    "decoder_width": np.random.choice([128, 256, 512, 1024, 2048, 4096]),
    "decoder_depth": np.random.choice([1, 2, 3, 4, 6, 8]),
    "use_batch_norm": np.random.choice([True, False]),
    "use_layer_norm": np.random.choice([True, False]),
    "step_size_lr": np.random.choice([45, 90, 180]),
    "attribute_dropout_rate": np.random.choice([0.05, 0.1, 0.25, 0.5, 0.75]),
    "cosine_scheduler": np.random.choice([True, False]),
    "scheduler_final_lr": np.random.choice([1e-3, 1e-4, 1e-5, 1e-6]),
    "n_latent": np.random.choice([16, 32, 64, 128, 256]),
    "n_latent_attribute_ordered": np.random.choice([128, 256, 512]),
    "n_latent_attribute_categorical": np.random.choice([2, 3, 4, 6, 8]),
    "reconstruction_penalty": np.random.choice([1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4]),
    "attribute_nn_width": np.random.choice([32, 64, 128, 256, 512]),
    "attribute_nn_depth": np.random.choice([1, 2, 3, 4, 6, 8]),
    "attribute_nn_lr": np.random.choice([1e-2, 1e-3, 1e-4]),
    "attribute_nn_wd": np.random.choice([1e-8, 4e-8, 1e-7]),
    "latent_lr": np.random.choice([1e-2, 1e-3, 1e-4]),
    "latent_wd": np.random.choice([1e-2, 1e-3, 1e-4]),
    "decoder_lr": np.random.choice([1e-2, 1e-3, 1e-4]),
    "decoder_wd": np.random.choice([1e-2, 1e-3, 1e-4]),
    "unknown_attribute_penalty": np.random.choice([1e-1, 1e0, 2e0, 5e0, 1e1, 2e1, 5e1, 1e2, 2e2])
}

