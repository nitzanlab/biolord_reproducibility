ood_drugs = [
    "Dacinostat",
    "Givinostat",
    "Belinostat",
    "Hesperadin",
    "Quisinostat",
    "Alvespimycin",
    "Tanespimycin",
    "TAK-901",
    "Flavopiridol",
]

validation_drugs = [
    "CUDC-907",
    "Luminespib",
    "Epothilone",
    "Raltitrexed",
    "CUDC-101",
    "Abexinostat",
    "Panobinostat",
    "AR-42",
    "Trichostatin",
    "M344",
    "Resminostat",
    "Pirarubicin", #ood
    "Mocetinostat", #no_ood
    "Pracinostat", #no_ood
    "Entinostat", #no_ood
    "Tucidinostat", #no_ood
    "Tacedinaline", #no_ood
    "Patupilone", #no_ood
    "GSK1070916", #no_ood
    "JNJ-26854165", #no_ood
    "Trametinib",  #no_ood
    "Dasatinib" #no_ood
    "YM155", #apoptosis
    "Barasertib", #cell cycle
    "Fulvestrant", #nuclear receptor
    "Nintedanib", #tyrosine
    "Rigosertib", #tyrosine
    "BMS-754807", #tyrosine
    "KW-2449", #tyrosine
    "Crizotinib", #tyrosin
    "ENMD-2076", #cell cycle
    "Alisertib", #cell cycle
    "(+)-JQ1", #epigenetic
]


groups = [
    "Epigenetic regulation",
    "Tyrosine kinase signaling",
    "DNA damage & DNA repair",
    "Cell cycle regulation"
]

groups_colors = {
    "Epigenetic regulation": "#FFE162",
    "Histone deacetylation" : "#FFE162",
    "Tyrosine kinase signaling": "#FF6464",
    "DNA damage & DNA repair":"#91C483",
    "Cell cycle regulation": "#7EC8E3",
    "NA":"#EEEEEE", 
}


decoder_width = 4096
decoder_depth = 4
latent_lr = 1e-4
latent_wd = 1e-4
decoder_lr = 1e-4
decoder_wd = 1e-4
attribute_dropout_rate = 0.1
attribute_nn_width = 2048
attribute_nn_depth = 2
attribute_nn_lr = 1e-2
attribute_nn_wd = 4e-8
unknown_attribute_noise_param = 2e+1
unknown_attribute_penalty = 1e-1
cosine_scheduler = True
train_classifiers = False

loss_ae = "gauss"
n_latent_attribute_ordered = 256
n_latent_attribute_categorical = 3
reconstruction_penalty = 1e+4
cosine_scheduler = True
scheduler_final_lr = 1e-5
step_size_lr = 45
use_batch_norm = False
use_layer_norm = False


module_params = {
    "decoder_width": decoder_width,
    "decoder_depth": decoder_depth,
    "attribute_nn_width":  attribute_nn_width,
    "attribute_nn_depth": attribute_nn_depth,
    "use_batch_norm": use_batch_norm,
    "use_layer_norm": use_layer_norm,
    "unknown_attribute_noise_param": unknown_attribute_noise_param,
    "seed": 42,
    "n_latent_attribute_ordered": n_latent_attribute_ordered,
    "n_latent_attribute_categorical": n_latent_attribute_categorical,
    "loss_ae": loss_ae,
    "reconstruction_penalty": reconstruction_penalty,
    "unknown_attribute_penalty": unknown_attribute_penalty,
    "attribute_dropout_rate": attribute_dropout_rate,
}

trainer_params = {
    "n_epochs_warmup": 0,
    "latent_lr": latent_lr,
    "latent_wd": latent_wd,
    "decoder_lr": decoder_lr,
    "decoder_wd": decoder_wd,
    "attribute_nn_lr": attribute_nn_lr,
    "attribute_nn_wd": attribute_nn_wd,
    "step_size_lr": step_size_lr,
    "cosine_scheduler": cosine_scheduler,
    "scheduler_final_lr": scheduler_final_lr,
}
