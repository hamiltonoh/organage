from importlib import resources
import pickle
import json
import dill
import pandas as pd
import numpy as np

# Class for OrganAge
class CreateOrganAgeObject:

    # init method or constructor
    def __init__(self,
                 path_scale_dict='v4_to_v4.1_scale_dict.json',
                 path_seed_dict='Bootstrap_and_permutation_500_seed_dict.json',
                 path_organ_plist_dict='tissue_pproteinlist_5k_dict_gtex_tissue_enriched_fc4_stable_proteins_seqid.json',
                 path_models_dir='data/ml_models/Covance/Zprot_stableps_perf95/'
                 ):

        self.data_and_model_paths = {"path_scale_dict": path_scale_dict,
                                     "path_seed_dict": path_seed_dict,
                                     "path_organ_plist_dict": path_organ_plist_dict,
                                     "path_models_dir": path_models_dir}
        self.load_data_and_models()

    def load_data_and_models(self):

        # Seqid:scale_factor dictionary
        scale_dict = json.load(resources.open_text("organage.data", self.data_and_model_paths["path_scale_dict"]))
        self.scale_dict = scale_dict

        # 500 bootstrap seed list
        all_seeds = json.load(resources.open_text("organage.data", self.data_and_model_paths["path_seed_dict"]))["BS_Seed"]
        self.all_seeds = all_seeds

        # organ:proteinlist dictionary
        organ_plist_dict = json.load(resources.open_text("organage.data", self.data_and_model_paths["path_organ_plist_dict"]))
        self.organ_plist_dict = organ_plist_dict

        # models
        models_dict = {}
        for organ in organ_plist_dict:
            models_dict[organ] = {}
            models_dict[organ]["aging_models"] = []

            # load protein scaler
            fn_protein_scaler = 'Covance_all_based_protein_zscore_scaler.pkl'
            loaded_model = pickle.loads(resources.read_binary('organage.data.ml_models.Covance.Zprot_stableps_perf95.' + organ, fn_protein_scaler))
            models_dict[organ]["prot_scaler"] = loaded_model

            # age gap scaler
            fn_agegap_scaler = 'Covance_all_Zprot_stableps_perf95_lasso_'+organ+'_agegap_zscore_scaler.pkl'
            loaded_model = pickle.loads(resources.read_binary('organage.data.ml_models.Covance.Zprot_stableps_perf95.' + organ, fn_agegap_scaler))
            models_dict[organ]["agegap_scaler"] = loaded_model

            # age prediction lowess
            fn_agepred_lowess = 'Covance_all_Zprot_stableps_perf95_lasso_' + organ + '_age_prediction_lowess.dill'
            loaded_model = dill.loads(resources.read_binary('organage.data.ml_models.Covance.Zprot_stableps_perf95.' + organ, fn_agepred_lowess))
            models_dict[organ]["age_prediction_lowess"] = loaded_model

            # load all aging models
            for seed in all_seeds:
                fn_aging_model = 'Covance_all_Zprot_stableps_perf95_lasso_'+organ+'_seed'+str(seed)+'_aging_model.pkl'
                loaded_model = pickle.loads(resources.read_binary('organage.data.ml_models.Covance.Zprot_stableps_perf95.'+ organ, fn_aging_model))
                models_dict[organ]["aging_models"].append(loaded_model)

        # save to object
        self.models_dict = models_dict