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
                 path_version_scale_factors='v4_to_v4.1_scale_dict.json',
                 path_organ_plist_dict='tissue_pproteinlist_5k_dict_gtex_tissue_enriched_fc4_stable_proteins_seqid.json',
                 path_bootstrap_seeds='Bootstrap_and_permutation_500_seed_dict.json',
                 path_models_dir='data/ml_models/Covance/Zprot_stableps_perf95/'
                 ):

        self.data_and_model_paths = {"path_version_scale_factors": path_version_scale_factors,
                                     "path_organ_plist_dict": path_organ_plist_dict,
                                     "path_bootstrap_seeds": path_bootstrap_seeds,
                                     "path_models_dir": path_models_dir
                                     }
        self.load_data_and_models()
        del self.data_and_model_paths


    def load_data_and_models(self):

        # Seqid:scale_factor dictionary
        version_scale_factors = json.load(resources.open_text("organage.data", self.data_and_model_paths["path_version_scale_factors"]))
        self.version_scale_factors = version_scale_factors

        # organ:proteinlist dictionary
        organ_plist_dict = json.load(resources.open_text("organage.data", self.data_and_model_paths["path_organ_plist_dict"]))
        self.organ_plist_dict = organ_plist_dict

        # 500 bootstrapped models
        bootstrap_seeds = json.load(resources.open_text("organage.data", self.data_and_model_paths["path_bootstrap_seeds"]))["BS_Seed"]
        models_dict = {}
        for organ in organ_plist_dict:
            models_dict[organ] = {}
            models_dict[organ]["aging_models"] = []

            # load protein zscore scaler
            fn_protein_scaler = 'Covance_all_based_'+organ+'_protein_zscore_scaler.pkl'
            loaded_model = pickle.loads(resources.read_binary('organage.data.ml_models.Covance.Zprot_stableps_perf95.' + organ, fn_protein_scaler))
            models_dict[organ]["prot_scaler"] = loaded_model

            # age gap zscore scaler
            fn_agegap_scaler = 'Covance_all_Zprot_stableps_perf95_lasso_'+organ+'_agegap_zscore_scaler.pkl'
            loaded_model = pickle.loads(resources.read_binary('organage.data.ml_models.Covance.Zprot_stableps_perf95.' + organ, fn_agegap_scaler))
            models_dict[organ]["agegap_scaler"] = loaded_model

            # age prediction lowess
            fn_agepred_lowess = 'Covance_all_Zprot_stableps_perf95_lasso_' + organ + '_age_prediction_lowess.dill'
            loaded_model = dill.loads(resources.read_binary('organage.data.ml_models.Covance.Zprot_stableps_perf95.' + organ, fn_agepred_lowess))
            models_dict[organ]["age_prediction_lowess"] = loaded_model

            # load all aging models
            for seed in bootstrap_seeds:
                fn_aging_model = 'Covance_all_Zprot_stableps_perf95_lasso_'+organ+'_seed'+str(seed)+'_aging_model.pkl'
                loaded_model = pickle.loads(resources.read_binary('organage.data.ml_models.Covance.Zprot_stableps_perf95.'+ organ, fn_aging_model))
                models_dict[organ]["aging_models"].append(loaded_model)

        # save to object
        self.models_dict = models_dict


    def add_data(self, md_hot, df_prot):
        # user inputs sample metadata and sample expression dataframes
        self.md_hot = md_hot
        self.df_prot = df_prot


    def normalize(self, assay_version="v4.1"):
        # normalizing protein levels
        df_prot_norm = self.df_prot.copy()
        if assay_version == "v4":
            for prot in df_prot_norm.columns:
                df_prot_norm[prot] = df_prot_norm[prot] * self.scale_dict[prot]
        # log
        df_prot_norm = np.log10(df_prot_norm)
        self.df_prot_norm = df_prot_norm