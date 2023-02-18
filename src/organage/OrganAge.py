from importlib import resources
import os
import pickle
import json
import dill
import pandas as pd
import numpy as np

# Class for OrganAge
class CreateOrganAgeObject:

    # init method or constructor
    def __init__(self,
                 path_scale_dict='data/v4_to_v4.1_scale_dict.json',
                 path_organ_plist_dict='data/tissue_pproteinlist_5k_dict_gtex_tissue_enriched_fc4_stable_proteins_seqid.json',
                 path_seed_dict='data/Bootstrap_and_permutation_500_seed_dict.json',
                 path_models_dir='data/ml_models/Covance/Zprot_stableps_perf95/'
                 ):

        self.model_fps = {"path_scale_dict":path_scale_dict,
                          "path_organ_plist_dict":path_organ_plist_dict,
                          "path_seed_dict":path_seed_dict,
                          "path_models_dir":path_models_dir}
        self.load_models()

    def load_models(self):

        # Seqid:scale_factor dictionary
        scale_dict = json.load(resources.open_text("organage.data", self.model_fps.model_fps["path_scale_dict"]))

        # organ:proteinlist dictionary
        scale_dict = json.load(resources.open_text("organage.data", self.model_fps.model_fps["path_organ_plist_dict"]))

        # 500 bootstrap seed list
        scale_dict = json.load(resources.open_text("organage.data", self.model_fps.model_fps["path_seed_dict"]))

