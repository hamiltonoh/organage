from organage.OrganAge import CreateOrganAgeObject
import pandas as pd

def test_OrganAge():
    data = CreateOrganAgeObject()
    md_hot = pd.read_csv("tests/md_hot.csv").set_index("ID")
    df_prot = pd.read_csv("tests/df_prot.csv").set_index("ID")

    # sample metadata data with Age and Sex_F
    data.add_data(md_hot, df_prot)
    data.normalize(assay_version="v4.1")

    dfres = data.estimate_organ_ages()
    return dfres
