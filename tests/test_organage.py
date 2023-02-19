from organage.OrganAge import CreateOrganAgeObject
import pandas as pd
import numpy as np

def test_add_data():
    data = CreateOrganAgeObject()

    md_hot = pd.read_csv("tests/md_hot.csv").set_index("Barcode")
    df_prot = pd.read_csv("tests/df_prot.csv").set_index("Barcode")

    # sample metadata data with Age and Sex_F
    data.add_data(md_hot, df_prot)
    assert len(data.md_hot) == 5, "user data loaded incorrectly"


