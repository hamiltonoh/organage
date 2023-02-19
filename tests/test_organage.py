from organage.OrganAge import CreateOrganAgeObject
import pandas as pd
import numpy as np

def test_add_data():
    data = CreateOrganAgeObject()

    # sample metadata data with Age and Sex_F
    rownames = list(range(5))
    md_hot = pd.DataFrame({"Age":[50,50,60,70,80], "Sex_F":[0,1,1,0,1]}, index=rownames)

    # sample expression data with Adipose proteins
    rows = np.arange(400,601,50)
    df_prot=  pd.DataFrame([rows+i*10 for i in range(5)],
                           columns=['8484-24', '15386-7', '6578-29', '3554-24', '18830-1'],
                           index=rownames)

    data.add_data(md_hot, df_prot)

    assert data.md_hot == md_hot, "sample metadata loaded incorrectly"
    assert data.df_prot == df_prot, "sample metadata loaded incorrectly"


