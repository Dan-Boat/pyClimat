# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 15:00:08 2023

@author: dboateng
Analyse the latitudinal and longitudinal gradient of MAP, MAT, and isotopes accross the Alps. 
Compare the results with proxies (check for new compilations or better rely on current one)
"""

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

np.random.seed(44)
# DATAFRAMES WITH TRIAL COLUMN ASSIGNED
df1 = pd.DataFrame(np.random.randn(5,5), columns=list(range(1,6))).assign(Trial=1)
df2 = pd.DataFrame(np.random.randn(5,5), columns=list(range(1,6))).assign(Trial=2)
df3 = pd.DataFrame(np.random.randn(5,5), columns=list(range(1,6))).assign(Trial=3)

cdf = pd.concat([df1, df2, df3])                                # CONCATENATE
mdf = pd.melt(cdf, id_vars=['Trial'], var_name=['Number'])      # MELT

print(mdf.head())
#    Trial Number     value
# 0      1      1 -0.750615
# 1      1      1 -1.715070
# 2      1      1 -0.963404
# 3      1      1  0.360856
# 4      1      1 -1.190504

ax = sns.boxplot(x="Trial", y="value", hue="Number", data=mdf)  # RUN PLOT   
plt.show()

plt.clf()
plt.close()