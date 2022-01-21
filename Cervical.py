#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.feature_selection import mutual_info_regression
from sklearn.decomposition import PCA

data = pd.read_csv('sobar-72.csv')
data.head()

print(" Este arquivo contém {0} linhas e {1} colunas.".format(data.shape[0], data.shape[1]))
data.isnull().any()

print("Dados duplicados: {0}".format(data.duplicated().sum()))

data.describe()
