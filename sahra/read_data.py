#%%
import pandas as pd
import numpy as np

# %%
# read in data
data_df = pd.read_csv("/home/ghalebik/projects/acic/adam/87-50.csv")

from sklearn.preprocessing import OneHotEncoder

z = data_df.z.values
y = data_df.y.values

# one hot encode categorical categories
object_cols = data_df.columns[(data_df.dtypes == "object")].values
x = data_df[
    [col for col in data_df.columns if "x_" in col and col not in object_cols]
].values
enc = OneHotEncoder(drop="first", sparse=False)
x_onehot = enc.fit_transform(data_df[object_cols])
x = np.concatenate((x, x_onehot), 1)

# scale data
from sklearn.preprocessing import StandardScaler

scaler = StandardScaler()
x = scaler.fit(x)

# %%
# compute true atet
atet = (data_df["y.1"] - data_df["y.0"])[z == 1].mean()
print(atet)

# %%
# estimate atet with regression without covariates
from sklearn.linear_model import LinearRegression

reg = LinearRegression().fit(z[:, None], y[:, None])
reg.score(z[:, None], y)
reg.coef_[0][0]

# %%
# estimate atet with regression with covariates
from sklearn.linear_model import LinearRegression

reg = LinearRegression().fit(np.concatenate((z[:, None], x), 1), y[:, None])
reg.predict(np.concatenate((z[:, None], x), 1), y)
reg.coef_[0][0]
