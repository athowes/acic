#%%
import pandas as pd
import numpy as np
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression

# %%
results_arr = np.zeros((5, 3, 5))
value_names = ["atet", "regress", "regress with covs", "g formula"]

for param_idx, param in enumerate([30, 40, 50]):
    for seed_idx, seed in enumerate(range(83, 88)):
        # read in data
        data_df = pd.read_csv(
            f"/home/ghalebik/projects/acic/adam/data/{seed}-{param}.csv"
        )

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
        scaler = StandardScaler()
        x = scaler.fit_transform(x)

        # compute true atet
        atet = (data_df["y.1"] - data_df["y.0"])[z == 1].mean()
        results_arr[seed_idx, param_idx, 0] = atet

        # %%
        # estimate atet with regression without covariates
        reg = LinearRegression().fit(z[:, None], y[:, None])
        results_arr[seed_idx, param_idx, 1] = reg.coef_[0][0]

        # %%
        # estimate atet with regression with covariates
        reg = LinearRegression().fit(np.concatenate((z[:, None], x), 1), y[:, None])
        results_arr[seed_idx, param_idx, 2] = reg.coef_[0][0]

        # %%
        # estimate ate with g-formula
        x_cols = [f"t_{i}" for i in range(x.shape[1])]
        post_data_df = pd.DataFrame(
            np.concatenate((y[:, None], z[:, None], x), 1), columns=["y", "z"] + x_cols
        )

        g = TimeFixedGFormula(
            post_data_df, exposure="z", outcome="y", outcome_type="normal"
        )
        g.outcome_model(model="z + " + " + ".join(x_cols))
        g.fit(treatment="none")
        r_none = g.marginal_outcome
        g.fit(treatment="all")
        r_all = g.marginal_outcome
        g.fit(treatment="(g['z']==0)")
        r_untr = g.marginal_outcome
        perc_treated = post_data_df["z"].mean()
        # results_arr[seed_idx, param_idx, 2] = ((r_all - r_none) - r_untr * (1 - perc_treated)) / perc_treated
        results_arr[seed_idx, param_idx, 3] = r_all - r_none

        # estimate ate with AIPTW
        sdr = AIPTW(post_data_df, exposure="z", outcome="y")
        sdr.exposure_model(" + ".join(x_cols))
        sdr.outcome_model("z + " + " + ".join(x_cols))
        sdr.fit()
        results_arr[seed_idx, param_idx, 4] = sdr.average_treatment_effect

        # estimate ate with IPW
        sdr = IPTW(post_data_df, exposure="z", outcome="y")
        sdr.exposure_model(" + ".join(x_cols))
        sdr.outcome_model("z + " + " + ".join(x_cols))
        sdr.fit()
        results_arr[seed_idx, param_idx, 4] = sdr.average_treatment_effect

# %%
# estimate atet with regression with covariates
import zepid as ze
from zepid.causal.gformula import TimeFixedGFormula
from zepid.causal.doublyrobust import AIPTW
from zepid.causal.ipw import IPTW
import supylearner as sl

from zepid import load_sample_data, spline, RiskDifference
from zepid.causal.gformula import TimeFixedGFormula, SurvivalGFormula
from zepid.causal.ipw import IPTW, IPMW
from zepid.causal.snm import GEstimationSNM

from zepid.causal.doublyrobust import TMLE


def SuPyFitter(X, y):
    svm = SVC(kernel="linear", probability=True, random_state=101)
    log1 = LogisticRegression(penalty="l1", random_state=201)
    log2 = LogisticRegression(penalty="l2", random_state=103)
    randf = RandomForestClassifier(random_state=141)
    adaboost = AdaBoostClassifier(random_state=505)
    bayes = GaussianNB()
    lib = [svm, log1, log2, randf, adaboost, bayes]
    libnames = ["SVM", "Log_L1", "Log_L2", "Random Forest", "AdaBoost", "Bayes"]
    sl = supylearner.SuperLearner(lib, libnames, loss="nloglik", K=10)
    sl.fit(X, y)
    sl.summarize()
    return sl


import numpy as np
import supylearner
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import (
    RandomForestClassifier,
    AdaBoostClassifier,
)  # Random Forest, AdaBoost
from sklearn.naive_bayes import GaussianNB


x_cols = [f"t_{i}" for i in range(x.shape[1])]
post_data_df = pd.DataFrame(
    np.concatenate((y[:, None], z[:, None], x), 1), columns=["y", "z"] + x_cols
)

g = TimeFixedGFormula(post_data_df, exposure="z", outcome="y", outcome_type="normal")
g.outcome_model(model="z + " + " + ".join(x_cols))
g.fit(treatment="none")
r_none = g.marginal_outcome
g.fit(treatment="all")
r_all = g.marginal_outcome

