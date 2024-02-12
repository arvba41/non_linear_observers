# %%
# %matplotlib
# %load_ext autoreload
# %autoreload 2

import numpy as np
import matplotlib.pyplot as plt
import json
import pandas as pd
from seaborn import despine


# %% Example functions
def fun(x):
    return 1 / (1 + np.exp(-x))  # sigmoid


def dfun(x):
    return np.exp(-x) / (1 + np.exp(-x)) ** 2  # derivative of sigmoid


# %% Read data
with open("result.json", "r") as f:
    results = json.load(f)

data = pd.DataFrame(results["data"])
spline = pd.DataFrame(results["spline"])


# %% Plot a little
_, ax = plt.subplots(num=10, clear=True)
ax.plot(spline.x, fun(spline.x), label="True sigmoid")
ax.plot(data.x, data.y, "ro", label="Interpolation points")
ax.plot(spline.x, dfun(spline.x), label="True derivative")
ax.legend(frameon=False)
despine(ax=ax)

_, ax = plt.subplots(num=20, clear=True)
ax.plot(spline.x, fun(spline.x), label="True sigmoid")
ax.plot(spline.x, spline.y, label="Interpolated sigmoid")
ax.legend(frameon=False)
despine(ax=ax)

_, ax = plt.subplots(num=30, clear=True)
ax.plot(spline.x, dfun(spline.x), label="True derivative")
ax.plot(spline.x, spline.dy, label="Interpolated derivative")
ax.legend(frameon=False)
despine(ax=ax)


# %%
plt.show()
