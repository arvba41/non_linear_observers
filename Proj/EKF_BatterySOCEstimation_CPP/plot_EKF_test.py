# %% import packages

import numpy as py
import matplotlib.pyplot as plt
import json 
import pandas as pd
from seaborn import despine

# %% Read data
with open("build/result.json", "r") as f:
    results = json.load(f)

data = pd.DataFrame(results["data"])
spline = pd.DataFrame(results["spline"])
simdata = pd.DataFrame(results["simdata"])
columb_counting = pd.DataFrame(results["columb_counting"])
EKF = pd.DataFrame(results["EKF"])

# # simulated data
# with open("simdata.json", "r") as f:
#     simresults = json.load(f)

# simdata = pd.DataFrame(simresults["Simdata"])

# %% Plot a littleimport pandas as pd
plt.figure(10)
plt.scatter(data.SOC / 24, data.OCV, marker= 'x', c='r', label="OCV curve")
plt.plot(spline.SOC / 24,spline.OCV, label="OCV interp1")
plt.xlabel("SOC [%]")
plt.ylabel("OCV [V]")
plt.legend()

# %% plotting simulated data
plt.figure(11)
plt.subplot(311)
plt.plot(simdata.Time,simdata.SOC)
plt.xlabel("t [s]")
plt.ylabel("SOC [-]")

plt.subplot(312)
plt.plot(simdata.Time,simdata.Current)
plt.xlabel("t [s]")
plt.ylabel("Current [A]")

plt.subplot(313)
plt.plot(simdata.Time,simdata.Voltage)
plt.xlabel("t [s]")
plt.ylabel("Voltage [V]")

# %% plotting columb counting data
plt.figure(12)
plt.plot(simdata.Time,simdata.SOC, label="Simulated")
plt.plot(simdata.Time,columb_counting.SOC / 24, label="columb counting")
plt.plot(simdata.Time,EKF.SOC / 24, label="EKF")
plt.xlabel("t [s]")
plt.ylabel("SOC [-]")
plt.legend()

#%%
plt.show()
# %%