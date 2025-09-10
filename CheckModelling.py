import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("Data_files/data_acc_veloc_0_heading_50_R_0_P_0_freq_100_turn_V_coo_gps.csv",  delimiter=";", index_col = False, header = None)

data.columns = ["Abx", "Aby", "Abz", "Ombx", "Omby", "Ombz", "biasAbx", "biasAby", "biasAbz", "biasOmbx", "biasOmby", "biasOmbz", "randAbx", "randAby", "randAbz", "randOmbx", "randOmby", "randOmbz", "Vgps_x", "Vgps_y", "randVgps_x", "randVgps_y", "phi_gps", "lambda_gps", "randPhi_gps", "randLambda_gps", "bla-bla"]

data.drop("bla-bla", axis = 1)

fig1, (fig1_ax1) = plt.subplots(ncols=1, nrows=1)
fig1.suptitle("Истинные скорости (GPS без погрешностей)")
# fig1_ax1.set_title("");
fig1_ax1.plot(data["Vgps_x"])
fig1_ax1.plot(data["Vgps_y"])
fig1_ax1.grid(True)

# fig3, (fig3_ax1) = plt.subplots(ncols = 1, nrows = 1)
# fig3.suptitle("Углы ориентации (Эйлера)")
# fig3_ax1.plot(np.rad2deg(data["heading"]))

fig5, (fig5_ax1, fig5_ax2, fig5_ax3) = plt.subplots(ncols=1, nrows=3)
fig5.suptitle("Показания ДУС")
fig5_ax1.set_title("$\Omega_{bx}$")
fig5_ax1.plot(np.rad2deg(data["Ombx"])*3600)
fig5_ax1.grid(True)
# 
fig5_ax2.set_title("$\Omega_{by}$")
fig5_ax2.plot(np.rad2deg(data["Omby"])*3600)
fig5_ax2.grid(True)
# 
fig5_ax3.set_title("$\Omega_{bz}$")
fig5_ax3.plot(np.rad2deg(data["Ombz"])*3600)
fig5_ax3.grid(True)

plt.show()