import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

path = sys.argv[1]
data = pd.read_csv(path,  delimiter=";", index_col = False, header = None)
# разбиваем путь к файлу по разделителю дерикторий и разбиваем только имя файла (так как оно по-любому должно идти в конце) на нижние подчеркивания
splitted_path = path.split("/")[-1].split("_")
num_freq = splitted_path.index("freq") # ищем порядковый номер слова "freq", следующий за ним будет значение частоты моделирования данных
freq = int(splitted_path[num_freq + 1])

data.columns = ["Abx", "Aby", "Abz", "Ombx", "Omby", "Ombz", "biasAbx", "biasAby", "biasAbz", "biasOmbx", "biasOmby", "biasOmbz", "randAbx", "randAby", "randAbz", "randOmbx", "randOmby", "randOmbz", "Vgps_x", "Vgps_y", "randVgps_x", "randVgps_y", "phi_gps", "lambda_gps", "randPhi_gps", "randLambda_gps", "bla-bla"]

data.drop("bla-bla", axis = 1)

time = np.linspace(0, len(data)/freq/60, len(data))

fig1, (fig1_ax1) = plt.subplots(ncols=1, nrows=1)
fig1.suptitle("Истинные скорости (GPS без погрешностей)")
# fig1_ax1.set_title("");
fig1_ax1.plot(time, data["Vgps_x"], label = "Vgps_x")
fig1_ax1.plot(time, data["Vgps_y"], label = "Vgps_y")
fig1_ax1.set_xlabel("мин")
fig1_ax1.set_ylabel("м/с")
fig1_ax1.legend(loc = 'best');
fig1_ax1.grid(True)

# fig3, (fig3_ax1) = plt.subplots(ncols = 1, nrows = 1)
# fig3.suptitle("Углы ориентации (Эйлера)")
# fig3_ax1.plot(np.rad2deg(data["heading"]))

fig5, (fig5_ax1, fig5_ax2, fig5_ax3) = plt.subplots(ncols=1, nrows=3)
fig5.suptitle("Показания ДУС")
fig5_ax1.set_title("$\Omega_{bx}$")
fig5_ax1.plot(time, np.rad2deg(data["Ombx"])*3600)
fig5_ax1.set_xlabel("мин")
fig5_ax1.set_ylabel("град/ч")
fig5_ax1.grid(True)
# 
fig5_ax2.set_title("$\Omega_{by}$")
fig5_ax2.plot(time, np.rad2deg(data["Omby"])*3600)
fig5_ax2.set_xlabel("мин")
fig5_ax2.set_ylabel("град/ч")
fig5_ax2.grid(True)
# 
fig5_ax3.set_title("$\Omega_{bz}$")
fig5_ax3.plot(time, np.rad2deg(data["Ombz"])*3600)
fig5_ax3.set_xlabel("мин")
fig5_ax3.set_ylabel("град/ч")
fig5_ax3.grid(True)


# Ускорения с акселерометров (идеальный вариант)
fig10, (fig10_ax1, fig10_ax2, fig10_ax3) = plt.subplots(ncols=1, nrows=3)
fig10.suptitle("Показания акслереометров")
fig10_ax1.set_title("$A_{bx}$")
fig10_ax1.plot(time, data["Abx"])
fig10_ax1.set_xlabel("мин")
# fig10_ax1.set_ylabel("$\dfrac{\\text{м}}{\\text{с}^2}$")
fig10_ax1.grid(True)
# 
fig10_ax2.set_title("$A_{by}$")
fig10_ax2.plot(time, data["Aby"])
fig10_ax2.set_xlabel("мин")
# fig10_ax2.set_ylabel("$\dfrac{\\text{м}}{\\text{с}^2}$")
fig10_ax2.grid(True)
# 
fig10_ax3.set_title("$A_{bz}$")
fig10_ax3.plot(time, data["Abz"])
fig10_ax3.set_xlabel("мин")
# fig10_ax3.set_ylabel("$\dfrac{\\text{м}}{\\text{с}^2}$")
fig10_ax3.grid(True)

plt.show()