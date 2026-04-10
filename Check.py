import pandas as pd
import sys

file = sys.argv[1]

data = pd.read_csv(file, sep = ";")
data.columns =  header=["Ax", "Ay", "Az", "Omx", "Omy", "Omz", "Abx", "Aby", "Abz", "Ombx", "Omby", "Ombz", 
                        "Arx", "Ary", "Arz", "Omrx", "Omry", "Omrz", "Ve", "Vn", "Vre", "Vrn", "Lat","Lon", "RLat", "RLon", None]

print(f"Дисперсия случайных ошибок акселерометров:");
print(f"Ax: {data["Arx"].var()}");
print(f"Ay: {data["Ary"].var()}");
print(f"Az: {data["Arz"].var()}");

print(f"Дисперсия случайных ошибок гироскопов:");
print(f"Omx: {data["Omrx"].var()}");
print(f"Omy: {data["Omry"].var()}");
print(f"Omz: {data["Omrz"].var()}");

print(f"Дисперсия случайных ошибок ГНСС:");
print(f"Vre: {data["Vre"].var()}");
print(f"Vrn: {data["Vrn"].var()}");
print(f"Lattitude: {data["RLat"].var()}");
print(f"Longitude: {data["RLon"].var()}");