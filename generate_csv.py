'''
Моделирует показания чувствительных элементов при движении с заданными кглами ориентации (углами Эйлера-Крылова)
и с заданными скоростями. Все сохраняет в csv файл, разделенный пробелами (так C++ лучше читает)
'''
import numpy as np
import struct
import matplotlib.pyplot as plt

def matrix_o_b(H, R, P):
    '''Функция формирует матрицу перехода (МНК) из опорной систему координат (географической) в связанную'''
    '''
    H - heading0 (курс)
    R - Roll (крен)
    P - Pitch (тангаж)
    '''
    C = np.array([
        [np.cos(R) * np.cos(H) + np.sin(R)*np.sin(H)* np.sin(P), np.sin(R)*np.cos(H)*np.sin(P) - np.cos(R)*np.sin(H), -np.sin(R) * np.cos(P)],
        [np.cos(P)*np.sin(H), np.cos(H)*np.cos(P), np.sin(P)],
        [np.sin(R)*np.cos(H) - np.cos(R)*np.sin(H)*np.sin(P), -np.sin(R)*np.sin(H) - np.cos(R)*np.cos(H)*np.sin(P), np.cos(R)*np.cos(P)]
    ])
    return C

def make_portable_speed(phi):
#     portable -- переносная. Формирование переносной угловой скорости (т.е. скорости Земли)
#     В осях опорной системы координат
    Om_p = np.array([0, U*np.cos(phi), U*np.sin(phi)]);
    return Om_p


def cross_sim(A):
    '''Кососиметрическая матрица для угловых скоростей'''
    return np.array([0,-A[2], A[1],A[2],0,-A[0], -A[1],A[0],0]).reshape(3,3)

U = np.deg2rad(15)/3600;
g = 9.81;
Re = 6400e3; # радиус Земли в метрах
a = 6378245;
b = 6356856;
e=np.sqrt(1 - b*b/a/a);

freq = 100; # частота (измерений за 1 секунду)
t_nav = 90*60 # в секундах
num_samples = t_nav*freq;
time_to_alignment = 5*60; # время выставки в секундах
'''
C_b_n = np.array([[],
                  [],
                  []])
'''

Vabs = 0 # модуль конечной линейной скорости, м/с
phi0 = np.deg2rad(55)
lmbda0 = np.deg2rad(33)

tics_to_turn = 1 #5*60*freq # время в тактах на разворот
tics_start_turn = (t_nav + time_to_alignment)*freq//2 - tics_to_turn // 2# время в тактах (от подачи питания, т.е начала выставки) начала поворота (time_to_alihnment+ (t_nav - time_to_alignemrnt)/2)
tics_stop_turn = (t_nav + time_to_alignment)*freq//2 + (tics_to_turn // 2 + 1)# время в тактах (от подачи питания, т.е начала выставки) окончания поворота 
angle_turn = np.deg2rad(90) # угол разворота, в рад
om_turn = - angle_turn / (tics_to_turn / freq ) #угловая скорость поворота в проекции на местную вертикаль (положительное мзменение курса по часовой, следовательно угловая скрость отрицательная), рад/с
Om_turn = np.zeros(3); Om_turn[2] = om_turn; # массив угловых скоростей разворота

'''тип файлы, бинарный или текстовый'''
extention_out_file = "csv" # bin (для бинарного) или csv (для текстового)

'''задаем ориентацию объекта'''
heading0 = np.deg2rad(0)
roll = np.deg2rad(0);
pitch = np.deg2rad(0);

# file_name = f"data_acc_veloc_{Vabs}_heading_{int(np.rad2deg(heading0))}_R_{int(np.rad2deg(roll))}_P_{int(np.rad2deg(pitch))}_freq_{freq}_turn_V_coo_gps.{extention_out_file}"
file_name = f"horizontal_gyr_bias_{Vabs}_heading_{int(np.rad2deg(heading0))}_R_{int(np.rad2deg(roll))}_P_{int(np.rad2deg(pitch))}_freq_{freq}_V_coo_gps_fast_turn.{extention_out_file}"
dest_dir = "/home/nikita/Документы/C_Cpp_progs/InertialNavigation/Data_files/" #дериктория назначения
C_n_b = matrix_o_b(heading0, roll, pitch)

'''Систематические дрейфы'''
bias_acc = np.array([1e-4 for _ in range(3)]);
bias_gyr = np.array([np.deg2rad(0.05)/3600, np.deg2rad(0.05)/3600, np.deg2rad(0.05)/3600]);
'''Случайные дрейфы'''
T_k_a = 1
beta_acc = 1/T_k_a
std_acc=0.1*g*1e-3
T_k_g = 2
beta_gyr = 1/T_k_g
std_gyr = np.deg2rad(0.02)/3600 #rad in sec 

mean = 0;

white_noise_acc= np.random.normal(mean, 1, size=(3, num_samples))
    
white_noise_gyr = np.random.normal(mean, 1, size=(3, num_samples))

colour_noise_acc = np.zeros((3, num_samples))
colour_noise_gyr = np.zeros((3, num_samples))


sko_gnss_n = 0.05 #in meter in seconds
sko_gnss_d = sko_gnss_n*np.sqrt(freq)
sko_gnss_pos_n = 0.2 #in meters
sko_gnss_pos_d = sko_gnss_pos_n/Re*np.sqrt(freq)

'''формирующий фильтр для экспоненциально косинусоидальной корреляционной функции'''
for i in range(1,num_samples):
    for jj in range(3):
        #акселерометры
        colour_noise_acc[jj][i] = colour_noise_acc[jj][i-1] * (1 - beta_acc/freq) + std_acc*np.sqrt(2*beta_acc/freq)*white_noise_acc[jj][i-1]
        #гироскопы  
        colour_noise_gyr[jj][i] = colour_noise_gyr[jj][i-1] * (1 - beta_gyr/freq) + std_gyr*np.sqrt(2*beta_gyr/freq)*white_noise_gyr[jj][i-1]

'''моделируем идеальный курс после начала навигации'''
# До поворота (включая время выставки) он постоянный и равный heading0, во время поворота он меняется с угловой скоростью om_turn; далее курс опять постоянный, равен heading0 + angle_turn
# heading = np.concatenate((np.array([heading0 for _ in range(0, tics_start_turn*freq)]), np.array([heading0 - i/freq*om_turn for i in range(0, (tics_stop_turn - tics_start_turn )*freq)]), np.array([heading0+angle_turn for _ in range(tics_stop_turn* freq, num_samples)]) )) #heading за все время моделирования (включая время выставки)
heading = np.array([heading0 for _ in range(0, tics_start_turn)] + [heading0 - i/freq*om_turn for i in range(0, (tics_stop_turn - tics_start_turn ))] + [heading0+angle_turn for _ in range(tics_stop_turn, num_samples)]) # так должно работать быстрее

'''Моделируем измерения в связанной системе координат и записываем в файл'''
# Ve = np.array([Vabs*np.sin(heading[iii]) for iii in range(time_to_alignment*freq, (t_nav)*freq)]) # np.linspace(Vabs*np.sin(heading0), Vabs*np.sin(heading0), (t_nav - time_to_alignment)*freq)
# Vn = np.array([Vabs*np.cos(heading[iii]) for iii in range(time_to_alignment*freq, (t_nav)*freq)]) # np.linspace(Vabs*np.cos(heading0), Vabs*np.cos(heading0), (t_nav - time_to_alignment)*freq)
# Ve = np.concatenate((np.array([0 for _ in range(time_to_alignment*freq)]), Ve))
# Vn = np.concatenate((np.array([0 for _ in range(time_to_alignment*freq)]), Vn))
'''так должно работать быстрее'''
Ve = [Vabs*np.sin(heading[iii]) for iii in range(time_to_alignment*freq, (t_nav)*freq)] # np.linspace(Vabs*np.sin(heading0), Vabs*np.sin(heading0), (t_nav - time_to_alignment)*freq)
Vn = [Vabs*np.cos(heading[iii]) for iii in range(time_to_alignment*freq, (t_nav)*freq)] # np.linspace(Vabs*np.cos(heading0), Vabs*np.cos(heading0), (t_nav - time_to_alignment)*freq)
Ve = np.array([0 for _ in range(time_to_alignment*freq)] + Ve)
Vn = np.array([0 for _ in range(time_to_alignment*freq)] + Vn)



allow_dA = True # ускорения как дифференциал скорости

# Изменение высоты полета
height_start = 0;
height_stop = 0;
# Height = np.linspace(0, 0, (t_nav - time_to_alignment)*freq)
# Height = np.concatenate((np.zeros(time_to_alignment*freq), Height))
Height = [height_start + (height_stop - height_start)/((t_nav - time_to_alignment)*freq - 1)*i for i in range((t_nav - time_to_alignment)*freq)] #изменение ввысоты от начального до конечного значения на всем протяжении полета
Height = [0 for _ in range(time_to_alignment*freq)] + Height # значение высоты вместе с выставкой

# Изменение широты
phi = [phi0] 
lmbda = [lmbda0]

for iii in range(1, num_samples): # цикл от 1 (практически с начала)
    phi.append(phi[-1] + Vn[iii]/(Re + Height[iii])/freq)
    lmbda.append(lmbda[-1] + Ve[iii]/((Re + Height[iii])*np.cos(phi[-1]))/freq)

round = 12; # количество знаков после запятой, с какой округлять и выводить в файл
if (extention_out_file == "bin"):
    type_open_file = "wb"
else:# в любом случае, чтобы создался файл и данные записались
    type_open_file = "wt" 
with (open(dest_dir + file_name, type_open_file) as file):
        '''
        последовательность данных
        ["Abx", "Aby", "Abz", "Ombx", "Omby", "Ombz", "biasAbx", "biasAby", "biasAbz", "biasOmbx", "biasOmby", "biasOmbz", "randAbx", "randAby", "randAbz", "randOmbx", "randOmby", "randOmbz", "Vgps_x", "Vgps_y", "randVgps_x", "randVgps_y", "phi_gps", "lambda_gps", "randPhi_gps", "randLambda_gps"]
        '''
        if (extention_out_file == "bin"):
            s = struct.Struct("<26d");
        for itr in range(0, num_samples-1):
            Rlambda = Re# / np.sqrt(1 - e **2 * np.sin(phi[itr]) **2);
            Rphi = Re# * (1 - e**2) / (np.sqrt(1 - e**2 * np.sin(phi[itr]) **2) * (1 - e**2 * np.sin(phi[itr]) **2));
            
            A_o = np.array([0.0, 0.0, g])

            Om_e = make_portable_speed(phi[itr]); # формирование скорости вращения Земли на широте

            # Дифференцирование скоростей и прибавление к ускорению
            if (itr>time_to_alignment*freq):
                dA = np.array([(Ve[itr + 1] - Ve[itr]), (Vn[itr + 1] - Vn[itr]),
                            (Height[itr + 1] - Height[itr])]) * freq
                # Изменение угловых скоростей от движения
                dOm_or = np.array([-Vn[itr] / (Rphi + Height[itr]), #
                                Ve[itr] / ( (Rlambda + Height[itr]) ), #* np.cos(phi[itr] )
                                Ve[itr] * np.tan(phi[itr]) / ((Rlambda + Height[itr]) ) ]) # r -- realtive (относительная)
                '''
                Ускорение Кориолиса равно удвоенному векторному произведению угловой 
                скорости переносного движения на относительную скорость точки (видимо имеется в виду линейную скорость)
                Только в случае с инерциальной навигацией "относительная" должна пониматься как относительная относительно
                инерциального пространства
                '''
                if ( (itr >= (tics_start_turn) ) and (itr < (tics_stop_turn) ) ): #во время движения. есть поворот
                    allow_turn = True; # добавление угловой скорости поворота
                    allow_centrifugal = True; # добавление центробедных сил 
                    C_n_b = matrix_o_b(heading[itr], roll, pitch);
                else: # во время движения. нет поворота
                    allow_turn = False; # добавление угловой скорости поворота
                    allow_centrifugal = False; # добавление центробедных сил 
                Coriolise = np.cross((2*Om_e + dOm_or + 0 * allow_turn * Om_turn), np.array([Ve[itr], Vn[itr], 0])); #+  np.cross(Om_e, np.array([Ve[itr], Vn[itr], 0])) 
            else: # время выставка, движения никакого нет, следовательно нет ни угловой скорости поворота, ни центробежных сил
                allow_turn = False; 
                allow_centrifugal = False; # добавление центробедных сил 
                dA = np.zeros(3)
                dOm_or = np.zeros(3) # r -- realtive (относительная)
                Coriolise = np.zeros(3);

            Om_b = C_n_b @ (Om_e + dOm_or + allow_turn * Om_turn);
            A_b = C_n_b @ (A_o + allow_dA * dA + Coriolise + ( - 0 * allow_centrifugal*np.array([0, om_turn*Vabs, 0])) ); 
            C_n_b = matrix_o_b(heading[itr], roll, pitch) # 
            if (extention_out_file == "csv"):
                '''Запись в файл в текстовом виде'''
                for i in range(3):
                    file.write(f"{np.round(A_b[i], round)};")
                for i in range(3):
                    file.write(f"{np.round(Om_b[i], round)};")
                
                for i in range(3):
                    file.write(f"{np.round(bias_acc[i], round)};") # Смещение нулей акселерометров
                for i in range(3):
                    file.write(f"{np.round(bias_gyr[i], round)};") # Смещение нулей гироскопов в рад/с
                for i in range(3):
                    file.write(f"{np.round(colour_noise_acc[i][itr], round)};") # Случайные погрешности акселерометров
                for i in range(3):
                    file.write(f"{np.round(colour_noise_gyr[i][itr], round)};") # Случайные погрешности гироскопов

                file.write(f"{np.round(Ve[itr], round)};") # Восточная скорость (идеал, будет воспринята как GPS)
                
                file.write(f"{np.round(Vn[itr], round)};") # Северная скорость (идеал, будет воспринята как GPS)
                
                file.write(f"{np.round(np.random.normal(loc = 0, scale = sko_gnss_d), round)};") # Один белый гауссовский шум для скорости GPS
                file.write(f"{np.round(np.random.normal(loc = 0, scale = sko_gnss_d), round)};") # Второй белый гауссовский шум для скорости GPS

                file.write(f"{np.round(phi[itr], round)};") # Широта в рад (идеал, будет воспринята как GPS)
                
                file.write(f"{np.round(lmbda[itr], round)};") # Долгота (идеал, будет воспринята как GPS)
                
                file.write(f"{np.round(np.random.normal(loc = 0, scale = sko_gnss_pos_d), round)};") # Один белый гауссовский шум для координат GPS
                file.write(f"{np.round(np.random.normal(loc = 0, scale = sko_gnss_pos_d), round)};") # Второй белый гауссовский шум для координат GPS
                
                file.write("\n")
                
            elif (extention_out_file == "bin"):
                
                data = s.pack(*A_b, *Om_b, bias_acc, bias_acc, bias_acc, bias_gyr, bias_gyr, bias_gyr, *colour_noise_acc[:, itr], *colour_noise_gyr[:, itr], Ve[itr], Vn[itr], np.random.normal(loc = 0, scale = sko_gnss_d), np.random.normal(loc = 0, scale = sko_gnss_d), phi[itr], lmbda[itr], np.random.normal(loc = 0, scale = sko_gnss_pos_d), np.random.normal(loc = 0, scale = sko_gnss_pos_d)) # формирование пакета (структуры)

                file.write(data); # запись этой структуры в файл
                
        
print(f"Generation done: {dest_dir + file_name}")
