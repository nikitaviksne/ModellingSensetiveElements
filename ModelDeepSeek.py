import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.spatial.transform import Rotation

class INSSimulatorENU:
    def __init__(self):
        # Параметры модели
        self.g = 9.81  # ускорение свободного падения, м/с²
        self.dt = 0.01  # шаг интегрирования, с
        self.t_end = 300  # время моделирования, с
        
        # Параметры Земли
        self.omega_earth = 7.292115e-5  # угловая скорость вращения Земли, рад/с
        self.R_earth = 6378137.0  # экваториальный радиус Земли, м
        self.e = 0.08181919  # эксцентриситет Земли
        
        # Начальные координаты (широта, долгота, высота)
        self.lat0 = np.radians(55.7558)  # Москва
        self.lon0 = np.radians(37.6173)
        self.h0 = 1000  # начальная высота, м
        
        # Параметры ЛА
        self.mass = 1000  # масса, кг
        self.Ixx = 1000   # момент инерции по оси X, кг*м²
        self.Iyy = 2000   # момент инерции по оси Y, кг*м²
        self.Izz = 1500   # момент инерции по оси Z, кг*м²
        
        # Истинные уходы гироскопов (рад/с)
        self.true_gyro_biases = np.array([0.01, -0.005, 0.007])
        
        # Истинные смещения акселерометров (м/с²)
        self.true_accel_biases = np.array([0.02, -0.01, 0.015])
        
        # Начальные условия в системе ENU
        self.init_state = np.zeros(12)
        self.init_state[0] = 0      # начальная позиция по востоку, м
        self.init_state[1] = 0      # начальная позиция по северу, м
        self.init_state[2] = 0      # начальная позиция по высоте, м
        self.init_state[3] = 100    # начальная скорость на восток, м/с
        self.init_state[4] = 0      # начальная скорость на север, м/с
        self.init_state[5] = 0      # начальная вертикальная скорость, м/с
        
        # Углы Эйлера (в радианах): крен, тангаж, рыскание
        self.init_state[6] = 0      # крен (вращение вокруг оси X - по правому крылу)
        self.init_state[7] = 0      # тангаж (вращение вокруг оси Y - продольной)
        self.init_state[8] = 0      # рыскание (вращение вокруг оси Z - нормальной)
        
        # Угловые скорости в связанной системе
        self.init_state[9] = 0      # p (вращение вокруг оси X)
        self.init_state[10] = 0     # q (вращение вокруг оси Y)
        self.init_state[11] = 0     # r (вращение вокруг оси Z)
        
    def calculate_gravity(self, lat, h):
        """Вычисление силы тяжести с учетом высоты и широты"""
        # Формула Соммильяны
        sin_lat = np.sin(lat)
        g0 = 9.780327 * (1 + 0.0053024 * sin_lat**2 - 0.0000058 * np.sin(2*lat)**2)
        return g0 * (1 - 2*h/self.R_earth)
    
    def calculate_transport_rate(self, v_n, v_e, lat, h):
        """Вычисление переносной угловой скорости"""
        # Радиус кривизны Земли в меридиональном направлении
        R_n = self.R_earth * (1 - self.e**2) / (1 - self.e**2 * np.sin(lat)**2)**(3/2)
        
        # Радиус кривизны Земли в первом вертикале
        R_e = self.R_earth / np.sqrt(1 - self.e**2 * np.sin(lat)**2)
        
        # Компоненты переносной угловой скорости
        omega_en_n = v_e / (R_e + h)
        omega_en_e = -v_n / (R_n + h)
        omega_en_u = -v_e * np.tan(lat) / (R_e + h)
        
        return np.array([omega_en_e, omega_en_n, omega_en_u])
    
    def aircraft_dynamics(self, t, state):
        """
        Динамика летательного аппарата в системе ENU с учетом Кориолиса
        state: [x_e, x_n, x_u, v_e, v_n, v_u, phi, theta, psi, p, q, r]
        """
        # Извлечение переменных состояния
        x_e, x_n, x_u, v_e, v_n, v_u, phi, theta, psi, p, q, r = state
        
        # Текущие координаты
        lat = self.lat0 + x_n / self.R_earth
        lon = self.lon0 + x_e / (self.R_earth * np.cos(self.lat0))
        h = self.h0 + x_u
        
        # Вычисление силы тяжести
        g = self.calculate_gravity(lat, h)
        
        # Вычисление переносной угловой скорости
        omega_en = self.calculate_transport_rate(v_n, v_e, lat, h)
        
        # Угловая скорость вращения Земли в системе ENU
        omega_ie_e = self.omega_earth * np.cos(lat)
        omega_ie_u = self.omega_earth * np.sin(lat)
        omega_ie_enu = np.array([omega_ie_e, 0, omega_ie_u])
        
        # Полная угловая скорость системы ENU
        omega_en_enu = omega_en + omega_ie_enu
        
        # Управляющие воздействия
        t_period = 60
        elevator = 0.5 * np.sin(2 * np.pi * t / t_period)
        aileron = 0.8 * np.sin(2 * np.pi * t / (t_period * 1.3))
        rudder = 0.3 * np.sin(2 * np.pi * t / (t_period * 1.7))
        throttle = 0.7 + 0.2 * np.sin(2 * np.pi * t / (t_period * 2))
        
        # Аэродинамические силы и моменты
        v_total = np.sqrt(v_e**2 + v_n**2 + v_u**2)
        lift = 0.5 * 1.225 * np.exp(-h/8500) * v_total**2 * 20 * (0.2 + elevator)
        drag = 0.5 * 1.225 * np.exp(-h/8500) * v_total**2 * 20 * 0.05
        side_force = 0.5 * 1.225 * np.exp(-h/8500) * v_total**2 * 5 * rudder
        
        L = 0.5 * 1.225 * np.exp(-h/8500) * v_total**2 * 20 * 5 * aileron
        M = 0.5 * 1.225 * np.exp(-h/8500) * v_total**2 * 20 * 2 * elevator
        N = 0.5 * 1.225 * np.exp(-h/8500) * v_total**2 * 20 * 5 * rudder
        
        thrust = throttle * 20000
        
        # Производные угловых скоростей
        p_dot = (L + (self.Iyy - self.Izz) * q * r) / self.Ixx
        q_dot = (M + (self.Izz - self.Ixx) * p * r) / self.Iyy
        r_dot = (N + (self.Ixx - self.Iyy) * p * q) / self.Izz
        
        # Матрица преобразования угловых скоростей в производные углов Эйлера
        euler_matrix = np.array([
            [1, np.sin(phi)*np.tan(theta), np.cos(phi)*np.tan(theta)],
            [0, np.cos(phi), -np.sin(phi)],
            [0, np.sin(phi)/np.cos(theta), np.cos(phi)/np.cos(theta)]
        ])
        
        euler_dot = euler_matrix @ np.array([p, q, r])
        
        # Преобразование сил из связанной системы в ENU
        rotation_body_to_enu = Rotation.from_euler('zyx', [psi, theta, phi]).as_matrix()
        
        forces_body = np.array([
            side_force,      # по оси X (правое крыло)
            thrust - drag,   # по оси Y (вперед)
            -lift            # по оси Z (вверх)
        ])
        
        forces_enu = rotation_body_to_enu @ forces_body
        
        # Ускорения в системе ENU с учетом Кориолиса и центробежных сил
        # a = f/m + g - (2ω_ie + ω_en) × v - ω_ie × (ω_ie × r)
        accel_enu = forces_enu / self.mass
        
        # Гравитация
        accel_enu[2] -= g  # гравитация действует вниз
        
        # Ускорение Кориолиса
        coriolis_accel = -2 * np.cross(omega_en_enu, [v_e, v_n, v_u])
        accel_enu += coriolis_accel
        
        # Центробежное ускорение (пренебрежимо мало для большинства применений)
        # centrifugal_accel = -np.cross(omega_ie_enu, np.cross(omega_ie_enu, [R_earth, 0, 0]))
        # accel_enu += centrifugal_accel
        
        # Производные состояния
        state_dot = np.zeros(12)
        state_dot[0:3] = [v_e, v_n, v_u]  # производные положения
        state_dot[3:6] = accel_enu  # производные скорости
        state_dot[6:9] = euler_dot  # производные углов Эйлера
        state_dot[9:12] = [p_dot, q_dot, r_dot]  # производные угловых скоростей
        
        return state_dot
    
    def simulate_ins(self):
        """Моделирование движения ЛА и генерация показаний ИНС"""
        t_eval = np.arange(0, self.t_end, self.dt)
        
        sol = solve_ivp(self.aircraft_dynamics, [0, self.t_end], self.init_state, 
                        t_eval=t_eval, method='RK45', rtol=1e-6)
        
        t = sol.t
        states = sol.y.T
        
        positions = states[:, 0:3]
        velocities = states[:, 3:6]
        euler_angles = states[:, 6:9]
        angular_velocities = states[:, 9:12]
        
        # Вычисление идеальных показаний акселерометров в связанной системе
        accelerations = np.zeros_like(velocities)
        
        for i in range(len(t)):
            # Текущие координаты
            lat = self.lat0 + positions[i, 1] / self.R_earth
            h = self.h0 + positions[i, 2]
            
            # Вычисление силы тяжести
            g = self.calculate_gravity(lat, h)
            
            # Вычисление переносной угловой скорости
            omega_en = self.calculate_transport_rate(velocities[i, 1], velocities[i, 0], lat, h)
            
            # Угловая скорость вращения Земли в системе ENU
            omega_ie_e = self.omega_earth * np.cos(lat)
            omega_ie_u = self.omega_earth * np.sin(lat)
            omega_ie_enu = np.array([omega_ie_e, 0, omega_ie_u])
            
            # Полная угловая скорость системы ENU
            omega_en_enu = omega_en + omega_ie_enu
            
            # Вычисление матрицы поворота из ENU в связанную систему
            phi, theta, psi = euler_angles[i]
            rotation_enu_to_body = Rotation.from_euler('zyx', [psi, theta, phi]).as_matrix().T
            
            # Ускорение в системе ENU (без гравитации)
            if i > 0:
                dt = t[i] - t[i-1]
                accel_enu = (velocities[i] - velocities[i-1]) / dt
            else:
                accel_enu = np.zeros(3)
            
            # Добавляем гравитацию
            accel_enu[2] += g
            
            # Добавляем ускорение Кориолиса
            coriolis_accel = 2 * np.cross(omega_en_enu, velocities[i])
            accel_enu += coriolis_accel
            
            # Преобразование в связанную систему
            accel_body = rotation_enu_to_body @ accel_enu
            
            # Добавление смещений акселерометров
            accel_body_with_bias = accel_body + self.true_accel_biases
            accelerations[i] = accel_body_with_bias
        
        # Добавление уходов к гироскопам
        gyro_with_biases = angular_velocities + self.true_gyro_biases
        
        return t, positions, velocities, euler_angles, accelerations, gyro_with_biases
    
    def plot_results(self, t, positions, velocities, euler_angles, accelerations, angular_velocities):
        """Визуализация результатов моделирования"""
        fig, axes = plt.subplots(3, 2, figsize=(15, 12))
        
        # Траектория в горизонтальной плоскости
        axes[0, 0].plot(positions[:, 0], positions[:, 1])
        axes[0, 0].set_xlabel('Восток (м)')
        axes[0, 0].set_ylabel('Север (м)')
        axes[0, 0].set_title('Траектория в горизонтальной плоскости (EN)')
        axes[0, 0].grid(True)
        axes[0, 0].axis('equal')
        
        # Высота
        axes[0, 1].plot(t, positions[:, 2])
        axes[0, 1].set_xlabel('Время (с)')
        axes[0, 1].set_ylabel('Высота (м)')
        axes[0, 1].set_title('Высота (Up)')
        axes[0, 1].grid(True)
        
        # Углы Эйлера
        axes[1, 0].plot(t, np.degrees(euler_angles[:, 0]), label='Крен (φ)')
        axes[1, 0].plot(t, np.degrees(euler_angles[:, 1]), label='Тангаж (θ)')
        axes[1, 0].plot(t, np.degrees(euler_angles[:, 2]), label='Рыскание (ψ)')
        axes[1, 0].set_xlabel('Время (с)')
        axes[1, 0].set_ylabel('Углы (градусы)')
        axes[1, 0].set_title('Углы Эйлера')
        axes[1, 0].legend()
        axes[1, 0].grid(True)
        
        # Угловые скорости в связанной системе
        axes[1, 1].plot(t, np.degrees(angular_velocities[:, 0]), label='p (ось X - правое крыло)')
        axes[1, 1].plot(t, np.degrees(angular_velocities[:, 1]), label='q (ось Y - продольная)')
        axes[1, 1].plot(t, np.degrees(angular_velocities[:, 2]), label='r (ось Z - нормальная)')
        axes[1, 1].set_xlabel('Время (с)')
        axes[1, 1].set_ylabel('Угл. скорость (град/с)')
        axes[1, 1].set_title('Угловые скорости в связанной системе')
        axes[1, 1].legend()
        axes[1, 1].grid(True)
        
        # Скорости в системе ENU
        axes[2, 0].plot(t, velocities[:, 0], label='V_e (восток)')
        axes[2, 0].plot(t, velocities[:, 1], label='V_n (север)')
        axes[2, 0].plot(t, velocities[:, 2], label='V_u (вверх)')
        axes[2, 0].set_xlabel('Время (с)')
        axes[2, 0].set_ylabel('Скорость (м/с)')
        axes[2, 0].set_title('Скорости в системе ENU')
        axes[2, 0].legend()
        axes[2, 0].grid(True)
        
        # Ускорения в связанной системе
        axes[2, 1].plot(t, accelerations[:, 0], label='A_x (правое крыло)')
        axes[2, 1].plot(t, accelerations[:, 1], label='A_y (продольная)')
        axes[2, 1].plot(t, accelerations[:, 2], label='A_z (нормальная)')
        axes[2, 1].set_xlabel('Время (с)')
        axes[2, 1].set_ylabel('Ускорение (м/с²)')
        axes[2, 1].set_title('Ускорения в связанной системе (с смещениями)')
        axes[2, 1].legend()
        axes[2, 1].grid(True)
        
        plt.tight_layout()
        plt.show()

# Запуск моделирования
if __name__ == "__main__":
    simulator = INSSimulatorENU()
    t, positions, velocities, euler_angles, accelerations, angular_velocities = simulator.simulate_ins()
    
    # Вывод истинных значений уходов
    print("Истинные уходы гироскопов (рад/с):", simulator.true_gyro_biases)
    print("Истинные смещения акселерометров (м/с²):", simulator.true_accel_biases)
    
    # Визуализация результатов
    simulator.plot_results(t, positions, velocities, euler_angles, accelerations, angular_velocities)
    
    # Сохранение данных для дальнейшего использования
    data = {
        'time': t,
        'position_enu': positions,
        'velocity_enu': velocities,
        'euler_angles': euler_angles,
        'acceleration_body': accelerations,
        'angular_velocity_body': angular_velocities,
        'gyro_biases': simulator.true_gyro_biases,
        'accel_biases': simulator.true_accel_biases
    }
    
    np.save('ins_simulation_enu_coriolis.npy', data)
    print("Данные сохранены в файл 'ins_simulation_enu_coriolis.npy'")