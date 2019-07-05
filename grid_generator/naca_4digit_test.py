# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt


class Naca_4_digit(object):
    def __init__(self, int_4, attack_angle_deg, resolution, auto_normalize = True, monotonization = False,
                 quasi_equidistant = True, length_adjust = False, closed_te = True, position="[0,1]",
                 from5digit = False):
        """
        
        :param int_4: NACA翼の4桁の型番
        :param attack_angle_deg: 迎角
        :param resolution: 返される点列の解像度
        :param auto_normalize: 長さの正規化の有無(NACA6192などは[0,1]に収まらないため)
        :param monotonization: 上側点・下側点がそれぞれx方向に単調増加になるような並べかえを実施(非推奨)
        :param quasi_equidistant: 近似的にx方向等間隔な点列を得る
        :param length_adjust:
        :param closed_te: NACA翼の後縁を閉じるオプション
        :param position:    [0,1]:[0,1]x[0,1]区間に配置，start_0:翼後縁にある[x,y]の始点を[0,0]に固定
        :param from5digit: NACA5桁翼を入力する場合の分岐
        """
        if from5digit == False:
            self.m = float(int_4[0]) / 100  # maximum camber
            self.p = float(int_4[1]) / 10  # position of the maximum camber
            self.t = float(int_4[2:4]) / 100  # maximum thickness
            self.load_setting(attack_angle_deg, resolution, quasi_equidistant, length_adjust)
            self.__y_c()
            self.__dyc_dx()

        self.closed = closed_te
        self.position = position
        self.__y_t()
        self.theta = np.arctan(self.dyc_dx)
        self.get_surface()
        self.get_center_and_length()
        
        if auto_normalize:
            self.normalize()
            
        if self.position == "start_0":
            self.set_zero_center()
        
        self.monotonized = monotonization
        if monotonization:
            self.monotonize()
        
        if quasi_equidistant == True:
            self.get_quasi_equidistant_line()
    
    def set_zero_center(self):
        self.x_u -= 1.0
        self.x_l -= 1.0
        self.y_u -= 0.5
        self.y_l -= 0.5
    
    def load_setting(self, attack_angle_deg, resolution, quasi_equidistant = True, length_adjust = False):
        self.use_quasi_equidistant = quasi_equidistant
        self.reshape = length_adjust
        if quasi_equidistant == True:
            self.resolution = 100 * resolution
        else:
            self.resolution = resolution
        self.new_resolution = resolution
        self.attack_angle = attack_angle_deg
        self.x = np.linspace(start = 0, stop = 1, num = self.resolution)
    
    def get_max_and_min(self):
        self.x_max = max(np.max(self.x_u), np.max(self.x_l))
        self.x_min = min(np.min(self.x_u), np.min(self.x_l))
        self.y_max = max(np.max(self.y_u), np.max(self.y_l))
        self.y_min = min(np.min(self.y_u), np.min(self.y_l))
    
    def get_center_and_length(self):
        self.get_max_and_min()
        self.center_x = 0.5 * (self.x_max + self.x_min)
        self.center_y = 0.5 * (self.y_max + self.y_min)
        self.length_x = (self.x_max - self.x_min)
        self.length_y = (self.y_max - self.y_min)
    
    def normalize(self):
        self.get_center_and_length()
        self.x_u = (self.x_u - self.center_x) / self.length_x
        self.x_l = (self.x_l - self.center_x) / self.length_x
        self.y_u = (self.y_u - self.center_y) / self.length_x + self.center_y
        self.y_l = (self.y_l - self.center_y) / self.length_x + self.center_y
        self.get_center_and_length()
        self.x_u -= self.x_min
        self.x_l -= self.x_min
    
    def monotonize(self):
        # x方向の座標変換を単調にするため，端部でx_u, x_lの一部の点を互いにやり取りする
        x_u = self.x_u.tolist()
        x_l = self.x_l.tolist()
        y_u = self.y_u.tolist()
        y_l = self.y_l.tolist()
        size_u = len(x_u)
        passing = np.zeros(size_u)
        for i in range(size_u - 1):
            if self.x_u[i] > self.x_u[i + 1]:
                passing[i] = 1
        
        k = 0
        for i in range(size_u - 1):
            if passing[i] == 1:
                x_l.insert(0, x_u.pop(k))
                y_l.insert(0, y_u.pop(k))
                k -= 1
            k += 1
        
        size_l = len(x_l)
        passing = np.zeros(size_l)
        for j in range(int(size_l / 2), size_l - 1):
            if x_l[j] > x_l[j + 1]:
                passing[j] = 1
        
        k = size_l - 1
        for j in range(size_l - 1, -1, -1):
            if passing[j] == 1:
                x_u.append(x_l.pop(k))
                y_u.append(y_l.pop(k))
                k -= 1
            k -= 1
        
        self.x_u = np.array(x_u)
        self.x_l = np.array(x_l)
        self.y_u = np.array(y_u)
        self.y_l = np.array(y_l)
    
    def __y_c(self):
        x_lt_p = lambda m, p, x: m / (p ** 2) * (2.0 * p * x - x ** 2)
        x_ge_p = lambda m, p, x: m / ((1 - p) ** 2) * ((1.0 - 2.0 * p) + 2.0 * p * x - x ** 2)
        
        m = self.m
        p = self.p
        x = self.x
        if ((p != 0) and (p != 1)):
            self.y_c = np.where(x < p, x_lt_p(m, p, x), x_ge_p(m, p, x))
        elif (p == 0):
            self.y_c = m * (1.0 - x ** 2)
        elif (p == 1):
            self.y_c = m * (2.0 * x - x ** 2)
    
    def __y_t(self):
        t = self.t
        x = self.x
        a0 = 0.2969
        a1 = -0.1260
        a2 = -0.3516
        a3 = 0.2843
        if self.closed:
            a4 = -0.1036
        else:
            a4 = -0.1015
        self.y_t = t / 0.2 * (a0 * np.sqrt(x) + a1 * x + a2 * x ** 2 + a3 * x ** 3 + a4 * x ** 4)
    
    def __dyc_dx(self):
        x_lt_p = lambda m, p, x: 2.0 * m / (p ** 2) * (p - x)
        x_ge_p = lambda m, p, x: 2.0 * m / ((1.0 - p) ** 2) * (p - x)
        
        m = self.m
        p = self.p
        x = self.x
        if ((p != 0) and (p != 1)):
            self.dyc_dx = np.where(x < p, x_lt_p(m, p, x), x_ge_p(m, p, x))
        elif (p == 0):
            self.dyc_dx = - 2.0 * m * x
        elif (p == 1):
            self.dyc_dx = 2.0 * m * (1.0 - x)
    
    def get_surface(self):
        # original NACA-4digit wings
        # upper
        vec_l = np.full((3, self.resolution), 1.0)
        vec_u = np.full((3, self.resolution), 1.0)
        
        vec_u[0] = self.x - self.y_t * np.sin(self.theta) - 0.5
        vec_u[1] = self.y_c + self.y_t * np.cos(self.theta)
        # lower
        vec_l[0] = self.x + self.y_t * np.sin(self.theta) - 0.5
        vec_l[1] = self.y_c - self.y_t * np.cos(self.theta)
        
        attack_angle = self.attack_angle / 180 * (np.pi)
        rotMat = np.array(
            [[np.cos(attack_angle), np.sin(attack_angle), 0], [- np.sin(attack_angle), np.cos(attack_angle), 0],
             [0, 0, 1]])
        
        rot_l = rotMat.dot(vec_l)
        rot_u = rotMat.dot(vec_u)
        
        if self.reshape == True:
            x_min = min(np.min(rot_l[0]), np.min(rot_u[0]))
            x_max = max(np.max(rot_l[0]), np.max(rot_u[0]))
            rate = 1.0 / (x_max - x_min)
            
            if rate != 1.0:
                expMat = np.array([[rate, 0, 0], [0, rate, 0], [0, 0, 1]])
                rot_l = expMat.dot(rot_l)
                rot_u = expMat.dot(rot_u)
        
        center = 0.5

        self.x_l = rot_l[0] + center
        self.y_l = rot_l[1] + center
        self.x_u = rot_u[0] + center
        self.y_u = rot_u[1] + center
    
    def plot(self):
        plt.xlim([0, 1])
        plt.ylim([0, 1])
        plt.plot(self.x_u, self.y_u)
        plt.plot(self.x_l, self.y_l)
        plt.show()
    
    def get_quasi_equidistant_line(self):
        new_resolution = self.new_resolution
        x_min = min(np.min(self.x_u), np.min(self.x_l))
        x_max = max(np.max(self.x_u), np.max(self.x_l))
        
        if self.reshape == False:
            self.equidistant_x = np.linspace(start = 0, stop = 1, num = new_resolution)
        else:
            self.equidistant_x = np.linspace(start = x_min, stop = x_max, num = new_resolution)
        self.equidistant_y_l = np.zeros(new_resolution)
        self.equidistant_y_u = np.zeros(new_resolution)
        for index in range(new_resolution):
            if ((x_min <= self.equidistant_x[index]) and (x_max >= self.equidistant_x[index])):
                self.equidistant_y_l[index] = self.y_l[np.argmin(np.abs(self.x_l - self.equidistant_x[index]))]
                self.equidistant_y_u[index] = self.y_u[np.argmin(np.abs(self.x_u - self.equidistant_x[index]))]
            else:
                self.equidistant_y_l[index] = -1.0  # 外れ値
                self.equidistant_y_u[index] = -1.0
    
    def plot_quasi_equidistant_shape(self):
        plt.xlim([0, 1])
        plt.ylim([0, 1])
        plt.plot(self.equidistant_x, self.equidistant_y_u, "o")
        plt.plot(self.equidistant_x, self.equidistant_y_l, "o")
        plt.show()
    
    def transform2complex(self):
        z_u_reverse = (self.x_u + 1j * self.y_u)[::-1]
        z_l = self.x_l + 1j * self.y_l
        if self.use_quasi_equidistant == True:
            return np.concatenate([z_u_reverse[::100], z_l[::100], z_u_reverse[0].reshape(-1)])
        else:
            if z_u_reverse[self.resolution - 1] == z_l[0]:
                return np.concatenate([z_u_reverse, z_l[1:], z_u_reverse[0].reshape(-1)])
            else:
                return np.concatenate([z_u_reverse, z_l, z_u_reverse[0].reshape(-1)])
    
    def one_stroke_ccw(self):
        if self.monotonized:
            print("one_stroke cannot use with monotonized")
            exit()
    
        self.ccw_x = np.concatenate([self.x_u[::-1], self.x_l[1:]])
        self.ccw_y = np.concatenate([self.y_u[::-1], self.y_l[1:]])

    def one_stroke_cw(self):
        self.one_stroke_ccw()
        self.cw_x = self.ccw_x[::-1]
        self.cw_y = self.ccw_y[::-1]


class Naca_5_digit(Naca_4_digit):
    def __init__(self, int_5, attack_angle_deg, resolution, auto_normalize = True, monotonization = False,
                 quasi_equidistant = True, length_adjust = False, closed_te = True, position="[0,1]",
                 from5digit = True):
        
        """

        :param int_5:   1st: 1~6, 2nd and 3rd [10,20,30,40,50,21,31,41,51], 4th and 5th [01~99]
        :param attack_angle_deg: recommended 0
        :param resolution:
        :param quasi_equidistant:
        :param length_adjust:
        :param from5digit:
        :param closed_te:
        """
        self.cl = float(int_5[0]) * (3.0 / 2.0) / 10  # designed lift_coefficient
        self.p = float(int_5[1]) / 2.0 / 10  # position of the maximum camber
        self.ref = int_5[2]  # enable / disable reflect
        self.t = float(int_5[3:5]) / 100.0  # maximum thickness
        
        self.camberline_plofile = int(int_5[0:3])
        self.camberline_plofile_table()
        self.load_setting(attack_angle_deg, resolution, quasi_equidistant, length_adjust)
        self.__y_c()
        self.__dyc_dx()
        self.closed = closed_te
        super(Naca_5_digit, self).__init__(int_5, attack_angle_deg, resolution, auto_normalize = auto_normalize,
                                           monotonization = monotonization, quasi_equidistant = quasi_equidistant,
                                           length_adjust = length_adjust, closed_te = closed_te,
                                           position = position, from5digit = True)
    
    def __y_c(self):
        x_lt_m_nr = lambda m, k1, x: k1 / 6.0 * (x ** 3 - 3.0 * m * x ** 2 + m ** 2 * (3.0 - m) * x)
        x_gt_m_nr = lambda m, k1, x: k1 / 6.0 * m ** 3 * (1.0 - x)
        
        x_lt_m_rf = lambda m, k1, k2_k1, x: k1 / 6.0 * ((x - m) ** 3 - k2_k1 * (1.0 - m) ** 3 * x - m ** 3 * x + m ** 3)
        x_gt_m_rf = lambda m, k1, k2_k1, x: k1 / 6.0 * (
                k2_k1 * (x - m) ** 3 - k2_k1 * (1.0 - m) ** 3 * x - m ** 3 * x + m ** 3)
        
        m = self.m
        k1 = self.k1
        x = self.x
        if int(self.ref) == 0:  # not reflected
            self.y_c = np.where(x < m, x_lt_m_nr(m, k1, x), x_gt_m_nr(m, k1, x))
        else:
            k2_k1 = self.k2byk1
            self.y_c = np.where(x < m, x_lt_m_rf(m, k1, k2_k1, x), x_gt_m_rf(m, k1, k2_k1, x))
    
    def __dyc_dx(self):
        x_lt_m_nr = lambda m, k1, x: k1 / 6.0 * (3.0 * x ** 2 - 6.0 * m * x + m ** 2 * (3.0 - m))
        x_gt_m_nr = lambda m, k1, x: - k1 / 6.0 * m ** 3
        
        x_lt_m_rf = lambda m, k1, k2_k1, x: k1 / 6.0 * (3.0 * (x - m) ** 2 - k2_k1 * (1.0 - m) ** 3 - m ** 3)
        x_gt_m_rf = lambda m, k1, k2_k1, x: k1 / 6.0 * (3.0 * k2_k1 * (x - m) ** 2 - k2_k1 * (1.0 - m) ** 3 - m ** 3)
        
        m = self.m
        k1 = self.k1
        x = self.x
        if int(self.ref) == 0:  # not reflected
            self.dyc_dx = np.where(x < m, x_lt_m_nr(m, k1, x), x_gt_m_nr(m, k1, x))
        
        else:
            k2_k1 = self.k2byk1
            self.dyc_dx = np.where(x < m, x_lt_m_rf(m, k1, k2_k1, x), x_gt_m_rf(m, k1, k2_k1, x))
    
    def camberline_plofile_table(self, fromfunc = True):
        if fromfunc:
            from scipy.optimize import newton
            xf = self.p
            mfunc = lambda m: m * (1 - np.sqrt(m / 3.0)) - xf
            self.m = m = newton(mfunc, 0.058)
            
            Qfunc = lambda m: (3 * m - 7 * m ** 2 + 8 * m ** 3 - 4 * m ** 4) / np.sqrt(m * (1 - m)) - 1.5 * (
                    1 - 2 * m) * (np.pi / 2.0 - np.arcsin(1 - 2 * m))
            self.k1 = 6.0 * self.cl / Qfunc(self.m)
            k2byk1 = lambda m, xf: (3 * (m - xf) ** 2 - m ** 3) / (1 - m) ** 3
            self.k2byk1 = k2byk1(self.m, xf)
        else:
            if self.camberline_plofile == 210:
                self.m = 0.058
                self.k1 = 361.4
            elif self.camberline_plofile == 220:
                self.m = 0.126
                self.k1 = 51.64
            elif self.camberline_plofile == 230:
                self.m = 0.2025
                self.k1 = 15.957
            elif self.camberline_plofile == 240:
                self.m = 0.29
                self.k1 = 6.643
            elif self.camberline_plofile == 250:
                self.m = 0.391
                self.k1 = 3.230
            
            elif self.camberline_plofile == 221:
                self.m = 0.130
                self.k1 = 51.990
                self.k2byk1 = 0.000764
            elif self.camberline_plofile == 231:
                self.m = 0.217
                self.k1 = 15.793
                self.k2byk1 = 0.00677
            elif self.camberline_plofile == 241:
                self.m = 0.318
                self.k1 = 6.520
                self.k2byk1 = 0.0303
            elif self.camberline_plofile == 251:
                self.m = 0.441
                self.k1 = 3.191
                self.k2byk1 = 0.1355
            else:
                print("this type wing is not defined")
                exit()


def main():
    deg = 0.0
    naca = Naca_4_digit(int_4 = "0012", attack_angle_deg = deg, resolution = 100, quasi_equidistant = True,
                        length_adjust = True)
    # naca.plot()
    # naca.plot_quasi_equidistant_shape()
    naca = Naca_5_digit(int_5 = "21001", attack_angle_deg = deg, resolution = 100, quasi_equidistant = True,
                        length_adjust = True, position = "start_0")

    #naca.plot()
    plt.plot(naca.x_u, naca.y_u)
    plt.plot(naca.x_l, naca.y_l)
    plt.show()
    #naca.plot_quasi_equidistant_shape()


if __name__ == '__main__':
    main()
