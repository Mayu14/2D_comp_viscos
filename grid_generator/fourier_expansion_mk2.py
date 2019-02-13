# coding: utf-8
# 配列と解像度を受け取ってFourier展開の係数，および必要に応じて誤差を返す
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from naca_4digit_test import Naca_4_digit
from numba import jit

class fourier_expansion(object):
    def __init__(self, x_u, y_u, x_l, y_l, n=128):
        # mirroring_array = lambda l: 2.0 * np.max(l) - l
        mirroring_array = lambda l, len2: len2 * np.max(l) - l
        convert_array = lambda u, l: np.concatenate((u, l[-1::-1]))

        z = np.concatenate([x_u[::-1], x_l]) + 1j * np.concatenate([y_u[::-1], y_l])
        x_argmin = np.argmin(np.real(z))
        x_argmax = np.argmax(np.real(z))
        if x_argmin < x_argmax:
            z_l = z[x_argmin:x_argmax]
            z_u = np.concatenate([z[x_argmax:], z[:x_argmin]])[::-1]
        else:    # x_argmax < x_argmin
            z_l = np.concatenate([z[x_argmin:], z[:x_argmax]])
            z_u = z[x_argmax:x_argmin][::-1]

        x_min = np.real(z[x_argmin])
        x_max = np.real(z[x_argmax])
        self.L = 2.0 * (x_max - x_min)

        # self.x = convert_array(x_u, mirroring_array(x_l))
        self.x = convert_array(np.real(z_u), mirroring_array(np.real(z_l), self.L))

        # self.y = convert_array(y_u - 0.5, y_l - 0.5)
        self.y = convert_array(np.imag(z_u) - 0.5, np.imag(z_l) - 0.5)
        self.resolution = x_u.shape[0]

        self.max_n = n
        self.length = 2.0
        self.fourier_expansion()
        
    def test_plot_raw_data(self):
        plt.xlim([0, 2])
        plt.ylim([-0.25, 0.25])
        plt.plot(self.x, self.y)
        plt.show()
        
    def test_plot_decryption_data(self):
        plot_resolution = 5000
        decryption = self.decryption(plot_resolution)
        x_star = np.linspace(start = 0, stop = 2, num = plot_resolution * 2)
        # plt.xlim([0, 2])
        # plt.ylim([-0.25, 0.25])
        plt.plot(x_star, decryption)
        plt.show()

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.set_title("Fourier Sr. Sampling")
        # ax.set_xlim(-0.01, 1.01)
        # ax.set_ylim(-0.2, 0.2)
        decryption_y_u = decryption[:plot_resolution]
        decryption_y_l = decryption[plot_resolution:][-1::-1]
        x_star = np.linspace(start = 0, stop = 1, num = plot_resolution)
        ax.plot(x_star, decryption_y_u, ".", color="darkred")
        ax.plot(x_star, decryption_y_l, ".", color="darkred")
        plt.show()

    @jit
    def fourier_expansion(self):
        L = self.length
        n_max = self.max_n + 1
        self.bn = np.zeros(n_max)
        for n in range(n_max):
            y_star = self.y * np.sin((n + 1) * np.pi * self.x / L)
            
            self.bn[n] = integrate.trapz(y=y_star, x=self.x, axis = -1)
    
    def decryption(self, plot_resolution):
        L = self.length
        n_max = self.max_n + 1
        x_star = np.linspace(start = 0, stop = 2, num = plot_resolution*2)
        decryption = np.zeros(plot_resolution * 2)
        for n in range(n_max):
            decryption += self.bn[n] * np.sin((n + 1) * np.pi * x_star / L)
        return decryption
    
def main():
    naca = Naca_4_digit(int_4 = "2612", attack_angle_deg = 0, resolution = 10000)
    fex = fourier_expansion(naca.x_u, naca.y_u, naca.x_l, naca.y_l, n=200)
    # fex.test_plot_raw_data()
    fex.test_plot_decryption_data()
    # resolution は100あれば十分っぽい
    error = np.sum(np.sqrt((fex.decryption(fex.resolution) - fex.y)**2)) / fex.resolution
    print(error)
    print(fex.bn)
    
if __name__ == '__main__':
    main()
