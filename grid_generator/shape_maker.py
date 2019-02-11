# coding: utf-8
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

class polygon(object):
    # 初期化
    def __init__(self, num_vertex, regular = False, star = False, inner_diameter=0.5, threshold= 10.0**(-9), resolution=1000):
        self.num_vertex = num_vertex
        self.inner_diameter = inner_diameter
        self.vertex = []

        self.threshold = threshold
        self.resolution = resolution

        self.radius = 0.5
        self.center = 0.5 * (1 + 1j)

        if regular:
            self.header = "polygon_regular_"
            self.set_regular_polygon()

        if star:
            self.header = "polygon_star_"
            self.set_star_polygon()

        self.make_edge()

    # coord_xyzは座標[x, y, z]のリスト
    def define_vertex(self, coord_xyz):
        self.vertex.append(np.array(coord_xyz))
        if len(self.vertex) == self.num_vertex:
            self.vertex.append(self.vertex[0])
            self.get_perimeter()
    
    # 周の長さを求める
    def get_perimeter(self):
        self.perimeter = 0.0
        for i in range(self.num_vertex):
            self.perimeter += np.sqrt(np.sum((self.vertex[i + 1] - self.vertex[i]) ** 2))

    def set_regular_polygon(self):
        self.vertex = self.radius * np.exp(1j * np.linspace(start = 0, stop = 2.0 * np.pi, num = self.num_vertex + 1))[1:] + self.center
        # self.sample_plot()

    def set_star_polygon(self):
        outer_vertex = self.radius * np.exp(1j * np.linspace(start=0, stop=2.0 * np.pi, num=self.num_vertex + 1)).reshape(1, -1) + self.center
        phase_differnce = 2.0 * np.pi / (2.0 * self.num_vertex)
        inner_vertex = self.radius * self.inner_diameter * np.exp(1j * np.linspace(start=0 - phase_differnce, stop=2.0 * np.pi - phase_differnce, num=self.num_vertex + 1)).reshape(1, -1) + self.center
        self.vertex = np.concatenate((inner_vertex, outer_vertex)).T.reshape(-1)[2:]
        self.num_vertex = self.vertex.shape[0]
        # self.sample_plot()

    def make_edge(self):
        trail = np.argmax(np.real(self.vertex)) # 後方
        lead = np.argmin(np.real(self.vertex))  # 前方

        # 上側と下側とを分離(ここまでの処理で反時計回りに格納されていることに注意)
        vertex_u = []
        vertex_l = []

        # 点番号を座標が小さい順に格納する
        tmp_vertex = []
        if lead < trail:    # upper -> lower -> upperの順
            for i in range(lead + 1):   # leadまで含める
                vertex_u.insert(0, i)
            for i in range(lead, trail+1):  # trailまで含める
                vertex_l.append(i)
            for i in range(trail, self.num_vertex):
                tmp_vertex.insert(0, i)
            vertex_u.extend(tmp_vertex)

        else:   # lower -> upper -> lowerの順
            for i in range(trail+1):    # trailまで含める
                tmp_vertex.append(i)
            for i in range(trail, lead+1):
                vertex_u.insert(0, i)
            for i in range(lead, self.num_vertex):
                vertex_l.append(i)
            vertex_l.extend(tmp_vertex)

        self.z_l = np.zeros(len(vertex_l), dtype=complex)
        self.z_u = np.zeros(len(vertex_l), dtype=complex)
        i = 0
        for p_num_u in vertex_u:
            self.z_u[i] = self.vertex[p_num_u]
            i += 1

        i = 0
        for p_num_l in vertex_l:
            self.z_l[i] = self.vertex[p_num_l]
            i += 1

        self.get_equidistant_points()

    def get_equidistant_points(self):
        x_eq = np.linspace(start=0, stop=1, num=self.resolution)
        f_u = interpolate.interp1d(np.real(self.z_u), np.imag(self.z_u), kind="linear")
        f_l = interpolate.interp1d(np.real(self.z_l), np.imag(self.z_l), kind="linear")

        plt.plot(x_eq, f_u(x_eq), "x")
        plt.plot(x_eq, f_l(x_eq), "o")
        self.sample_plot()




    def sample_plot(self):
        plt.plot(np.real(self.vertex), np.imag(self.vertex))
        plt.show()


def main():
    pol = polygon(6, star = True, inner_diameter=0.7)


if __name__ == '__main__':
    main()
