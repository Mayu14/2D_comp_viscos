# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt

class polygon(object):
    # 初期化
    def __init__(self, num_vertex, regular = False, star = False, inner_diameter=0.5):
        self.num_vertex = num_vertex
        self.inner_diameter = inner_diameter
        self.vertex = []
        
        if regular:
            self.header = "polygon_regular_"
            self.set_regular_polygon()

        if star:
            self.header = "polygon_star_"
            self.set_star_polygon()
    
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
        self.vertex = np.exp(1j * np.linspace(start = 0, stop = 2.0 * np.pi, num = self.num_vertex + 1))
        self.sample_plot()

    def set_star_polygon(self):
        outer_vertex = np.exp(1j * np.linspace(start=0, stop=2.0 * np.pi, num=self.num_vertex + 1)).reshape(1, -1)
        phase_differnce = 2.0 * np.pi / (2.0 * self.num_vertex)
        inner_vertex = self.inner_diameter * np.exp(1j * np.linspace(start=0 - phase_differnce, stop=2.0 * np.pi - phase_differnce, num=self.num_vertex + 1)).reshape(1, -1)
        self.vertex = np.concatenate((inner_vertex, outer_vertex)).T.reshape(-1)
        self.sample_plot()

    def sample_plot(self):
        plt.plot(np.real(self.vertex), np.imag(self.vertex))
        plt.show()


def main():
    pol = polygon(6, star = True, inner_diameter=0.7)


if __name__ == '__main__':
    main()
