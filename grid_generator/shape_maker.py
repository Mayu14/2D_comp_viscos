# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt

class polygon(object):
    # 初期化
    def __init__(self, num_vertex, regular = False, star = False, inner_diameter=0.5, threshold= 10.0**(-9), resolution=1000):
        self.num_vertex = num_vertex
        self.inner_diameter = inner_diameter
        self.vertex = []

        self.threshold = threshold
        self.resolution = resolution

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
        self.vertex = np.exp(1j * np.linspace(start = 0, stop = 2.0 * np.pi, num = self.num_vertex + 1))[1:]
        # self.sample_plot()

    def set_star_polygon(self):
        outer_vertex = np.exp(1j * np.linspace(start=0, stop=2.0 * np.pi, num=self.num_vertex + 1)).reshape(1, -1)
        phase_differnce = 2.0 * np.pi / (2.0 * self.num_vertex)
        inner_vertex = self.inner_diameter * np.exp(1j * np.linspace(start=0 - phase_differnce, stop=2.0 * np.pi - phase_differnce, num=self.num_vertex + 1)).reshape(1, -1)
        self.vertex = np.concatenate((inner_vertex, outer_vertex)).T.reshape(-1)[2:]
        # self.sample_plot()

    def make_edge(self):
        # 辺のy両端の複素座標から傾きと切片求める
        def get_a_and_b(z1, z2):
            if abs(np.real(z2) - np.real(z1)) < self.threshold:
                a = np.inf  # 直線がy軸と並行
            else:
                a = (np.imag(z2) - np.imag(z1)) / (np.real(z2) - np.real(z1))

            if a != np.inf:
                b = np.imag(z1) - a * np.real(z1)
            else:
                b = 0
            return a, b

        self.num_edge = self.vertex.shape[0]    # 植木算
        self.edge_param = np.zeros((self.num_edge, 2))    # 傾き + 切片

        self.edge_a = np.zeros(self.num_edge)
        self.edge_b = np.zeros(self.num_edge)

        for i in range(self.num_edge - 1):
            self.edge_a[i], self.edge_b[i] = get_a_and_b(self.vertex[i], self.vertex[i + 1])

        i = self.num_edge - 1
        self.edge_a[i], self.edge_b[i] = get_a_and_b(self.vertex[i], self.vertex[0])

        # 上側と下側とを分離
        self.edge_u = np.zeros(self.num_edge)
        trail = np.argmax(np.real(self.vertex))
        lead = np.argmin(np.real(self.vertex))

        # ここ，リストに変更する
        """
        # 辺番号を座標が小さい順に格納する
        upper -> lower -? upperのときは
        upper:逆順に格納，2回目のupperは
        lower:順番に格納
        """
        if lead < trail:    # upper -> lower -> upperの順
            for i in range(lead):
                self.edge_u[i] = 1
            for i in range(trail, self.num_edge):
                self.edge_u[i] = 1

        else:   # lower -> upper -> lowerの順
            for i in range(trail, lead):
                self.edge_u[i] = 1

    def get_equidistant_points(self):
        y_u = np.zeros(self.resolution)
        y_l = np.zeros(self.resolution)



    def sample_plot(self):
        plt.plot(np.real(self.vertex), np.imag(self.vertex))
        plt.show()


def main():
    pol = polygon(6, star = True, inner_diameter=0.7)


if __name__ == '__main__':
    main()
