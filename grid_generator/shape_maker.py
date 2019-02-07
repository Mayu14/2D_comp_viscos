# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt


class polygon(object):
    # 初期化
    def __init__(self, num_vertex, regular = False, star = False):
        self.num_vertex = num_vertex
        self.vertex = []
        
        if regular:
            self.set_regular_polygon(num_vertex)
    
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
    
    def set_regular_polygon(self, num_vertex):
        circle = np.exp(1j * np.linspace(start = 0, stop = 2.0 * np.pi, num = num_vertex + 1))


def main():
    pol = polygon(8, regular = True)


if __name__ == '__main__':
    main()
