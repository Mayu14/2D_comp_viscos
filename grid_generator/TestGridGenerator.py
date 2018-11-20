# coding: utf-8
import numpy as np

def main(xi_max, eta_max):
    def output_vtk_tri(fname, path):
        fname = path + "tri_" + fname + ".vtk"
        with open(fname, 'w') as f:
            point_number = str((xi_max) * eta_max)
            cell_number = str(2 * (xi_max - 1) * (eta_max - 1))
            cell_vertex_number = str(2 * 4 * (xi_max - 1) * (eta_max - 1))
            pid = lambda i, j: i + xi_max * j

            # 四角形セルの対角線長さを計算し，左下→右上の対角線(diag1)が，左上→右下の対角線(diag2)より短いときTrueを返す
            def set_relatively_not_flat_tri(i, ip1, j):
                diag1 = (grid_x[ip1, j + 1] - grid_x[i, j]) ** 2 + (grid_y[ip1, j + 1] - grid_y[i, j]) ** 2
                diag2 = (grid_x[i, j + 1] - grid_x[ip1, j]) ** 2 + (grid_y[i, j + 1] - grid_y[ip1, j]) ** 2

                if diag1 < diag2:
                    return True
                else:
                    return False

            def cell_structure(i, j):
                # 短い対角線を新しい辺として四角形を三角形に分割する
                if i != xi_max - 1:
                    ip1 = i + 1
                else:
                    ip1 = 0

                flag = set_relatively_not_flat_tri(i, ip1, j)
                # flag == True:左下と右上を結ぶ対角線を新たな辺として，左上の三角形と右下の三角形に分割
                if flag:
                    tri1 = "3 " + str(pid(i, j)) + " "  + str(pid(ip1, j + 1)) + " " + str(pid(i, j + 1))
                    tri2 = "3 " + str(pid(i, j)) + " " + str(pid(ip1, j)) + " " + str(pid(ip1, j + 1))
                else:
                    tri1 = "3 " + str(pid(i, j)) + " " + str(pid(ip1, j)) + " " + str(pid(i, j + 1))
                    tri2 = "3 " + str(pid(ip1, j)) + " " + str(pid(ip1, j + 1)) + " " + str(pid(i, j + 1))
                return tri1 + "\n" + tri2

            # header
            f.write("# vtk DataFile Version 3.0\n")
            f.write("Unstructured Grid tri example\n")
            f.write("ASCII\nDATASET UNSTRUCTURED_GRID\n")
            f.write("POINTS " + point_number + " double\n")
            # point coordinates
            for j in range(eta_max):
                for i in range(xi_max):
                    f.write(str(grid_x[i, j]) + " " + str(grid_y[i, j]) + " 0.0\n")

            # cell structure
            f.write("CELLS " + cell_number + " " + cell_vertex_number + "\n")
            for j in range(eta_max - 1):
                for i in range(xi_max - 1):
                    f.write(cell_structure(i, j) + "\n")

            # cell types
            f.write("CELL_TYPES " + cell_number + "\n")
            for j in range(eta_max - 1):
                for i in range(xi_max - 1):
                    f.write("5\n5\n")

    grid_x = np.zeros((xi_max, eta_max))
    grid_y = np.zeros((xi_max, eta_max))
    
    dx = 1.0 / (xi_max - 1)
    dy = 1.0 / (eta_max - 1)
    
    for i in range(xi_max):
        for j in range(eta_max):
            grid_x[i, j] = i * dx
            grid_y[i, j] = j * dy
    
    output_vtk_tri(fname = "SquareGrid", path = "")
    
    
if __name__ == '__main__':
    Nx = 100
    Ny = 100
    main(Nx, Ny)
    