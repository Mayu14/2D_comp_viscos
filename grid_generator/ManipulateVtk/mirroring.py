# -- coding: utf-8 --
import numpy as np

class tri_grid(object):
    def __init__(self, fname, center, line_vector, threshold = 10.0**(-9), easy_mode=0):
        self.fname = fname
        self.new_name = "mirror_" + fname

        self.center = center    # 反転の中心点
        self.line_vector = line_vector  # 反転の基準方向
        self.threshold = threshold

        self.header = ""
        self.read_data()

        if easy_mode == 0  :    # 直線x=0にて反転処理
            self.mirror_x0()

        self.output_vtk()

    def read_data(self):
        with open(self.fname, "r") as f:
            phase = 0
            p_count = 0
            c_count = 0
            for line in f:
                if ((phase == 2) and (c_count < self.c_num)):
                    self.c_strct[c_count] = np.array(line.split()[1:], dtype=int)
                    c_count += 1

                if ((phase == 1) and (p_count == self.p_num)):
                    self.c_num = int(line.split()[1])
                    self.c_strct = np.zeros((self.c_num, 3), dtype=int)
                    phase = 2

                if ((phase == 1) and (p_count < self.p_num)):
                    self.p_coord[p_count] = np.array(line.split(), dtype=float)
                    p_count += 1

                if (phase == 0):
                    if line.find("POINTS") == 0:
                        self.p_num = int(line.split()[1])
                        self.p_coord = np.zeros((self.p_num, 3), dtype=float)
                        phase = 1
                    else:
                        self.header += line

    def mirror_x0(self):
        # 新旧番号対応表作成
        self.number_map1 = np.arange(self.p_num)    # 旧番号から新番号を検索する
        self.add_p_num = 0
        # 除外要素の特定
        for i in range(self.p_num):
            if abs(self.p_coord[i, 0]) < self.threshold:
                self.number_map1[i] = i
            else:
                self.number_map1[i] = self.p_num + self.add_p_num
                self.add_p_num += 1

        # 反転後座標の作成
        self.add_p_coord = np.zeros((self.add_p_num, 3), dtype=float)
        p_count = 0
        self.number_map2 = np.arange(self.p_num + self.add_p_num)   # 新番号から旧番号を検索する
        for i in range(self.p_num):
            if (self.number_map1[i] >= self.p_num):
                self.number_map2[self.p_num + p_count] = i
                self.add_p_coord[p_count] = np.array([-self.p_coord[i, 0], self.p_coord[i, 1], self.p_coord[i, 2]])
                p_count += 1

        self.new_p_coord = np.concatenate([self.p_coord, self.add_p_coord])
        self.new_p_num = self.new_p_coord.shape[0]

        # 反転後セル構造の作成
        self.add_c_strct = np.zeros((self.c_num, 3), dtype=int)
        for i in range(self.c_num):
            self.add_c_strct[i] = np.array([self.number_map1[self.c_strct[i,0]], self.number_map1[self.c_strct[i,1]], self.number_map1[self.c_strct[i,2]]])

        self.new_c_strct = np.concatenate([self.c_strct, self.add_c_strct])
        self.new_c_num = self.new_c_strct.shape[0]

    def output_vtk(self):
        with open(self.new_name, "w") as mf:
            mf.write(self.header)

            mf.write("POINTS " + str(self.new_p_num) + " double\n")
            [mf.write(str(self.new_p_coord[i, 0]) + " " + str(self.new_p_coord[i, 1]) + " " + str(self.new_p_coord[i, 2]) + " "  + "\n") for i in range(self.new_p_num)]

            mf.write("CELLS " + str(self.new_c_num) + " " + str(4*self.new_c_num) + "\n")
            [mf.write("3 " + str(self.new_c_strct[i, 0]) + " " + str(self.new_c_strct[i, 1]) + " " + str(self.new_c_strct[i, 2]) + " " + "\n") for i in range(self.new_c_num)]

            mf.write("CELL_TYPES " + str(self.new_c_num) + "\n")
            [mf.write("5\n") for i in range(self.new_c_num)]




def set_path():
    return input("please input vtkname")



def main():
    fname = "Square_Half_fine.vtk"

    center = [0.0, 0.0]
    line_vector = [0.0, 1.0]

    threshold = 10.0**(-9)

    tri = tri_grid(fname, center, line_vector, easy_mode=0)





if __name__ == '__main__':
    main()