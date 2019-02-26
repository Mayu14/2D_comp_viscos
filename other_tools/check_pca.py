# coding: utf-8
import numpy as np
from read_training_data_viscos import read_csv_type3
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

def pca(source, fname_lift_train, fname_lift_test, fname_shape_train, fname_shape_test):
    
    # r_rate = [1, 2, 4, 8]
    r_rate = [1]
    # s_rate = [2, 4, 8]
    s_rate = [8]
    # s_skiptype = [True, False]
    s_skiptype = True
    # r_rate = [1, 2]
    # r_rate = [4, 8]
    # r_rate = [16, 32]
    # r_rate = [64, 160]

    for sr in s_rate:
        for rr in r_rate:
            if rr == 1:
                s_odd = 0  # 全部読みだす
            elif fname_shape_train.find("fourier") != -1:
                s_odd = 3  # 前方から読み出す(fourier用)
            else:
                s_odd = 4  # 全体にわたって等間隔に読み出す(equidistant, dense用)

            X_train, y_train, scalar = read_csv_type3(source, fname_lift_train, fname_shape_train, shape_odd = s_odd, read_rate = rr, skip_rate=sr, return_scalar = True)

            """
            X_train, y_test = read_csv_type3(source, fname_lift_test, fname_shape_test, shape_odd = s_odd,
                                                      read_rate = rr, skip_rate = sr, scalar = scalar)
            """
            pca = PCA(n_components = 8)
            pca.fit(X_train)
            # transformed = pca.fit_transform(X_train)
            print(pca.explained_variance_ratio_)
            print(sum(pca.explained_variance_ratio_))
            """
            exit()

            for i in range(len(transformed)):
                plt.scatter(transformed[i,0], transformed[i,1])
            plt.title("(iv)")
            plt.xlabel("pc1")
            plt.ylabel("pc2")
            plt.show()
            """
            

def cos_

if __name__ == '__main__':
    source = "G:\\Toyota\\Data\\Compressible_Invicid\\training_data\\"
    fname_lift_train = "NACA4\\s1122_e9988_s4_a014.csv"
    fname_lift_test = "NACA5\\s21011_e25190_s1_a014.csv"
    fname_shape_train = "NACA4\\shape_fourier_1112_9988_s04.csv"
    fname_shape_test = "NACA5\\shape_fourier_21011_25190_s1.csv"

    pca(source, fname_lift_train, fname_lift_test, fname_shape_train, fname_shape_test)
    