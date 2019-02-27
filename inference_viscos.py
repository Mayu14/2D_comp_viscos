# -- coding: utf-8 --
import numpy as np
import os
import glob
from keras.models import model_from_json
from read_training_data_viscos import read_csv_type3
from scatter_plot_viscos import make_scatter_plot
from sklearn.metrics import r2_score, mean_squared_error
import pandas as pd
import matplotlib.pyplot as plt

def inference(source, x_test, y_test, case_name, scatter=True, anglerplot=False, return_r2rms=False):
    json_name = "learned\\" + case_name + "_mlp_model_.json"
    weight_name = "learned\\" + case_name + "_mlp_weight.h5"

    # あとでこの辺を自由に変更できるようにする
    # fname_lift_train = "NACA4\\s0000_e5000_a040_odd.csv"
    # fname_shape_train = "NACA4\\shape_fourier_5000_odd.csv"
    # X_test, y_test = read_csv_type3(source, fname_lift_train, fname_shape_train, shape_odd=0, read_rate=1)

    model = model_from_json(open(source + json_name).read())
    model.load_weights(source + weight_name)

    model.summary()

    model.compile(loss="mean_squared_error",
                  optimizer='Adam')

    # score = model.evaluate()
    # print('test loss :', score[0])
    # print('test accuracy :', score[1])
    
    y_predict = model.predict(x_test)
    r2 = r2_score(y_test, y_predict)
    rms = np.sqrt(mean_squared_error(y_test, y_predict))
    if scatter:
        case_name += "r2_" + str(r2) + "_rms_" + str(rms)
        make_scatter_plot(y_test[:, 1], y_predict[:, 1], "CL(Exact)", "CL(Predict)",
                          path = "G:\\Toyota\\Data\\Compressible_Invicid\\fig_post\\", fname = "CL_" + case_name)
        make_scatter_plot(y_test[:, 0], y_predict[:, 0], "CD(Exact)", "CD(Predict)",
                          path = "G:\\Toyota\\Data\\Compressible_Invicid\\fig_post\\", fname = "CD_" + case_name)
        make_scatter_plot(y_test[:,1], y_predict[:,1], "CL(Exact)", "CL(Predict)", path="G:\\Toyota\\Data\\Compressible_Invicid\\fig_post\\", fname="CL_"+case_name, fix_scale = True)
        make_scatter_plot(y_test[:,0], y_predict[:,0], "CD(Exact)", "CD(Predict)",
                          path = "G:\\Toyota\\Data\\Compressible_Invicid\\fig_post\\", fname = "CD_" + case_name, fix_scale = True)

    if anglerplot:
        tekito = (99 + 13) * 40  # 22012
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(x_test[tekito:tekito + 10, 0], y_test[tekito:tekito + 10], marker="o", markersize=10, color="dodgerblue", label = "Exact")
        ax.plot(x_test[tekito:tekito + 10, 0], y_predict[tekito:tekito + 10], marker="o", markersize=10, color="mediumvioletred", label = "Estimate")
        ax.plot(x_test[tekito:tekito + 10, 0], y_test[tekito:tekito + 10], color="dodgerblue", linestyle="dashdot", linewidth=1.0)
        ax.plot(x_test[tekito:tekito + 10, 0], y_predict[tekito:tekito + 10], color="mediumvioletred", linestyle="dashdot", linewidth=1.0)
        ax.legend(bbox_to_anchor=(0, 1), loc="upper left", borderaxespad=1, fontsize=12)
        ax.set_xlabel("Angle of Attack [deg]")
        ax.set_ylabel("$\it{C_{L}}$")
        ax.set_title("NACA22012 Wing $\it{C_{L}}$ , $\it{C_{D}}$ distribution")
        ax.grid(True)

        fig.savefig(source + case_name + "_NACA22012.png")

    if return_r2rms:
        return r2, rms

# 保存先を検索し，ありそうなファイル名を検索，発見したらリストに追加
def case_name_list_generator(source, fname_lift_test, some_case_test=False, some_case = [], scatter=True, anglerplot=False):
    casename_list = []
    top = "learned\\"
    head_list = ["fourierSr_", "concertrate_", "equidistant_"]
    mid1_list = [str(200000), str(100000), str(50000), str(25000)]
    mid2_list = ["_less_angle_", "_less_shape_"]
    tail_rate = [1, 2, 4, 8]
    tail_total = 200
    bottom1 = "_mlp_model_.json"
    bottom2 = "_mlp_weight.h5"
    fname_shape_list = ["NACA5\\shape_fourier_21011_25190_s1.csv"]
    some_case_index = 0
    some_case_total = len(some_case)
    for i in range(3):
        head = head_list[i]
        fname_shape_test = fname_shape_list[i]
        for mid1 in mid1_list:
            for mid2 in mid2_list:
                for rr in tail_rate:
                    tail = str(int(tail_total / rr))
                    fname0 = head + mid1 + mid2 + tail
                    fname1 = source + top + fname0
                    if os.path.exists(fname1 + bottom1):
                        if os.path.exists(fname1 + bottom2):
                            if rr == 1:
                                s_odd = 0  # 全部読みだす
                            elif head == head_list[0]:
                                s_odd = 3  # 前方から読み出す(fourier用)
                                
                            else:
                                s_odd = 4  # 全体にわたって等間隔に読み出す(equidistant, dense用)
                            if some_case_test:
                                if some_case[some_case_index] == fname0:
                                    some_case_index += 1
                                    x_test, y_test = read_csv_type3(source, fname_lift_test, fname_shape_test,
                                                                    total_data = 560 * 14, shape_odd = s_odd, read_rate = rr)
                                    inference(source, x_test, y_test, fname0, scatter, anglerplot)
                            else:
                                x_test, y_test = read_csv_type3(source, fname_lift_test, fname_shape_test,
                                                                total_data = 560 * 14, shape_odd = s_odd, read_rate = rr)
                                inference(source, x_test, y_test, fname0)
                            
                            if (some_case_total != 0) and (some_case_index == some_case_total):
                                print("all process finished")
                                exit()
                                
    return casename_list

def some_case_test(source, fname_lift_test, fname_shape_test, fname_lift_train, fname_shape_train):
    # head_list = ["fourierSr_", "concertrate_", "equidistant_"]
    # mid1_list = [str(200000), str(100000), str(50000), str(25000)]
    # mid2_list = ["_less_angle_", "_less_shape_"]
    # tail_rate = [1, 2, 4, 8]
    # tail_total = 200
    # some_case = ["fourierSr_100000_less_shape_200", "fourierSr_50000_less_angle_200", "fourierSr_25000_less_angle_200"]
    
    # case_name_list_generator(source, fname_lift_test, some_case_test = True, some_case = some_case, scatter = True, anglerplot = True)
    #"""
    # nlist = ["25000", "50000", "100000", "200000"]
    nlist = ["200000"]
    vlist = ["200"]
    """
    dens_list = [[1024],[512,512], [16,32,64],[64,128,256],[128,128,128],[256,256,256]]
    dens_list = [[1024, 1024], [1024,1024,1024], [1024,1024,1024,1024], [1024,1024,1024,1024,1024],
                 [512, 512, 512,512,512], [32,64,128,256,512,1024],[2,16],[2,16,32]]
    """
    # dens_list = [[2048],[2048,2048], [2048,2048,2048],[2048,2048,2048,2048],[2048,2048,2048,2048,2048],[2048,2048,2048,2048,2048,2048],[2048,2048,2048,2048,2048,2048,2048]]
    dens_list = [[1024, 1024, 1024, 1024, 1024],
                 [1024, 1024, 1024, 1024, 1024, 1024], [1024, 1024, 1024, 1024, 1024, 1024,1024]]
    dr = [[12, 24, 48, 96, 192, 384]]
    
    for i in range(5):
        for j in range(7):
            dr.append([2 ** (i + 7)] * (j + 1))
    
        dr.append([2 ** (i + 5), 2 ** (i + 6), 2 ** (i + 7)])
    dens_list = dr
    dens_name = []
    
    for dense in dens_list:
        name = str(len(dense)) + "L"
        for i in range(len(dense)):
            name += "_" + str(dense[i])
        dens_name.append(name)
    some_case = []
    for num in nlist:
        for dname in dens_name:
            for vec in vlist:
                some_case.append("fourierSr_" + num + "_less_angle_" + dname + "_" + vec)

    X_train, y_train, scalar = read_csv_type3(source, fname_lift_train, fname_shape_train, shape_odd = 0, read_rate = 1,
                                      total_data = 0, return_scalar = True)
    x_test, y_test = read_csv_type3(source, fname_lift_test, fname_shape_test,
                                    total_data = 0, shape_odd = 0, read_rate = 1, scalar = scalar)

    mid = []
    for fname0 in some_case:
        r2_test, rms_test = inference(source, x_test, y_test, fname0, scatter=True, anglerplot=True, return_r2rms=True)
        r2_train, rms_train = inference(source, X_train, y_train, fname0, scatter=False, anglerplot=False, return_r2rms=True)
        mid.append([fname0[27:], r2_train, rms_train, r2_test, rms_test])

    columns = ["case", "train_r2", "train_rms", "test_r2", "test_rms"]
    dtypes = {'case': 'object', 'train_r2': 'float64', 'train_rms': 'float64', 'test_r2': 'float64',
              'test_rms': 'float64'}
    df = pd.DataFrame(mid, columns=columns, dtype=dtypes)
    i = 0
    find = True
    while find:
        fname = source + "inference" + str(i).zfill(4) + ".csv"
        find = os.path.exists(fname)
        i += 1

    df.to_csv(fname)

def all_case_test(source, fname_lift_test, fname_shape_test, fname_lift_train, fname_shape_train):
    dir = source + "learned\\"
    fname_head = "fourierSr_200000_less_angle_"
    fname_tail = "_200_mlp.json"

    file_list = []
    file_list = glob.glob(source + fname_head + "\*" + fname_tail)

    for fname in file_list:
        dens_name = fname[len(fname_head):].replace(fname_tail, "")

    print(dens_name)


if __name__ == '__main__':
    source = "G:\\Toyota\\Data\\Compressible_Invicid\\training_data\\"
    fname_lift_test = "NACA5\\s21011_e25190_s1_a014.csv"
    fname_lift_train = "NACA4\\s1122_e9988_s4_a014.csv"
    fname_shape_train = "NACA4\\shape_fourier_1112_9988_s04.csv"
    fname_shape_test = "NACA5\\shape_fourier_21011_25190_s1.csv"

    some_case_test(source, fname_lift_test)
    # case_name_list_generator(source, fname_lift_test)

    
