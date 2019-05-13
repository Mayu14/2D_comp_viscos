# -*- coding: UTF-8 -*-
# 単純な関数のデータを与えて元のデータを予測させてみる
# 学習用データが(TesraK80の)メモリに乗らないため,Generatorを使ってバッチごとにデータをロードさせる感じで
# 20万件のデータを200件ずつ取り出す
from keras.models import Sequential, Model
from keras.layers.core import Dense, Dropout, Activation
from keras.layers.normalization import BatchNormalization
from keras.layers import LeakyReLU, PReLU, Input
from keras.callbacks import EarlyStopping, TensorBoard
import keras.backend.tensorflow_backend as KTF
import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt
import os
from read_training_data_viscos import read_csv_type3
#from scatter_plot import make_scatter_plot
from other_tools.dataset_reduction import data_reduction
from other_tools.complex_layer import MLPLayer
from math import floor



def batch_iter(data, labels, batch_size, shuffle=True):
    num_batches_per_epoch = int((len(data) - 1) / batch_size) + 1

    def data_generator():
        data_size = len(data)
        while True:
            # Shuffle the data at each epoch
            if shuffle:
                shuffle_indices = np.random.permutation(np.arange(data_size))
                shuffled_data = data[shuffle_indices]
                shuffled_labels = labels[shuffle_indices]
            else:
                shuffled_data = data
                shuffled_labels = labels

            for batch_num in range(num_batches_per_epoch):
                start_index = batch_num * batch_size
                end_index = min((batch_num + 1) * batch_size, data_size)
                X = shuffled_data[start_index: end_index]
                y = shuffled_labels[start_index: end_index]
                yield X, y

    return num_batches_per_epoch, data_generator()

def save_my_log(source, case_number, fname_lift_train, fname_shape_train, model_sum):
    with open(source + str(case_number).zfill(4) + "_log.txt", "w") as f:
        f.write("case number :" + str(case_number).zfill((3)) + "\n")
        f.write("training_data of Lift :" + fname_lift_train + "\n")
        f.write("training_data of Shape :" + fname_shape_train + "\n")
        f.write("model summary" + "\n")
        f.write(str(model_sum) + "\n")

def get_case_number(source, env, case_number):
    flag = 0
    source = source + "learned\\"
    if env == "Colab":
        source = source.replace("\\", "/")
        case_number += 10000
    while flag == 0:
        if os.path.exists(source + str(case_number).zfill(5) + "_mlp_model_.json"):
            case_number += 1
        else:
            flag = 1
    return str(case_number).zfill(5)

# case_numberから何のデータだったか思い出せない問題が起きたのでファイル名の命名規則を変更する
# (形状)_(データ数)とする
def get_case_number_beta(case_number, dense_list, rr, sr, skiptype, cluster, preprocess="",
                         criteria_method="nearest_centroid", shape_data=200, total_data=200000, resnet=False,
                         highway=False, densenet=False, useBN = True, useDrop = True, dor = 0.3):
    if int(case_number) / 1000 == 0:
        head = "fourierSr"
    elif int(case_number) / 1000 == 1:
        head = "equidistant"
    elif int(case_number) / 1000 == 2:
        head = "concertrate"
    else:
        "case number error"
        exit()
    mid1 = str(int(total_data / sr))
    if skiptype:
        mid2 = "less_angle"
    else:
        mid2 = "less_shape"
    tail = str(int(shape_data / rr))
    
    tail2 = str(cluster).zfill(5)
    
    mid3 = str(len(dense_list)) + "L"
    
    for i in range(len(dense_list)):
        mid3 += "_" + str(dense_list[i])
        if i == 10:
            mid3 += "..._"
            break
    mid2 += "_new"
    cm = ""
    if criteria_method == "farthest_from_center":
        cm += "_FFC"

    if resnet:
        cm += "_resnet"
    if highway:
        cm += "_highway"
    if densenet:
        cm += "_denseblock"
    elif (resnet==False) and (highway==False):
        cm += "_FC"
    if useBN:
        cm += "_BN"
    if useDrop:
        cm += "_DR" + str(dor)

    return head + "_" + mid1 + "_" + mid2 + "_" + mid3 + "_" + tail + "_" + tail2 + preprocess + cm

def main(fname_lift_train, fname_shape_train, fname_lift_test, fname_shape_test, case_number, case_type=3, env="Lab", validate=True, gpu_mem_usage=0.5):
    # r_rate = [1, 2, 4, 8]
    # s_rate = [1, 2, 4, 8]
    # s_skiptype = [True, False]
    s_skiptype = True
    r_rate = [1]
    s_rate = [1]
    # r_rate = [1, 2]
    # r_rate = [4, 8]
    # r_rate = [16, 32]
    # r_rate = [64, 160]
    sr = 1
    rr = 1
    # dr = [[12, 24, 48, 96, 192, 384]]
    block_total = 8
    dr = [[128]*block_total]
    weight_layer_list = [[3]*block_total]
    """
    dr = [[100] * 100]
    dr = []
    dr.append(202)
    for i in range(1, 18):
        dr.append(int(dr[i - 1] / 5.0 * 4.0))

    dr = [dr]
    print(dr)
    exit()
    """
    # for i in range(5):
    data_reduction_test = False
    # i = 1
    # for j in range(6, 20):
        # dr.append([1024]*(j+1))
    criteria_method = "farthest_from_center"
    useBN_list = [True, False]
    before_activate = True
    dor_list = [0.0, 0.1, 0.2, 0.3, 0.4]
    resnet = False
    high_way = False
    dense_net = False
    # preprocesses = ["None"]
    preprocesses = [""]#, "rbf", "poly", "linear", "cosine", "sigmoid", "PCA"]
    dense_list = dr[0]
    # for dense_list in dr:
    # for reduct in range(40):
    for dor in dor_list:
        # cluster = 500 * (reduct + 1)
        cluster = 22680
        # for rr in r_rate:
        i = 0
        for useBN in useBN_list:
            preprocess = ""
            if i == 0:
                dense_net = True
                high_way = False
            elif i == 2:
                resnet = True
                dense_net = False
            elif i == 4:
                high_way = True
                resnet = False
            else:
                high_way = False

            i += 1
            
            if rr == 1:
                s_odd = 0   # 全部読みだす
            elif fname_shape_train.find("fourier") != -1:
                s_odd = 3   # 前方から読み出す(fourier用)
            else:
                s_odd = 4   # 全体にわたって等間隔に読み出す(equidistant, dense用)

            config = tf.ConfigProto()
            config.gpu_options.per_process_gpu_memory_fraction = gpu_mem_usage
            KTF.set_session(tf.Session(config = config))
            old_session = KTF.get_session()

            with tf.Graph().as_default():
                source = "Compressible_Invicid\\training_data\\"
                if env == "Lab":
                    source = "G:\\Toyota\\Data\\" + source
                    # case_num = get_case_number(source, env, case_number)
                    case_num = get_case_number_beta(case_number, dense_list, rr, sr, s_skiptype, cluster, preprocess, criteria_method, resnet=resnet, highway=high_way, densenet=dense_net, useBN = useBN, useDrop = useDrop, dor = dor)
                    log_name = "learned\\" + case_num + "_tb_log.hdf5"
                    json_name = "learned\\" + case_num + "_mlp_model_.json"
                    weight_name = "learned\\" + case_num + "_mlp_weight.h5"
                elif env == "Colab":
                    source = "/content/drive/Colab Notebooks/" + source.replace("\\", "/")
                    case_num = get_case_number(source, env, case_number)
                    log_name = "learned/" + case_num + "_log.hdf5"
                    json_name = "learned/" + case_num + "_mlp_model_.json"
                    weight_name = "learned/" + case_num + "_mlp_weight.h5"
                print(case_num)


                session = tf.Session('')
                KTF.set_session(session)
                KTF.set_learning_phase(1)
                
                # model = Sequential()
                if case_type == 3:
                    # ここ書き換えポイント
                    X_train, y_train, scalar = read_csv_type3(source, fname_lift_train, fname_shape_train, shape_odd = s_odd, read_rate = rr, skip_rate=sr, total_data = 0, return_scalar = True)
                    
                    if data_reduction_test:
                        X_train, y_train = data_reduction(X_train, y_train, reduction_target = cluster, output_csv = False, preprocess = preprocess, criteria_method = criteria_method)
                        
                    if validate:
                        x_test, y_test = read_csv_type3(source, fname_lift_test, fname_shape_test, total_data = 0, shape_odd=s_odd, read_rate = rr, scalar = scalar)

                
                        
                input_vector_dim = X_train.shape[1]

                def simple_network(dense_list, inputs, resnet = False, highwaynet = False, dense_net = False, units_list = None):
                    # leaky_relu = LeakyReLU()
                    # input layer
                    MLP = MLPLayer(activate_before_fc=before_activate, batch_normalization=useBN, dropout=useDrop,
                                   dropout_rate=dor, dropout_timing = "final", gate_bias=-3, growth_rate = 32)

                    x = Dense(units = dense_list[0])(inputs)
                    # mid layer
                    for i in range(1, len(dense_list)):
                        leaky_relu = LeakyReLU()
                        if dense_net:
                            x = MLP.denseblock(inputs=x, Activator=leaky_relu,
                                               weight_layer_number=weight_layer_list[i])
                            x = Dense(units=dense_list[i])(x)
                        else:
                            # units_list = [dense_list[i], int(dense_list[i]/2), dense_list[i]]
                            units_list = [dense_list[i]*weight_layer_list]
                            if resnet:
                                x = MLP.residual(inputs = x, Activator = leaky_relu,
                                                 weight_layer_number=weight_layer_list[i], units_list=units_list)
                            else:
                                if highwaynet:
                                    x = MLP.highway(inputs = x, Activator = leaky_relu,
                                                    weight_layer_number=weight_layer_list[i], units_list=units_list)
                                else:
                                    x = MLP.fully_connected(units=dense_list[i], inputs=x, Activator=leaky_relu)

                    return x
                
                with tf.name_scope("inference") as scope:
                    inputs = Input(shape = (input_vector_dim,))
                    
                    x = simple_network(dense_list, inputs, resnet = resnet, highwaynet = high_way, dense_net=dense_net, units_list=units_list)
                    
                    """
                    x = Dense(units = dense_net_list[0])(inputs)
                    # x = Activation(LeakyReLU())(x)
                    x = LeakyReLU()(x)
                    """
                    """
                    for i in range(1, len(dense_net_list)):
                        dense_block =
                        x = Dense(units = dense_list[0])(x)
                        leaky_relu = LeakyReLU()
                        # x = residual(inputs=x, Activator=leaky_relu, batch_normalization=True, dropout=True)
                        x = highway(inputs=x, Activator=leaky_relu, batch_normalization=True, dropout=True)
                    """
                    
                    # output layer
                    predictions = Dense(units = 2, activation = None)(x)
                
                model = Model(inputs = inputs, outputs = predictions)

                save_my_log(source, case_number, fname_lift_train, fname_shape_train, model.summary())
                # es_cb = EarlyStopping(monitor='val_loss', patience=0, verbose=0, mode='auto')
                tb_cb = TensorBoard(log_dir=source + log_name, histogram_freq=0, write_grads=True)

                model.compile(loss="mean_squared_error",
                              optimizer='Adam')

                batch_size = y_train.shape[0]
                threshold = 30000
                if batch_size > threshold:
                    split = floor(float(batch_size) / threshold)
                    batch_size = floor(float(batch_size) / split)
                    
                train_steps, train_batches = batch_iter(X_train, y_train, batch_size)
                if validate:
                    valid_steps, valid_batches = batch_iter(x_test, y_test, batch_size)
                #"""
                model.fit(x=X_train, y=y_train,
                          batch_size=batch_size, nb_epoch=1000,
                          validation_split=0.1, callbacks=[tb_cb])
                #"""
                """
                model.fit_generator(train_batches, train_steps,
                                    epochs=1000,
                                    validation_data=valid_batches,
                                    validation_steps=valid_steps,
                                    callbacks=[tb_cb])
                """
                # X_train: [number, angle, shape001, shape002, ..., shapeMAX]
                # y_train: [number, lift]
                # 適当に中央付近の翼を抜き出しての-40-38degreeをプロットさせてみる
                """
                tekito = 1306 * 40  # NACA2613 or NACA2615
                plt.figure()
                plt.plot(X_train[tekito:tekito+40, 0], y_train[tekito:tekito+40])
                plt.plot(X_train[tekito:tekito+40, 0], model.predict(X_train)[tekito:tekito+40])
                plt.savefig(source + case_num + "_train.png")

                y_predict = model.predict(x_test)
                tekito = (99 + 13) * 40 # 22012
                plt.figure()
                plt.plot(x_test[tekito:tekito+40, 0], y_test[tekito:tekito+40])
                plt.plot(x_test[tekito:tekito+40, 0], y_predict[tekito:tekito+40])
                plt.savefig(source + case_num + "_test.png")
                
                make_scatter_plot(y_test, y_predict, "CL(Exact)", "CL(Predict)", path="G:\\Toyota\\Data\\Incompressible_Invicid\\fig\\", fname=case_num)
                """

            json_string = model.to_json()
            open(source + json_name, 'w').write(json_string)
            model.save_weights(source + weight_name)
            KTF.set_session(old_session)


if __name__ == '__main__':
    # env_in = input("Please set envirionment: 0:Lab, 1:Colab")
    env_in = os.name
    if env_in == "nt":
        env = "Lab"
    else:
        env = "Colab"

    # shape_type = input("please set shape_type: 0:fourier, 1:equidistant, 2:dense")
    # for i in range(1,3):
    i = 0
    shape_type = str(i)
    fname_lift_train = "NACA4\\s1122_e9988_s4_a014.csv"
    fname_lift_test = "NACA5\\s21011_e25190_s1_a014.csv"

    if shape_type == str(0):
        fname_shape_train = "NACA4\\shape_fourier_1112_9988_s04.csv"
        fname_shape_test = "NACA5\\shape_fourier_21011_25190_s1.csv"
        case_number = 0
        #"""
    elif shape_type == str(1):
        fname_shape_train = "NACA4\\shape_equidistant_1112_9988_s04.csv"
        fname_shape_test = "NACA5\\shape_equidistant_21011_25190_s1.csv"
        case_number = 1000

    elif shape_type == str(2):
        fname_shape_train = "NACA4\\shape_crowd_0.1_0.15_30_50_20_1112_9988_d4.csv"
        fname_shape_test = "NACA5\\shape_crowd_0.1_0.15_30_50_20_560_new.csv"
        case_number = 2000
        #"""
    else:
        print("shape_type error")
        exit()

    if env == "Colab":
        fname_lift_train = fname_lift_train.replace("\\", "/")
        fname_shape_train = fname_shape_train.replace("\\", "/")
        fname_lift_test = fname_lift_test.replace("\\", "/")
        fname_shape_test = fname_shape_test.replace("\\", "/")

    main(fname_lift_train, fname_shape_train, fname_lift_test, fname_shape_test, case_number, case_type=3, env=env, validate=False)