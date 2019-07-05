# coding: utf-8
import os
import numpy as np
from sklearn.model_selection import train_test_split
from keras.models import Model, model_from_json
from keras.layers import Dense, Dropout, Flatten, Input, Activation, Concatenate, Lambda
from keras.layers import Conv2D, MaxPooling2D
from keras.optimizers import Adadelta
import tensorflow as tf
from keras import backend as K
from keras.callbacks import EarlyStopping, ModelCheckpoint, TensorBoard
from keras.utils import plot_model
from sklearn.metrics import r2_score
from other_tools.ResNet_builder import ResnetBuilder

import matplotlib.pyplot as plt
import csv
from itertools import product

def load_data(path, size, type=4):
    if type == 4:
        header = "train_"
    elif type == 5:
        header = "test_"
    fname = path + header + str(size).zfill(3) + "_" + str(size).zfill(3) + ".npz"
    with np.load(fname, allow_pickle = True) as f:
        x_train_img, x_train_param, y_train = f["x_train_img"], f["x_train_param"], f["y_train"]

    x_train_img = x_train_img.astype(np.float32).reshape(x_train_img.shape[0], size, size, 1)
    x_train_img /= 255.0
    x_train_param[:, 0] /= 39.0   # angle
    x_train_param[:, 1] /= 100.0  # aspect
    return x_train_img, x_train_param, y_train

def main(mach=False, reynolds=False):
    prm_num = 2
    if mach:
        prm_num += 1
    if reynolds:
        prm_num += 1

    source = "G:\\Toyota\\Data\\Compressible_Invicid\\training_data\\"
    temporary = "D:\\Toyota\\Data\\temp\\"
    npz_path_train = source + "NACA4\\"
    npz_path_test = source + "NACA5\\"
    fig_dir = source + "cnn\\fig\\"
    source += "cnn\\"
    source_learned = source + "learned\\"
    size = 256
    x_train_img, x_train_param, y_train = load_data(path=npz_path_train, size=size, type=4)
    x_test_img, x_test_param, y_test = load_data(path=npz_path_test, size=size, type=5)
    
    x_train_img, x_valid_img, x_train_param, x_valid_param, y_train, y_valid = train_test_split(x_train_img, x_train_param, y_train, test_size = 0.175)

    config = tf.ConfigProto()
    config.gpu_options.per_process_gpu_memory_fraction = 0.75

    filter_size = [2,4]

    filter_lists = []
    for layer in range(10, 15):
        p = product(filter_size, repeat=layer)
        for v in p:
            v = list(v)
            flag = True
            for j in range(layer - 1):
                if v[0] != 2:
                    flag = False

                if v[j] > v[j + 1]:
                    flag = False

            if flag:
                filter_lists.append(v)

    filter_lists = [[2, 2, 2]]    # fixed

    drop_rates = [0.0, 0.25, 0.5]
    dropList = []
    p = product(drop_rates, repeat=2)
    for v in p:
        dropList.append(list(v))
    dropList = [[0.5, 0.25]]

    dense_unit = [128, 256, 512, 1024, 2048]
    unit_lists = []
    p = product(dense_unit, repeat=2)
    for v in p:
        unit_lists.append(list(v))
    unit_lists = [[256, 128]]
    # for filter_list in filter_lists:
    filter_list = filter_lists[0]
    drop = dropList[0]
    # for drop in dropList:
    for unit in unit_lists:
        name = "size_" + str(size)
        """
        for i, filter in enumerate(filter_list):
            name += "_f" + str(i) + "_" + str(filter)
        """
        name += "ResNet-18"
        name += "_u1_" + str(unit[0]) + "_u2_" + str(unit[1])
        # name += "_d1_" + str(drop[0]) + "_d2_" + str(drop[1])
        
        
        with tf.Session(config=config) as sess:
            K.set_session(sess)

            model_json = source_learned + name + ".json"
            weight_hdf5 = source_learned + name + ".hdf5"
            log_name = source_learned + name + "_log.hdf5"
            baseSaveDir = temporary + "checkpoints\\"
            chkpt = baseSaveDir + 'MLP_.{epoch:02d}-{val_loss:.2f}.hdf5'

            if not os.path.exists(model_json):
                input_shape = (size, size, 1)
                num_classes = unit[0] - prm_num
    
                # model = ResnetBuilder.build_resnet_18(input_shape, num_classes)
                inputs_img, x = ResnetBuilder.build_resnet_18_continues(input_shape, num_classes)

                inputs_prm = Input(shape = (prm_num,), name = "input_prm")
                """
                inputs_img = Input(shape = input_shape, name="input_img")
                x = Conv2D(filter_list[0], kernel_size = (3, 3), input_shape = (size, size, 1))(inputs_img)
                x = Activation("relu")(x)
                for i, n in enumerate(filter_list):
                    if i != 0:
                        x = Conv2D(n, (3, 3))(x)
                        x = Activation("relu")(x)
                x = MaxPooling2D(pool_size = (2, 2))(x)
                x = Dropout(drop[0])(x)
                x = Flatten()(x)
                x = Dropout(drop[1])(x)
                """
                x_prm = Lambda(lambda x: x, input_shape = (prm_num,), output_shape=(prm_num,))(inputs_prm)
                x = Concatenate()([x, x_prm])
                x = Dense(unit[0])(x)
                x = Activation("relu")(x)
                x = Dense(unit[1])(x)
                x = Activation("relu")(x)
                x = Dense(2)(x)
                outputs = Activation("linear")(x)

                model = Model(inputs=[inputs_img, inputs_prm], outputs=outputs)
                epoch = 50
            else:
                model = model_from_json(open(model_json).read())
                model.load_weights(weight_hdf5)
                epoch = 1000

            model.compile(loss="mean_squared_error",
                          optimizer=Adadelta(),
                          metrics=["mae"])

            cp_cb = ModelCheckpoint(filepath=chkpt, monitor="val_loss", verbose=0, save_best_only=True, mode='auto')
            tb_cb = TensorBoard(log_dir=log_name, histogram_freq=0, write_grads=True)
            es_cb = EarlyStopping(monitor='val_loss', patience=50, verbose=0, mode='auto')

            fname_model = fig_dir + name + "_model.png"
            fname_hist = fig_dir + name  + "_history.png"
            plot_model(model, to_file=fname_model, show_shapes=True)

            batch_size = 16
            # """
            print(name)
            history = model.fit([x_train_img, x_train_param],
                                y=y_train,
                                verbose=0,
                                batch_size=batch_size,
                                epochs=epoch,
                                validation_split=0.1,
                                validation_data = ([x_valid_img, x_valid_param], y_valid),
                                callbacks=[cp_cb, tb_cb, es_cb])

            json_string = model.to_json()
            open(model_json, 'w').write(json_string)
            model.save_weights(weight_hdf5)
            print(history.history)
            fig = plt.figure()
            """
            ax = fig.add_subplot(1, 2, 1)
            ax.plot(history.history["acc"], marker = ".", label = "acc")
            ax.plot(history.history["val_acc"], marker = ".", label = "val_acc")
            ax.set_title("model accuracy")
            ax.grid(True)
            ax.set_xlabel("epoch")
            ax.set_ylabel("accuracy")
            ax.legend(loc = "best")
            """
            #ax = fig.add_subplot(1, 2, 2)
            ax = fig.add_subplot(1, 1, 1)
            ax.plot(history.history['loss'], marker = '.', label = 'loss')
            ax.plot(history.history['val_loss'], marker = '.', label = 'val_loss')
            ax.set_title('model loss')
            ax.grid(True)
            ax.set_xlabel('epoch')
            ax.set_ylabel('loss')
            ax.legend(loc = 'best')

            #plt.show()
            plt.savefig(fname_hist)
            plt.close()
            score_valid = model.evaluate([x_valid_img, x_valid_param], y_valid, verbose = 0)
            score_test = model.evaluate([x_test_img, x_test_param], y_test, verbose=0)
            print(score_valid)
            print(score_test)
            y_pred = model.predict([x_valid_img, x_valid_param])
            r2_valid = r2_score(y_valid, y_pred)
            y_pred_test = model.predict([x_test_img, x_test_param])
            r2_test = r2_score(y_test, y_pred_test)
            print(r2_valid, r2_test)
            log = [r2_valid, r2_test, score_valid[0], score_valid[1], score_test[0], score_test[1]]
            log.extend(name.split("_"))
            with open(source + "log.csv", "a") as f:
                writer = csv.writer(f, lineterminator='\n')
                writer.writerow(log)

        K.clear_session()

if __name__ == '__main__':
    main()
    