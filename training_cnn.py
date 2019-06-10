# coding: utf-8
import numpy as np
from sklearn.model_selection import train_test_split
from keras.models import Sequential, Model
from keras.layers import Dense, Dropout, Flatten, Input, Activation, Concatenate, Lambda
from keras.layers import Conv2D, MaxPooling2D
from keras.layers.normalization import BatchNormalization
from keras.layers import LeakyReLU, PReLU
from keras.callbacks import EarlyStopping, TensorBoard
import keras.backend.tensorflow_backend as KTF
import tensorflow as tf

import matplotlib.pyplot as plt
import os
from read_training_data_viscos import read_csv_type3
from other_tools.dataset_reduction import data_reduction
from other_tools.complex_layer import MLPLayer
from math import floor

def load_data(path, size):
    fname = path + "train_" + str(size).zfill(3) + "_" + str(size).zfill(3) + ".npz"
    with np.load(fname, allow_pickle = True) as f:
        x_train_img, x_train_param, y_train = f["arr_0"], f["arr_1"], f["arr_2"]
    
    return x_train_img, x_train_param, y_train

def main(mach=False, reynolds=False):
    prm_num = 2
    if mach:
        prm_num += 1
    if reynolds:
        prm_num += 1

    source = "G:\\Toyota\\Data\\Compressible_Invicid\\training_data\\"
    npz_path = source + "NACA4\\"
    size = 512
    x_train_img, x_train_param, y_train = load_data(path=npz_path, size=size)
    
    x_train_img, x_valid_img, x_train_param, x_valid_param, y_train, y_valid = train_test_split(x_train_img, x_train_param, y_train, test_size = 0.175)

    inputs_img = Input(shape = (512, 512, 1), name="input_img")
    inputs_prm = Input(shape = (prm_num, ), name="input_prm")
    x = Conv2D(32, kernel_size = (3, 3), input_shape = (512, 512, 1))(inputs_img)
    x = Activation("relu")(x)
    x = Conv2D(64, (3, 3))(x)
    x = Activation("relu")(x)
    x = MaxPooling2D(pool_size = (2, 2))(x)
    x = Dropout(0.25)(x)
    x = Flatten()(x)
    x = Dropout(0.5)(x)
    x_prm = Lambda(lambda x: x, output_shape=(prm_num,))(inputs_prm)
    x = Concatenate()([x, x_prm])
    x = Dense(256)(x)
    x = Activation("relu")(x)
    x = Dense(128)(x)
    outputs = Activation("linear")(x)

    model = Model(inputs=[inputs_img, inputs_prm], outputs=outputs)

    model.compile(loss="mean_squared_error",
                  optimizer='Adam')


    batch_size = 512
    # """
    model.fit({"input_img": x_train_img, "input_prm": x_train_param},
              y=y_train,
              batch_size=batch_size, epochs=50,
              validation_split=0.1)

    json_string = model.to_json()
    open(source + "testcnn.json", 'w').write(json_string)
    model.save_weights(source + "testweight.h5")

if __name__ == '__main__':
    main()
    