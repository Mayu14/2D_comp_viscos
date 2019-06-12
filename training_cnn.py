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
from keras import callbacks
from sklearn.metrics import r2_score

import matplotlib.pyplot as plt

def load_data(path, size):
    fname = path + "train_" + str(size).zfill(3) + "_" + str(size).zfill(3) + ".npz"
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
    npz_path = source + "NACA4\\"
    size = 128
    x_train_img, x_train_param, y_train = load_data(path=npz_path, size=size)
    
    x_train_img, x_valid_img, x_train_param, x_valid_param, y_train, y_valid = train_test_split(x_train_img, x_train_param, y_train, test_size = 0.175)

    config = tf.ConfigProto()
    config.gpu_options.per_process_gpu_memory_fraction = 0.5
    sess = tf.Session(config=config)
    K.set_session(sess)

    model_json = source + "testcnn" + str(size).zfill(3) + ".json"
    weight_hdf5 = source + "testweight" + str(size).zfill(3) + ".hdf5"
    if not os.path.exists(model_json):
        inputs_img = Input(shape = (size, size, 1), name="input_img")
        inputs_prm = Input(shape = (prm_num, ), name="input_prm")
        x = Conv2D(32, kernel_size = (3, 3), input_shape = (size, size, 1))(inputs_img)
        x = Activation("relu")(x)
        x = Conv2D(64, (3, 3))(x)
        x = Activation("relu")(x)
        x = MaxPooling2D(pool_size = (2, 2))(x)
        x = Dropout(0.25)(x)
        x = Flatten()(x)
        x = Dropout(0.5)(x)
        x_prm = Lambda(lambda x: x, input_shape = (prm_num,), output_shape=(prm_num,))(inputs_prm)
        x = Concatenate()([x, x_prm])
        x = Dense(256)(x)
        x = Activation("relu")(x)
        x = Dense(128)(x)
        x = Activation("relu")(x)
        x = Dense(2)(x)
        outputs = Activation("linear")(x)

        model = Model(inputs=[inputs_img, inputs_prm], outputs=outputs)
    else:
        model = model_from_json(open(model_json).read())
        model.load_weights(weight_hdf5)

    model.compile(loss="mean_squared_error",
                  optimizer=Adadelta(),
                  metrics=["mae"])

    es_cb = callbacks.EarlyStopping(monitor='val_loss', patience=0, verbose=0, mode='auto')

    batch_size = 1
    # """

    history = model.fit({"input_img": x_train_img, "input_prm": x_train_param},
                        y=y_train,
                        batch_size=batch_size,
                        epochs=50,
                        validation_split=0.1,
                        validation_data = ([x_valid_img, x_valid_param], y_valid),
                        callbacks=[es_cb])

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
    ax = fig.add_subplot(1, 2, 2)
    ax.plot(history.history['loss'], marker = '.', label = 'loss')
    ax.plot(history.history['val_loss'], marker = '.', label = 'val_loss')
    ax.set_title('model loss')
    ax.grid(True)
    ax.set_xlabel('epoch')
    ax.set_ylabel('loss')
    ax.legend(loc = 'best')

    plt.show()

    score = model.evaluate([x_valid_img, x_valid_param], y_valid, verbose = 0)
    print(score)
    y_pred = model.predict([x_valid_img, x_valid_param])
    r2 = r2_score(y_valid, y_pred)
    print(r2)

if __name__ == '__main__':
    main()
    