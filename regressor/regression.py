# coding: utf-8
import os
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from sklearn.preprocessing import MinMaxScaler
import tensorflow as tf

from keras.models import Model, model_from_json
from keras.layers import Input, Dense, Conv2D, Conv2DTranspose, Flatten, Reshape, Concatenate, Dropout, Lambda, Subtract
from keras.layers import BatchNormalization, PReLU, Activation, LeakyReLU, Embedding
from keras.callbacks import EarlyStopping, ModelCheckpoint, TensorBoard
from keras.optimizers import Adadelta, Adam
from keras.losses import categorical_crossentropy, mean_squared_error, mean_absolute_error
from keras.utils import plot_model
from keras.initializers import RandomNormal
from keras import regularizers
from keras.backend import tensorflow_backend as KTF
from keras import backend as K

import optuna
from optuna.integration import KerasPruningCallback

from training.read_training_data import load_mixed
from other_tools.scatter_plot import make_scatter_plot

def __errors(y_pred, y_true):
    return y_pred.flatten() - y_true.flatten()

def max_error(y_pred, y_true):
    err = __errors(y_pred, y_true)
    return err[np.argmax(np.abs(err))]

def min_error(y_pred, y_true):
    err = __errors(y_pred, y_true)
    return err[np.argmin(np.abs(err))]

def rmsq_error(y_pred, y_true):
    err = __errors(y_pred, y_true)
    return np.sqrt(np.sum(err**2) / err.shape[0])

def abs_median_error(y_pred, y_true):
    err = __errors(y_pred, y_true)
    return np.median(np.abs(err))

def error_print(y_pred, y_true, csv=False):
    r2 = r2_score(y_pred, y_true)
    rmsE = rmsq_error(y_pred, y_true)
    medE = abs_median_error(y_pred, y_true)
    maxE = max_error(y_pred, y_true)
    minE = min_error(y_pred, y_true)
    if not csv:
        return "r2:{4}, rmse:{0}_median:{3}_max:{1}_min:{2}\n".format(rmsE, maxE, minE, medE, r2)
    else:
        return "{4},{0},{1},{2},{3}\n".format(rmsE, medE, maxE, minE,r2)

def fnameGen(title, r2 = "", num=None, old=False, noExt=False, name=None):
    def set_new_name(basename = "log", extension = ".txt", zfill = 5):
        new_name = lambda i: basename + "_" + str(i).zfill(zfill) + extension
        i = 0
        while os.path.exists(new_name(i)):
            i += 1
        return new_name(i)

    if not noExt:
        string = title.split(".")
        title = string[0]
        ext = "." + string[1]
    else:
        title = ""
        ext = ""
    if name is None:
        if old:
            if not r2 == "":
                r2 = "_R2_{0:.5f}".format(r2)
            basename = CODE_NAME + "_" + shapeType.replace("Concat", "") + reductTarget + "_{0}{1}".format(title, r2)
            if not num is None:
                return Path(basename + "_" + str(num).zfill(5) + ext)
            return Path(set_new_name(basename, ext))
        else:
            arch = r2
            basename = CODE_NAME + "_" + shapeType.replace("Concat", "") + reductTarget + "_{0}_{1}".format(title,arch)
            return Path(basename)
    else:
        basename = CODE_NAME + "_" + name.replace("_", ",") + "_{0}".format(title)
        return Path(basename + ext)

# COMMON
CODE_NAME = "CV_mixed"
shapeType = "EquidistantConcat"#""FourierConcat"
reductTarget = "kmeans_norm_pca_9000"#"NACA_half_4"#"kmeans_norm_pca_4500"#"NACA_similar"#"NACA_half"#"AoA"

epoch = 5000
batch_size = 500

baseSaveDir = Path("G:\\Toyota\\Data\\Compressible_Viscos\\rhoSimpleFoam\\training_data\\learned\\" + shapeType.replace("Concat", ""))
modelDir = baseSaveDir / Path("model")
figsDir = baseSaveDir / Path("figs")

n_classes = 2   # lift, drag

load_model = False
epoch_overwrite = 1

load_shapeType = shapeType.replace("Concat", "")
r2_rec = ""#0.99348
num = 0

if load_model:
    epoch = epoch_overwrite
# for MLP
input_dim = 203 # s:200, AoA, Ma, Re

# for CNN
msdf_res = 128   # pixel
img_shape = (msdf_res, msdf_res, 3)
n_labeles = 2   # Ma, Re

def regressor(droprate=0.25, opt=False):
    if opt:
        model = model_from_json(open('keras_model_opt2_{0}.json'.format(shapeType)).read())
        return model
    else:
        inputs = Input(shape = (input_dim,), name = "input_layer")
        reg = regularizers.l2(0.00001)
        rand_init = RandomNormal(stddev = 0.02)
    
        x = PReLU(shared_axes = [1])(inputs)
        
        x = Dense(units = 151, activation = None, kernel_initializer = rand_init,
                  kernel_regularizer = reg, bias_regularizer = reg, name = "dense_01")(x)
        x = BatchNormalization(axis = -1)(x)
        x = PReLU(shared_axes = [1])(x)
        
        x = Dense(units = 101, activation = None, kernel_initializer = rand_init,
                  kernel_regularizer = reg, bias_regularizer = reg, name = "dense_02")(x)
        x = BatchNormalization(axis = -1)(x)
        x = PReLU(shared_axes = [1])(x)
        # x = Dropout(0.125)(x)
        x = Dense(units = 51, activation = None, kernel_initializer = rand_init,
                    kernel_regularizer = reg, bias_regularizer = reg, name = "dense_03")(x)
        x = BatchNormalization(axis = -1)(x)
        x = PReLU(shared_axes = [1])(x)
        if droprate != 0.0:
            x = Dropout(droprate)(x)
    
        prediction = Dense(units = n_classes, activation = None, kernel_initializer = RandomNormal(stddev = 0.02),
                    kernel_regularizer = reg, bias_regularizer = reg, name = "output_layer")(x)
        return Model(inputs = [inputs], outputs = [prediction])


def regressorOPT(trial):
    inputs = Input(shape = (input_dim,), name = "input_layer")
    alpha = trial.suggest_loguniform(f"alpha", 1e-6, 1e-1)
    reg = regularizers.l2(alpha)
    rand_init = RandomNormal(stddev = 0.02)
    
    x = PReLU(shared_axes = [1])(inputs)
    n_layers = trial.suggest_int("n_layers", 1, 5)
    for i in range(n_layers):
        num_hidden = int(trial.suggest_int(f"n_units_l{i}", 2, 251))
        dropout = trial.suggest_uniform(f"dropout_l{i}", 0.0, 0.5)
        x = Dense(units = num_hidden, activation = None, kernel_initializer = rand_init,
                  kernel_regularizer = reg, bias_regularizer = reg)(x)
        x = BatchNormalization(axis = -1)(x)
        x = PReLU(shared_axes = [1])(x)
        x = Dropout(dropout)(x)
        
    prediction = Dense(units = n_classes, activation = None, kernel_initializer = RandomNormal(stddev = 0.02),
                       kernel_regularizer = reg, bias_regularizer = reg, name = "output_layer")(x)
    return Model(inputs = [inputs], outputs = [prediction])


def regressorOPT_PClowd_noaoa(trial, layers_pts=None, layers_prms=None, aoa_dims=0, set_vars=3, max_points=100, pts_dim=2, commonliftdrag=False, pclowdBN=False, paramBN=True):
    alpha = trial.suggest_loguniform(f"alpha", 1e-6, 1e-1)
    reg = regularizers.l2(alpha)
    rand_init = RandomNormal(stddev = 0.02)
    
    def dense_block(num_points, layers, bn, lastShape = None, j=None):
        inputs = Input(shape = (num_points,), name = "input_layer")
        x = PReLU(shared_axes = [1])(inputs)
        
        n_layers = trial.suggest_int(f"n_layers{j}", 1, 5)
        for i in range(n_layers):
            num_hidden = int(trial.suggest_int(f"n_units_l{i}{j}", 2, 251))
            dropout = trial.suggest_uniform(f"dropout_l{i}{j}", 0.0, 0.5)
            x = Dense(units = num_hidden, activation = None, kernel_initializer = rand_init,
                      kernel_regularizer = reg, bias_regularizer = reg)(x)
            if bn:
                x = BatchNormalization(axis = -1)(x)
            x = PReLU(shared_axes = [1])(x)
            x = Dropout(dropout)(x)
        
        if lastShape is None:
            outputs = Activation("linear")(x)
        else:
            outputs = Dense(units = lastShape, activation = None, kernel_initializer = rand_init,
                            kernel_regularizer = reg, bias_regularizer = reg)(x)
        return Model([inputs], [outputs])
    
    def sets_block(num_points, layers, bn, j=None):
        input_pts = Input(shape = (num_points,))
        pts = dense_block(num_points, layers, bn, j=j)(input_pts)
        Adder = Lambda(lambda x: K.sum(x, axis = 1), name = "Adder", output_shape = (1,))
        pts = Adder(pts)
        output = Reshape((1,))(pts)
        return Model(inputs = [input_pts], outputs = [output])
    
    def aero_block(lastShape = 1):
        input_pts = Input(shape = (max_points * pts_dim,), name = "input_set")
        input_params = Input(shape = (2 + aoa_dims,), name = "input_params")
        
        x = sets_block(max_points * pts_dim, layers_pts, bn = False, j=0)(input_pts)  # Sxy
        
        if set_vars == 3:
            x_pts = Lambda(lambda x: x[:, :max_points], name = "x_pts", output_shape = (max_points,))(input_pts)  # sx
            y_pts = Lambda(lambda x: x[:, max_points:], name = "y_pts", output_shape = (max_points,))(input_pts)  # sy
            x_pts = sets_block(max_points, layers_pts, bn = pclowdBN, j=1)(x_pts)
            y_pts = sets_block(max_points, layers_pts, bn = pclowdBN, j=2)(y_pts)
            x = Concatenate(name = "HiddenConcat")([x_pts, y_pts, x, input_params])
        else:
            x = Concatenate()([x, input_params])
        
        aero = dense_block(set_vars + aoa_dims + 2, layers_prms, bn = paramBN, lastShape = lastShape, j=3)(x)
        return Model([input_pts, input_params], [aero])
    
    input_pts = Input(shape = (max_points * pts_dim,), name = "input_set")
    input_params = Input(shape = (2 + aoa_dims,), name = "input_params")
    if not commonliftdrag:
        lift = aero_block()([input_pts, input_params])
        drag = aero_block()([input_pts, input_params])
        predicted = Concatenate()([lift, drag])
    else:
        predicted = aero_block(lastShape = 2)([input_pts, input_params])
    
    return Model(inputs = [input_pts, input_params], outputs = [predicted])


def objective(trial):
    K.clear_session()
    # model = regressorOPT(trial)
    model = regressorOPT_PClowd_noaoa(trial)
    optimizer = "adam"#trial.suggest_categorical("optimizer", ["sgd", "adam", "rmsprop"])
    model.compile(optimizer = optimizer,
                  loss="mean_squared_error",
                  metrics = ["mae"])
    history = model.fit(x_train, y_train,
                        batch_size = BATCHSIZE,
                        callbacks = [KerasPruningCallback(trial, "val_acc")],
                        epochs = EPOCHS,
                        validation_data = (x_test, y_test),
                        verbose = 2)

    # モデルの評価
    model_json = model.to_json()
    with open('keras_model_optpc_{0}.json'.format(shapeType), 'w') as f_model:
        f_model.write(model_json)
    model.save_weights('keras_model_optpc_{0}.hdf5'.format(shapeType))

    # 最小値探索なので
    # return np.amin(history.history['val_mean_absolute_error'])
    return -r2_score(y_test, model.predict(x_test))

def regressor_ii():
    from keras.models import Sequential
    model = Sequential()
    model.add(Dense(units=2, input_dim=203))
    model.add(LeakyReLU())
    model.add(Dense(units=16))
    model.add(LeakyReLU())
    model.add(Dense(units=2))
    return model

def regressorMLP(layers, activator, bn=True, dr=0.25):
    inputs = Input(shape=(203,), name="input_layer")
    reg = regularizers.l2(0.00001)
    rand_init = RandomNormal(stddev=0.02)

    if activator == PReLU:
        x = PReLU(shared_axes=[1])(inputs)
    else:
        x = activator()(inputs)

    for layer in layers:
        x = Dense(units=layer, activation=None, kernel_initializer=rand_init,
                  kernel_regularizer=reg, bias_regularizer=reg)(x)
        if bn:
            x = BatchNormalization(axis=-1)(x)
        if activator == PReLU:
            x = PReLU(shared_axes=[1])(x)
        else:
            x = activator()(x)

    if dr != 0.0:
        x = Dropout(dr)(x)
    prediction = Dense(units=2, activation=None,
                       kernel_initializer=RandomNormal(stddev=0.02),
                       kernel_regularizer=reg, bias_regularizer=reg, name="output_layer")(x)
    return Model(inputs=[inputs], outputs=[prediction])

def regressor_PClowd(layers_pts, layers_prms, aoa_dims=1, set_vars=3, max_points=100, pts_dim=2, commonliftdrag=False, pclowdBN=False, paramBN=True, opt=False):
    if opt:
        model = model_from_json(open('keras_model_opt2_{0}.json'.format(shapeType)).read())
        return model
    else:
        reg = regularizers.l2(0.00001)
        rand_init = RandomNormal(stddev = 0.02)
        def dense_block(num_points, layers, bn, lastShape = None):
            inputs = Input(shape = (num_points,), name = "input_layer")
            x = PReLU(shared_axes = [1])(inputs)
        
            for layer in layers:
                x = Dense(units = layer, activation = None, kernel_initializer = rand_init,
                          kernel_regularizer = reg, bias_regularizer = reg)(x)
                if bn:
                    x = BatchNormalization(axis = -1)(x)
                x = PReLU(shared_axes = [1])(x)
            
            if lastShape is None:
                outputs = Activation("linear")(x)
            else:
                outputs = Dense(units = lastShape, activation = None, kernel_initializer = rand_init,
                                kernel_regularizer = reg, bias_regularizer = reg)(x)
            return Model([inputs], [outputs])
            
        def sets_block(num_points, layers, bn):
            input_pts = Input(shape=(num_points,))
            pts = dense_block(num_points, layers, bn)(input_pts)
            Adder = Lambda(lambda x: K.sum(x, axis = 1), name="Adder", output_shape = (1,))
            pts = Adder(pts)
            output = Reshape((1,))(pts)
            return Model(inputs = [input_pts], outputs = [output])
    
        def aero_block(lastShape=1):
            input_pts = Input(shape = (max_points * pts_dim,), name = "input_set")
            input_params = Input(shape = (2 + aoa_dims,), name = "input_params")
        
            x = sets_block(max_points*pts_dim, layers_pts, bn=False)(input_pts) # Sxy
            
            if set_vars==3:
                x_pts = Lambda(lambda x: x[:, :max_points], name = "x_pts", output_shape = (max_points,))(input_pts)    # sx
                y_pts = Lambda(lambda x: x[:, max_points:], name = "y_pts", output_shape = (max_points,))(input_pts)    # sy
                x_pts = sets_block(max_points, layers_pts, bn = pclowdBN)(x_pts)
                y_pts = sets_block(max_points, layers_pts, bn = pclowdBN)(y_pts)
                x = Concatenate(name="HiddenConcat")([x_pts, y_pts, x, input_params])
            else:
                x = Concatenate()([x, input_params])
    
            aero = dense_block(set_vars+aoa_dims+2, layers_prms, bn=paramBN, lastShape = lastShape)(x)
            return Model([input_pts, input_params], [aero])
    
        input_pts = Input(shape = (max_points * pts_dim,), name = "input_set")
        input_params = Input(shape = (2 + aoa_dims,), name = "input_params")
        if not commonliftdrag:
            lift = aero_block()([input_pts, input_params])
            drag = aero_block()([input_pts, input_params])
            predicted = Concatenate()([lift, drag])
        else:
            predicted = aero_block(lastShape=2)([input_pts, input_params])
        
        return Model(inputs = [input_pts, input_params], outputs = [predicted])

def regressor_cnn(resnet=True):
    __s = Input(shape = img_shape, name = "input_img")
    __l = Input(shape = (n_labeles,), name = "input_label")
    reg = regularizers.l2(0.00001)
    rand_init = RandomNormal(stddev = 0.02)
    from regressor.res_net_builder import ResnetBuilder
    if resnet:
        model = ResnetBuilder.build_resnet_18(img_shape, num_outputs = 201)
        __h = model(__s)
    else:
        # 1st conv
        __h = Conv2D(filters = 64, kernel_size = 3, strides = 2, padding = 'same',
                     activation = None, kernel_initializer = rand_init,
                     kernel_regularizer = reg, bias_regularizer = reg)(__s)
        __h = PReLU(shared_axes = [1, 2, 3])(__h)
        # 2nd conv
        __h = Conv2D(filters = 128, kernel_size = 3, strides = 2, padding = 'same',
                     activation = None, kernel_initializer = rand_init,
                     kernel_regularizer = reg, bias_regularizer = reg)(__h)
        __h = BatchNormalization(axis = -1)(__h)
        __h = PReLU(shared_axes = [1, 2, 3])(__h)
        # 3nd conv
        __h = Conv2D(filters = 256, kernel_size = 4, strides = 1, padding = 'valid',
                     activation = None, kernel_initializer = rand_init,
                     kernel_regularizer = reg, bias_regularizer = reg)(__h)
        __h = BatchNormalization(axis = -1)(__h)
        __h = PReLU(shared_axes = [1, 2, 3])(__h)
        # 4th conv
        __h = Conv2D(filters = 512, kernel_size = 3, strides = 2, padding = 'same',
                     activation = None, kernel_initializer = rand_init,
                     kernel_regularizer = reg, bias_regularizer = reg)(__h)
        __h = BatchNormalization(axis = -1)(__h)
        __h = PReLU(shared_axes = [1, 2, 3])(__h)
        __h = Dropout(0.25)(__h)
        # 5th fc
        __h = Flatten()(__h)

    __h = BatchNormalization(axis = -1)(__h)
    __h = Concatenate()([__h, __l])
    if resnet:
        model2 = regressor(droprate=0.125)
        __y = model2(__h)
    else:
        __h = Dense(units = 256, activation = None, kernel_initializer = RandomNormal(stddev = 0.02),
                    kernel_regularizer = reg, bias_regularizer = reg)(__h)
        # __h = Activation("relu")(__h)
        __h = PReLU(shared_axes=[1])(__h)
        __h = Dropout(0.5)(__h)
        __y = Dense(units = n_classes, activation = None, kernel_initializer = RandomNormal(stddev = 0.02),
                    kernel_regularizer = reg, bias_regularizer = reg)(__h)

    return Model([__s, __l], __y, name = "regressor")


def getNewestModel(model, dirname, name=""):
    from glob import glob
    target = os.path.join(dirname, '*{0}*'.format(name))
    files = [(f, os.path.getmtime(f)) for f in glob(target)]
    if len(files) == 0:
        return model
    else:
        newestModel = sorted(files, key=lambda files: files[1])[-1]
        model.load_weights(newestModel[0])
        return model

def relationalError(y_true, y_pred, percentage=True, absoluteValue=False):
    def check_dimension(y_true, y_pred):
        if y_true.shape != y_pred.shape:
            raise ValueError
        return
    check_dimension(y_true, y_pred)
    n_feat = y_true.shape[1]
    input_true = Input(shape = (n_feat,))
    input_pred = Input(shape = (n_feat,))
    # error = Subtract()([input_pred, input_true])
    # outputs = Lambda(lambda x: x[0]/x[1])([error, input_true])
    if absoluteValue:
        outputs = Lambda(lambda x: (x[1] - x[0]))([input_true, input_pred])
    else:
        outputs = Lambda(lambda x: (x[1] - x[0]) / x[0])([input_true, input_pred])
    model = Model(inputs =[input_true, input_pred], outputs=outputs)
    relationalError = model.predict([y_true, y_pred])
    if percentage and (not absoluteValue):
        return relationalError * 100
    else:
        return relationalError

def debug(x_train, x_test, y_train, y_test):
    def exchange(x_train, x_test, y_train, y_test, idx):
        xMax = x_test[idx]
        yMax = y_test[idx]
        x_test = np.concatenate([x_test[:idx], x_test[idx+1:]])
        y_test = np.concatenate([y_test[:idx], y_test[idx+1:]])
        
        x_train = np.concatenate([x_train, xMax.reshape(1, -1)])
        y_train = np.concatenate([y_train, yMax.reshape(1, -1)])
        return x_train, x_test, y_train, y_test
    
    if np.max(y_train[:, 0]) < np.max(y_test[:, 0]):
        idx = np.argmax(y_test[:, 0])
        x_train, x_test, y_train, y_test = exchange(x_train, x_test, y_train, y_test, idx)
    
    if np.max(y_train[:, 1]) < np.max(y_test[:, 1]):
        idx = np.argmax(y_test[:, 1])
        x_train, x_test, y_train, y_test = exchange(x_train, x_test, y_train, y_test, idx)
    
    if np.min(y_train[:, 0]) > np.min(y_test[:, 0]):
        idx = np.argmin(y_test[:, 0])
        x_train, x_test, y_train, y_test = exchange(x_train, x_test, y_train, y_test, idx)

    if np.min(y_train[:, 1]) > np.min(y_test[:, 1]):
        idx = np.argmin(y_test[:, 1])
        x_train, x_test, y_train, y_test = exchange(x_train, x_test, y_train, y_test, idx)
    
    return x_train, x_test, y_train, y_test

def split_data(x_data, rescale_x = None, updateRescaler = False):
    dataNum = x_data.shape[0]
    if rescale_x is None:
        objS = x_data[:, 1:201]
        return [objS, np.concatenate([x_data[:, 0].reshape(-1, 1), x_data[:, 201:]], axis = 1)]
    else:
        from sklearn.preprocessing import StandardScaler
        x_data = rescale_x(x_data)
        s_data = x_data[:, 1:201].reshape(dataNum, 2, -1)
        a_data = x_data[:, 0] * np.pi / 180.0
        rot = np.array([[np.cos(a_data), -np.sin(a_data)], [np.sin(a_data), np.cos(a_data)]])
        
        s_data = np.einsum("ijk,jli->ilk", s_data, rot).reshape(dataNum, -1)
        scalar = StandardScaler()
        x_data = np.concatenate([s_data, x_data[:, 201:]], axis = 1)
        x_data = scalar.fit_transform(x_data)
        if updateRescaler:
            return [x_data[:, :200], x_data[:, 200:]], scalar.inverse_transform
        else:
            return [x_data[:, :200], x_data[:, 200:]], rescale_x
        
def nnr_main(reductTarget, shapeType=shapeType, iterator=None):
    config = tf.ConfigProto(gpu_options = tf.GPUOptions(allow_growth = True))
    session = tf.Session(config = config)
    KTF.set_session(session)
    
    if not "msdf" in shapeType.lower():
        x_train, x_test, y_train, y_test, rescale_x, rescale_y3, p_train, p_test = load_mixed(mode=shapeType, n_samples = 50000, test_size=0.01, reductTarget = reductTarget, split45 = False, with_param = True)
        # x_train, x_test, y_train, y_test = debug(x_train, x_test, y_train, y_test)

        # scalar = MinMaxScaler()
        # scalar.fit(y_train)
        rescale_y = rescale_y3
        # rescale_y2 = lambda y_data: rescale_y3(scalar.inverse_transform(y_data))
        # y_train = scalar.transform(y_train)
        # y_test = scalar.transform(y_test)
        # rescale_y = rescale_y2
        if "pclowd" in shapeType.lower():
            if "noaoa" in shapeType.lower():
                inputs_test, rescale_x = split_data(x_test, rescale_x, updateRescaler = False)
                inputs_train, rescale_x = split_data(x_train, rescale_x, updateRescaler = True)
            else:
                inputs_test = split_data(x_test)
                inputs_train = split_data(x_train)
            
        else:
            inputs_train = x_train
            inputs_test = x_test

    else:
        x_img_train, x_img_test, x_label_train, x_label_test, y_train, y_test, rescale_x_img, rescale_x_label, rescale_y = load_mixed(
            mode = shapeType, n_samples = 50000, msdf_res = msdf_res, reductTarget="")
        inputs_train = [x_img_train, x_label_train]
        inputs_test = [x_img_test, x_label_test]

    name = None
    optModel = True
    load_path_json = modelDir / Path(fnameGen("model.json", r2=r2_rec, num=0))
    load_path_weight = modelDir / Path(fnameGen("weight.hdf5", r2=r2_rec, num=0))
    if load_path_json.exists() and load_path_weight.exists():
        model = model_from_json(open(load_path_json).read())
        model.load_weights(load_path_weight)
        print("learned model loaded")
    else:
        if not "msdf" in shapeType.lower():
            if not "test" in shapeType.lower():
                if not "pclowd" in shapeType.lower():
                    model = regressor(opt=optModel)

                else:
                    layers_pts = [101,101,151]#[25,151]
                    layers_prms = [51, 151]
                    if "noaoa" in shapeType.lower():
                        aoa_dims = 0
                    else:
                        aoa_dims = 1
                    
                    model = regressor_PClowd(layers_pts, layers_prms, aoa_dims=aoa_dims, set_vars=3)
            else:
                if iterator is None:
                    layers_pts = [25, 151]
                    layers_prms = [51, 151]

                if "noaoa" in shapeType.lower():
                    aoa_dims = 0
                else:
                    aoa_dims = 1
                if "incomp_model" in shapeType.lower():
                    model = regressor_ii()
                elif "pclowd" in shapeType.lower():
                    if "normal" in shapeType.lower():
                        model = regressor_PClowd(layers_pts, layers_prms, aoa_dims=aoa_dims, set_vars=3, commonliftdrag=False)
                    elif "common" in shapeType.lower():
                        model = regressor_PClowd(layers_pts, layers_prms, aoa_dims=aoa_dims, set_vars=3, commonliftdrag=True)
                    elif "s1only" in shapeType.lower():
                        model = regressor_PClowd(layers_pts, layers_prms, aoa_dims=aoa_dims, set_vars=1)
                    elif "c1only" in shapeType.lower():
                        model = regressor_PClowd(layers_pts, layers_prms, aoa_dims = aoa_dims, set_vars = 1, commonliftdrag = True)
                    elif "pbn" in shapeType.lower():
                        model = regressor_PClowd(layers_pts, layers_prms, aoa_dims=aoa_dims, set_vars=3, pclowdBN=True)
                    elif "gspclowd" in shapeType.lower():
                        layer1, layer2, activator, bn1, bn2, dr, name = iterator.__next__()
                        print(name)
                        model = regressor_PClowd(layer1, layer2, aoa_dims=aoa_dims, set_vars=3, commonliftdrag=False, pclowdBN=bn1, paramBN=bn2)
                    else:
                        raise ValueError
                elif "nodrop" in shapeType.lower():
                    model = regressor(droprate=0.0)
                elif "drop0125" in shapeType.lower():
                    model = regressor(droprate=0.125)
                elif "drop0375" in shapeType.lower():
                    model = regressor(droprate=0.375)
                elif "drop0500" in shapeType.lower():
                    model = regressor(droprate=0.5)
                elif "gsdense" in shapeType.lower():
                    if iterator is None:
                        raise ValueError
                    layer, activator, bn, dr, name = iterator.__next__()
                    print(name)
                    model = regressorMLP(layer, activator, bn, dr)
                else:
                    raise ValueError
        else:
            model = regressor_cnn()
    
        model.compile(loss=mean_squared_error,
                      optimizer=Adam(lr = 0.0002, beta_1 = 0.5),
                      metrics=["mae"])

        chkptDir = Path("checkpoint2")
        cname = fnameGen("dMLP_.tmp", name = name).name
        chkpt = chkptDir / Path(cname + '.{epoch:02d}-{val_loss:.2f}.hdf5')

        # es_cb = EarlyStopping(monitor='val_loss', patience=500, verbose=0, mode='auto')
        cp_cb = ModelCheckpoint(filepath=str(chkpt), monitor="val_loss", verbose=0, save_best_only=True, mode='auto')

        model.summary()

        history = model.fit(x=inputs_train,
                            y=y_train,
                            verbose=1,
                            batch_size=batch_size,
                            epochs=epoch,
                            validation_split=0.1,
                            callbacks = [cp_cb]
                            )
        model = getNewestModel(model, str(chkptDir), name=cname)

        hidden_save = False
        if hidden_save:
            if not "c1only" in shapeType and not "common" in shapeType:
                if "s1only" in shapeType:
                    hidden_layer_model_cl = Model(inputs=model.get_layer("model_8").get_input_at(0), outputs=model.get_layer("model_8").get_layer("HiddenConcat").output)
                    hidden_layer_model_cd = Model(inputs=model.get_layer("model_16").get_input_at(0), outputs=model.get_layer("model_16").get_layer("HiddenConcat").output)
                else:
                    hidden_layer_model_cl = Model(inputs = model.get_layer("model_4").get_input_at(0), outputs = model.get_layer("model_4").get_layer("concatenate_1").output)
                    hidden_layer_model_cd = Model(inputs = model.get_layer("model_8").get_input_at(0), outputs = model.get_layer("model_8").get_layer("concatenate_2").output)
                    
                hidden_output_cl = hidden_layer_model_cl.predict(inputs_test)
                hidden_output_cd = hidden_layer_model_cd.predict(inputs_test)
            
                np.savez("hOutCV3.npz", hidden_output_cl = hidden_output_cl, hidden_output_cd = hidden_output_cd, inputs_train_0 = inputs_train[0], inputs_train_1=inputs_train[1],
                                inputs_test_0 = inputs_test[0], inputs_test_1=inputs_test[1], y_train = y_train, y_test = y_test, p_train = p_train, p_test = p_test)
            else:
                if "common" in shapeType:
                    hidden_layer_model = Model(inputs = model.get_layer("model_8").get_input_at(0),
                                                  outputs = model.get_layer("model_8").get_layer("HiddenConcat").output)
                else:
                    hidden_layer_model = Model(inputs=model.get_layer("model_4").get_input_at(0), outputs=model.get_layer("model_4").get_layer("concatenate_1").output)
                hidden_layer_output = hidden_layer_model.predict(inputs_test)
                np.savez("hOutCV3.npz", hidden_output = hidden_layer_output,
                         inputs_train_0 = inputs_train[0], inputs_train_1 = inputs_train[1],
                         inputs_test_0 = inputs_test[0], inputs_test_1 = inputs_test[1], y_train = y_train,
                         y_test = y_test, p_train = p_train, p_test = p_test)

        path_history = figsDir / fnameGen(title="history.png", r2="", name=name)
        # history plot
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(history.history['loss'], marker=".", label='loss')
        ax.plot(history.history['val_loss'], marker='.', label='val_loss')
        ax.set_title('model loss history')
        ax.grid(True)
        ax.set_xlabel('epoch')
        ax.set_ylabel('loss')
        ax.legend(loc='best')
        plt.savefig(path_history)
        plt.close()

    model.summary()
    y_pred_t = rescale_y(model.predict(inputs_train))
    y_pred = rescale_y(model.predict(inputs_test))
    y_train = rescale_y(y_train)
    y_test = rescale_y(y_test)

    r2_train = r2_score(y_train, y_pred_t)
    r2_pred = r2_score(y_test, y_pred)
    r2_train_cl = r2_score(y_train[:, 0], y_pred_t[:, 0])
    r2_pred_cl = r2_score(y_test[:, 0], y_pred[:, 0])
    r2_train_cd = r2_score(y_train[:, 1], y_pred_t[:, 1])
    r2_pred_cd = r2_score(y_test[:, 1], y_pred[:, 1])
    
    rError = relationalError(y_true = y_test, y_pred = y_pred, absoluteValue = True)
    print(relationalError(y_true = y_test, y_pred = y_pred, absoluteValue = True))
    np.save(str(fnameGen(title = "relationalError.npz", name=name)), rError)
    print(rError.shape)
    cle = rError[:, 0]
    cde = rError[:, 1]

    print(r2_train)
    print(r2_pred)
    print(r2_train_cl, r2_pred_cl, np.mean(abs(cle)), np.median(np.abs(cle)))
    print(r2_train_cd, r2_pred_cd, np.mean(abs(cde)), np.median(np.abs(cde)))
    # """
    plt.rcParams["font.size"] = 18
    fig = plt.figure()
    ax = fig.add_subplot(111)
    from matplotlib.colors import LogNorm
    H = ax.hist2d(cle, cde, bins=[np.linspace(-0.5,0.5,51),np.linspace(-0.5,0.5,51)],
                  norm=LogNorm())
    ax.set_title("Error 2D-Histgram :{0}".format(shapeType.replace("Concat", "")), fontsize=22)
    ax.set_xlabel("$C_L$ Error")
    ax.set_ylabel("$C_D$ Error")
    
    fig.colorbar(H[3], ax=ax)
    # plt.show()
    plt.savefig(str(figsDir / fnameGen("errHist2D.png", r2 = "", name=name)))
    plt.close()
    # check(cle)
    # check(cde)
    # """
    r2_pred=""
    path_json = modelDir / fnameGen(title="model.json", r2 = r2_pred, name=name)
    path_weight = modelDir / fnameGen(title="weight.hdf5", r2 = r2_pred, name=name)

    if not load_model:
        with open(path_json, 'w') as fout:
            model_json_str = model.to_json(indent=4)
            fout.write(model_json_str)
            ## model weights
        model.save_weights(path_weight)
        if not "incomp_model" in shapeType.lower():
            try:
                plot_model(model, to_file=modelDir / fnameGen(title = "arch.png", r2 = r2_pred, name=name))
            except:
                pass

    
    make_scatter_plot(data_a = y_train[:, 0], data_b = y_pred_t[:, 0], label_a = "observed $C_L$ (train data)",
                      label_b = "predicted $C_L$",
                      fname = str(figsDir / fnameGen("yy_train_cl.png", r2 = r2_pred, name=name)), dash_line = "no", max_value=2.0, min_value=-1.0,
                      title = "yyplot {0} $C_L$ (train)".format(shapeType.replace("Concat", "")))
    make_scatter_plot(data_a = y_train[:, 1], data_b = y_pred_t[:, 1], label_a = "observed $C_D$ (train data)",
                      label_b = "predicted $C_D$",
                      fname = str(figsDir / fnameGen("yy_train_cd.png", r2 = r2_pred, name=name)), dash_line = "no", max_value=2.0, min_value=0.0,
                      title = "yyplot {0} $C_D$ (train)".format(shapeType.replace("Concat", "")))
    make_scatter_plot(data_a = y_test[:, 0], data_b = y_pred[:, 0], label_a = "observed $C_L$ (test data)",
                      label_b = "predicted $C_L$",
                      fname = str(figsDir / fnameGen("yy_test_cl.png", r2 = r2_pred, name=name)), dash_line = "no", max_value=2.0, min_value=-1.0,
                      title = "yyplot {0} $C_L$ (test)".format(shapeType.replace("Concat", "")))
    make_scatter_plot(data_a = y_test[:, 1], data_b = y_pred[:, 1], label_a = "observed $C_D$ (test data)",
                      label_b = "predicted $C_D$",
                      fname = str(figsDir / fnameGen("yy_test_cd.png", r2 = r2_pred, name=name)), dash_line = "no", max_value=2.0, min_value=0.0,
                      title = "yyplot {0} $C_D$ (test)".format(shapeType.replace("Concat", "")))
    # log0 = str(history.history)

    log = str(r2_train) + " " + str(r2_pred) + "\n"
    print("NNR:" + log)
    # log += str(log0) + "\n\n"
    r2_train_cl = r2_score(y_train[:, 0], y_pred_t[:, 0])
    r2_pred_cl = r2_score(y_test[:, 0], y_pred[:, 0])
    r2_train_cd = r2_score(y_train[:, 1], y_pred_t[:, 1])
    r2_pred_cd = r2_score(y_test[:, 1], y_pred[:, 1])
    if not "msdf" in shapeType.lower():
        with open("log_new_test651.csv", "a") as f:
            f.write("{0},{1},{2},{3},{4},mix\n".format(x_train.shape[0], r2_train, r2_pred, reductTarget, fnameGen("log.csv", name=name)))
            f.write("{0},{1},{2},{3},{4},{5},cl\n".format(x_train.shape[0], r2_train_cl, r2_pred_cl, reductTarget, np.mean(np.abs(cle)),fnameGen("log.csv", name=name)))
            f.write("{0},{1},{2},{3},{4},{5},cd\n".format(x_train.shape[0], r2_train_cd, r2_pred_cd, reductTarget, np.mean(np.abs(cde)), fnameGen("log.csv", name=name)))

    log2 = error_print(y_pred[:, 0], y_test[:, 0], csv=True).replace("\n", ",")
    log2 += error_print(y_pred[:, 1], y_test[:, 1], csv=True).replace("\n", ",")
    log2 += str(fnameGen("", noExt=True, name=name)).replace("_", ",") + "\n"

    # with open("newLog634{0}.csv".format(shapeType), "a") as f:
    with open("newLog651.csv", "a") as f:
        log2.replace("\n", ",") + shapeType + ","
        f.write(log2)
    print(log2)
    return log2

from itertools import product
def permutationClowd():
    activators = [PReLU]
    # bns = [True, False]
    bns1 = [False]#,True]
    bns2 = [True]
    drs = [0.0]#125, 0.25, 0.375, 0.5]
    layer_num = [3]  # [2, 3, 4]
    units_num = [25, 51, 101, 151]
    layers = []
    for layer in layer_num:
        layers.extend(list(product(units_num, repeat = layer)))
    layers2 = [(51,151)]
    # layers = [(101,151,101)]
    yield len(list(product(layers, layers2, activators, bns1, bns2, drs)))
    flag = False
    # target = [(101, 101, 51), (51, 151), PReLU, False, True, 0.0]
    target = None
    for layer1, layer2, activator, bn1, bn2, dr in product(layers, layers2, activators, bns1, bns2, drs):
        name = "{0}_{1}_{2}_{3}_{4}_{5}".format(layer1, layer2, activator, bn1, bn2, dr).replace("<", "").replace(">", "").replace(".",
                                                                                                            "").replace(
            "(", "").replace(")", "").replace("'", "")
        if target is not None:
            if flag:
                yield layer1, layer2, activator, bn1, bn2, dr, name
            if [layer1, layer2, activator, bn1, bn2, dr] == target:
                flag = True
        else:
            yield [layer1, layer2, activator, bn1, bn2, dr, name]

def permutationDense():
    activators = [PReLU]
    bns = [True, False]
    drs = [0.0, 0.125, 0.25, 0.375]
    layer_num = [3]  # [2, 3, 4]
    units_num = [25, 51, 101, 151]
    layers = []
    for layer in layer_num:
        layers.extend(list(product(units_num, repeat=layer)))

    yield len(list(product(layers, activators, bns, drs)))
    # target = None
    target = [(151, 101, 51), PReLU, False, 0.125]
    flag = False
    for layer, activator, bn, dr in product(layers, activators, bns, drs):
        name = "{0}_{1}_{2}_{3}".format(layer, activator, bn, dr).replace("<","").replace(">","").replace(".","").replace("(","").replace(")","").replace("'","")
        if target is not None:
            if flag:
                yield layer, activator, bn, dr, name
            if [layer, activator, bn, dr] == target:
                flag = True
        else:
            yield layer, activator, bn, dr, name

def optimize_network():
    study_name = "study2_{0}".format(shapeType)
    study = optuna.create_study(
        study_name=study_name,
        storage = "sqlite:///{0}.db".format(study_name),
        load_if_exists=True,
        direction = "minimize",
        pruner = optuna.pruners.MedianPruner()
    )
    study.optimize(objective, n_trials = 1000)
    pruned_trials = [t for t in study.trials if t.state == optuna.structs.TrialState.PRUNED]
    complete_trials = [t for t in study.trials if t.state == optuna.structs.TrialState.COMPLETE]
    print("Study statistics: ")
    print("  Number of finished trials: ", len(study.trials))
    print("  Number of pruned trials: ", len(pruned_trials))
    print("  Number of complete trials: ", len(complete_trials))

    # 結果の表示
    print("Best trial:")
    trial = study.best_trial
    print("  Value: ", trial.value)
    print("  Params: ")
    for key, value in trial.params.items():
        print(f"    {key}: {value}")

    
if __name__ == '__main__':
    """check
    shapeType = "CircumConcat"
    model = regressor(opt=True)
    model.summary()
    exit()
    """
    np.random.seed(1)
    tf.random.set_random_seed(1)
    """OPTUNA
    BATCHSIZE = 5000
    EPOCHS = 150
    shapeTypes = ["CircumConcat_pclowd_noaoa"]#["CrowdConcat", "CircumConcat"]
    for shapeType in shapeTypes:
        reductTarget = "kmeans_norm_pca_9000"
        x_train, x_test, y_train, y_test, rescale_x, rescale_y3, p_train, p_test = load_mixed(mode = shapeType,
                                                                                              n_samples = 50000,
                                                                                              test_size = 0.01,
                                                                                              reductTarget = reductTarget,
                                                                                              split45 = False,
                                                                                              with_param = True)
        x_test, rescale_x = split_data(x_test, rescale_x, updateRescaler = False)
        x_train, rescale_x = split_data(x_train, rescale_x, updateRescaler = True)
        optimize_network()
    exit()
    #"""
    """permutation
    head = "kmeans_norm_pca_"
    shapeTypeTail = "gspclowd_noaoa"#"gsdense"#
    shapeTypes = ["test_CircumConcat"]
    reductTarget = head + str(5000)
    # iterator = permutationDense()
    iterator = permutationClowd()
    stop = iterator.__next__()
    print(stop)
    for i in range(stop):
        shapeType = shapeTypes[0] + shapeTypeTail
        print(shapeType)
        nnr_main(reductTarget, shapeType, iterator=iterator)
    exit()
    #"""
    """single
    head = "kmeans_norm_pca_"
    # shapeTypeHead = ["test_CircumConcat"]#PClowdNoAoA"]#["CrowdConcatPClowd", "EquidistantConcatPClowd"]#, "CrowdConcat", "EquidistantConcat"]#["FourierConcat"]#["CrowdConcat", "EquidistantConcat"]#["FourierConcat", "CrowdConcat", "EquidistantConcat"]
    # ["incomp_model", "pclowd_common", "pclowd_s1only", "pclowd_pbn", "nodrop", ]
    shapeTypeTail = ""#"pclowds_c1only_noaoa"
    shapeTypes = ["MSDFResNew4"]#["CircumConcat"]  # "drop0125", "drop0375", "drop0500"]
    # reductTarget = head + str(5000)
    # iterator = None
    #
    for res in range(1):
        msdf_res = 2**(res + 4)  # pixel
        img_shape = (msdf_res, msdf_res, 3)
        shapeType = shapeTypes[0] + str(msdf_res)  # "gspclowd_noaoa"#"gsdense"#
        batch_size = int(3000)
        print(shapeType)
        nnr_main(reductTarget, shapeType)
        exit()
    exit()
    #"""
    #"""Continuous Test
    """
    shapeTypes = ["FourierConcat", "EquidistantConcat"]  # ["test_CircumConcat"]#"drop0125", "drop0375", "drop0500"]
    reductTarget = head + str(5000)
    iterator = None
    for shapeType in shapeTypes:
        print(shapeType)
        nnr_main(reductTarget, shapeType, iterator = iterator)
    """
    shapeType = "CircumConcatPClowdNoAoA"#"FourierConcat"
    heads = ["kmeans_norm_pca_"]
    for head in heads:
        # for i in range(54, 235):
        for i in range(1):
            
            print(reductTarget)
            reductTarget = head + str(5000)
            nnr_main(reductTarget, shapeType)
    #"""