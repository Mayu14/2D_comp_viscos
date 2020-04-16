# coding: utf-8
from glob import glob
import os
from math import ceil

import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import tensorflow as tf

from keras.models import Model, model_from_json
from keras.layers import Input, Dense, Conv2D, Conv2DTranspose, Flatten, Reshape, Concatenate, Dropout, Lambda, Subtract
from keras.layers import BatchNormalization, PReLU, Activation, LeakyReLU
from keras.callbacks import EarlyStopping, ModelCheckpoint, TensorBoard
from keras.optimizers import Adadelta, Adam
from keras.losses import categorical_crossentropy, mean_squared_error, mean_absolute_error
from keras.utils import plot_model
from keras.initializers import RandomNormal
from keras import regularizers
from keras.backend import tensorflow_backend as KTF

from sklearn.metrics import r2_score

shapeType = "FourierConcat"
reductTarget = ""

class AutoEncoder(object):
    batch_size = 15000
    epoch = 1000
    l2_penalty = 0.00001
    rand_stddev = 0.02
    baseSaveDir = Path("G:\\Toyota\\Data\\Compressible_Viscos\\rhoSimpleFoam\\training_data\\learned\\" + shapeType.replace("Concat",
                                                                                                             ""))
    modelDir = baseSaveDir / Path("model")
    
    def __init__(self, input_dim, mid_dim, layer_num=0, mode="Fourier", deep=True, x_train=None, x_test=None, train_now=True):
        self.input_dim = input_dim
        self.mid_dim = mid_dim
        self.layer_num = layer_num
        # self.deepNet = deep
        self.mode = mode
        self.x_train = x_train
        self.x_test = x_test
        
        self.name = "m{3}_i{0}_m{1}_d{2}".format(str(self.input_dim[0]), str(self.mid_dim[0]), str(self.layer_num), self.mode.lower())
        self.name_ae = "ae_" + self.name
        self.name_encoder = "encoder_" + self.name
        self.name_decoder = "decoder_" + self.name
        model_tail = "model.json"
        weight_tail = "weight.json"
        self.load_path_json_ae = self.modelDir / Path(self.name_ae + model_tail)
        self.load_path_weight_ae = self.modelDir / Path(self.name_ae + weight_tail)
        self.load_path_json_encoder = self.modelDir / Path(self.name_encoder + model_tail)
        self.load_path_weight_encoder = self.modelDir / Path(self.name_encoder + weight_tail)
        self.load_path_json_decoder = self.modelDir / Path(self.name_decoder + model_tail)
        self.load_path_weight_decoder = self.modelDir / Path(self.name_decoder + weight_tail)

        if not self.load_model():
            print("start new AE training")
            self.buildALL()
            if train_now:
                self.train()
        else:
            print("AE load completed.")

    def build_encoder(self):
        units_diff = ceil((self.input_dim[0] - self.mid_dim[0]) / (self.layer_num + 1))
        
        reg = regularizers.l2(self.l2_penalty)
        rand_init = RandomNormal(stddev = self.rand_stddev)
        
        inputs = Input(shape = self.input_dim, name = "input_of_encoder")
        __x = Lambda(lambda x: x)(inputs)
        units = self.input_dim[0]
        for i in range(self.layer_num + 1):
            units = max(units - units_diff, self.mid_dim[0])
            __x = PReLU(shared_axes = [1])(__x)
            __x = Dense(units = units, activation = None, kernel_initializer = rand_init,
                        kernel_regularizer = reg, bias_regularizer = reg)(__x)
            __x = BatchNormalization(axis = -1)(__x)
            
        __x = PReLU(shared_axes = [1])(__x)
        outputs = Dense(units = self.mid_dim[0], activation = None, kernel_initializer = rand_init,
                        kernel_regularizer = reg, bias_regularizer = reg,
                        name = "output_of_encoder")(__x)
        return Model([inputs], [outputs], name = "encoder")
        

    def build_decoder(self):
        units_diff = ceil((self.input_dim[0] - self.mid_dim[0]) / (self.layer_num + 1))
        
        reg = regularizers.l2(self.l2_penalty)
        rand_init = RandomNormal(stddev = self.rand_stddev)
    
        inputs = Input(shape = self.mid_dim, name = "input_of_decoder")
        __x = Lambda(lambda x: x)(inputs)
        units = self.mid_dim[0]
        for i in range(self.layer_num + 1):
            units = min(units + units_diff, self.input_dim[0])
            __x = PReLU(shared_axes = [1])(__x)
            __x = Dense(units = units, activation = None, kernel_initializer = rand_init,
                        kernel_regularizer = reg, bias_regularizer = reg)(__x)
            __x = BatchNormalization(axis = -1)(__x)
            
        __x = PReLU(shared_axes = [1])(__x)
        outputs = Dense(units = self.input_dim[0], activation = None, kernel_initializer = rand_init,
                        kernel_regularizer = reg, bias_regularizer = reg,
                        name = "output_of_decoder")(__x)
        return Model([inputs], [outputs], name = "decoder")
    
    def buildALL(self, returnModel=False):
        x_origin = Input(shape = self.input_dim, name = "original")
        
        encoder = self.build_encoder()
        x_encoded = encoder(x_origin)
        
        
        decoder = self.build_decoder()
        x_decoded = decoder(x_encoded)
        
        auto_encoder = Model([x_origin], [x_decoded], name="auto_encoder")
        auto_encoder_optimizer = Adam(lr = 0.0002, beta_1 = 0.5)
        auto_encoder.compile(optimizer = auto_encoder_optimizer,
                             loss = mean_squared_error,
                             metrics = ["mae"])

        if not returnModel:
            self.auto_encoder, self.encoder, self.decoder = auto_encoder, encoder, decoder
        else:
            return auto_encoder, encoder, decoder

    def load_model(self, x_train=None, threshold=0.9):
        flag = False
        if (self.load_path_json_ae.exists() and self.load_path_weight_ae.exists()) and (
                self.load_path_json_encoder.exists() and self.load_path_weight_encoder.exists()) and (
                self.load_path_weight_decoder.exists() and self.load_path_json_decoder.exists()):
            self.auto_encoder = model_from_json(open(self.load_path_json_ae).read())
            self.auto_encoder.load_weights(self.load_path_weight_ae)

            self.encoder = model_from_json(open(self.load_path_json_encoder).read())
            self.encoder.load_weights(self.load_path_weight_encoder)

            self.decoder = model_from_json(open(self.load_path_json_decoder).read())
            self.decoder.load_weights(self.load_path_weight_decoder)

            x_train = self.check_x_input(x_train)
            r2 = r2_score(x_train, self.auto_encoder.predict(x_train))
            if r2 > threshold:
                flag = True # load completely
        return flag

    def check_x_input(self, x_train):
        flag = False
        use_input = False
        if not (self.x_train is None):
            flag = True
        if not (x_train is None):
            flag = True
            use_input = True

        if not flag:
            raise ValueError

        if not use_input:
            return self.x_train
        else:
            return x_train

    def train(self, x_train=None):
        def getNewestModel(model, dirname):
            target = os.path.join(dirname, '*')
            files = [(f, os.path.getmtime(f)) for f in glob(target)]
            if len(files) == 0:
                return model
            else:
                newestModel = sorted(files, key = lambda files: files[1])[-1]
                model.load_weights(newestModel[0])
                return model

        x_train = self.check_x_input(x_train)

        self.chkptDir = Path("checkpoint")
        if not self.chkptDir.exists():
            self.chkptDir.mkdir()
        self.chkpt = self.chkptDir / Path('AE_' + str(self.mid_dim[0]) + '.{epoch:02d}-{val_loss:.2f}.hdf5')
        self.cp_cb = ModelCheckpoint(filepath=str(self.chkpt), monitor="val_loss", verbose=0, save_best_only=True,
                                     mode='auto')

        config = tf.ConfigProto(gpu_options = tf.GPUOptions(allow_growth = True))
        session = tf.Session(config = config)
        KTF.set_session(session)
        self.encoder.summary()
        self.decoder.summary()
        history = self.auto_encoder.fit(x=x_train,
                                        y = x_train,
                                        verbose = 0,
                                        batch_size = self.batch_size,
                                        epochs = self.epoch,
                                        validation_split = 0.1,
                                        callbacks = [self.cp_cb])
        
        self.auto_encoder = getNewestModel(self.auto_encoder, str(self.chkptDir))

        with open(self.load_path_json_ae, 'w') as fout:
            model_json_str = self.auto_encoder.to_json(indent=4)
            fout.write(model_json_str)
        self.auto_encoder.save_weights(self.load_path_weight_ae)

        with open(self.load_path_json_encoder, 'w') as fout:
            model_json_str = self.encoder.to_json(indent=4)
            fout.write(model_json_str)
        self.encoder.save_weights(self.load_path_weight_encoder)

        with open(self.load_path_json_decoder, 'w') as fout:
            model_json_str = self.decoder.to_json(indent=4)
            fout.write(model_json_str)
        self.decoder.save_weights(self.load_path_weight_decoder)

        x_pred = self.auto_encoder.predict(x_train)
        print(r2_score(x_train, x_pred))
    
    def predict(self, x_test=None):
        x_test = self.check_x_input(x_test)
        return self.auto_encoder.predict(x_test)

    def transform(self, x_original):
        return self.encoder.predict(x_original)

    def inverse_transform(self, x_encoded):
        return self.decoder.predict(x_encoded)
        
def rmse(x_true, x_pred):
    return np.sqrt(np.sum((x_true - x_pred)**2, axis = 1))

def checkStat(x_true, x_pred):
    error = rmse(x_true, x_pred)
    eAverage = np.mean(error)
    eVariance = np.var(error)
    eMax = np.max(error)
    eMin = np.min(error)
    return np.array([eAverage, eMax, eMin, eVariance])

def linearMethod(start=2, stop=203):
    from sklearn.decomposition import PCA
    from training.read_training_data import load_mixed
    x_train, x_test, y_train, y_test, rescale_x, rescale_y3 = load_mixed(mode = shapeType, n_samples = 50000,
                                                                         test_size = 0.01, reductTarget = reductTarget,
                                                                         split45 = False)

    s_train = x_train[:,1:201]
    s_test = x_test[:,1:201]

    names = ["# of units", "average", "max", "min", "variance"]
    data = np.zeros((stop - start + 1, 5))
    for i in range(stop - start + 1):
        n_components = i + start
        pca = PCA(n_components = n_components)
        pca.fit(s_train)
        s_encode = pca.transform(s_train)
        s_decode = pca.inverse_transform(s_encode)
        data[i, 0] = n_components
        data[i, 1:] = checkStat(s_train, s_decode)
    
    fig = plt.figure()
    plt.title("RMSE vs. # of components (by PCA)")
    ax = fig.add_subplot(111)
    for i in range(1,4):
        ax.plot(data[:, 0], data[:, i], ".", label=names[i])
    ax.set_xlabel("Number of Principle Components")
    ax.set_ylabel("Root Mean Squared Error (log scale)")
    ax.legend()
    ax.grid()
    ax.set_yscale("log")
    plt.show()

def nonlinearMethod(start=1, stop=101):
    input_dim = (200,)
    from training.read_training_data import load_mixed
    x_train, x_test, y_train, y_test, rescale_x, rescale_y3 = load_mixed(mode = shapeType, n_samples = 50000,
                                                                         test_size = 0.01, reductTarget = reductTarget,
                                                                         split45 = False, __dimReduct=None)
    s_train = x_train[:,1:201]
    s_test = x_test[:,1:201]
    s_train = np.concatenate([s_train, s_test])
    names = ["# of units", "average", "max", "min", "variance"]
    data = np.zeros((stop - start + 1, 5))
    data_test = np.zeros((stop - start + 1, 5))
    for i in range(stop - start + 1):
        mid_dim = (i + start, )
        ae = AutoEncoder(input_dim, mid_dim, layer_num = 3, deep=True, x_train=s_train, x_test=s_test)
        # ae.train()
        s_decoded_train = ae.predict(s_train)
        # x_decoded_test = ae.predict()
        data[i, 0] = mid_dim[0]
        data[i, 1:] = checkStat(s_train, s_decoded_train)
    np.save("AEdeep2.npz", data)
    fig = plt.figure()
    plt.title("RMSE vs. # of mid-layer units (by AE)")
    ax = fig.add_subplot(111)
    for i in range(1, 4):
        ax.plot(data[:, 0], data[:, i], ".", label = names[i])
    ax.set_xlabel("Number of Mid-Layer's units")
    ax.set_ylabel("Root Mean Squared Error (log scale)")
    ax.legend()
    ax.grid()
    ax.set_yscale("log")
    plt.show()
    
def main():
    pass

if __name__ == '__main__':
    # linearMethod()
    nonlinearMethod()
