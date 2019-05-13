# coding: utf-8
import keras.backend as K
from keras.layers import Dense, Activation, Dropout, Multiply, Add, Lambda, Concatenate
from keras.layers.normalization import BatchNormalization
import keras.initializers

# inputs: (tensor) 入力テンソル
# Activator: (function) 活性化関数(種類指定ではなく，関数で与える)例：Activation=Activation("sigmoid")
# batch_normalization: (bool) BNを有効にするかどうか
# dropout: (bool) Dropoutを入れるかどうか
# dropout_rate: (float) Dropoutさせる割合(0.0 ~ 1.0)
# weight_layer_number: (int) skipしない側の層数
# units_list: (list[int, int, ...]) weight_layerの素子数
def residual(inputs, Activator, batch_normalization=False, dropout=False, dropout_rate=0.3,
             weight_layer_number=1, units_list=None):

    if type(units_list) == type(None):
        units_list = [K.int_shape(inputs)[-1]] * weight_layer_number

    fx = Lambda(lambda x: x, output_shape=(units_list[0],))(inputs)  # fx = x
    # F(x)
    for i in range(weight_layer_number):
        if batch_normalization:
            fx = BatchNormalization()(fx)

        fx = Activator(fx)
        if (dropout) and (i + 1 == weight_layer_number):    # only Final Layer
            fx = Dropout(rate=dropout_rate)(fx)

        fx = Dense(units_list[i], activation=None, use_bias=True, kernel_initializer='he_normal',
                           bias_initializer='zeros', kernel_regularizer=None, bias_regularizer=None,
                           activity_regularizer=None, kernel_constraint=None, bias_constraint=None)(fx)

    fx = Add()([fx, inputs])
    return Activator(fx)    # F(x) + x

# inputs: (tensor) 入力テンソル
# Activator: (function) 活性化関数(種類指定ではなく，関数で与える)例：Activation=Activation("sigmoid")
# gate_bias: (float) ゲートバイアス(シグモイド関数側)の初期化設定(デフォルトが推奨値)
# batch_normalization: (bool) BNを有効にするかどうか
# dropout: (bool) Dropoutを入れるかどうか
# dropout_rate: (float) Dropoutさせる割合(0.0 ~ 1.0)
# weight_layer_number: (int) skipしない側の層数
# units_list: (list[int, int, ...]) weight_layerの素子数
# referred by https://gist.github.com/iskandr/a874e4cf358697037d14a17020304535
def highway(inputs, Activator=Activation("tanh"), gate_bias=-3,
            batch_normalization=False, dropout=False, dropout_rate=0.3,
            weight_layer_number=1, units_list=None):

    if type(units_list) == type(None):
        units_list = [K.int_shape(inputs)[-1]] * weight_layer_number

    gate_bias_initializer = keras.initializers.Constant(gate_bias)
    gate = Dense(units=units_list[0], bias_initializer=gate_bias_initializer)(inputs)
    gate = Activation("sigmoid")(gate)  # s

    negated_gate = Lambda(
        lambda x: 1.0 - x,
        output_shape=(units_list[0],))(gate)  # (1 - s)

    # F(x)
    fx = Lambda(lambda x: x, output_shape=(units_list[0],))(inputs)  # fx = x
    for i in range(weight_layer_number):
        if batch_normalization:
            fx = BatchNormalization()(fx)

        fx = Activator(fx)

        if (dropout) and (i + 1 == weight_layer_number):  # only Final Layer
            fx = Dropout(rate=dropout_rate)(fx)

        fx = Dense(units_list[i], activation=None, use_bias=True, kernel_initializer='glorot_uniform',
                   bias_initializer='zeros', kernel_regularizer=None, bias_regularizer=None,
                   activity_regularizer=None, kernel_constraint=None, bias_constraint=None)(fx)

    fx_gated = Multiply()([gate, fx]) # sF(x)
    identity_gated = Multiply()([negated_gate, inputs]) # (1-s)x

    fx = Add()([fx_gated, identity_gated])   # sF(x) + (1-s)x
    return Activator(fx)

# refered by https://github.com/flyyufelix/DenseNet-Keras/blob/master/densenet121.py
# dense_blockはDense + Activation + Dropoutをまとめた関数として与える．(inputを入れれば全結合層の出力が得られる形にする)
def densenet(inputs, Activator, growth_rate = 32, batch_normalization=False, dropout=False, dropout_rate=0.3,
                   weight_layer_number=1):

    units = K.int_shape(inputs)[-1]
    concat_feat = Lambda(lambda x: x, output_shape=(units,))(inputs)   # concat_feat = x
    for i in range(1, weight_layer_number + 1):
        fx = Lambda(lambda x: x, output_shape=(units,))(concat_feat)
        if batch_normalization:
            fx = BatchNormalization()(fx)

        fx = Activator(fx)
        if (dropout) and (i == weight_layer_number):    # only Final Layer
            fx = Dropout(rate=dropout_rate)(fx)

        fx = Dense(growth_rate, activation=None, use_bias=True, kernel_initializer='he_normal',
                      bias_initializer='zeros', kernel_regularizer=None, bias_regularizer=None,
                      activity_regularizer=None, kernel_constraint=None, bias_constraint=None)(fx)

        concat_feat = Concatenate([concat_feat, fx])

    return Activator(concat_feat)    # F(x) + x

def rec23(x):
    print(x)
    if x < 3:
        exit()
    return rec23(int(x/3*2))

if __name__ == '__main__':
    # import keras
    # print(keras.__version__)
    # residual_dense()
    rec23(204)
    exit()