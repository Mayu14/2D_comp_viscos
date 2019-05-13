# coding: utf-8
import keras.backend as K
from keras.layers import Dense, Activation, Dropout, Multiply, Add, Lambda
from keras.layers.normalization import BatchNormalization
import keras.initializers

# inputs: 入力テンソル
# Activation: 活性化関数(種類指定ではなく，関数で与える)例：Activation=Activation("sigmoid")
def residual(inputs, Activator, batch_normalization=False, dropout=False, dropout_rate=0.3):
    units = K.int_shape(inputs)[-1]

    fx = Lambda(lambda x: x, output_shape=(units,))(inputs)  # fx = x
    if batch_normalization:
        fx = BatchNormalization()(fx)

    fx = Activator(fx)
    if dropout:
        fx = Dropout(rate=dropout_rate)(fx)
    # F(x)
    fx = Dense(units, activation=None, use_bias=True, kernel_initializer='he_normal',
                       bias_initializer='zeros', kernel_regularizer=None, bias_regularizer=None,
                       activity_regularizer=None, kernel_constraint=None, bias_constraint=None)(fx)

    fx = Add()([fx, inputs])
    return Activator(fx)    # F(x) + x

# referred by https://gist.github.com/iskandr/a874e4cf358697037d14a17020304535
def highway(inputs, Activator=Activation("tanh"), gate_bias=-3,
            batch_normalization=False, dropout=False, dropout_rate=0.3):

    units = K.int_shape(inputs)[-1]

    gate_bias_initializer = keras.initializers.Constant(gate_bias)
    gate = Dense(units=units, bias_initializer=gate_bias_initializer)(inputs)
    gate = Activation("sigmoid")(gate)  # s

    negated_gate = Lambda(
        lambda x: 1.0 - x,
        output_shape=(units,))(gate)  # (1 - s)

    fx = Lambda(lambda x: x, output_shape=(units,))(inputs)  # fx = x
    if batch_normalization:
        fx = BatchNormalization()(fx)

    fx = Activator(fx)  # F(x)

    if dropout:
        fx = Dropout(rate=dropout_rate)(fx)

    fx = Dense(units, activation=None, use_bias=True, kernel_initializer='glorot_uniform',
                       bias_initializer='zeros', kernel_regularizer=None, bias_regularizer=None,
                       activity_regularizer=None, kernel_constraint=None, bias_constraint=None)(fx)

    fx_gated = Multiply()([gate, fx]) # sF(x)
    identity_gated = Multiply()([negated_gate, inputs]) # (1-s)x

    fx = Add()([fx_gated, identity_gated])   # sF(x) + (1-s)x
    return Activator(fx)

# dense_blockはDense + Activation + Dropoutをまとめた関数として与える．(inputを入れれば全結合層の出力が得られる形にする)
def residual_dense(dense_block, inputs, *other_inputs):
    fx_N = dense_block(inputs)    # F(x_N)

    sumlist = [fx_N, inputs]

    other_inputs = list(other_inputs)
    if len(other_inputs) != 0:
        sumlist.extend(other_inputs)

    return Add()(sumlist)

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