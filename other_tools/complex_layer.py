# coding: utf-8
import keras.backend as K
from keras.layers import Dense, Activation, Dropout, Multiply, Add, Lambda, Concatenate
from keras.layers.normalization import BatchNormalization
import keras.initializers

class MLPLayer():
    def __init__(self, activate_before_fc=True, batch_normalization=False, dropout=False,
                 dropout_rate=0.3,
                 dropout_timing = "final", gate_bias=-3, growth_rate = 32):
        """
        :param activate_before_fc: (bool) 活性化関数を全結合層の前に入れるかどうか
        :param batch_normalization: (bool) BNを有効にするかどうか
        :param dropout: (bool) Dropoutを入れるかどうか
        :param dropout_rate: (float) Dropoutさせる割合(0.0 ~ 1.0)
        :param dropout_timing: (character) Dropoutを入れるタイミング(every, even, finalから選ぶ)
        :param gate_bias:
        :param growth_rate:
        """
        self.activate_before_fc = activate_before_fc
        self.batch_normalization = batch_normalization

        self.dropout = dropout
        self.dropout_default = dropout
        self.dropout_rate = dropout_rate
        self.dropout_timing = "final"   # "even"

        self.gate_bias = gate_bias
        self.growth_rate = growth_rate

    def dropout_switch(self, layer_number, max_layer_number):
        if self.dropout_timing == "final":
            if layer_number == max_layer_number:
                self.dropout = self.dropout_default
            else:
                self.dropout = False

        elif self.dropout_timing == "even":
            if layer_number%2 == 0:
                self.dropout = self.dropout_default
            else:
                self.dropout = False
        elif self.dropout_timing == "every":
            self.dropout = self.dropout_default
        else:
            print('Invalid value is assigned to "self.dropout_timing"')

    def fully_connected(self, units, inputs, Activator):
        """
        :param units: 出力の素子数
        :param inputs: (tensor) 入力テンソル
        :param Activator: (function) 活性化関数(種類指定ではなく，関数で与える)例：Activation=Activation("sigmoid")
        :return:
        """
        input_units = K.int_shape(inputs)[-1]
        fx = Lambda(lambda x: x, output_shape=(input_units,))(inputs)  # fx = x
        # F(x)
        if self.batch_normalization:
            fx = BatchNormalization()(fx)

        if self.activate_before_fc:
            fx = Activator()(fx)

        if self.dropout:
            fx = Dropout(rate=self.dropout_rate)(fx)

        fx = Dense(units, activation=None, use_bias=True, kernel_initializer='he_normal',
                   bias_initializer='zeros', kernel_regularizer=None, bias_regularizer=None,
                   activity_regularizer=None, kernel_constraint=None, bias_constraint=None)(fx)

        if self.activate_before_fc == False:
            fx = Activator()(fx)

        return fx

    def residual(self, inputs, Activator, weight_layer_number=1, units_list=None):
        """
        :param inputs: (tensor) 入力テンソル
        :param Activator: (function) 活性化関数(種類指定ではなく，関数で与える)例：Activation=Activation("sigmoid")
        :param activate_before_fc: (bool) 活性化関数を全結合層の前に入れるかどうか
        :param weight_layer_number: (int) skipしない側の層数
        :param units_list: (list[int, int, ...]) weight_layerの素子数
        :return: (tensor)
        """

        # 指定が無ければ同じサイズの
        input_units = K.int_shape(inputs)[-1]
        if type(units_list) == type(None):
            units_list = [input_units] * weight_layer_number


        fx = Lambda(lambda x: x, output_shape=(units_list[0],))(inputs)
        for i in range(weight_layer_number):
            self.dropout_switch(i+1, weight_layer_number)
            fx = self.fully_connected(units=units_list[i], inputs=fx, Activator=Activator())

        # Addレイヤーの直前ではデータのサイズが一致していなければならない
        if units_list[-1] != input_units:
            fx = Dense(input_units)(fx)

        fx = Add()([fx, inputs])
        return Activator()(fx)    # F(x) + x


    def highway(self, inputs, Activator, weight_layer_number=1, units_list=None):
        """
        # referred by https://gist.github.com/iskandr/a874e4cf358697037d14a17020304535
        :param inputs: (tensor) 入力テンソル
        :param Activator: (function) 活性化関数(種類指定ではなく，関数で与える)例：Activation=Activation("sigmoid")
        :param weight_layer_number: (int) skipしない側の層数
        :param units_list: (list[int, int, ...]) weight_layerの素子数
        :return: (tensor)
        """

        input_units = K.int_shape(inputs)[-1]
        if type(units_list) == type(None):
            units_list = [input_units] * weight_layer_number

        gate_bias_initializer = keras.initializers.Constant(self.gate_bias)
        gate = Dense(units=units_list[0], bias_initializer=gate_bias_initializer)(inputs)
        gate = Activation("sigmoid")(gate)  # s

        negated_gate = Lambda(lambda x: 1.0 - x, output_shape=(units_list[0],))(gate)  # (1 - s)

        # F(x)
        fx = Lambda(lambda x: x, output_shape=(units_list[0],))(inputs)
        for i in range(weight_layer_number):
            self.dropout_switch(i + 1, weight_layer_number)
            fx = self.fully_connected(units=units_list[i], inputs=fx, Activator=Activator())

        if units_list[-1] != input_units:
            fx = Dense(input_units)(fx)

        fx_gated = Multiply()([gate, fx]) # sF(x)
        identity_gated = Multiply()([negated_gate, inputs]) # (1-s)x

        fx = Add()([fx_gated, identity_gated])   # sF(x) + (1-s)x
        return Activator()(fx)

    # refered by https://github.com/flyyufelix/DenseNet-Keras/blob/master/densenet121.py
    def denseblock(self, inputs, Activator, weight_layer_number=1):
        """
        :param inputs: (tensor) 入力テンソル
        :param Activator: (function) 活性化関数(種類指定ではなく，関数で与える)例：Activation=Activation("sigmoid")
        :param weight_layer_number: (int) skipしない側の層数
        :return: (tensor)
        """

        input_units = K.int_shape(inputs)[-1]
        concat_feat = Lambda(lambda x: x, output_shape=(input_units,))(inputs)   # concat_feat = x

        for i in range(weight_layer_number):
            fx = self.fully_connected(units=self.growth_rate, inputs=concat_feat, Activator=Activator())
            concat_feat = Concatenate()([concat_feat, fx])

        return Activator()(concat_feat)    # F(x) + x


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