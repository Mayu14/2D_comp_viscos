# coding: utf-8
from functools import reduce

import keras.backend as K
from keras.layers import (Activation, Add, GlobalAveragePooling2D,
                          BatchNormalization, Conv2D, Dense, Flatten, Input,
                          MaxPooling2D)

from keras.models import Model
from keras.regularizers import l2


# refered by http://pynote.hatenablog.com/entry/keras-resnet-implementation
def compose(*funcs):
    """
    複数の層を結合するためのヘルパー関数
    :param funcs: f(x), g(x), h(x)
    :return: h(g(f(x))) = compose(f, g, h)
    """
    if funcs:
        return reduce(lambda f, g: lambda *args, **kwargs: g(f(*args, **kwargs)), funcs)
    else:
        raise ValueError('Composition of empty sequence not supported.')


def ResNetConv2D(*args, **kwargs):
    """conv layerを作成する
    """
    conv_kwargs = {
        'strides': (1, 1),
        'padding': 'same',
        'kernel_initializer': 'he_normal',
        'kernel_regularizer': l2(1.e-4)
    }
    conv_kwargs.update(kwargs)
    return Conv2D(*args, **conv_kwargs)


def bn_relu_conv(*args, **kwargs):
    """ batch normalization -> ReLU -> conv を作成する
    """
    return compose(
        BatchNormalization(),
        Activation('relu'),
        ResNetConv2D(*args, **kwargs)
    )


def shortcut(x, residual):
    """ shortcut connectionを作成
    """
    x_shape = K.int_shape(x)
    residual_shape = K.int_shape(residual)
    
    if x_shape == residual_shape:
        # x と residual の形状が同じ場合何もしない
        shortcut = x
    else:
        # x と residual の形状が異なる場合，線形変換を行い，形状を一致させる
        stride_w = int(round(x_shape[1] / residual_shape[1]))
        stride_h = int(round(x_shape[2] / residual_shape[2]))
        
        shortcut = Conv2D(filters = residual_shape[3],
                          kernel_size = (1, 1),
                          strides = (stride_w, stride_h),
                          kernel_initializer = 'he_normal',
                          kernel_regularizer = l2(1.e-4))(x)
    
    return Add()([shortcut, residual])


def basic_block(filters, first_strides, is_first_block_of_first_layer):
    """
    building blockの作成
    :param filters: (integer) フィルター数
    :param first_strides: (tuple) 最初の畳み込み層のストライド
    :param is_first_block_of_first_layer: (logical) max pooling直後のresidual blockか
    :return: building block layers
    """
    
    def f(x):
        if is_first_block_of_first_layer:
            # conv1でbatch normalization -> reluを適用しているため
            # max poolingの直後のresidual blockは畳み込みから
            conv1 = ResNetConv2D(filters = filters, kernel_size = (3, 3))(x)
        else:
            conv1 = bn_relu_conv(filters = filters, kernel_size = (3, 3),
                                 strides = first_strides)(x)
        
        conv2 = bn_relu_conv(filters = filters, kernel_size = (3, 3))(conv1)
        
        return shortcut(x, conv2)
    
    return f


def bottleneck_block(filters, first_strides, is_first_block_of_first_layer):
    """
    building blockの作成
    :param filters: (integer) フィルター数
    :param first_strides: (tuple) 最初の畳み込み層のストライド
    :param is_first_block_of_first_layer: (logical) max pooling直後のresidual blockか
    :return: building block layers
    """
    
    def f(x):
        if is_first_block_of_first_layer:
            # conv1でbatch normalization -> reluを適用しているため
            # max poolingの直後のresidual blockは畳み込みから
            conv1 = ResNetConv2D(filters = filters, kernel_size = (3, 3))(x)
        else:
            conv1 = bn_relu_conv(filters = filters, kernel_size = (1, 1),
                                 strides = first_strides)(x)
        
        conv2 = bn_relu_conv(filters = filters, kernel_size = (3, 3))(conv1)
        conv3 = bn_relu_conv(filters = filters * 4, kernel_size = (1, 1))(conv2)
        
        return shortcut(x, conv3)
    
    return f


def residual_blocks(block_function, filters, repetitions, is_first_layer):
    """
    residual blockを反復する構造を作成
    :param block_function: residual_blockを作成する関数
    :param filters: フィルター数
    :param repetitions: 何回繰り返すか
    :param is_first_layer: max pooling直後かどうか
    :return:
    """
    
    def f(x):
        for i in range(repetitions):
            # conv3_x, conv4_x, conv5_xの最初の畳み込みはプーリング目的　-> strides (2,2)
            # ただし，conv2_xの最初の畳み込みは直前のmax poolingでプーリング済み -> strides (1, 1)
            first_strides = (2, 2) if i == 0 and not is_first_layer else (1, 1)
            
            x = block_function(filters = filters, first_strides = first_strides,
                               is_first_block_of_first_layer = (i == 0 and is_first_layer))(x)
        
        return x
    
    return f


class ResnetBuilder():
    @staticmethod
    def build(input_shape, num_outputs, block_type, repetitions, continues=False):
        """
        ResNet モデルを作成するFactoryクラス
        :param input_shape: 入力形状
        :param num_outputs: ネットワークの出力数
        :param block_type: residual blockの種類 ('basic' or 'bottleneck')
        :param repetitions: 同じresidual blockの反復回数
        :return:
        """
        
        # block_typeに応じてresidual blockを生成する関数を変更
        if block_type == 'basic':
            block_fn = basic_block
        elif block_type == 'bottleneck':
            block_fn = bottleneck_block
        
        # モデルを作成
        input = Input(shape = input_shape)
        
        # conv1 (bn -> relu -> conv)
        conv1 = compose(ResNetConv2D(filters = 64, kernel_size = (7, 7), strides = (2, 2)),
                        BatchNormalization(),
                        Activation('relu'))(input)
        
        # pool
        pool1 = MaxPooling2D(
            pool_size = (3, 3), strides = (2, 2), padding = 'same'
        )(conv1)
        
        # conv2_x, conv3_x, conv4_x, conv5_x
        block = pool1
        filters = 64
        for i, r in enumerate(repetitions):
            block = residual_blocks(block_function = block_fn, filters = filters,
                                    repetitions = r, is_first_layer = (i == 0))(block)
            filters *= 2
        
        # batch normalization -> ReLU
        block = compose(BatchNormalization(),
                        Activation('relu'))(block)
        
        # global average pooling
        pool2 = GlobalAveragePooling2D()(block)
        
        if not continues:
            # dense
            fc1 = Dense(units = num_outputs,
                        kernel_initializer = 'he_normal',
                        activation = "softmax")(pool2)
            return Model(inputs = input, outputs = fc1)
        else:
            fc1 = Dense(units = num_outputs,
                        kernel_initializer = 'he_normal')(pool2)
            return input, fc1
    
    @staticmethod
    def build_resnet_18(input_shape, num_outputs):
        return ResnetBuilder.build(
            input_shape, num_outputs, 'basic', [2, 2, 2, 2])

    @staticmethod
    def build_resnet_18_continues(input_shape, num_outputs):
        return ResnetBuilder.build(
            input_shape, num_outputs, 'basic', [2, 2, 2, 2], continues = True)
    
    @staticmethod
    def build_resnet_34(input_shape, num_outputs):
        return ResnetBuilder.build(
            input_shape, num_outputs, 'basic', [3, 4, 6, 3])
    
    @staticmethod
    def build_resnet_50(input_shape, num_outputs):
        return ResnetBuilder.build(
            input_shape, num_outputs, 'bottleneck', [3, 4, 6, 3])
    
    @staticmethod
    def build_resnet_101(input_shape, num_outputs):
        return ResnetBuilder.build(
            input_shape, num_outputs, 'bottleneck', [3, 4, 23, 3])
    
    @staticmethod
    def build_resnet_152(input_shape, num_outputs):
        return ResnetBuilder.build(
            input_shape, num_outputs, 'bottleneck', [3, 8, 36, 3])


if __name__ == '__main__':
    input_shape = (128, 128, 1)
    num_classes = 200
    
    # model = ResnetBuilder.build_resnet_18(input_shape, num_classes)
    input, fc1 = ResnetBuilder.build_resnet_18_continues(input_shape, num_classes)
    model = Model(inputs=input, outputs=fc1)
    
    model.summary()
    model.count_params()
    
    from keras.utils import plot_model
    
    plot_model(model, to_file = 'resnet-model.png',
               show_shapes = True, show_layer_names = True)
