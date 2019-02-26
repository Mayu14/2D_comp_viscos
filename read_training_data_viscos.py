# -- coding: utf-8 --
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
# source:データの置いてあるディレクトリのパス(絶対or相対)
# fpath_lift:揚力係数の入ったcsvデータのsourceからの相対パス
# fpath_shape:形状データの入ったcsvデータのsourceからの相対パス
# shape_odd:物体形状ベクトルの次元の奇偶等（読み飛ばしに関する変数）
# read_rate:物体形状ベクトルの次元をread_rateで割った値に変更する
# skip_rate:揚力係数データ(教師データ)数をskip_rateで割った数に減らす
# regularize:正規化の有無
# scalar:正規化する際に、読み込んだデータを基準にz_scoreを求める場合"None"，別データ基準に正規化する場合はsklearnのStandardScalarを渡す
# return_scalar:scalarも一緒に返す場合True
def read_csv_type3(source, fpath_lift, fpath_shape, shape_odd = 0, read_rate = 1, total_data=1620*14, skip_rate = 1,
                   regularize=True, scalar="None", return_scalar=False):
    name = ("naca4", "angle", "lift_coef", "drag_coef")
    # data_typeとusecolsを書き換えることで情報量の増加に対応
    data_type = {"naca4": int, "angle":float, "lift_coef":float, "drag_coef":float}
    skip_row = np.arange(start=0, stop=total_data, step=skip_rate, dtype=int).tolist()

    df_l = pd.read_csv(source + fpath_lift, header=None,
                     usecols=[1, 3, 4, 5], names=name, dtype=data_type,
                       skiprows=skip_row)
    
    def make_use_cols_for_shape(data, shape_odd, rate):
        # dataは，shapeを何項まで計算したのかによって変更(今回は200固定でよさげ)
        # shape_oddによって全部読む・奇数のみ読む・偶数のみ読む，前半だけ読む、みたいな順番変更
        # dataのうち1/rateだけ読む
        odd_num = lambda index: (2 * index) + 1
        even_num = lambda index: (index + 1) * 2
        original_num = lambda index: index + 1
        skip_n_num = lambda index: (rate * index) + 1

        if((shape_odd != 0) and (rate == 1) or (rate == 0)):
            print("rate error!")
            exit()
        num = int(data / rate)
        if shape_odd == 1:    # 奇数番のみを抽出
            next = odd_num
        elif shape_odd == 2:  # 偶数番のみを抽出
            next = even_num
        elif shape_odd == 3:  # 並び順は同じだけど，数は半分
            next = original_num
        elif shape_odd == 4:
            next = skip_n_num   # 最後まで読むけどrateごとにスキップ
        else:   # 元々のまま読み出す
            next = original_num
            num = data

        col = [0]
        name = ["naca4"]
        data_type = {"naca4":int}
        for index in range(num):
            col.append(next(index))
            name.append("t" + str(index).zfill(3))
            data_type["t" + str(index).zfill(3)] = float
        return col, name, data_type

    col, name, data_type = make_use_cols_for_shape(data=200, shape_odd=shape_odd, rate=read_rate)

    df_s = pd.read_csv(source + fpath_shape, header=None,
                     usecols=col, names=name, dtype=data_type
                     )# .set_index("naca4")

    """このままだとクッソ重いので名前を被らせてメモリ節約する
    本当に書きたい処理はこれ
    df_xy = pd.merge(df_l, df_s, on="naca4")
    y_train = df_xy["lift_coef"].values
    X_train = df_xy.drop("lift_coef", axis=1).drop("naca4", axis=1).values
    """
    X_train = pd.merge(df_l, df_s, on="naca4", how="inner")
    # y_train_l = X_train["lift_coef"].values.reshape(-1, 1)
    # y_train_d = X_train["drag_coef"].values.reshape(-1, 1)
    # y_train = np.concatenate((y_train_l.T, y_train_d.T)).T
    y_train = np.concatenate((X_train["lift_coef"].values.reshape(-1, 1).T, X_train["drag_coef"].values.reshape(-1, 1).T)).T
    X_train = X_train.drop("lift_coef", axis=1).drop("naca4", axis=1).values

    if regularize:
        if scalar == "None":
            scalar = StandardScaler()
            scalar.fit(X_train)
        X_train = scalar.transform(X_train)
        
        if return_scalar:
            return X_train, y_train, scalar
    
    return X_train, y_train

if __name__ == '__main__':
    # 自宅で作成したのでLaboratory用に書き換える
    source = "D:\\Dropbox\\shareTH\\program\\keras_training\\"
    # source = "D:\\Toyota\\Dropbox\\shareTH\\program\\laplace\\"
    # ちなみにNACA4とNACA5は基本動作がほぼ一緒だから使い回せるtype3とtype4は同じやり方でok
    fpath_lift = "s1122_e9988_s4_a013.csv"
    fpath_shape = "shape_fourier_1112_9988_s04.csv"
    shape_odd = 0
    read_rate = 1
    read_csv_type3(source, fpath_lift, fpath_shape, shape_odd, read_rate, total_data=0)
