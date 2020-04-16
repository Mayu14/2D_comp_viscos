# coding: utf-8
import math
import pandas as pd
import numpy as np
import re
from pathlib import Path
from sklearn.model_selection import train_test_split
from sklearn.cluster import KMeans, MiniBatchKMeans
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA, KernelPCA
from functools import partial
from itertools import product
from other_tools.setup_foam import load_patternList
from auto_encoder.encodeAndDecode import AutoEncoder

# parameters
baseDir = Path("G:\\Toyota\\Data\\Compressible_Viscos\\rhoSimpleFoam\\training_data\\")
aeroData = baseDir / Path("aerocharacteristics\\aerodynamicForce_mk2_CV0.csv")
aeroData2 = baseDir / Path("aerocharacteristics\\aerodynamicForce_mk2_CV_append0.csv")

c_pkl = baseDir / Path("shape_vectors\\crowd_CV240.pkl")
shapeDataC4 = baseDir / Path("shape_vectors\\shape_crowd_0.1_0.15_30_50_20_CV_240_xy.csv")
shapeDataC5 = baseDir / Path("shape_vectors\\shape_crowd_all_mk2_xy.csv")
shapeDataC = [shapeDataC4, shapeDataC5, c_pkl]

e_pkl = baseDir / Path("shape_vectors\\equidist_CV240.pkl")
shapeDataE4 = baseDir / Path("shape_vectors\\shape_equidistant_CV_240.csv")
shapeDataE5 = baseDir / Path("shape_vectors\\shape_equidistant_all_mk2.csv")
shapeDataE = [shapeDataE4, shapeDataE5, e_pkl]

f_pkl = baseDir / Path("shape_vectors\\fourier_CV240.pkl")
shapeDataF4 = baseDir / Path("shape_vectors\\shape_modified_fourier_CV_240.csv")
shapeDataF5 = baseDir / Path("shape_vectors\\shape_modified_fourier_all_mk2.csv")
shapeDataF = [shapeDataF4, shapeDataF5, f_pkl]

cf_pkl = baseDir / Path("shape_vectors\\circum_CV240.pkl")
shapeDataCF4 = baseDir / Path("shape_vectors\\shape_circumferential_CV_240.csv")
shapeDataCF5 = baseDir / Path("shape_vectors\\shape_circumferential_all_mk2.csv")
shapeDataCF = [shapeDataCF4, shapeDataCF5, cf_pkl]

np.random.seed(1)
shapeN = 200

threshold = 3000

DEBUG = False

def make_df(shapeData, concat=False, caseName = "CV_mixed", reductTarget="", split45=False):
    def read_pkl():
        df = pd.read_pickle(fname)
        if not reductTarget == "":
            df_r = pd.read_pickle(fname_r)
        else:
            df_r = pd.DataFrame()
        return df, df_r
        
    def make_pkl():
        names = ["NACA", "AoA", "Ma", "Re", "Cm", "Cd", "Cl", "Cl(r)", "Cl(f)"]
        aeroDf = pd.read_csv(aeroData, names = names)
    
        names2 = ["NACA"]
        [names2.append("b" + str(i + 1).zfill(3)) for i in range(shapeN)]
        shapeDf4 = pd.read_csv(shapeData[0], names = names2)
        shapeDf5 = pd.read_csv(shapeData[1], names = names2)
        shapeDf = pd.concat([shapeDf4, shapeDf5])

        df = pd.merge(aeroDf, shapeDf, on = "NACA", how = "left")

        if concat:
            df = df[df['Ma'] != 0.5]
            df = df[df['Re'] != 100000]
            df = df[df['Re'] != 10000]
            aeroDf = pd.read_csv(aeroData2, names = names)
            df2 = pd.merge(aeroDf, shapeDf, on = "NACA", how = "left")
            df = pd.concat([df, df2])

        if split45:
            df5 = df[df["NACA"] >= 10000]
            df = df[df["NACA"] < 10000]
        else:
            df5 = pd.DataFrame()
            
        df_r = pd.DataFrame()
        reductList = generateReductList(caseName, reductTarget, df)
        if reductList:  # reductList is not empty
            if not "kmeans" in reductTarget.lower():
                for keys in reductList:
                    df_r = pd.concat([df[df[keys[0]] == keys[1]], df_r])
                    df = df[df[keys[0]] != keys[1]]
            else:
                reductList.sort()
                reductList.reverse()
                reductList = np.array(reductList, dtype = int).reshape(-1)
                df = df.reset_index(drop = True)
                df_r = df.copy()
                df = df_r.T[reductList].T
                df_r = df_r.drop(df_r.index[reductList])
        
            df_r = df_r.sort_values(by = ["NACA", "Ma", "Re", "AoA"], ascending = True)
            df_r["NACA"] = df_r["NACA"].astype(float)
            if (not split45) and (not DEBUG):
                df_r.to_pickle(fname_r)
    
        if not df.isnull().values.sum() == 0:
            print("DataFrame has NaN value")
            raise ValueError
        else:
            df = df.sort_values(by = ["NACA", "Ma", "Re", "AoA"], ascending = True)
            if (not split45) and (not DEBUG):
                df.to_pickle(fname)
        return df, df_r, df5
    ###
    fdir = shapeData[2].parent
    fname_str = shapeData[2].name
    if not reductTarget == "":
        fname_str = reductTarget + "_less_" + fname_str
    if concat:
        fname_str = "concat" + fname_str
    fname = fdir / Path(fname_str)
    fname_r = fdir / Path("r_" + fname_str)
    
    if ((fname.exists() and reductTarget == "") or (fname.exists() and fname_r.exists() and (not reductTarget == ""))) and (not split45) and (not DEBUG):
        try:
            df, df_r = read_pkl()
        except ModuleNotFoundError as e:
            print("catch ModuleNotFoundError:", e)
            df, df_r, df5 = make_pkl()
    else:
        df, df_r, df5 = make_pkl()
    def correctMach(df):
        gamma = 1.4
        machs = [0.3, 0.4]
        labels = ["Cm", "Cd", "Cl", "Cl(r)", "Cl(f)"]
        for mach, label in product(machs, labels):
            coef = gamma * mach ** 2
            df.loc[df["Ma"] == mach, label] = df[label].loc[df["Ma"] == mach] * coef
            # df[labels].loc[df["Ma"] == mach] *= coef
        return df

    df = correctMach(df)

    if not df_r.empty:
        df_r = correctMach(df_r)
    df["NACA"] = df["NACA"].astype(float)
    df = [df, df_r]
    if not split45:
        return df
    else:
        df5 = correctMach(df5)
        def gen_df(df, df5):
            yield df
            df5 = df5.sort_values(by = ["NACA", "Ma", "Re", "AoA"], ascending = True)
            df = [df5, pd.DataFrame()]
            yield df
        return gen_df(df, df5)
        
        
def generateReductList(caseName = "CV_mixed", reductTarget = "", df=None):
    def genListWithKey(strkey, arrValues):
        lists = []
        for val in arrValues:
            lists.append([strkey, val])
        return lists

    if (not "aoa" in reductTarget.lower()) and (not "naca" in reductTarget.lower()) and (not "kmeans" in reductTarget.lower()):
        reductTarget = "NoReduct"
    
    reductList = []
    # rewrite every time
    # AoA reuction for "CV_mixed"
    num = re.sub(r'\D', '', reductTarget)
    if num != "":
        rate = int(num)
    else:
        rate = None
    if "aoa" in reductTarget.lower():
        if not rate is None:
            angleList = [-7.5, -5.0, -2.5, 0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0,
                         32.5, 35.0]
            n = len(angleList)
            red_angleList = []
            total = 0
            while n - 1 >= total:
                red_angleList.append(angleList[round(total)])
                total += float(n / rate)
            reductList.extend(genListWithKey("AoA", red_angleList))
        else:
            raise ValueError

    # wing reduction for "CV_mixed"
    if "naca" in reductTarget.lower():
        if (rate is None) and ("half" in reductTarget.lower()):
            raise ValueError

        for digit in [True, False]:
            reList, maList, aoaList, i1_list, i2_list, i3_list, i2_len = load_patternList(naca4 = digit, case = caseName)
            i4_list = []
            for i1, i2, i3 in product(i1_list, i2_list, i3_list):
                if (((i1 == 0) and (i2 == 0)) or (i1 != 0)):
                    i4_list.append(int(str(i1) + str(i2).zfill(i2_len) + str(i3)))
            n = len(i4_list)
            red_i4List = []
            total = 0
            while n - 1 >= total:
                red_i4List.append(i4_list[round(total)])
                total += float(n / rate)
            reductList.extend(genListWithKey("NACA", red_i4List))
            
    
    if "kmeans" in reductTarget.lower() and not df is None:
        print(reductTarget)
        k_cluster = int(re.sub(r'\D', '', reductTarget))
        midFile = f_pkl.parent / Path(reductTarget + ".npz")
        if midFile.exists():
            loaded_array = np.load(midFile)
            #  __x = loaded_array["__x"]
            # centroids = loaded_array["centroids"]
            nearest_indices = loaded_array["nearest_indices"]
        else:
            __x, _ = split_column(df, y_keys = [])
            
            if "norm" in reductTarget.lower():
                scalar = StandardScaler()
                __x = scalar.fit_transform(__x)
            
            if "pca" in reductTarget.lower():
                if "kpca" in reductTarget.lower():
                    pca = KernelPCA(n_components = 20, kernel = "rbf", n_jobs = -1, random_state = 1)
                else:
                    pca = PCA(n_components=20, random_state = 1)
                __x = pca.fit_transform(__x)
    
            debug = False
            if debug:
                samples = __x.shape[0]
                nearest_indices = np.random.choice(samples, k_cluster, replace = False)
            else:
                if k_cluster > threshold:
                    split = math.floor(float(k_cluster) / threshold)
                    batch_size = math.floor(float(k_cluster) / split)
                    kmeans = MiniBatchKMeans(n_clusters = k_cluster,  # クラスタ数
                                             batch_size = batch_size,
                                             init = "k-means++",  # セントロイド初期化アルゴリズム
                                             n_init = 10,  #
                                             init_size = k_cluster + batch_size,
                                             max_iter = 300,  # 最大反復回数
                                             tol = 1e-04,  # 収束判定条件
                                             random_state = 1)
                else:
                    kmeans = KMeans(n_clusters=k_cluster,  # クラスタ数
                                    init="k-means++",  # セントロイド初期化アルゴリズム
                                    n_init=10,  #
                                    max_iter=300,  # 最大反復回数
                                    tol=1e-04,  # 収束判定条件
                                    n_jobs=-1,  # 並列実行数(-1 = 環境限界まで利用)
                                    random_state=1)
                kmeans.fit(__x)
                centroids =  kmeans.cluster_centers_
                nn = NearestNeighbors(metric = "minkowski")
                nn.fit(__x)
                distances, nearest_indices = nn.kneighbors(centroids, n_neighbors = 1)
                np.savez_compressed(str(midFile), __x = __x, centroids = centroids, nearest_indices = nearest_indices)
            
        reductList = nearest_indices.reshape(-1).tolist()
    return reductList

def load_df_without_shape(concat=False, keyLists=None):
    dfs = make_df(shapeDataF, concat)
    df = dfs[0]
    defaultKeys = ["NACA", "Ma", "Re", "AoA"]
    if keyLists is None:
        keyLists = ["Cm", "Cd", "Cl", "Cl(r)", "Cl(f)"]
    defaultKeys.extend(keyLists)
    return df[defaultKeys]
    
def sample_plot(df, caseName="CV"):
    import matplotlib.pyplot as plt

    digits = [True, False]
    for digit in digits:
        reList, maList, aoaList, i1_list, i2_list, i3_list, i2_len = load_patternList(naca4 = digit, case = caseName)
        row = int(len(reList) / 2)
        for i1, i2, i3 in product(i1_list, i2_list, i3_list):
            if (((i1 == 0) and (i2 == 0)) or (i1 != 0)):
                naca_digit = str(i1) + str(i2).zfill(i2_len) + str(i3)
                sub_df = df[df["NACA"] == float(naca_digit)]

                titleSize = 18
                axisLabelSize = 14
                fig = plt.figure(figsize=(16, 7*row), dpi=200)
                fig.suptitle("NACA" + naca_digit + ": Aerodynamic Forces", fontsize=titleSize + 4)
                colors = ["red", "blue", "red"]
                colorsD = ["darkorange", "dodgerblue", "red"]
                for i, re in enumerate(reList):
                    subsub = sub_df[sub_df["Re"] == re]
                    ax = fig.add_subplot(row,2,i+1)
                    ax.set_title("Re: $10^" + str(int(np.log10(re))) + "$", fontsize = titleSize)
                    ax.set_ylim([-1.0, 2.0])
                    ax.set_xlabel("AoA [degree]", fontsize = axisLabelSize)
                    ax.set_ylabel("$C_L$ & $C_D$", fontsize = axisLabelSize)
                    for j, ma in enumerate(maList):
                        sss_df = subsub[subsub["Ma"] == ma]
                        ax.plot(sss_df["AoA"], sss_df["Cl"], marker="o", linestyle="dashed", color=colors[j], label = "$C_L$ Ma:" + str(ma), markersize=10, linewidth=4)
                        ax.plot(sss_df["AoA"], sss_df["Cd"], marker="^", linestyle="dotted", color=colorsD[j], label = "$C_D$ Ma:" + str(ma), markersize=10, linewidth=4)
                    ax.legend()

                fpath = baseDir / Path("overview_5\\" + caseName + "_NACA" + naca_digit + ".png")
                plt.savefig(fpath)
                plt.close()
    
def logarithmic_Tr(df, key, base=10):
    logarithm = lambda x: np.log10(x) / np.log10(base)
    df[key] = df[key].apply(logarithm)
    return df

def standardization(df, default_method="unbiased", dict=None):
    def transform_interface(df, method, key=None):
        def transform(sr, method):
            def regularize(sr):
                if type(sr) is pd.core.frame.DataFrame:
                    return (sr - sr.apply(np.min)) / (sr.apply(np.max) - sr.apply(np.min))
                elif type(sr) is pd.core.series.Series:
                    return (sr - sr.min()) / (sr.max() - sr.min())
                
            def standardize(sr, ddof):
                std = partial(np.std, ddof=ddof)
                if type(sr) is pd.core.frame.DataFrame:
                    return (sr - sr.apply(np.average)) / sr.apply(std)
                elif type(sr) is pd.core.series.Series:
                    return (sr - sr.mean()) / sr.std(ddof=ddof)
            
            if method == "minmax":
                return regularize(sr)
            elif method == "std":
                return standardize(sr, 0)
            elif method == "unbiased":
                return standardize(sr, 1)
            elif method is None:
                return sr
            else:
                raise ValueError
        
        if key is None:
            std_df = transform(sr=df, method = method)
        else:
            std_df = transform(sr = df[key], method = method)
        return std_df

    std_df = transform_interface(df, default_method)
    
    if dict is not None:
        for key, method in dict.items():
            std_df[key] = transform_interface(df[key], method)
    
    return std_df

def p_process(p_data):
    logRe = p_data[:, 2].reshape(-1, 1)
    re = 10**logRe
    naca = np.array(p_data[:, 3], dtype=int).reshape(-1, 1)
    thickness = (naca%100)/100.0
    max_camber_pos = np.where(naca < 10 ** 5, ((naca - (naca // 1000) * 1000) // 100) / 10,
                              ((naca - (naca // 10000) * 10000) // 1000) / 20)
    p_data_2 = np.concatenate([re, thickness, max_camber_pos], axis=1)
    # AoA, Ma, LogRe, NACA, Cm, Cl, Cm, Re, Thickness, MaxCamberPosition
    p_data = np.concatenate([p_data, p_data_2], axis=1)
    # max_camber = np.where(naca < 10**5, naca//1000, -1)
    return p_data

def split_column(df, y_keys, returnParam=False):
    def extractParam(df):
        p_data = []
        p_data.append(df["AoA"].values.reshape(-1, 1).T)
        p_data.append(df["Ma"].values.reshape(-1, 1).T)
        p_data.append(df["Re"].values.reshape(-1, 1).T)
        p_data.append(df["NACA"].values.reshape(-1,1).T)
        p_data.append(df["Cm"].values.reshape(-1, 1).T)
        p_data.append(df["Cl"].values.reshape(-1, 1).T)
        p_data.append(df["Cd"].values.reshape(-1, 1).T)
        return np.concatenate(p_data).T.astype(np.float64)
    
    if df.empty:
        return None, None
    else:
        if returnParam:
            return extractParam(df)
        
        all_y_keys = ["NACA", "Cm", "Cd", "Cl", "Cl(r)", "Cl(f)"]
        y_data = []
        for key in y_keys:
            y_data.append(df[key].values.reshape(-1, 1).T)
        if not y_keys:
            y_data = None
        else:
            y_data = np.concatenate(y_data).T.astype(np.float64)
    
        for key in all_y_keys:
            df = df.drop(key, axis = 1)
        x_data = df.values.astype(np.float64)
        # AoA, Ma, log(Re), b1, b2, ..., b200
        feature_index = np.concatenate([np.array([0]), np.arange(shapeN) + 3, np.array([1,2])])
        x_data = x_data[:, feature_index]
        # AoA, b1, b2, ..., b200, Ma, log(Re)
        return x_data, y_data

def __load_msdf(logarithmic=["Re"], msdf_res = 16, caseName="CV_mixed", reductTarget=""):
    if not reductTarget == "":
        reductTarget += "_less_"

    shapeDataM4 = baseDir / Path("input_vectors\\NACA4\\{0}{1}_c_v_4.npz".format(msdf_res, reductTarget))
    shapeDataM5 = baseDir / Path("input_vectors\\NACA5\\{0}{1}_c_v_5.npz".format(msdf_res, reductTarget))
    if not shapeDataM4.exists():
        print(shapeDataM4)
        print("reduction npz file not exists")
        exit()
    loaded_array = np.load(shapeDataM4)
    X_train_img = loaded_array["x_img"]
    X_train_label = loaded_array["x_label"]
    y_train = loaded_array["y_train"]
    
    loaded_array = np.load(shapeDataM5)
    X_test_img = loaded_array["x_img"]
    X_test_label = loaded_array["x_label"]
    y_test = loaded_array["y_train"]
    del loaded_array
    
    X_train_img = np.concatenate([X_train_img, X_test_img])
    X_train_label = np.concatenate([X_train_label, X_test_label])
    y_train = np.concatenate([y_train, y_test])
    
    if logarithmic == ["Re"]:
        X_train_label[:, 1] = np.log10(X_train_label[:, 1])
    else:
        raise ValueError
    
    return [X_train_img, X_train_label], y_train, None, None

def main(shapeType="Fourier", logarithmic=["Re"],
         x_keys=["AoA", "Re", "Ma", "b"], y_keys=["Cl", "Cd"], msdf_res=16, caseName = "CV_mixed", reductTarget="", split45=False, with_param=False):
    """
    
    :param shapeType:
    :param logarithmic:
    :param x_keys:
    :param y_keys:
    :param standardize_method:
    :param standardized_dict: sample_dict = {"NACA": None, "AoA": "minmax", "Ma": "minmax", "Re": "minmax"}
    :return:
    """
    # shapeType = "Crowd"#"Equidistant"#"Fourier"
    concat = "Concat" in shapeType
    if "fourier" in shapeType.lower():
        dfs = make_df(shapeDataF, concat, caseName, reductTarget, split45)
    elif "equidistant" in shapeType.lower():
        reductTarget += "equidistant"
        dfs = make_df(shapeDataE, concat, caseName, reductTarget, split45)
    elif "crowd" in shapeType.lower():
        reductTarget += "crowd"
        dfs = make_df(shapeDataC, concat, caseName, reductTarget, split45)
    elif "circum" in shapeType.lower():
        reductTarget += "circum"
        dfs = make_df(shapeDataCF, concat, caseName, reductTarget, split45)
    elif "msdf" in shapeType.lower():
        # reductTarget += "msdf"
        return __load_msdf(logarithmic = logarithmic, msdf_res = msdf_res, caseName = caseName, reductTarget = reductTarget)
    else:
        raise ValueError

    # sample_plot(dfs[0],  caseName = "CV_mixed")
    if not split45:
        for df in dfs:
            if df.empty:
                pass
            else:
                for key in logarithmic:
                    df = logarithmic_Tr(df, key)
                """
                if not standardize_method is None:
                    std_df = standardization(df, default_method = standardize_method, dict = standardized_dict)
                    # std_df = standardization(df, default_method = None, dict = None)
                else:
                    std_df = df
                """
        x__, y__ = split_column(dfs[0], y_keys)
        x_r, y_r = split_column(dfs[1], y_keys)
        if not with_param:
            return x__, y__, x_r, y_r
        else:
            p__ = split_column(dfs[0], y_keys = None, returnParam = True)
            p_r = split_column(dfs[1], y_keys = None, returnParam = True)
            return x__, y__, x_r, y_r, p__, p_r
    else:
        def gen_xy(dfs, logarithmic, y_keys):
            for df in dfs:
                for key in logarithmic:
                    df[0] = logarithmic_Tr(df[0], key)
                x__, y__ = split_column(df[0], y_keys)
                if not with_param:
                    yield x__, y__, None, None
                else:
                    p__ = split_column(df[0], y_keys = None, returnParam = True)
                    yield x__, y__, p__, None
        return gen_xy(dfs, logarithmic, y_keys)
    
def mean_and_std(x_mixed):
    x_mean = x_mixed.mean(axis = 0)
    x_std = x_mixed.std(axis = 0, ddof = 1)
    return x_mean, x_std

def reductDimension(x_data, scale_x, rescale_x, n_dim_with_method="pca10", excludeFeat=[0,1,2], mode="FourierConcat"):
   originDim = x_data.shape[1]
   targetFeat = np.ones(originDim, dtype = bool)
   targetFeat[excludeFeat] = False
   num_excluded = len(excludeFeat)

   def splitFeat(x_data):
       return x_data[:, targetFeat], x_data[:, excludeFeat].reshape(-1, num_excluded)  # reshape: for len(excludeFeat) == 1
   concatFeat = lambda x_target, x_excluded: np.concatenate([x_excluded, x_target], axis = 1)
   
   x_target, x_excluded = splitFeat(x_data)

   n_feat = int(re.sub(r'\D', '', n_dim_with_method))
   if "pca" in n_dim_with_method.lower():
       if "kpca" in n_dim_with_method.lower():
           pca = KernelPCA(n_components = n_feat, kernel = "rbf", n_jobs = -1, random_state = 1)
       else:
           pca = PCA(n_components = n_feat, random_state = 1)
           
       x_target = pca.fit_transform(x_target)

   elif "ae" in n_dim_with_method.lower():
       pca = AutoEncoder(input_dim=(x_target.shape[1],), mid_dim=(n_feat,), x_train=x_target, mode=mode.lower().replace("concat", ""))
       x_target = pca.transform(x_target)
   else:
       raise ValueError

   x_data = concatFeat(x_target, x_excluded)

   scale_x1, rescale_x1 = scale_x, rescale_x
   x_mean, x_std = x_data.mean(axis = 0), x_data.std(axis = 0, ddof = 1)
   scale_x2 = lambda x_data: (x_data - x_mean) / x_std
   rescale_x2 = lambda x_data_in: x_data_in * x_std + x_mean
   if np.any(x_std == 0.0):
       scale_x2, rescale_x2, _, _ = overwriteScaling(x_data, x_data[:2])
   x_data = scale_x2(x_data)
   
   def scale_x12(x_data):
       x_target, x_excluded = splitFeat(x_data)
       x_data = concatFeat(pca.transform(x_target), x_excluded)
       return scale_x2(x_data)
   
   def rescale_x12(x_data):
       x_data = rescale_x2(x_data)  # pca直後に戻す
       s_data = pca.inverse_transform(x_data[:, num_excluded:])  # shapeのみpca逆変換(次元復元)
       x_data = np.concatenate([x_data[:, :num_excluded].reshape(-1, num_excluded), s_data], axis = 1)  # angle & shape concat
       return rescale_x1(x_data)  # 元のスケールに戻す

   return x_data, scale_x12, rescale_x12

def sortFeatures(x_data):
    if x_data is None:
        return None
    else:
        nFeat = x_data.shape[1]
        idx = np.arange(nFeat, dtype=int)
        maRe = idx[-2:]
        shape = idx[1:-2]
        idx = np.concatenate([np.array([0]), maRe, shape])
        return x_data[:, idx]

def overwriteScaling(x_data, y_data):
    def scaleAndRescale(data):
        scalar = StandardScaler()
        scalar.fit(data)
        return scalar.transform, scalar.inverse_transform
    scale_x, rescale_x = scaleAndRescale(x_data)
    scale_y, rescale_y = scaleAndRescale(y_data)
    return scale_x, rescale_x, scale_y, rescale_y
    

def load_mixed(mode="FourierConcat", n_samples=20000, test_size = 0.33, forGan=False, msdf_res=16, caseName = "CV_mixed", reductTarget="", split45=False, __dimReduct=None, with_param=False, shuffle=True):
    # AoA, b1, b2, ..., b200, Ma, log(Re)
    if not split45:
        if not with_param:
            x_mixed, y_mixed, x_mixed_r, y_mixed_r = main(shapeType = mode, msdf_res = msdf_res, caseName = caseName, reductTarget=reductTarget, split45=split45)
        else:
            x_mixed, y_mixed, x_mixed_r, y_mixed_r, p_mixed, p_mixed_r = main(shapeType = mode, msdf_res = msdf_res,
                                                                              caseName = caseName,
                                                                              reductTarget = reductTarget,
                                                                              split45 = split45,
                                                                              with_param = with_param)
    else:
        reductTarget += "Split"
        generator = main(shapeType = mode, msdf_res = msdf_res, caseName = caseName, reductTarget = reductTarget, split45 = split45)
        if not with_param:
            x_mixed, y_mixed, x_mixed_r, y_mixed_r = generator.__next__()
        else:
            x_mixed_r, y_mixed_r = None, None
            x_mixed, y_mixed, p_mixed, _ = generator.__next__()

    if "msdf" in mode.lower():
        return __load_mixed_msdf(x_mixed[0], x_mixed[1], y_mixed, n_samples, test_size)

    if not __dimReduct is None:
        x_mixed, x_mixed_r = sortFeatures(x_mixed), sortFeatures(x_mixed_r)
        
    x_mean, x_std = mean_and_std(x_mixed)
    y_mean, y_std = mean_and_std(y_mixed)
    
    scale_x = lambda x_data: (x_data - x_mean) / x_std
    rescale_x = lambda x_data: x_data * x_std + x_mean
    scale_y = lambda y_data: (y_data - y_mean) / y_std
    rescale_y = lambda y_data: y_data * y_std + y_mean

    if np.any(x_std == 0.0) or np.any(y_std == 0.0):
        scale_x, rescale_x, scale_y, rescale_y = overwriteScaling(x_mixed, y_mixed)

    # if shuffle:
    samples = x_mixed.shape[0]
    if n_samples > samples:
        n_samples = samples
    
    if shuffle:
        index = np.random.choice(samples, n_samples, replace = False)  # 重複なし並び替え
    else:
        index = np.arange(n_samples)
        
    x_mixed = scale_x(x_mixed)[index]
    y_mixed = scale_y(y_mixed)[index]
    if with_param:
        p_mixed = p_mixed[index]

    if forGan:
        if not with_param:
            return x_mixed, x_mean, x_std, rescale_x, y_mixed, rescale_y
        else:
            return x_mixed, x_mean, x_std, rescale_x, y_mixed, rescale_y, p_process(p_mixed), p_process(p_mixed_r)
        
    if not split45:
        if not with_param:
            x_train, x_test, y_train, y_test = train_test_split(x_mixed, y_mixed, test_size = test_size, random_state = 42)
        else:
            x_train, x_test, y_train, y_test, p_train, p_test = train_test_split(x_mixed, y_mixed, p_mixed, test_size = test_size,
                                                                random_state = 42)
        if not x_mixed_r is None:
            x_mixed_r = scale_x(x_mixed_r)
            y_mixed_r = scale_y(y_mixed_r)
            x_test = np.concatenate([x_test, x_mixed_r])
            y_test = np.concatenate([y_test, y_mixed_r])
            if with_param:
                p_test = np.concatenate([p_test, p_mixed_r])
        if not __dimReduct is None:
            x_train, scale_x2, rescale_x = reductDimension(x_train, scale_x, rescale_x, n_dim_with_method = __dimReduct, mode=mode)
            x_test = scale_x2(x_test)
        
        if not with_param:
            return x_train, x_test, y_train, y_test, rescale_x, rescale_y
        else:
            return x_train, x_test, y_train, y_test, rescale_x, rescale_y,p_process(p_train), p_process(p_test)

    else:
        x_train, y_train = x_mixed, y_mixed
        x_test, y_test, p_test, _ = generator.__next__()
        x_test = scale_x(x_test)
        if not __dimReduct is None:
            x_train, scale_x2, rescale_x = reductDimension(x_train, scale_x, rescale_x, n_dim_with_method = __dimReduct, mode=mode)
            x_test = scale_x2(x_test)
        if not with_param:
            return x_train, x_test, y_train, scale_y(y_test), rescale_x, rescale_y
        else:
            return x_train, x_test, y_train, scale_y(y_test), rescale_x, rescale_y, p_process(p_mixed), p_process(p_test)
        
def __load_mixed_msdf(x_img, x_label, y_mixed, n_samples, test_size):
    x_mean, x_std = mean_and_std(x_label)
    y_mean, y_std = mean_and_std(y_mixed)

    scale_x_img = lambda x_data: (x_data / 255.0)  # - 127.5) / 127.5
    rescale_x_img = lambda x_data: (x_data * 255).astype("int64")  # 127.5 + 127.5)
    scale_x_label = lambda x_data: (x_data - x_mean) / x_std
    rescale_x_label = lambda x_data: x_data * x_std + x_mean
    scale_y = lambda y_data: (y_data - y_mean) / y_std
    rescale_y = lambda y_data: y_data * y_std + y_mean

    samples = y_mixed.shape[0]
    if n_samples > samples:
        n_samples = samples
    index = np.random.choice(samples, n_samples, replace = False)  # 重複なし並び替え

    x_img = scale_x_img(x_img)[index]
    x_label = scale_x_label(x_label)[index]
    y_mixed = scale_y(y_mixed)[index]

    x_img_train, x_img_test, x_label_train, x_label_test, y_train, y_test = train_test_split(x_img, x_label, y_mixed, test_size = test_size, random_state = 1)
    return x_img_train, x_img_test, x_label_train, x_label_test, y_train, y_test, rescale_x_img, rescale_x_label, rescale_y
    
def load4gan(mode="FourierConcat", n_samples=20000, test_size=0.33, shuffle=False):
    x_mixed, x_mean, x_std, rescale_x, y_mixed, rescale_y = load_mixed(mode, n_samples, forGan = True, shuffle = shuffle)

    x_as = x_mixed[:, :-2]  # AoA and Shape
    x_as_mean = x_mean[:-2]
    x_as_std = x_std[:-2]

    x_mr = x_mixed[:, -2:]  # Ma, Re
    x_mr_mean = x_mean[-2:]
    x_mr_std = x_std[-2:]
    
    rescale_x_as = lambda x_data: x_data * x_as_std + x_as_mean
    def rescale_x_mr(x_data, invLog=False):
        x_tmp = x_data * x_mr_std + x_mr_mean
        if invLog:
            x_tmp[:, 1] = np.log10(x_tmp[:, 1])
        return x_tmp

    x_train, x_test, x_as_train, x_as_test, x_mr_train, x_mr_test, y_train, y_test = train_test_split(x_mixed, x_as, x_mr, y_mixed, test_size = test_size, random_state = 1)
    return x_train, x_test, x_as_train, x_as_test, x_mr_train, x_mr_test, y_train, y_test, rescale_x, rescale_x_as, rescale_x_mr, rescale_y

def reduceAngle(x_data, rescale_x, updateScaler = False):
    from sklearn.preprocessing import StandardScaler
    x_data = rescale_x(x_data)
    n_sample = x_data.shape[0]
    s_data = x_data[:, 1:].reshape(n_sample, 2, -1)
    a_data = x_data[:, 0] * np.pi / 180.0
    rot = np.array([[np.cos(a_data), -np.sin(a_data)], [np.sin(a_data), np.cos(a_data)]])
    s_data = np.einsum("ijk,jli->ilk", s_data, rot).reshape(n_sample, -1)
    scalar = StandardScaler()
    s_data = scalar.fit_transform(s_data)
    if updateScaler:
        return s_data, scalar.inverse_transform
    else:
        return s_data, rescale_x

def create_optional_information(s_data = None, rescale_s = None, x_data = None, rescale_x = None, train = True, n_feat=2, notUseS=False):
    from sympy import Polygon
    def get_aspect(s_data):
        params = np.zeros((total, n_feat))
        for i in range(total):
            print(i, total)
            r_data = s_data[i, :].reshape(2, -1).T
            pca = PCA(n_components = 2)
            pca.fit(r_data)
            r_data = pca.transform(r_data)
            length = (np.max(r_data, axis = 0) - np.min(r_data, axis = 0))
            params[i, 0] = length[1] / length[0]
            params[i, 1] = np.sum(
                np.sqrt(np.sum(np.diff(r_data, append = r_data[0].reshape(-1, 2), axis = 0) ** 2, axis = 1)))
            if n_feat == 3:
                params[i, 2] = Polygon(*(r_data.tolist())).area
        return params
    
    if s_data is None and x_data is None:
        raise ValueError
    elif x_data is not None:
        total = x_data.shape[0]
    else:
        total = s_data.shape[0]
    
    # n_feat = 3
    params = np.zeros((total, n_feat))
    from pathlib import Path
    case_name = Path("prmF_{0}{1}{2}.npy".format(n_feat, total, train))

    if s_data is None and not notUseS:
        s_data, rescale_s = reduceAngle(x_data=x_data, rescale_x=rescale_x, updateScaler=True)
        s_data = rescale_s(s_data)
        params[:, 0:n_feat] = get_aspect(s_data)
        np.save(str(case_name), params)
        return params, s_data, rescale_s

    if case_name.exists():
        print("calculated data is exists.")
        params = np.load(str(case_name))
    else:
        print("calculated data is not exists.")
        params[:, 0:n_feat] = get_aspect(rescale_s(s_data))
        np.save(str(case_name), params)
    return params

def load4ganMk2(mode = "FourierConcat", n_samples = 50000, test_size = 0.33, n_feat=2):
    x_train, x_test, y_train, y_test, rescale_x, rescale_y = load_mixed(mode = "CircumConcat", n_samples = n_samples,
                                                                        forGan = False, msdf_res = 16,
                                                                        caseName = "CV_mixed",
                                                                        reductTarget = "", split45 = False, shuffle = False)

    p_train = create_optional_information(s_data = None, rescale_s = None, x_data = x_train, rescale_x = rescale_x, train = True, n_feat=n_feat, notUseS=True)
    p_test = create_optional_information(s_data = None, rescale_s = None, x_data = x_test, rescale_x = rescale_x, train = False, n_feat=n_feat, notUseS=True)
    scalar = StandardScaler()
    p_train = scalar.fit_transform(p_train)
    p_test = scalar.transform(p_test)
    rescale_p = scalar.inverse_transform

    x_train, x_test, x_as_train, x_as_test, x_mr_train, x_mr_test, y_train, y_test, rescale_x, rescale_x_as, rescale_x_mr, rescale_y = load4gan(
        mode = mode, n_samples = n_samples, test_size = test_size, shuffle = False)

    def generate_label(x_mr_data, p_data, y_data, rescale_x_mr, rescale_p, rescale_y):
        def rescale_l(l_data):
            l_data[:, 0:2] = rescale_y(l_data[:, 0:2])
            # l_data[:, 2:4] = rescale_y(l_data[:, 2:4])
            l_data[:, 2:4] = rescale_x_mr(l_data[:, 2:4])
            l_data[:, 4:4+n_feat] = rescale_p(l_data[:, 4:4+n_feat])
            return l_data
        
        l_data = np.concatenate([y_data, x_mr_data, p_data], axis=1)
        return l_data, rescale_l
    
    l_train, rescale_l = generate_label(x_mr_train, p_train, y_train, rescale_x_mr, rescale_p, rescale_y)
    l_test, rescale_l = generate_label(x_mr_test,p_test, y_test, rescale_x_mr, rescale_p, rescale_y)
    
    return x_as_train, x_as_test, x_mr_train, x_mr_test, l_train, l_test, rescale_x_as, rescale_x_mr, rescale_l

if __name__ == '__main__':
    # load4gan()
    reductTarget = ""
    x_train, x_test, y_train, y_test, rescale_x, rescale_y = load_mixed(mode = "CircumConcat", n_samples = 20000,
                                                                        forGan = False, msdf_res = 16,
                                                                        caseName = "CV_mixed",
                                                                        reductTarget = reductTarget, split45 = False)
    print(x_train.shape)
    exit()
    for i in range(28):
        target = 500 + 500 * i
        reductTarget = "kmeans_norm_kpca_" + str(target)
        x_train, x_test, y_train, y_test, rescale_x, rescale_y = load_mixed(mode = "FourierConcat", n_samples = 20000, forGan = False, msdf_res = 16, caseName = "CV_mixed",
                   reductTarget = reductTarget, split45 = True)

    print(x_train.shape)
    print(x_test.shape)
    exit()
    x,y, xr, yr = main(shapeType = "FourierConcat", reductTarget = "AoA")
    print(x)
    print(x.shape)
    print(y.shape)
    