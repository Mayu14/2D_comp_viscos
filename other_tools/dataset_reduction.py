# coding: utf-8
import numpy as np
import os
from sklearn.cluster import KMeans, MiniBatchKMeans
from sklearn.neighbors import NearestNeighbors
from math import floor

# X_data:[n_samples, m_features]
# y_data:[n_samples, k_features]
# reduction_target: (integer) Number of data after deletion
# output_csv: (bool)
# preprocess: (charcter) Types of preprocessing to be performed before clustering processing (ex."PCA")
# update_xy: (bool) Return preprocessed data
# criteria_method:(character) Selection criteria for extracting data from each cluster
# fname: (character) name of csv file
# path:  (character) name of path
def data_reduction(X_data, y_data, reduction_target = 10000, output_csv = False, preprocess = "None", update_xy = False,
                   criteria_method = "nearest_centroid",
                    fname = "reduct", path = "G:\\Toyota\\Data\\Compressible_Invicid\\training_data\\NACA4\\reduct_csv\\"):
    # X_data:[N_samples, M_features]
    # y_data:[N_samples, K_features]
    def get_nearest_indices(X_data, k_cluster, criteria_method):
        threshold = 3000
        if k_cluster < threshold:
            kmeans = KMeans(n_clusters = k_cluster,  # クラスタ数
                            init = "k-means++",  # セントロイド初期化アルゴリズム
                            n_init = 10,  #
                            max_iter = 300,  # 最大反復回数
                            tol = 1e-04,  # 収束判定条件
                            n_jobs = -1,  # 並列実行数(-1 = 環境限界まで利用)
                            random_state = 0)
        else:
            split = floor(float(k_cluster) / threshold)
            batch_size = floor(float(k_cluster) / split)
            kmeans = MiniBatchKMeans(n_clusters = k_cluster,  # クラスタ数
                                     batch_size = batch_size,
                                     init = "k-means++",  # セントロイド初期化アルゴリズム
                                     n_init = 10,  #
                                     max_iter = 300,  # 最大反復回数
                                     tol = 1e-04,  # 収束判定条件
                                     random_state = 0)
            
        
        y_cluster = kmeans.fit_predict(X_data)
        
        if criteria_method == "nearest_centroid":
            nn = NearestNeighbors(metric = "minkowski")
            nn.fit(X_data)
            distances, nearest_indices = nn.kneighbors(kmeans.cluster_centers_, n_neighbors = 1)
        elif criteria_method == "farthest_from_center":
            center = np.sum(X_train, axis = 0) / X_train.shape[0]
            nearest_indices = np.zeros(k_cluster, dtype = int)
            for i in range(k_cluster):
                # dif = X_train[y_cluster == i, :] - center
                # norm = np.linalg.norm(dif, axis = 1)
                # max_index = np.argmax(norm)
                nearest_indices[i] = np.argmax(np.linalg.norm((X_train[y_cluster == i, :] - center), axis = 1))
        elif criteria_method == "triangle_only":
            exit()
        elif criteria_method == "triangle_with_center":
            exit()
        
        return nearest_indices
    
    # output_indicesは常時Trueに変更予定
    def reduction2k_datas(X_data, y_data, k_cluster, indices = None, output_indices = True, criteria_method = "nearest_centroid"):
        nearest_indices = indices
        if type(indices) == type(None):
            nearest_indices = get_nearest_indices(X_data, k_cluster, criteria_method)
        X_data = X_data[nearest_indices, :]
        y_data = y_data[nearest_indices, :]
        if output_indices:
            return X_data, y_data, nearest_indices
        else:
            return X_data, y_data, None
    
    def indices2csv(csvname, nearest_indices):
        np.savetxt(csvname, nearest_indices, delimiter = ",")
    
    def csv2indices(csvname):
        if os.path.exists(csvname):
            return np.loadtxt(csvname, delimiter = ",", dtype = "int")
        else:
            return None

    def preprocessing(X_data, preprocess):
        def pre_pca(X_data):
            from sklearn.decomposition import PCA
            pca = PCA(n_components = X_data.shape[1])   # 全成分残す
            pca.fit(X_data)
            return pca.fit_transform(X_data)
        
        def pre_kpca(X_data, kernel="rbf", gamma = 20.0):
            from sklearn.decomposition import KernelPCA
            kpca = KernelPCA(n_components = X_data.shape[1], kernel = kernel, gamma = gamma)
            kpca.fit(X_data)
            return kpca.fit_transform(X_data)

        if preprocess == "None":
            pass
        elif preprocess == "PCA":
            X_data = pre_pca(X_data)
        elif (preprocess == "rbf") or (preprocess == "poly") or (preprocess == "linear") or (
                preprocess == "sigmoid") or (preprocess == "cosine"):
            X_data = pre_kpca(X_data, preprocess)
        else:
            print("Error! input priprocess")
            exit()
        return X_data
    
    def get_fname(preprocess, criteria_method):
        if preprocess == "None":
            mid = ""
        elif preprocess == "PCA":
            mid = "_PCA_"
        elif (preprocess == "rbf") or (preprocess == "poly") or (preprocess == "linear") or (
                preprocess == "sigmoid") or (preprocess == "cosine"):
            mid = "_" + preprocess + "_"
        else:
            print("Error! input preprocess")
            exit()
        if criteria_method == "nearest_centroid":
            pass
        elif criteria_method == "farthest_from_center":
            mid += "FFC_"
        elif criteria_method == "triangle_only":
            mid += "Tri_"
        elif criteria_method == "triangle_with_center":
            mid += "TwC_"
        else:
            print('please check value of "criteria_method"')
            exit()
        return mid
        
    header = path + fname
    footer = str(reduction_target).zfill(5) + ".csv"
    
    if ((update_xy == False) and (preprocess != "None")):
        X_origin = X_data
        y_origin = y_data
    
    mid = get_fname(preprocess, criteria_method)
    csvname = header + mid + footer
    nearest_indices = csv2indices(csvname)
    if type(nearest_indices) == type(None):
        X_data = preprocessing(X_data, preprocess)
    
    X_data, y_data, nearest_indices = reduction2k_datas(X_data = X_data, y_data = y_data, k_cluster = reduction_target, indices=nearest_indices)
    
    if output_csv:
        indices2csv(csvname, nearest_indices)
    
    if update_xy:
        return X_data, y_data
    else:
        return X_origin[nearest_indices, :], y_origin[nearest_indices, :]
    
if __name__ == '__main__':
    source = "G:\\Toyota\\Data\\Compressible_Invicid\\training_data\\"
    fname_lift_train = "NACA4\\s1122_e9988_s4_a014.csv"
    fname_shape_train = "NACA4\\shape_fourier_1112_9988_s04.csv"
    from read_training_data_viscos import read_csv_type3

    X_train, y_train, scalar = read_csv_type3(source, fname_lift_train, fname_shape_train, shape_odd = 0, read_rate = 1,
                                              total_data = 0, return_scalar = True)

    for i in range(40):
        cluster = 500 * (i + 1)
        print(cluster)
        data_reduction(X_train, y_train, preprocess = "rbf" ,reduction_target = cluster, output_csv = False, criteria_method = "farthest_from_center")
        
    
    
    
    
    