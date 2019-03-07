# coding: utf-8
import numpy as np
import os
from sklearn.cluster import KMeans
from sklearn.neighbors import NearestNeighbors


def data_reduction(X_data, y_data, reduction_target = 10000, output_csv = False,
                    fname = "reduct", path = "G:\\Toyota\\Data\\Compressible_Invicid\\training_data\\NACA4\\reduct_csv\\"):
    # X_data:[N_samples, M_features]
    # y_data:[N_samples, K_features]
    def get_nearest_indices(X_data, k_cluster):
        kmeans = KMeans(n_clusters = k_cluster,  # クラスタ数
                        init = "k-means++",  # セントロイド初期化アルゴリズム
                        n_init = 10,  #
                        max_iter = 300,  # 最大反復回数
                        tol = 1e-04,  # 収束判定条件
                        n_jobs = -1,  # 並列実行数(-1 = 環境限界まで利用)
                        random_state = 0)
        
        y_cluster = kmeans.fit_predict(X_data)
        
        nn = NearestNeighbors(metric = "minkowski")
        nn.fit(X_data)
        distances, nearest_indices = nn.kneighbors(kmeans.cluster_centers_, n_neighbors = 1)
        
        return nearest_indices
    
    def reduction2k_datas(X_data, y_data, k_cluster, indices = None, output_indices = False):
        nearest_indices = indices
        if type(indices) == type(None):
            nearest_indices = get_nearest_indices(X_data, k_cluster)
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
    
    csvname = path + fname + str(reduction_target).zfill(5) + ".csv"
    nearest_indices = csv2indices(csvname)
    X_data, y_data, nearest_indices = reduction2k_datas(X_data = X_data, y_data = y_data, k_cluster = reduction_target, indices=nearest_indices, output_indices=output_csv)
    
    if output_csv:
        indices2csv(csvname, nearest_indices)
    
    return X_data, y_data
    
if __name__ == '__main__':
    source = "G:\\Toyota\\Data\\Compressible_Invicid\\training_data\\"
    fname_lift_train = "NACA4\\s1122_e9988_s4_a014.csv"
    fname_shape_train = "NACA4\\shape_fourier_1112_9988_s04.csv"
    from read_training_data_viscos import read_csv_type3

    X_train, y_train, scalar = read_csv_type3(source, fname_lift_train, fname_shape_train, shape_odd = 0, read_rate = 1,
                                              total_data = 0, return_scalar = True)
    
    for i in range(40):
        cluster = 500 * (i + 1)
        data_reduction(X_train, y_train, reduction_target = cluster, output_csv = True)
    
    
    
    
    
    