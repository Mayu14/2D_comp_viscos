# coding: utf-8
import numpy as np
import os
from sklearn.cluster import KMeans, MiniBatchKMeans
from sklearn.neighbors import NearestNeighbors
from sklearn.mixture import GaussianMixture as GMM
from sklearn.metrics import explained_variance_score
from sklearn.decomposition import KernelPCA
# import hdbscan
from math import floor
from copy import deepcopy
import matplotlib.pyplot as plt

# X_data:[n_samples, m_features]
# y_data:[n_samples, k_features]
# reduction_target: (integer) Number of data after deletion for kmeans, hdbscan cannnot be specified
# output_csv: (bool)
# preprocess: (charcter) Types of preprocessing to be performed before clustering processing (ex."PCA")
# update_xy: (bool) Return preprocessed data
# criteria_method:(character) Selection criteria for extracting data from each cluster
# main_process: (character) set clustering algorithm
# cv_types: (character) define Covariance Matrix type (only GMM)
# check_plot : (bool) plot dataset on 2D-field after clustering
# force_overwrite: (bool) overwrite existed csv data
# fname: (character) name of csv file
# path:  (character) name of path
def data_reduction(X_data, y_data, reduction_target = 10000, output_csv = False, preprocess = "None", update_xy = False,
                   criteria_method = "nearest_centroid", main_process = "kmeans++", cv_types="all", check_plot = False, force_overwrite = False,
                    fname = "reduct", path = "G:\\Toyota\\Data\\Compressible_Invicid\\training_data\\NACA4\\reduct_csv\\"):
    # X_data:[N_samples, M_features]
    # y_data:[N_samples, K_features]
    def get_nearest_indices(X_data, k_cluster, criteria_method, main_process):
        # main_process
        def main_kmeans(k_cluster, X_data):
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
                                         init_size = k_cluster + batch_size,
                                         max_iter = 300,  # 最大反復回数
                                         tol = 1e-04,  # 収束判定条件
                                         random_state = 0)

            return kmeans.fit_predict(X_data), kmeans.cluster_centers_
            
        # cvにはcv_typesのうちいくつかをリストで入力する．
        def main_gmm(k_cluster, X_data, cv = "all"):
            bic = []
            lowest_bic = np.infty
            if cv == "all":
                cv_types = ['spherical', 'tied', 'diag', 'full']
            else:
                cv_types = [cv]
            for cv_type in cv_types:
                gmm = GMM(n_components = k_cluster, covariance_type = cv_type, warm_start=True)
                gmm.fit(X_data)
                bic.append(gmm.bic(X_data))
                if bic[-1] < lowest_bic:
                    lowest_bic = bic[-1]
                    best_gmm = gmm
            
            print(bic)
            return best_gmm.predict(X_data), best_gmm.means_
        
        if main_process == "kmeans++":
            y_cluster, centroids = main_kmeans(k_cluster, X_data)
        
        elif main_process == "gmm":
            y_cluster, centroids = main_gmm(k_cluster, X_data, cv_types)
        
        else:
            print("process name error!, please check variable 'main process'")
            exit()
        print("main processing is finished.")
        
        # post_process
        def nearest_centroid(X_data, centroids, metric="minkowski"):
            nn = NearestNeighbors(metric = metric)
            nn.fit(X_data)
            distances, nearest_indices = nn.kneighbors(centroids, n_neighbors = 1)
            return nearest_indices
        
        def farthest_from_center(X_data, y_cluster, center):
            nearest_indices = np.zeros(k_cluster, dtype = int)
            print(y_cluster, y_cluster.shape, k_cluster)
            
            for i in range(k_cluster):
                dif = X_data[y_cluster == i, :] - center
                norm = np.linalg.norm(dif, axis = 1)
                max_index = np.argmax(norm)
                nearest_indices[i] = max_index
                # nearest_indices[i] = np.argmax(np.linalg.norm((X_data[y_cluster == i, :] - center), axis = 1))
            return nearest_indices
        
        print(preprocess, main_process, criteria_method)
        if criteria_method == "nearest_centroid":
            nearest_indices = nearest_centroid(X_data, centroids)
                
        elif criteria_method == "farthest_from_center":
            center = np.sum(X_data, axis = 0) / X_data.shape[0]
            nearest_indices = farthest_from_center(X_data, y_cluster, center)
            
        elif criteria_method == "triangle_only":
            exit()
        elif criteria_method == "triangle_with_center":
            exit()
        
        return nearest_indices
    
    # output_indicesは常時Trueに変更予定
    def reduction2k_datas(X_data, y_data, k_cluster, indices = None, output_indices = True, criteria_method = "nearest_centroid", main_process="kmeans++"):
        nearest_indices = indices
        if type(indices) == type(None):
            nearest_indices = get_nearest_indices(X_data, k_cluster, criteria_method, main_process)
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

    def preprocessing(X_data, preprocess, force_overwrite=False):
        def pre_pca(X_data):
            from sklearn.decomposition import PCA
            original_dimension = X_data.shape[1]
            pca = PCA(n_components = original_dimension)   # 全成分残す
            pca.fit(X_data)
            cumulative_contribution_ratio = 0.0
            new_dimension = 1
            for i in range(original_dimension):
                cumulative_contribution_ratio += pca.explained_variance_ratio_[i]
                if cumulative_contribution_ratio > 0.99:
                    new_dimension = i + 1
                    break
            pca = PCA(n_components = new_dimension)
            """
            pca.fit(X_data)
            print(pca.explained_variance_ratio_)
            X_pca = pca.transform(X_data)
            print(np.sum(pca.explained_variance_ratio_))
            print(explained_variance_score(X_data, pca.inverse_transform(X_pca)))
            exit()
            """
            return pca.fit_transform(X_data)
            
        
        def pre_kpca(X_data, kernel="rbf", gamma = 20.0):
            # for i in range(2, X_data.shape[1], 3):
            components = X_data.shape[1]
            if kernel=="rbf":
                gamma = 2.0 ** (-6)
            elif kernel=="sigmoid":
                gamma = 0.013
            elif kernel=="poly":
                gamma = 0.053
                # degree = 3    # default
            elif kernel == "linear":
                components -= 1
                # kpca = KernelPCA(n_components = X_data.shape[1], kernel = kernel, gamma = gamma, fit_inverse_transform = True)
                # kpca = KernelPCA(n_components = X_data.shape[1], kernel = kernel, gamma = gamma, fit_inverse_transform = True, n_jobs = -1)
            kpca = KernelPCA(n_components = components, kernel = kernel, gamma = gamma, n_jobs = -1)
                # kpca.fit(X_data)
                # X_kpca = kpca.transform(X_data)
            # kpca = KernelPCA(n_components = X_data.shape[1], kernel = kernel, gamma = gamma)
            return kpca.fit_transform(X_data)

        save_path = "G:\\Toyota\\Data\\Compressible_Invicid\\training_data\\NACA4\\reduct_csv\\"
        csvname = save_path + preprocess + ".csv"
        if os.path.exists(csvname) and (force_overwrite == False):
            X_data = np.loadtxt(csvname, delimiter = ",", dtype = "float")
        else:
            if preprocess == "None":
                pass
            elif preprocess == "PCA":
                X_data = pre_pca(X_data)
            elif (preprocess == "rbf") or (preprocess == "poly") or (preprocess == "linear") or (
                    preprocess == "sigmoid") or (preprocess == "cosine"):
                gamma = 1.0 / X_data.shape[0]
                X_data = pre_kpca(X_data, preprocess, gamma)
            else:
                print("Error! input priprocess")
                exit()
            np.savetxt(csvname, X_data, delimiter = ",")
        return X_data
    
    def get_fname(preprocess, criteria_method, main_process, cv_types):
        mid = "Pre"
        if (preprocess == "None") or (preprocess == "PCA") or (preprocess == "rbf") or (preprocess == "poly") or (
                preprocess == "linear") or (preprocess == "sigmoid") or (preprocess == "cosine"):
            mid += "_" + preprocess + "_"
        else:
            print("Error! input preprocess")
            exit()

        mid += "Main_"
        
        if (main_process == "kmeans++"):
            mid += main_process
        elif (main_process == "gmm"):
            mid += main_process
            mid += "cv_" + cv_types
        else:
            print('please check value of "main_process"')
            exit()

        mid += "_Post_"
        if criteria_method == "nearest_centroid":
            mid += "NeC_"
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
    
    # if ((update_xy == False) or (preprocess != "None")):
    X_origin = deepcopy(X_data)
    y_origin = deepcopy(y_data)

    mid = get_fname(preprocess, criteria_method, main_process, cv_types)
    print(mid + footer)
    csvname = header + mid + footer
    
    if force_overwrite:
        nearest_indices = None
    else:
        nearest_indices = csv2indices(csvname)
        
    if type(nearest_indices) == type(None):
        X_data = preprocessing(X_data, preprocess)
        print("preprocessing is finished.")
        print(X_data)
    X_data, y_data, nearest_indices = reduction2k_datas(X_data = X_data, y_data = y_data, k_cluster = reduction_target, criteria_method = criteria_method, indices=nearest_indices, main_process=main_process)


    if output_csv:
        indices2csv(csvname, nearest_indices)

    def eval_of_variance(X):
        """
        :param X: ndarray(n_components, n_features)
        :return: trace of covariance matrix
        """
        return np.trace(np.cov(X.T))

    #print(X_data)
    #print(X_data.shape)
    #print(X_origin.shape)

    new_var = eval_of_variance(X_data)
    org_var = eval_of_variance(X_origin)
    print(reduction_target, new_var, org_var, new_var / org_var)

    def make_plot(nearest_indices, X_origin, y_origin, preprocess, criteria_method, main_process, fname):
        from sklearn.decomposition import PCA
        import matplotlib.pyplot as plt

        plot_pca = PCA(n_components = 2)
        tr_origin = plot_pca.fit_transform(X_origin)
        X_data4plot, y_data4plot, indices = reduction2k_datas(X_data = X_origin, y_data = y_origin, k_cluster = reduction_target, output_indices = False, criteria_method = criteria_method, indices=nearest_indices, main_process=main_process)
        transformed = plot_pca.transform(X_data4plot)

        for i in range(2):
            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            ax.set_xlim(np.min(tr_origin[:, 0]), np.max(tr_origin[:, 0]))
            ax.set_xlabel("1st principle component")
            ax.set_ylim(np.min(tr_origin[:, 1]), np.max(tr_origin[:, 1]))
            ax.set_ylabel("2nd principle component")
            ax.set_title("reduced by " + preprocess + " + " + criteria_method + " + " + main_process)
    
            ax.plot(transformed[:, 0], transformed[:, 1], color = "limegreen", marker = ".", linestyle = "None", label = "reduced")
            tail = ""
            if i == 1:
                ax.plot(tr_origin[:, 0], tr_origin[:, 1], color = "dodgerblue", marker = ".", linestyle = "None", alpha = 0.1,
                        label = "origin")
                tail += "with_original"
                
            ax.legend()
            
            plt.savefig(fname + tail + ".png")
            plt.close()

    if check_plot:
        make_plot(nearest_indices, X_origin, y_origin, preprocess, criteria_method, main_process, csvname)

    if update_xy:
        return X_data, y_data
    else:
        return X_origin[nearest_indices, :], y_origin[nearest_indices, :]

# def evaluate_gmm():


if __name__ == '__main__':
    source = "G:\\Toyota\\Data\\Compressible_Invicid\\training_data\\"
    fname_lift_train = "NACA4\\s1122_e9988_s4_a014.csv"
    fname_shape_train = "NACA4\\shape_fourier_1112_9988_s04.csv"
    from read_training_data_viscos import read_csv_type3

    X_train, y_train, scalar = read_csv_type3(source, fname_lift_train, fname_shape_train, shape_odd = 0, read_rate = 1,
                                              total_data = 0, return_scalar = True)

    main_processes = ["gmm"]#kmeans++"]#, "kmeans++"]
    # main_process = "kmeans++"
    preprocesses = ["linear"]#, "PCA", "rbf", "poly", "linear", "cosine", "sigmoid"]
    # preprocesses = ["None"]
    postprocesses = ["nearest_centroid"]#, "farthest_from_center"]
    
    for preproc in preprocesses:
        for mainproc in main_processes:
            for postproc in postprocesses:
                for i in range(40, 0, -1):
                    cluster = 500 * (i + 1)
                    data_reduction(X_train, y_train, preprocess = preproc ,reduction_target = cluster, output_csv = True, main_process = mainproc, criteria_method = postproc, cv_types = "tied", check_plot = False, force_overwrite=False)
                
    
    
    
    
    