# coding: utf-8
import numpy as np
from training.read_training_data import load_mixed
from sklearn.decomposition import PCA, KernelPCA
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors
from matplotlib import pyplot as plt

np.random.seed(1)
shapeType = "FourierConcat"
reductTarget = "kmeans_norm_pca_9000"
use_kpca = False
not_y = False

def plotxy(x, xlabel, xcmap, y=None, ylabel=None, ycmap=None, yalpha=1.0, kpca=True):
    plt.rcParams["font.size"] = 18
    if kpca:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1, aspect="equal")
    else:
        fig = plt.figure(figsize = (10, 10))
        ax = fig.add_subplot(1, 1, 1)
    if (not not_y) and (not y is None):
        ax.scatter(y[:, 0], y[:, 1], marker="x", c=ycmap, s=10, label=ylabel, alpha=yalpha)
    ax.scatter(x[:, 0], x[:, 1], c=xcmap, s=2, label=xlabel)
    
    if kpca:
        ax.set_title('2D-Projection of Dataset by KPCA (Gauss kernel)', fontsize=22)
        ax.set_xlabel('1st principal component in space induced by $\phi$')
    else:
        ax.set_title('2D-Projection of Dataset by PCA', fontsize = 22)
        ax.set_xlabel('1st principal component')
    ax.set_ylabel('2nd component')
    ax.legend(loc = 'best')
    plt.show()

def getShapeName(shapeType):
    if "fourier" in shapeType.lower():
        return "fourier"
    elif "equidistant" in shapeType.lower():
        return "equidistant"
    elif "crowd" in shapeType.lower():
        return "concentrate"
    else:
        raise ValueError

def main(shapeType="FourierConcat", onlyALL=False, use_kpca=True, split45=True):
    if split45:
        x4, x5, _, _, _, _ = load_mixed(mode = shapeType, n_samples = 50000, reductTarget = reductTarget, split45 = split45)
    else:
        x4, x5, _, _, _, _ = load_mixed(mode = shapeType, n_samples = 50000, reductTarget = reductTarget, split45 = split45)
    if not onlyALL:
        x_AoA_4, _, _, _, _, _ = load_mixed(mode = shapeType, reductTarget = "AoA_14", split45 = True)
        x_NACA4, _, _, _, _, _ = load_mixed(mode = shapeType, reductTarget = "NACA_150", split45 = True)

    if use_kpca:
        kpca = KernelPCA(n_components = 2, kernel = "rbf", n_jobs = -1, random_state = 1)
    else:
        pca = PCA(n_components = 2, random_state = 1)
        ss = StandardScaler()
        kpca = make_pipeline(pca, ss)
    
    shapeName = getShapeName(shapeType)
    # p_x4 = pca.fit_transform(x4)
    kp_x4 = kpca.fit_transform(x4)
    kp_x5 = kpca.transform(x5)
    print(kp_x5)
    cmaps = ["blue", "lightcoral", "green", "purple", "pink", "lime"]
    if split45:
        plotxy(x=kp_x4, xlabel = "NACA4digit_ALL:{0}".format(shapeName), xcmap = cmaps[0], y=kp_x5, ylabel = "Test Data", ycmap = cmaps[1], kpca = use_kpca)
    else:
        plotxy(x = kp_x4, xlabel = "Training Data:{0}".format(shapeName), xcmap = cmaps[0], y = kp_x5,
               ylabel = "Test Data", ycmap = cmaps[1], kpca = use_kpca)
    #exit()
    if not onlyALL:
        # kp_x5 = kpca.transform(x5)
        kp_x_AoA4 = kpca.transform(x_AoA_4)
        plotxy(x = kp_x_AoA4, xlabel = "Reduce_AoA_N=3140:{0}".format(shapeName), xcmap = cmaps[2], y = kp_x4, ylabel = "NACA4digit_ALL", ycmap = "dodgerblue", yalpha=0.01, kpca = use_kpca)
        
        kp_x_NACA4 = kpca.transform(x_NACA4)
        plotxy(x = kp_x_NACA4, xlabel = "Reduce_Shape_N=3585:{0}".format(shapeName), xcmap = cmaps[3], y = kp_x4, ylabel = "NACA4digit_ALL",
               ycmap = "dodgerblue", yalpha=0.01, kpca = use_kpca)
        # x_c7000, _, _, _, _, _ = load_mixed(mode = shapeType, reductTarget = "kmeans_norm_pca_7000", split45 = True)
        x_c3500, _, _, _, _, _ = load_mixed(mode = shapeType, reductTarget = "kmeans_norm_pca_3500", split45 = True)
        # kp_x_c7000 = kpca.transform(x_c7000)
        kp_x_c3500 = kpca.transform(x_c3500)
        plotxy(x = kp_x_c3500, xlabel = "Clustering_N=3500:{0}".format(shapeName), xcmap = cmaps[5], y = kp_x4, ylabel = "NACA4digit_ALL",
               ycmap = "dodgerblue", yalpha = 0.01, kpca = use_kpca)


def cluster_img(shapeType="FourierConcat", split45=False):
    x4, x5, _, _, _, _ = load_mixed(mode = shapeType, n_samples = 50000, reductTarget = "", split45 = split45)
    print(x4.shape, x5.shape)
    xAll = np.concatenate([x4, x5], axis=0)
    print(xAll.shape)
    from sklearn.decomposition import PCA, KernelPCA
    from sklearn.cluster import KMeans
    from sklearn.preprocessing import StandardScaler
    
    xAll = xAll
    
    scalar = StandardScaler()
    pca = PCA(n_components = 2)
    pca2 = PCA(n_components = 2)
    __x = pca.fit_transform(xAll)
    __x = scalar.fit_transform(__x)
    
    kmeans = KMeans(n_clusters = 10,  # クラスタ数
                    init = "k-means++",  # セントロイド初期化アルゴリズム
                    n_init = 10,  #
                    max_iter = 300,  # 最大反復回数
                    tol = 1e-04,  # 収束判定条件
                    n_jobs = -1,  # 並列実行数(-1 = 環境限界まで利用)
                    random_state = 1)
    kmeans.fit(__x)
    y_pred = kmeans.predict(__x)
    
    x2 = pca2.fit_transform(xAll)
    x2 = scalar.fit_transform(x2)

    centroids = kmeans.cluster_centers_
    nearest_neighbours = True
    centroid_plt = False
    if nearest_neighbours:
        nn = NearestNeighbors(metric = "minkowski")
        nn.fit(__x)
        distances, nearest_indices = nn.kneighbors(centroids, n_neighbors = 1)
        x2 = x2[nearest_indices, :].reshape(-1, 2)
        y_pred = kmeans.predict(x2)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(x2[:, 0], x2[:, 1])
    ax.set_title("Data-points", fontsize=22)
    ax.scatter(x2[:, 0], x2[:, 1], c = y_pred)
    # ax.set_title("Data-points split by 10 Clusters")
    ax.set_title("Data-points reducted")
    if centroid_plt and not nearest_neighbours:
        plt.scatter(centroids[:, 0], centroids[:, 1],
                    marker = 'x', s = 169, linewidths = 3,
                    color = 'r', zorder = 10)
    ax.set_xlabel("1st principle component (standardized)")
    ax.set_ylabel("2nd principle component (standardized)")
    plt.show()
    
if __name__ == '__main__':
    cluster_img()
    exit()
    
    shapeTypes = ["FourierConcat", "EquidistantConcat", "CrowdConcat"]
    for shape in shapeTypes:
        main(shapeType = shape, onlyALL = True, use_kpca = False, split45 = False)