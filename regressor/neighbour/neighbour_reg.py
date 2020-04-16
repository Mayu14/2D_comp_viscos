# coding: utf-8
# sample from https://github.com/tksmd/hannari-python-2/blob/master/demo.ipynb
import numpy as np
# from sklearn.preprocessing import PolynomialFeatures
# from sklearn.linear_model import Lasso, Ridge, ElasticNet
from sklearn.neighbors import NearestNeighbors
# from sklearn.pipeline import make_pipeline
# from sklearn.model_selection import GridSearchCV
from sklearn.metrics import r2_score
# from sklearn.decomposition import KernelPCA
# from sklearn.preprocessing import StandardScaler

from training.read_training_data import load_mixed
import pickle
from pathlib import Path

#shape_expression = "FourierConcat"
shape_expression = "CircumConcat"
n_samples_total = 300000
n_feat_dim = 30
# poly_dim = 2
regressor = "neighbour"
dimReductMethod = "PCA"
reductTarget = "kmeans_norm_pca_3600"


def main(reductTarget = "kmeans_norm_pca_300", n_neighbors=16):
    dimReduce = dimReductMethod + str(n_feat_dim)
    X_train, X_test, y_train, y_test, rescale_x, rescale_y = load_mixed(mode = shape_expression,
                                                                        n_samples = n_samples_total,
                                                                        test_size = 0.01,
                                                                        reductTarget = reductTarget,
                                                                        caseName = "CV_mixed",
                                                                        split45 = False,
                                                                        __dimReduct = dimReduce)
    
    x = X_train
    y = y_train
    t = X_test
    # print(poly_dim)
    print(x.shape)
    print(t.shape)
    
    alpha = 0.002  # 0.002
    print(alpha)

    if n_neighbors != 5:
        name = "Nbrs{0}".format(n_neighbors)
    else:
        name = "Neighbour"
    baseName = 'm{0}_f{1}_tr{2}_p{3}_r{4}_a{5}.pkl'.format(name, n_feat_dim, x.shape[0], dimReductMethod,
                                                           reductTarget, alpha)
    
    if not "fourier" in shape_expression.lower():
        baseName = baseName.replace(".pkl", "") + shape_expression + ".pkl"
    print(baseName)
    fname_model = Path("model_{0}".format(baseName))
    
    if fname_model.exists():
        print("load {0}".format(fname_model))
        model = pickle.load(open(fname_model, "rb"))

    else:
        model = NearestNeighbors(n_neighbors=n_neighbors)
        model.fit(x)

        pickle.dump(model, open(fname_model, "wb"))
        
    def predict(x_test):
        weight = 9
        distances, indices = model.kneighbors(x_test)
        test_sample = distances.shape[0]
        test_feat = y_train.shape[1]
        y_test = np.zeros((test_sample, test_feat))
        i = 0
        import warnings
        with warnings.catch_warnings():
            warnings.filterwarnings("error")
            for index, distance in zip(indices, distances):
                try:
                    weight = (1.0 / distance)**(weight)
                except:
                    weight = np.zeros(distance.shape[0])
                    weight[np.argmin(distance)] = 1
                summation = np.dot(weight, y_train[index])
                
                y_test[i, :] = summation / (np.sum(weight))
                i += 1
        return y_test

    
    y_predicted = predict(x)
    t_predicted = predict(t)
    print(t_predicted)
        # r2_train = r2_score(y_train, y_predicted)
        # r2_test = r2_score(y_test, t_predicted)
        
        # print(r2_train, r2_test)
        
    
    
    
    r2_train = r2_score(y_train, y_predicted)
    r2_train_cl = r2_score(y_train[:, 0], y_predicted[:, 0])
    r2_train_cd = r2_score(y_train[:, 1], y_predicted[:, 1])
    r2_test = r2_score(y_test, t_predicted)
    r2_test_cl = r2_score(y_test[:, 0], t_predicted[:, 0])
    r2_test_cd = r2_score(y_test[:, 1], t_predicted[:, 1])
    
    print(r2_train, r2_test)
    with open("log_CF500_p2_500_{0}.txt".format(n_neighbors), "a") as f:
        f.write("{0},{1},{2},{3},{4},mix\n".format(x.shape[0], r2_train, r2_test, reductTarget, baseName))
        f.write("{0},{1},{2},{3},{4},cl\n".format(x.shape[0], r2_train_cl, r2_test_cl, reductTarget, baseName))
        f.write("{0},{1},{2},{3},{4},cd\n".format(x.shape[0], r2_train_cd, r2_test_cd, reductTarget, baseName))

    from regressor.regression import error_print
    with open("knnN9W16Circum9000.csv", "a") as f:
        f.write(error_print(t_predicted[:, 0], y_test[:, 0], csv = True).replace("\n", ","))
        f.write(error_print(t_predicted[:, 1], y_test[:, 1], csv = True).replace("\n", ","))
        f.write("knn,{0}\n".format(x.shape[0]))

    
    if plotFigure:
        y_pred = t_predicted
        from other_tools.scatter_plot import make_scatter_plot
        make_scatter_plot(data_a = y_test[:, 0], data_b = y_pred[:, 0], label_a = "observed $C_L$ (test data)",
                          label_b = "predicted $C_L$",
                          fname = "yyplotkNN", dash_line = "no",
                          max_value = 2.0, min_value = -1.0,
                          title = "yyplot k-NN $C_L$ (train)")
    
        from regressor.regression import relationalError
        rError = relationalError(y_true = y_test, y_pred = y_pred, absoluteValue = True)
        cle = rError[:, 0]
        cde = rError[:, 1]
        import matplotlib.pyplot as plt
        plt.rcParams["font.size"] = 18
        fig = plt.figure()
        ax = fig.add_subplot(111)
        from matplotlib.colors import LogNorm
        H = ax.hist2d(cle, cde, bins = [np.linspace(-0.5, 0.5, 51), np.linspace(-0.5, 0.5, 51)],
                      norm = LogNorm())
        ax.set_title("Error 2D-Histgram :k-NN", fontsize = 22)
        ax.set_xlabel("$C_L$ Error")
        ax.set_ylabel("$C_D$ Error")
    
        fig.colorbar(H[3], ax = ax)
        # plt.show()
        plt.savefig("errHist2D.png")
        plt.close()
        
if __name__ == '__main__':
    np.random.seed(1)
    np.seterr(divide = "raise", invalid="raise")
    plotFigure = True
    head = "kmeans_norm_pca_"
    for i in range(1):
        data = (i+1)*9000#30 + 2 * (i) + 2
    # for i in range(28):
        # for j in range(100):
        #data = 500 * (i + 1)
        # n_neighbors = j + 1
        # if n_neighbors != 5:
        reductTarget = head + str(data)
        main(reductTarget)
