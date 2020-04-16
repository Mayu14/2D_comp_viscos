# -*- coding: utf-8 -*- %reset -f
"""
@author: Hiromasa Kaneko
Fast optimization of SVR hyperparameters with Gaussian kernel
"""

import time

import matplotlib.figure as figure
import matplotlib.pyplot as plt
import numpy as np
# import pandas as pd
from sklearn import model_selection, svm, datasets
from sklearn.model_selection import train_test_split, GridSearchCV
# from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.decomposition import PCA
import pickle
import os

from other_tools.scatter_plot import make_scatter_plot

def svr_main(Xtrain, ytrain, Xtest, ytest, fname, n = 2):
    if n != 0:
        pca = PCA(n_components = n)
        Xtrain = np.concatenate([Xtrain[:, :3], pca.fit_transform(Xtrain[:, 3:])], axis = 1)
        Xtest = np.concatenate([Xtest[:, :3], pca.transform(Xtest[:, 3:])], axis = 1)
    # Standarize X and y    # 不偏分散
    autoscaled_Xtrain = (Xtrain - Xtrain.mean(axis = 0)) / Xtrain.std(axis = 0, ddof = 1)
    autoscaled_ytrain = (ytrain - ytrain.mean()) / ytrain.std(ddof = 1)
    autoscaled_Xtest = (Xtest - Xtrain.mean(axis = 0)) / Xtrain.std(axis = 0, ddof = 1)
    # autoscaled_Xtrain = pca(autoscaled_Xtrain, n)
    # autoscaled_Xtest = pca(autoscaled_Xtest, n)
    
    # Measure time in hyperparameter optimization
    start_time = time.time()
    
    fold_number = 5  # "fold_number"-fold cross-validation
    if os.path.exists(fname + "model.sav"):
        regression_model = pickle.load(open(fname + "model.sav", "rb"))
    else:
        # Settings
        svr_cs = 2 ** np.arange(-5, 11, dtype = float)  # Candidates of C
        svr_epsilons = 2 ** np.arange(-10, 1, dtype = float)  # Candidates of epsilon
        svr_gammas = 2 ** np.arange(-20, 11, dtype = float)  # Candidates of gamma
        
        """
        number_of_training_samples = 1000
        number_of_test_samples = 1000

        # Generate samples for demonstration
        X, y = datasets.make_regression(n_samples = number_of_training_samples + number_of_test_samples, n_features = 100,
                                        n_informative = 100, noise = 100, random_state = 0)

        # Divide samples into training samples and test samples
        Xtrain, Xtest, ytrain, ytest = train_test_split(X, y, test_size = number_of_test_samples, random_state = 0)
        """
        
        # Optimize gamma by maximizing variance in Gram matrix
        numpy_autoscaled_Xtrain = np.array(autoscaled_Xtrain)
        iter = 10
        estimated_gamma = []
        size = 1000
        if autoscaled_Xtrain.shape[0] < size:
            size = autoscaled_Xtrain.shape[0]
        for i in range(iter):
            sample_X = numpy_autoscaled_Xtrain[
                       np.random.choice(numpy_autoscaled_Xtrain.shape[0], size = size, replace = False), :]
            variance_of_gram_matrix = list()
            for svr_gamma in svr_gammas:
                """
                gram_matrix = np.exp(
                    -svr_gamma * ((numpy_autoscaled_Xtrain[:, np.newaxis] - numpy_autoscaled_Xtrain) ** 2).sum(axis=2))
                """
                gram_matrix = np.exp(
                    -svr_gamma * ((sample_X[:, np.newaxis] - sample_X) ** 2).sum(axis = 2))
                variance_of_gram_matrix.append(gram_matrix.var(ddof = 1))
            estimated_gamma.append(
                svr_gammas[np.where(variance_of_gram_matrix == np.max(variance_of_gram_matrix))[0][0]])
        optimal_svr_gamma = np.average(np.array(estimated_gamma))
        # print(estimated_gamma)
        # print(optimal_svr_gamma)
        
        # Optimize epsilon with cross-validation
        svr_model_in_cv = GridSearchCV(svm.SVR(kernel = 'rbf', C = 3, gamma = optimal_svr_gamma),
                                       {'epsilon': svr_epsilons},
                                       cv = fold_number)
        svr_model_in_cv.fit(autoscaled_Xtrain, autoscaled_ytrain)
        optimal_svr_epsilon = svr_model_in_cv.best_params_['epsilon']
        
        # Optimize C with cross-validation
        svr_model_in_cv = GridSearchCV(
            svm.SVR(kernel = 'rbf', epsilon = optimal_svr_epsilon, gamma = optimal_svr_gamma),
            {'C': svr_cs}, cv = fold_number)
        svr_model_in_cv.fit(autoscaled_Xtrain, autoscaled_ytrain)
        optimal_svr_c = svr_model_in_cv.best_params_['C']
        
        # Optimize gamma with cross-validation (optional)
        svr_model_in_cv = GridSearchCV(svm.SVR(kernel = 'rbf', epsilon = optimal_svr_epsilon, C = optimal_svr_c),
                                       {'gamma': svr_gammas}, cv = fold_number)
        svr_model_in_cv.fit(autoscaled_Xtrain, autoscaled_ytrain)
        optimal_svr_gamma = svr_model_in_cv.best_params_['gamma']
        
        # Check time in hyperparameter optimization
        elapsed_time = time.time() - start_time
        print("Elapsed time in hyperparameter optimization: {0} [sec]".format(elapsed_time))
        
        # Check optimized hyperparameters
        print("C: {0}, Epsion: {1}, Gamma: {2}".format(optimal_svr_c, optimal_svr_epsilon, optimal_svr_gamma))
        
        # Construct SVR model
        regression_model = svm.SVR(kernel = 'rbf', C = optimal_svr_c, epsilon = optimal_svr_epsilon,
                                   gamma = optimal_svr_gamma)
        regression_model.fit(autoscaled_Xtrain, autoscaled_ytrain)
        
        pickle.dump(regression_model, open(fname + "model.sav", "wb"))
    
    # Calculate y of training dataset
    calculated_ytrain = np.ndarray.flatten(regression_model.predict(autoscaled_Xtrain))
    calculated_ytrain = calculated_ytrain * ytrain.std(ddof = 1) + ytrain.mean()
    # r2, RMSE, MAE
    r2 = float(1 - sum((ytrain - calculated_ytrain) ** 2) / sum((ytrain - ytrain.mean()) ** 2))
    print("r2: {0}".format(float(1 - sum((ytrain - calculated_ytrain) ** 2) / sum((ytrain - ytrain.mean()) ** 2))))
    print("RMSE: {0}".format(float((sum((ytrain - calculated_ytrain) ** 2) / len(ytrain)) ** 0.5)))
    print("MAE: {0}".format(float(sum(abs(ytrain - calculated_ytrain)) / len(ytrain))))
    # yy-plot
    plt.figure(figsize = figure.figaspect(1))
    plt.scatter(ytrain, calculated_ytrain)
    YMax = np.max(np.array([np.array(ytrain), calculated_ytrain]))
    YMin = np.min(np.array([np.array(ytrain), calculated_ytrain]))
    plt.plot([YMin - 0.05 * (YMax - YMin), YMax + 0.05 * (YMax - YMin)],
             [YMin - 0.05 * (YMax - YMin), YMax + 0.05 * (YMax - YMin)], 'k-')
    plt.ylim(YMin - 0.05 * (YMax - YMin), YMax + 0.05 * (YMax - YMin))
    plt.xlim(YMin - 0.05 * (YMax - YMin), YMax + 0.05 * (YMax - YMin))
    plt.xlabel("Actual Y")
    plt.ylabel("Calculated Y")
    # plt.show()
    plt.savefig(fname + "_yy_train.png")
    plt.close()
    
    # Estimate y in cross-validation
    estimated_y_in_cv = np.ndarray.flatten(
        model_selection.cross_val_predict(regression_model, autoscaled_Xtrain, autoscaled_ytrain, cv = fold_number))
    estimated_y_in_cv = estimated_y_in_cv * ytrain.std(ddof = 1) + ytrain.mean()
    # r2cv, RMSEcv, MAEcv
    r2cv = float(1 - sum((ytrain - estimated_y_in_cv) ** 2) / sum((ytrain - ytrain.mean()) ** 2))
    print("r2cv: {0}".format(float(1 - sum((ytrain - estimated_y_in_cv) ** 2) / sum((ytrain - ytrain.mean()) ** 2))))
    print("RMSEcv: {0}".format(float((sum((ytrain - estimated_y_in_cv) ** 2) / len(ytrain)) ** 0.5)))
    print("MAEcv: {0}".format(float(sum(abs(ytrain - estimated_y_in_cv)) / len(ytrain))))
    # yy-plot
    plt.figure(figsize = figure.figaspect(1))
    plt.scatter(ytrain, estimated_y_in_cv)
    YMax = np.max(np.array([np.array(ytrain), estimated_y_in_cv]))
    YMin = np.min(np.array([np.array(ytrain), estimated_y_in_cv]))
    plt.plot([YMin - 0.05 * (YMax - YMin), YMax + 0.05 * (YMax - YMin)],
             [YMin - 0.05 * (YMax - YMin), YMax + 0.05 * (YMax - YMin)], 'k-')
    plt.ylim(YMin - 0.05 * (YMax - YMin), YMax + 0.05 * (YMax - YMin))
    plt.xlim(YMin - 0.05 * (YMax - YMin), YMax + 0.05 * (YMax - YMin))
    plt.xlabel("Actual Y")
    plt.ylabel("Estimated Y in CV")
    plt.savefig(fname + "_yy_cv.png")
    plt.close()
    
    # Estimate y of test dataset
    predicted_ytest = np.ndarray.flatten(regression_model.predict(autoscaled_Xtest))
    predicted_ytest = predicted_ytest * ytrain.std(ddof = 1) + ytrain.mean()
    # r2p, RMSEp, MAEp
    r2p = float(1 - sum((ytest - predicted_ytest) ** 2) / sum((ytest - ytest.mean()) ** 2))
    print("r2p: {0}".format(float(1 - sum((ytest - predicted_ytest) ** 2) / sum((ytest - ytest.mean()) ** 2))))
    print("RMSEp: {0}".format(float((sum((ytest - predicted_ytest) ** 2) / len(ytest)) ** 0.5)))
    print("MAEp: {0}".format(float(sum(abs(ytest - predicted_ytest)) / len(ytest))))
    # yy-plot
    plt.figure(figsize = figure.figaspect(1))
    plt.scatter(ytest, predicted_ytest)
    YMax = np.max(np.array([np.array(ytest), predicted_ytest]))
    YMin = np.min(np.array([np.array(ytest), predicted_ytest]))
    plt.plot([YMin - 0.05 * (YMax - YMin), YMax + 0.05 * (YMax - YMin)],
             [YMin - 0.05 * (YMax - YMin), YMax + 0.05 * (YMax - YMin)], 'k-')
    plt.ylim(YMin - 0.05 * (YMax - YMin), YMax + 0.05 * (YMax - YMin))
    plt.xlim(YMin - 0.05 * (YMax - YMin), YMax + 0.05 * (YMax - YMin))
    plt.xlabel("Actual Y")
    plt.ylabel("Estimated Y")
    plt.savefig(fname + "_yy_test.png")
    plt.close()
    make_scatter_plot(data_a = ytrain, data_b = calculated_ytrain, label_a = "observed (train data)",
                      label_b = "predicted",
                      fname = fname + "train_svr")
    make_scatter_plot(data_a = ytest, data_b = predicted_ytest, label_a = "observed (test data)", label_b = "predicted",
                      fname = fname + "test_svr")
    with open("log_data_fm50_kpca_feat.csv", "a") as f:
        f.write("{0},{1},{2},{3}\n".format(Xtrain.shape[0], r2, r2cv, r2p))


def pca(x, n = 2):
    pca = PCA(n_components = n)
    pca.fit(x)
    return pca.fit_transform(x)


def main():
    from training.read_training_data import load_mixed
    split45 = False
    shape_expressions = ["CircumConcat"]#["FourierConcat"]#, "EquidistantConcat", "CrowdConcat"]
    n_samples_total = 300000
    n_feat_dim = 20#12#30
    poly_dim = 5

    dimReductMethod = "kpca"  # "PCA"
    reductTargetHeader = "kmeans_norm_pca_"
    print(dimReductMethod)
    for n_feat_dim in range(5,30):
        dimReduce = dimReductMethod + str(n_feat_dim)
    # for i in range(235):
        i = 9
        for shape_expression in shape_expressions:
            # reductTarget = reductTargetHeader + str(300 * (i+1))
            # reductTarget = reductTargetHeader + str(30 + 2 * i + 2)
            reductTarget = reductTargetHeader + str(50 * (i + 1))
            X_train, X_test, y_train, y_test, rescale_x, rescale_y = load_mixed(mode = shape_expression,
                                                                                n_samples = n_samples_total,
                                                                                test_size = 0.001,
                                                                                reductTarget = reductTarget,
                                                                                caseName = "CV_mixed",
                                                                                split45 = split45,
                                                                                __dimReduct = dimReduce)
    
            print(X_train.shape, X_test.shape)
    
        #for n in range(1,60):
            n = 0   # use dimReduce
        
            print("\n\nn={0}".format(n))
            case_name = "_{0}_n{1}".format(str(300 * (i+1)).zfill(4), n_feat_dim)
            if not "fourier" in shape_expression:
                case_name = shape_expression + case_name
            if split45:
                case_name += "Split"
            
            svr_main(Xtrain = X_train, ytrain = y_train[:, 0].flatten(), Xtest = X_test, ytest = y_test[:, 0].flatten(),
                     fname = "cv_CL{0}".format(case_name), n = n)
            svr_main(Xtrain=X_train, ytrain=y_train[:, 1].flatten(), Xtest=X_test, ytest=y_test[:, 1].flatten(),
                     fname="cv_CD{0}".format(case_name), n=n)


if __name__ == '__main__':
    main()