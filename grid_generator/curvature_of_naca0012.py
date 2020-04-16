# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

def y_t(x, t=0.12):
    return t / 0.2 * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2 + 0.2843 * x ** 3 - 0.1015 * x ** 4)
    
def y_t_prime(x, t=0.12):
    return t / 0.2 * (0.2969 * (0.5 / np.sqrt(x)) - 0.1260 - 0.3516 * (2 * x) + 0.2843 * (3 * x ** 2) - 0.1015 * (4 * x ** 3))

def y_t_doubleprime(x, t=0.12):
    return t / 0.2 * (0.2969 * (-0.25 / np.sqrt(x)**3) - 0.3516 * (2) + 0.2843 * (6 * x) - 0.1015 * (12 * x ** 2))

def curvature_radius(x, t=0.12):
    return (1 + y_t_prime(x, t))**(3/2) / (y_t_doubleprime(x, t))

def curvature(x, t=0.12, abs = True):
    if abs:
        return np.abs(1.0 / curvature_radius(x, t))
    else:
        return 1.0 / curvature_radius(x, t)

def get_mean_square_error(num_invR, ex_invR):
    return np.sqrt((num_invR - ex_invR)**2)

def get_error(num_invR, ex_invR):    # ex(R) = num(R) + eps(R)
    return ex_invR - num_invR

def fitting(num_invR, ex_invR):
    x_data = num_invR
    y_data = ex_invR
    func1 = lambda x, a, b, c, d, e, f: a + b*x
    func2 = lambda x, a, b, c, d, e, f: a + b*x + c*x**2
    func3 = lambda x, a, b, c, d, e, f: a + b * x + c * x ** 2 + d * x**3
    func4 = lambda x, a, b, c, d, e, f: a + b * x + c * x ** 2 + d * x ** 3 + e*x**4
    func5 = lambda x, a, b, c, d, e, f: a + b * x + c * x ** 2 + d * x ** 3 + e*x**4 + f*np.sqrt(x)
    func6 = lambda x, a, b, c, d, e, f: a + b * x + c * x ** 2 + d * np.exp(x) + e * np.log(x) + f * np.sqrt(x)
    func7 = lambda x, a, b, c, d, e, f: a + b * x + c * x ** 2 + d * np.exp(e*x)
    func8 = lambda x, a, b, c, d, e, f: a + b * x + c * x ** 2 + d * x ** 3 + e*x**4 + f*x**5
    func9 = lambda x, a, b, c, d, e, f: a + b * x + c * x ** 2 + d * np.exp(x)
    
    func = func7
    param_init = np.zeros(6)
    param_opt, cov = scipy.optimize.curve_fit(func, x_data, y_data, p0=param_init)
    print(param_opt)
    # y = func2(x_data, param_opt[0], param_opt[1], param_opt[2])
    # y = func3(x_data, param_opt[0], param_opt[1], param_opt[2], param_opt[3])
    y = func(x_data, param_opt[0], param_opt[1], param_opt[2], param_opt[3], param_opt[4], param_opt[5])
    return y, param_opt
    
def main():
    # x = np.linspace(start = 0.000001, stop = 1.00, num = 1000)
    # R = curvature_radius(x)
    # y = y_t(x)
    # plt.plot(x, y)
    # plt.plot(x, 1/R)
    # plt.show()
    num_x, num_invR = readnumericalR()
    ex_invR = curvature(num_x)
    total = num_invR.shape[0]
    min_r = int(0.05*total)
    max_r = int(0.95*total)
    mmNX = num_x[min_r:max_r]
    mmNR = num_invR[min_r:max_r]
    mmER = ex_invR[min_r:max_r]
    # error = get_mean_square_error(num_invR, ex_invR)
    error = get_error(num_invR, ex_invR)
    # plt.plot(num_x, num_invR)
    # plt.plot(num_x, ex_invR)
    # plt.plot(ex_invR, error, "x")
    # plt.plot(num_invR[min_r:max_r], error[min_r:max_r], "x")
    
    mm_new_NR, params = fitting(mmNR, mmER)
    # plt.plot(mmNX, mm_new_NR, "x")
    new_error = get_error(mm_new_NR, mmER)
    sum = np.sum(np.sqrt(new_error**2))
    print(sum)
    plt.plot(mmNR, new_error, "x")
    plt.show()
    
def readnumericalR():
    path = "G:\\Toyota\\validQ\\centrifugal\\"
    fname = "curvature.txt"
    with open(path + fname, "r") as f:
        num = int(f.readline())
        x = np.zeros(num)
        y = np.zeros(num)
        invR = np.zeros(num)
        i = 0
        for line in f:
            x[i] = line.split()[0]
            y[i] = float(line.split()[1]) - 0.5
            invR[i] = line.split()[2]
            i += 1
    
        """
        plt.plot(x, y, "x")
        plt.plot(x, invR, "x")
        plt.show()
        exit()
        """
    return x, invR

if __name__ == '__main__':
    main()
    # readnumericalR()
