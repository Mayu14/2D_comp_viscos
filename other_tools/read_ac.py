# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os

def common(cp=False):
    #path = "D:\\Toyota\\github\\2D_comp_viscos\\flow_solver\\EulerSolver2_2018\\\ResultC\\"
    fname_head = "NACA0012_01_VALIDATE_M080_A01_25"
    path = "D:\\Toyota\\Downloads\\mk12\\ResultC\\"
    # fname_head = "NACA0012_02_VALIDATE_M050_A02_00"

    fname_tail = "th_AC.dat"
    if cp:
        fname_head = "CP_" + fname_head
        fname_tail = "th.dat"
    return path, fname_head, fname_tail

def dat2csv(end = 101):
    path, fname_head, fname_tail = common()

    start = 1

    ac = np.zeros((end - start, 2))

    for i in range(end - start):
        fname_mid = str(start + i)
        fname = path + fname_head + fname_mid + fname_tail

        with open(fname, "r") as f:
            ac[i, :] = np.array(f.readline().split(), float)

    csvname = path + "new_" + fname_head + str(end - start).zfill(3) + ".csv"
    np.savetxt(csvname, ac, delimiter=",")
    exit()

def dat2csv_cp(end=100, gif=False):
    path, fname_head, fname_tail = common(cp=True)
    step = end
    if gif:
        step = 1
        fig = plt.figure()
        ims = []

    for i in range(step, end + 1):
        cp = []

        fname = path + fname_head + str(i) + fname_tail
        with open(fname, "r") as f:
            for line in f:
                cp.append(line.split())
            cp.pop()
        cp = np.array(cp, dtype=float)
        print(cp)

        if gif:
            plt.grid(which="major", color="black", linestyle="-")
            plt.xlim([0, 1])
            plt.ylim([1.5, -1.0])
            line, = plt.plot(cp[:, 0], cp[:, 1], ".", color="blue")
            ims.append([line])
        else:
            csvname = path + "new_" + fname_head + str(end).zfill(3) + ".csv"
            np.savetxt(csvname, cp, delimiter=",")

    if gif:
        ani = animation.ArtistAnimation(fig, ims)
        ani.save("anim.gif", writer="imagemagick")
        ani.save("anim.mp4", writer="ffmpeg")
        plt.show()


def test():
    import sklearn.mixture
    
    
    
if __name__ == '__main__':
    a = 0.2864
    b = 2.1 * 10e-7
    c = 0.0019
    d = 0.24715
    # test()
    end = 50
    dat2csv_cp(end = end, gif = True)
    dat2csv(end=end)
    
    
