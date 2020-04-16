# coding: utf-8
import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.animation as animation
import os
import time, datetime
# import pandas as pd
import csv

def sort_mtime(rootdir, return_latest=False):
    xs = []
    ts = []
    for root, dir, files in os.walk(rootdir):
        for f in files:
            path = os.path.join(root, f)
            xs.append(path)
            ts.append(int(path.split("/")[-2]))

    if not return_latest:
        for mtime, path in sorted(xs):
            name = os.path.basename(path)
            t = datetime.datetime.fromtimestamp(mtime)
            print(t, name)
    else:
        # if "NACA2122_aoa17.5_ma0.4_re100" in rootdir:
            # print(xs)
            # print(xs[0][-1])
            # exit()
        # return xs[0][-1]
        # return xs[ts.index(max(ts))]
        return xs

def read_fCoeffs(path, developed_time = 1.0, converge=True):
    def append_coeff(coeff, log):
        time_step = float(log[0])
        coeff[0] += float(log[1])
        coeff[1] += float(log[2])
        coeff[2] += float(log[3])
        coeff[3] += float(log[4])
        coeff[4] += float(log[5])
        return time_step, coeff
    
    internal_dir = os.path.join("postProcessing", "forces")
    baseDir = os.path.join(path, internal_dir)
    # fname = 'forceCoeffs_0.0123741.dat'
    datas = sort_mtime(baseDir, return_latest = True)
    label = ['time', 'Cm', 'Cd', 'Cl', 'Cl(f)', 'Cl(r)']
    
    coeffs = []
    times = []
    for data in datas:
        coeff = np.zeros(5)
        if os.path.exists(data):
            with open(data, "r") as f:
                if not converge:
                    text = f.readlines()
                    total_data = 0
                    time_step = 0
                    for line in text:
                        log = line.split()
                        if not (log[0] == '#'):
                            if float(log[0]) > developed_time:
                                total_data += 1
                                time_step, coeff = append_coeff(coeff, log)
                    coeff /= total_data
                else:
                    log = f.readlines()[-1].split()
                    time_step, coeff = append_coeff(coeff, log)
                times.append(time_step)
                coeffs.append(coeff)
            # print(coeff)
    return coeffs[times.index(max(times))]

def main(type=0):
    if type == 0:
        converge = True
    else:
        converge = False
    run_case = "CV"
    fname = str(type) + "_" + run_case + ".txt"
    root = "/" + os.path.join("work", "A", "FMa", "FMa037", "openFoam_" + run_case)
    print(root)
    with open(fname, "r") as f:
        lines = f.readlines()
    
    data = []
    for line in lines[1:]:
        case_name = line.split()[-1]
        if not os.path.exists(os.path.join(root, case_name, "100000")): # except no converge case
            coeff = read_fCoeffs(path = os.path.join(root, case_name), converge = converge)
            case_sp = case_name.split("_")
            dlist = [case_sp[0].replace("NACA", ""), case_sp[1].replace("aoa", ""), case_sp[2].replace("ma", ""),
                     case_sp[3].replace("re", ""), coeff[0], coeff[1], coeff[2], coeff[3], coeff[4]]
            data.append(dlist)
            """
            naca4.append(case_sp[0].replace("NACA", ""))
            aoa.append(case_sp[1].replace("aoa", ""))
            ma.append(case_sp[2].replace("ma", ""))
            re.append(case_sp[3].replace("re", ""))
            cm.append(coeff[0])
            cd.append(coeff[1])
            cl.append(coeff[2])
            clf.append(coeff[3])
            clr.append(coeff[4])
            """
    """
    df = pd.DataFrame(
        data={"NACA":naca4, "AoA":aoa, "Ma":ma, "Re":re, "Cm":cm, "Cd":cd, "Cl":cl, "Cl(f)":clf, "Cl(r)":clr},
        columns = ["NACA", "AoA", "Ma", "Re", 'Cm', 'Cd', 'Cl', 'Cl(f)', 'Cl(r)']
    )
    
    df.to_csv("aerodynamicForce_CV" + str(type) + ".csv", sep = ",")
    """
    with open("aerodynamicForce_mk2_" + run_case + str(type) + ".csv", "w") as f:
        writer = csv.writer(f)
        writer.writerows(data)
    
    exit()

if __name__ == '__main__':
    main(type=0)
    path = "D:\\Toyota\\work2\\mod_21v\\"#postProcessing\\forces"
    # print(sort_mtime(path, return_latest = True))
    # exit()
    developed_time = 90.0
    read_fCoeffs(path, developed_time)


