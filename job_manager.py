# -- coding: utf-8 -
import os
from auto_throwing_job import job_throwing, run_unix

def get_fnamelist(path="", sortby_size=True):
    """
    フォルダ内のファイル名のリストを返す
    :param path:    調べるフォルダ
    :param sortby_size: ファイルサイズ順に並べるときTrue, Falseでタイムスタンプ順
    :return: ファイル名のリスト
    """
    if sortby_size:
        cmd = "ls -lS " + path
    else:
        cmd = "ls -lt " + path
    flist = run_unix(cmd).split("\n")
    fnamelist = []
    for line in flist:
        line = line.split()
        if len(line) > 2:
            fnamelist.append(line[8])

    return fnamelist

def auto_rename(program_name = "ES2"):
    i = 1
    while i > 0:
        if os.path.exists(program_name + str(i).zfill(6)):
            i += 1
        else:
            program_name += str(i).zfill(6)
            i = 0
    return program_name
    
    
def make_case_list(parallel, workingpath="/work/A/FMa/FMa037/Case5/", gridpath="/work/A/FMa/FMa037/mayu_grid/",
                   first=False, case = "Case5",
                   aoa_min=0.0, aoa_max=39.0, aoa_pattern=14, ma_min=0.15, ma_max=1.80, ma_pattern=12,
                   re_min=10000, re_max=1000000,re_pattern=5):
    """フォルダ構成のイメージ
    /workingpath/ES2
    /workingpath/Resume/*
    /workingpath/ResultU/*
    /workingpath/ResultR/*
    /workingpath/ResultC/*
    /workingpath/Running/*
    /workingpath/Complete/*
    """

    def get_case_name(naca, mach, reynolds, angle, case, gridpath, workingpath):
        """
        再計算時に使用する文字列を返す
        :param naca: (int) naca翼の4桁or5桁の整数
        :param mach: (float) マッハ数
        :param reynolds: (float) レイノルズ数
        :param angle: (float) 角度
        :param case: (character) ケース名
        :param gridpath: (character) 格子保存先
        :param savepath: (character) 計算データ保存先
        :return: (list) case_list
        """
        grid = gridpath + "NACA" + naca + ".mayu"
        file = "NACA" + str(naca) + "__AoA" + str(angle) + "__Ma" + str(mach) + "__Re" + str(reynolds) + "__" + case
        aoa = str(angle) + "d0"
        resfile = workingpath + "Running/" + file + "__Resume.vtk"
        return [grid, file, case, aoa, resfile]

    resumepath = workingpath + "Resume/"
    runningpath = workingpath + "Running/"

    if first:   # 初回のみ
        del_aoa = (aoa_max - aoa_min) / aoa_pattern
        del_mach = (ma_max - ma_min) / ma_pattern
        if re_pattern != 1:
            rat_re = (re_max / re_min) ** (1.0 / (re_pattern - 1))
        else:
            rat_re = 1

        length = 0
        case_list = []
        fname_list = get_fnamelist(path=gridpath)
        for fname in fname_list:
            naca = fname.split(".")[0][4:]  # "NACAxxxxx.mayu"拡張子と頭文字外す
            for i in range(re_pattern):
                if type(re_min) == type("infty"):   # 文字の場合は数式使わない
                    reynolds = re_min
                else:
                    reynolds = re_min * (rat_re ** (i - 1))
                for j in range(ma_pattern):
                    mach = ma_min + del_mach * (j - 1)
                    for k in range(aoa_pattern):
                        angle = aoa_min + del_aoa * (k - 1)
                        case_list.append(get_case_name(naca, mach, reynolds, angle, case, gridpath, workingpath))
                        length += 1
                        if length == parallel:
                            # f90ファイル書き換え
                            # makefile書き換え&make
                            # qsub書き換え・ジョブ投入
                            program_name = auto_rename(program_name = "ES2_start")
                            jobname = "sN" + naca
                            job_throwing(parallel, case_list, jobname = jobname, program_name = program_name, first = True)

    else:
        import shutil
        case_list = []
        fname_list = get_fnamelist(path=resumepath)
        
        """ test
        fname_list = ["NACA0012__AoA30.0__Ma0.80__Re50000__case5__128th.vtk",
                     "NACA21092__AoA23.5__Ma0.30__Re4000__case5__232th.vtk"]
        """
        for i in range(parallel):
            fname = (fname_list[i].split("__"))
            naca = fname[0][4:]
            angle = fname[1][3:]
            mach = fname[2][2:]
            reynolds = fname[3][2:]
            case = fname[4]
            case_list.append(get_case_name(naca, mach, reynolds, angle, case, gridpath, workingpath))
            # resumeファイルをrunning中に変更する
            resfile = resumepath + case_list[i][1] + "__Resume.vtk"
            shutil.move(resfile, runningpath)
            
    return case_list


def count_fnumber(path=""):
    """
    フォルダ内のファイル数をカウントする
    :param path: (character) 絶対or相対パス(無ければ現在のディレクトリ)
    :return: (int) ファイル数
    """
    cmd = "ls -U " + path + "| wc -l"
    return run_unix(cmd)



def main(parallel = 280, first=False, debug=False):
    #"""
    if not debug:
        workingpath = "/work/A/FMa/FMa037/Case5/"
        resumepath = "/work/A/FMa/FMa037/Case5/Resume/"
        gridpath = "/work/A/FMa/FMa037/mayu_grid/"
    else:
        parallel = 40
        workingpath = "/work/A/FMa/FMa037/SD_test/"
        resumepath = "/work/A/FMa/FMa037/SD_test/Resume/"
        gridpath = "/work/A/FMa/FMa037/mayu_SD/"
        aoa_min = 0.0
        aoa_max = 50.0 - 1.25
        aoa_pattern = 40
        ma_min = 0.80
        ma_max = ma_min
        ma_pattern = 1
        re_min = "infinity"
        re_max = re_min
        re_pattern = 1
        case_list = make_case_list(parallel, workingpath, gridpath, first=first, case="SD_test", aoa_min=aoa_min,
                                   aoa_max=aoa_max, aoa_pattern=aoa_pattern, ma_min=ma_min, ma_max=ma_max,
                                   ma_pattern=ma_pattern, re_min=re_min, re_max=re_max, re_pattern=re_pattern)
        print(case_list)
        exit()
    # workingpath = "/mnt/d/Toyota/github/2D_comp_viscos/"
    # resumepath = workingpath + "flow_solver/EulerSolver2_2018"
    # gridpath = "/mnt/g/Toyota/Data/grid_vtk/valid/mayu/"
    if first:
        aoa_min = 0.0
        aoa_max = 39.0
        aoa_pattern = 14
        ma_min = 0.75
        ma_max = ma_min
        ma_pattern = 1
        re_min = "infinity"
        re_max = re_min
        re_pattern = 1
        case_list = make_case_list(parallel, workingpath, gridpath, first=first, case="Case5", aoa_min=aoa_min,
                                   aoa_max=aoa_max, aoa_pattern=aoa_pattern, ma_min=ma_min, ma_max=ma_max,
                                   ma_pattern=ma_pattern, re_min=re_min, re_max=re_max, re_pattern=re_pattern)

    else:
        fnumber = int(count_fnumber(path=resumepath))
        if fnumber >= parallel:  # Resumeファイルの数が並列数より多いときジョブ投げる準備
            case_list = make_case_list(parallel, workingpath=workingpath, gridpath=gridpath)
            # f90ファイル書き換え
            # makefile書き換え&make
            # qsub書き換え・ジョブ投入
            jobname = "re_" + str(fnumber)
            program_name = auto_rename()
            job_throwing(parallel, case_list, jobname = jobname, program_name = program_name)

if __name__ == '__main__':
    main(first=True, debug=True)