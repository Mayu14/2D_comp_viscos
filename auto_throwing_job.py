# coding: utf-8
import os
import subprocess
def run_unix(cmd):
    """
    |で区切られたコマンドに対応可能なように変更した
    :param cmd: (character) linuxのシェルコマンド
    :return: コマンドの画面出力
    """
    cmd_list = cmd.split("|")
    proc_list = []
    i = 0
    for cmd in cmd_list:
        if i == 0:
            stdin = None
        else:
            stdin = proc_list[i - 1].stdout
        proc_list.append(
            subprocess.Popen(cmd.split(), stdin=stdin,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE))
        i += 1

    out, err = proc_list[i - 1].communicate()
    if err:
        print(err.strip().decode("utf8"))
        exit()

    return out.strip().decode("utf8")

def run_unix_one_liner(cmd):
    proc = subprocess.Popen(cmd, shell = True,
                            stdout = subprocess.PIPE,
                            stderr = subprocess.PIPE)
    out, err = proc.communicate()
    return out


def prepare_directory(workingpath):
    dir_list = ["ResultU", "ResultC", "ResultR", "Complete/vtk", "Complete/dat", "Resume", "dummy"]
    for dir in dir_list:
        if os.path.exists(workingpath + os.sep + dir) == False:
            os.makedirs(workingpath + dir)


def update_JPCalcCaseAutoFill_f90(case_list, first=False,
                                  path = "/home/FMa/FMa037/2D_comp_viscos/flow_solver/EulerSolver2_2018/source/Routine4JobParallel/"):
    """
    fortranプログラムを書き換える(JPCalcCaseAutoFill_originというファイルを作っておく)
    :param case_list:   計算条件に関するリスト
    :param path:        プログラムのあるフォルダ
    :return:            None
    """
    # path = "/mnt/d/Toyota/github/2D_comp_viscos/flow_solver/EulerSolver2_2018/source/Routine4JobParallel/"  # debug
    
    fname = path + "JPCalcCaseAutoFill"
    with open(fname + "_origin", "r") as f:
        main_program = f.read()

    s1 = "    "
    s2 = s1 + s1
    s3 = s1 + s2
    s4 = s2 + s2
    
    with open(fname + ".f90", "w") as f:
        f.write(main_program)
        f.write(
            s1 + "subroutine grid_change(UConf)\n" + s2 + "implicit none\n" + s2 + "type(Configulation), intent(inout) :: UConf\n")
        f.write(s3 + 'UConf%UseResume = 1\n' + s3 + 'UConf%UseMUSCL = 0\n')
        f.write(s3 + 'allocate(UConf%ResumeInFlowVars(5), UConf%ResumeOutFlowVars(5))\n')
        f.write(s3 + 'UConf%ResumeInFlowVars = 0.0d0\n' + s3 + 'UConf%ResumeInFlowVars(1) = 1.0d0\n')
        f.write(s3 + 'UConf%ResumeInFlowVars(5) = 1.0d0\n')
        for i in range(len(case_list)):
            mach = case_list[i][1].split("__")[2][2:]
            f.write(s3)
            if i != 0:
                f.write("else ")
            f.write("if(Uconf%my_rank == " + str(i) + ") then\n")
            f.write(s4 + 'Uconf%cGridName = trim(adjustl("' + case_list[i][0] + '"))\n')
            f.write(s4 + 'Uconf%cFileName = trim(adjustl("' + case_list[i][1] + '"))\n')
            f.write(s4 + 'UConf%cCaseName = trim(adjustl("' + case_list[i][2] + '"))\n')
            f.write(s4 + 'UConf%dAttackAngle = dPi / 180.0d0 * ' + case_list[i][3] + '\n')
            if first:
                f.write(s4 + 'UConf%UseResume = 0\n')
            else:
                f.write(s4 + 'UConf%ResumeFileName = trim(adjustl("' + case_list[i][4] + '"))\n')
            f.write(s4 + 'UConf%ResumeInFlowVars(2) = ' + mach + 'd0\n')
        f.write(s3 + "end if\n" + s3 + 'UConf%ResumeOutFlowVars = UConf%ResumeInFlowVars\n')
        f.write(s2 + "return\n" + s1 + "end subroutine grid_change\nend subroutine JPCalcCaseAutoFill\n")

def update_makefile(program_name, makefile_path = "/home/FMa/FMa037/2D_comp_viscos/flow_solver/EulerSolver2_2018/",
                    debug = False, mpi=True, compiler = "mpiifort", program_path="/work/A/FMa/FMa037/Case5/"):
    """
    Makefileを書き換える
    :param program_name:    出力されるプログラムの名称
    :param path:            Makefileのおいてあるパス
    :param debug:           デバッグ用のコンパイルオプションにするかどうか
    :param compiler:        どのコンパイラを使うか(mpiif90とかに書き換えるかも)
    :return: None
    """
    # makefile_path = "/mnt/d/Toyota/github/2D_comp_viscos/flow_solver/EulerSolver2_2018/"  # debug
    body = "EXE_DIR = ../../../../../.." + program_path + "\n"
    body += "EXE = " + program_name + "\n"
    body += "FC = " + compiler + "\nLD =  " + compiler + "\nIDIR =\n"
    if debug:
        body += "CFLAGS = -check uninit -check pointers -check bounds -qopenmp -stand f90 -fp-stack-check -O0 -traceback -warn all -ftrapuv -debug full\n"
    else:
        body += "CFLAGS = -xhost -qopenmp -O3\n"

    
    fname = makefile_path + "Makefile"
    
    with open(fname + "_header", "r") as f:
        header = f.read()
    
    with open(fname + "_footer", "r") as f:
        footer = f.read()
    
    with open(fname, "w") as f:
        f.write(header)
        f.write(body)
        f.write(footer)
    
    # clock skew detected 対策
    run_unix("sleep 5")

def generate_qsub(fname, jobname, parallel, program, mpi = True, comargs = "", sync = "n", program_path="flow_solver/EulerSolver2_2018/bin/Debug/"):
    """
    qsubスクリプトを生成する
    :param fname: qsubスクリプト名
    :param jobname: jobの名前
    :param parallel: 並列数
    :param program: プログラムの名前
    :param mpi: mpi並列を使うかどうか
    :param comargs: コマンドライン引数
    :param sync: ジョブ終了までフォアグラウンドで動かすかどうか
    :return:
    """
    
    def set_jobclass_and_time(parallel, mpi = False):
        if mpi:
            time = "20:00:00"
            if parallel < 961:
                jobclass = "dma.LS"
            elif parallel < 2881:
                jobclass = "dma.LM"
            else:
                print("use-mpi and parallel > 2880")
                print("please check the parallel number and parallel type")
                exit()
        
        else:  # openmp
            time = "99:00:00"
            if parallel < 41:
                jobclass = "dma.LS"
            elif parallel < 81:
                jobclass = "smb.A"
            elif parallel < 161:
                jobclass = "smb.B"
            else:
                print("use-openmp and parallel > 160")
                print("please check the parallel number and parallel type")
                exit()
        return jobclass, time
    
    jobclass, time = set_jobclass_and_time(parallel, mpi)
    header = "#!/bin/bash\n###UGEオプション####\n#$ -P GR06APR18\n#$ -cwd\n"
    if mpi:
        pe = " impi_pslots "
    else:
        pe = " OpenMP "
    body = "#$ -jc " + jobclass + "\n#$ -N " + jobname + "\n#$ -l h_rt=" + time + "\n#$ -pe" + pe + str(
        parallel) + "\n#$ -sync " + sync + "\n"
    footer = ". /etc/profile.d/modules.sh\nmodule load intel/2018.2.046\n"
    run = ""
    if mpi:
        run = "mpirun "
    last = run + "." + program_path + program + " " + comargs
    with open(fname, "w") as f:
        f.write(header)
        f.write(body)
        f.write(footer)
        f.write(last)


def job_throwing(parallel, case_list, qsubname="auto_qsub.sh", jobname="auto_gen_job", program_name = "EulerSolver2", make_script="make.sh", first=False, program_path=""):
    """
    ジョブを自動で投げるためのプログラム
    :param parallel: (int) スパコン上のジョブ並列数
    :param case_list: (list) [grid's full path, filename (without extension), casename, AoA [deg], ResumeFile's full path]
    :param qsubname: job script's name
    :param jobname: job's name
    :param program_name: program's name
    :param make_script: filename that run make
    :return: None
    """
    # f90ファイル書き換え
    update_JPCalcCaseAutoFill_f90(case_list, first)
    
    # makefile書き換え&make
    update_makefile(program_name, debug = False, program_path=program_path)
    
    run_unix("./" + make_script) # cd /home/FMa/FMa037/2D_comp_viscos/flow_solver/EulerSover2_2018/; make
    # qsub書き換え・ジョブ投入
    generate_qsub(qsubname, jobname, parallel, program_name, program_path=program_path)
    run_unix("qsub " + qsubname)


if __name__ == '__main__':
    qsubname = "auto_qsub.sh"
    jobname = "auto_job"
    parallel = 10
    case_list = [
        ['/mnt/g/Toyota/Data/grid_vtk/valid/mayu/NACA0012.mayu', 'NACA0012__AoA30.0__Ma0.80__Re50000__case5', 'case5',
         '30.0d0', '/mnt/d/Toyota/github/2D_comp_viscos/Running/NACA0012__AoA30.0__Ma0.80__Re50000__case5__Resume.vtk'],
        ['/mnt/g/Toyota/Data/grid_vtk/valid/mayu/NACA21092.mayu', 'NACA21092__AoA23.5__Ma0.30__Re4000__case5', 'case5',
         '23.5d0', '/mnt/d/Toyota/github/2D_comp_viscos/Running/NACA21092__AoA23.5__Ma0.30__Re4000__case5__Resume.vtk']]
    job_throwing(parallel, case_list, qsubname, jobname, program_name = "ES2", first = True, program_path = "")