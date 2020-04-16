# coding: utf-8
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

def gen_casename(shape, angle, mach, reynolds):
    if not mach == "incomp":
        ma = mach
    else:
        ma = "Inf"
    if reynolds == "inviscid":
        re = "Inf"
    else:
        re = reynolds
    return shape + "_aoa" + str(angle) + "_ma" + str(ma) + "_re" + str(re)
    

def gen_if_block(casename, i, num_iter, max_iter):
    block = ""
    if num_iter == 1:
        block += "if "
        if not (i == 0):
            block = 'else ' + block
        block += '(my_rank == ' + str(i) + ') then\n'
        block += "\tif "
    else:
        block += "\telse if "
    block += '(iter == ' + str(num_iter) + ') then\n'
    block += "\t\tdirectory = '" + casename + "'\n"
    if num_iter == max_iter:
        block += "\tend if\n"
    return block

def gen_if_block_all(offset, split, max_iter, old_style=False, incomp="", from_csv=False, csvName="no_finish.csv", reverse=False):
    shapeHeader = 'NACA'
    ifblock = ""
    i = 0
    if not from_csv:
        angleList = [-5.0, 0.0, 5.0, 10.0, 15.0, 20.0]
        if incomp == "incomp":
            reynoldsList = [100, 1000, 10000, 100000]
            machList = ["incomp"]
        elif incomp == "CV_append":
            reynoldsList = [100, 1000]
            machList = [0.3, 0.4]
            angleList = [-7.5, -2.5, 2.5, 7.5, 12.5, 17.5, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0]
        elif incomp == "validate":
            reynoldsList = [1000]
            machList = [0.1]
            angleList = [-5.0, -3.0, -2.0, 0.0, 2.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 14.0, 16.0]
        elif incomp == "invMk2":
            reynoldsList = ["inviscid"]
            machList = [0.5, 0.8, 1.1, 1.4]
            angleList = [-10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0]
        else:
            reynoldsList = [100, 1000, 100000, 10000]
            machList = [0.3, 0.4, 0.5]

            

        if old_style:
            for i1 in range(1, 10, 3):
                for i2 in range(1, 10, 3):
                    for i34 in range(12, 100, 20):
                        int_4 = str(i1) + str(i2) + str(i34)
                        shape = shapeHeader + int_4
                        for angle in angleList:
                            for reynolds in reynoldsList:
                                for mach in machList:
                                    if ((i >= offset) and (i < split + offset)):
                                        case = gen_casename(shape, angle, mach, reynolds)
                                        total_num = i - offset
                                        my_rank, num_iter = divmod(total_num, max_iter)
                                        ifblock += gen_if_block(case, my_rank, num_iter + 1, max_iter)
                                    i += 1
    
        else:
            if reverse:
                dir = -1
            else:
                dir = 1
            i3_list = [12, 22, 32, 42, 52, 62, 72, 82]
            naca45 = [True, False]
    
            for naca4 in naca45[::dir]:
                if naca4:
                    i2_len = 1
                    i1_list = [0, 2, 5, 7, 9]
                    i2_list = [0, 1, 3, 4, 6, 7]
                else:
                    i2_len = 2
                    i1_list = [1, 2, 3, 4, 5]
                    i2_list = [10, 21, 30, 41, 50]
                if incomp == "validate":
                    i3_list = [2]
                    i2_len = 1
                    i1_list = [4]
                    i2_list = [4]

                for mach in machList[::dir]:
                    for reynolds in reynoldsList[::dir]:
                        for angle in angleList[::dir]:
                            for i1 in i1_list[::dir]:
                                for i2 in i2_list[::dir]:
                                    if ((i1 == 0) and (i2 == 0)) or (i1 != 0):
                                        for i34 in i3_list[::dir]:
                                            int_4 = str(i1) + str(i2).zfill(i2_len) + str(i34).zfill(2)
                                            shape = shapeHeader + int_4

                                            if ((i >= offset) and (i < split + offset)):
                                                case = gen_casename(shape, angle, mach, reynolds)
                                                total_num = i - offset
                                                my_rank, num_iter = divmod(total_num, max_iter)
                                                ifblock += gen_if_block(case, my_rank, num_iter + 1, max_iter)
                                            i += 1
    else:
        import csv
        with open(csvName) as f:
            reader = csv.reader(f)
            for row in reader:
                for case in row:
                    if ((i >= offset) and (i < split + offset)):
                        total_num = i - offset
                        my_rank, num_iter = divmod(total_num, max_iter)
                        print(total_num, max_iter, my_rank, num_iter)
                        ifblock += gen_if_block(case, my_rank, num_iter + 1, max_iter)
                    i += 1
                
    ifblock += '\tend if\n'
    return ifblock

def gen_fsource(offset=960, split=960, use_offset=False, path='', mach=1.4, gen_if_block_all=gen_if_block_all, iter=1, app="rhoSimpleFoam"):
    runOF_1 = """
    program runOF
    use mpi
    !$ use omp_lib
    implicit none

    integer :: ierr, my_rank, PETOT
    integer :: iLoop, """
    runOF_2 = "iter = " + str(iter) + "\n"
    runOF_3 = """
    double precision :: time_remaining, maxCo
    
    time_remaining = 70000.0d0
    maxCo = 0.9d0

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, PETOT, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

    do iLoop = 1, iter
        call main_proc(my_rank, iLoop, time_remaining)
    end do

    call MPI_FINALIZE(ierr)

    stop
contains
    """
    runOF = runOF_1 + runOF_2 + runOF_3

    sub_main_proc_1 = """
    subroutine main_proc(my_rank, iter, time_remaining)
        implicit none
        integer, intent(in) :: my_rank, iter
        double precision, intent(inout) :: time_remaining
        
        character(len=1000) :: chd, bmd, shm, exm, crp, rcf, rcf_timeout, command, directory
        character(len=1000) :: cmdmsg, cFileName
        integer :: my_rank_2, cmdstat, exitstat, iUnit, access, offset = """
    if use_offset:
        sub_main_proc_2 = str(offset) + '\n'
    else:
        sub_main_proc_2 = str(0) + '\n'

    sub_main_proc_3 = """
        double precision :: time_start, time_elapsed
        
        time_start = MPI_WTIME()
        my_rank_2 = my_rank + offset
        directory = set_directory(my_rank_2, iter)
        iUnit = 2 * (my_rank_2 + 10)
        chd = "cd "//trim(adjustl(directory))//"/; "
        bmd = "blockMesh > log.blockMesh 2>&1; "
        shm = "snappyHexMesh -overwrite > log.snappyHexMesh 2>&1; "
        exm = "extrudeMesh > log.extrudeMesh 2>&1; "
        crp = "createPatch -overwrite > log.createPatch 2>&1; "
        rcf = """

    if mach == "incomp":
        sub_main_proc_4 = '"simpleFoam > log.simpleFoam 2>&1; " \n'
    else:
        sub_main_proc_4 = '"{0} > log.{0} 2>&1; " \n'.format(app)

    sub_main_proc_5 = """
        command = trim(adjustl(chd))//trim(adjustl(bmd))//trim(adjustl(shm))//trim(adjustl(exm))//trim(adjustl(crp))
        if (access(trim(adjustl(directory))//"/log.extrudeMesh", " ") /= 0) then
            call execute_command_line(trim(adjustl(command)), wait=.true., exitstat=exitstat, cmdstat=cmdstat, cmdmsg=cmdmsg)
        end if
        time_elapsed = MPI_WTIME()
        time_remaining = time_remaining - (time_elapsed - time_start)

        rcf_timeout = set_timeout(rcf, time_remaining)
        command = trim(adjustl(chd))//trim(adjustl(rcf_timeout))
        !if (access(trim(adjustl(directory))//"/100000/U", " ") /= 0) then
        if (access(trim(adjustl("0"))//"/"//trim(adjustl(directory)), " ") /= 0) then
            call execute_command_line(trim(adjustl(command)), wait=.true., exitstat=exitstat, cmdstat=cmdstat, cmdmsg=cmdmsg)
            exitstat = re_run(exitstat, chd, rcf, directory, time_remaining, my_rank_2, maxCo)
        else
            exitstat = 1
        end if

        write(cFileName, '("", i1)') exitstat
        cFileName = trim(adjustl(cFileName))//"/"//trim(adjustl(directory))

        open(unit = iUnit, file = trim(adjustl(cFileName)), status = 'unknown')
            write(iUnit, *) ""
        close(iUnit)

        return
    end subroutine main_proc
    """
    sub_main_proc = sub_main_proc_1 + sub_main_proc_2 + sub_main_proc_3 + sub_main_proc_4 + sub_main_proc_5

    sub_re_run_1 = """
    recursive integer function re_run(exitstat, chd, rcf, directory, time_remaining, my_rank, maxCo) result(exit_status)
        implicit none
        integer, intent(in) :: exitstat
        character(len=1000), intent(in) :: chd, rcf, directory
        double precision, intent(in) :: time_remaining
        integer, intent(in) :: my_rank
        double precision, intent(in) :: maxCo
        character(len=1000) :: command, rcf_timeout
        integer :: cmdstat
        character(len=1000) :: cmdmsg
        double precision :: time_start, time_elapsed, time_remaining_copy
        double precision :: dumping_ratio = 0.9d0, epsilon = 0.01d0

        time_start = MPI_WTIME()
        ! exitstat :: rhoCentralFoamのreturn statusを受け取って...
        if (time_remaining < 0.0d0) then
            exit_status  = 1
        else if (dumping_ratio * maxCo < epsilon) then	!異常終了した場合(発散など)
            exit_status  = 2
        else if (exitstat /= 0) then	!異常終了した場合(発散など)
            exit_status  = 2
        """
    sub_re_run_2 = ""
    if not (mach == "incomp"):
        sub_re_run_2 = """
            ! ここにcontrolDictの書き直し部分を追加すべし (directoryはこのための引数)
            call reWriteMaxCo(directory, my_rank, dumping_ratio * maxCo)
            ! rhoCentralFoamのみ再走させて
            rcf_timeout = set_timeout(rcf, time_remaining)
            command = trim(adjustl(chd))//trim(adjustl(rcf_timeout))
            write(6,*) trim(adjustl(command))
            call execute_command_line(trim(adjustl(command)), wait=.true., exitstat=exit_status, cmdstat=cmdstat, cmdmsg=cmdmsg)
            time_elapsed = MPI_WTIME()
            time_remaining_copy = time_remaining - (time_elapsed - time_start)
                ! 異常終了してないか確認
            write(6,*) exitstat, time_remaining_copy
            exit_status = re_run(exit_status, chd, rcf, directory, time_remaining_copy, my_rank, dumping_ratio * maxCo)
            """
        
    sub_re_run_3 = """
        else
            exit_status = 0	!正常終了時
        end if

        return
    end function re_run
    """
    sub_re_run = sub_re_run_1 + sub_re_run_2 + sub_re_run_3

    sub_submodules = """
    function set_timeout(command, remaining_time) result(command_with_timeout)
        character(len=1000), intent(in) :: command
        double precision, intent(in) :: remaining_time
        character(len=1000) :: command_with_timeout, remaining_sec

        write(remaining_sec, '("timeout -k 3s ", i5)') int(remaining_time)
        command_with_timeout = trim(adjustl(remaining_sec))//"s "//trim(adjustl(command))

        return
    end function

    function make_command(my_rank) result (command)
        integer, intent(in) :: my_rank
        character(len=1000) :: command

        write(command, "('python hoge.py ', i7.7)") my_rank ! my_rank is args
        return
    end function make_command

    subroutine reWriteMaxCo(directory, my_rank, maxCo)
        character(len=1000), intent(in) :: directory
        integer, intent(in) :: my_rank
        double precision, intent(in) :: maxCo
        character(len=1000) :: cDictName, cDictName_copy, command, cmdmsg
        character(len=2000) :: cLine, new_txt
        integer :: iLoop, access, iUnit, iUnit1, exitstat, cmdstat
        double precision :: dumping_ratio = 0.9d0
        cDictName = trim(adjustl(directory))//"/system/controlDict"
        cDictName_copy = trim(adjustl(cDictName))//"old"

        if (access(cDictName_copy, " ") == 0) then
            open(unit=iUnit, file=trim(adjustl(cDictName_copy)), status='unknown')
            close(iUnit, status='delete')
        end if
        command = "mv "//trim(adjustl(cDictName))//" "//trim(adjustl(cDictName_copy))
        call execute_command_line(trim(adjustl(command)), wait=.true., exitstat=exitstat, cmdstat=cmdstat, cmdmsg=cmdmsg)

        iUnit = 2 * (my_rank + 10)
        iUnit1 = iUnit + 1
        open(unit = iUnit, file = trim(adjustl(cDictName_copy)), status = 'unknown')
        open(unit = iUnit1, file = trim(adjustl(cDictName)), status = 'unknown')

            do iLoop = 1, 24
                read(iUnit, "(a)") cLine
                write(iUnit1, *) trim(adjustl(cLine))
            end do
            read(iUnit, "(a)") cLine
            write(new_txt, "('maxCo ', f7.5)") maxCo * dumping_ratio
            write(iUnit1,*) trim(adjustl(new_txt))//" ;"
            do iLoop = 26, 52
                read(iUnit, "(a)") cLine
                write(iUnit1, *) trim(adjustl(cLine))
            end do
        close(iUnit)
        close(iUnit1)
        return
    end subroutine reWriteMaxCo
    """

    sub_set_dir_1 = """
    function set_directory(my_rank, iter) result(directory)
        integer, intent(in) :: my_rank, iter
        character(len=1000) :: directory
    """
    sub_set_dir_2 = gen_if_block_all(offset, split, incomp=mach, max_iter=iter)

    sub_set_dir_3 ="""
        return
    end function set_directory
end program runOF
    """

    sub_set_dir = sub_set_dir_1 + sub_set_dir_2 + sub_set_dir_3

    program = runOF + sub_main_proc + sub_re_run + sub_submodules + sub_set_dir
    fname = path + get_program_name(split, offset, iter) + '.f90'

    with open(fname, 'w') as f:
        f.write(program)

def throw_qsub(split, offset, iter):
    qsub_head = """
#!/bin/bash
#$ -P GR06APR18
#$ -cwd
"""
    qsub_mid = """
#$ -l h_rt=20:00:00
#$ -sync n
"""
    thread_num = int(split / iter)
    if  thread_num ==  1:
        jobclass = "#$ -jc single"
        para = ""
    elif thread_num <= 960:
        jobclass = "#$ -jc dma.LS"
        para = '#$ -pe impi_pslots ' + str(thread_num) + "\n"
    elif thread_num < 2880:
        jobclass = "#$ -jc dma.LM"
        para = '#$ -pe impi_pslots ' + str(thread_num) + "\n"
    else:
        raise ValueError
        
    name = '#$ -N OF' + str(offset).zfill(5) + "\n"
    qsub_tail = """
. /etc/profile.d/modules.sh
module load intel/2018.2.046
module load openfoam/6_impi
"""
    fname = "auto_" + get_program_name(split, offset, iter) + '.sh'
    run = 'mpirun ./' + get_program_name(split, offset, iter)
    qsub = qsub_head + jobclass + qsub_mid + name + para + qsub_tail + run

    with open(fname, 'w') as f:
        f.write(qsub)
    cmd = 'qsub ' + fname
    run_unix(cmd)

def get_program_name(split, offset, iter):
    name = 'runf90_' + str(split) + "_" + str(offset).zfill(5) + "_" + str(iter).zfill(2)
    if split > 1:
        return name
    else:
        return "sr_" + 'runf90_' + str(split) + "_" + str(offset).zfill(5)

def compile(split, offset, iter, path='', run=False):
    p_name = get_program_name(split, offset, iter)
    chd = ''
    if not path == '':
        chd = 'cd ' + path + ';'

    cmd = chd + 'mpiifort -xhost -qopenmp ' + p_name + '.f90 -qopenmp -o ' + p_name
    if run:
        cmd += '; ./' + p_name
    run_unix(cmd)

def main(total=19200, parallel = 960, iter = 5, incomp=False, reverse=True, app="rhoSimpleFoam"):
    """
    処理対象はoffset -> split + offsetまでのsplit個
    parallel_num = split / iter
    :param total:
    :param incomp:
    :return:
    """
    mach = "CV_append"#0.5
    dir = 1
    if incomp:
        mach = "incomp"
    if reverse:
        dir = -1
    split = parallel * iter
    split_num = int(total / split)
    for i in range(split_num)[::dir]:
        offset = i * split

        gen_fsource(offset, split, mach = mach, iter=iter, app=app)
        compile(split, offset, iter)
        throw_qsub(split, offset, iter)
        # exit()

def restart(total=1588, parallel=1588, iter = 1, incomp=False, reverse=False, csvName="res.csv"):
    from functools import partial
    func = partial(gen_if_block_all, from_csv=True, csvName=csvName)
    mach = 0.5
    dir = 1
    if incomp:
        mach = "incomp"
    if reverse:
        dir = -1
    split = parallel * iter
    split_num = int(total / split)
    for i in range(split_num)[::dir]:
        if (i == split_num - 1) and (not total == parallel):
            split = total - split * split_num
        else:
            split = parallel * iter
        offset = i * split
        print(offset, split, iter)
        gen_fsource(offset, split, mach=mach, gen_if_block_all=func, iter=iter)
        compile(split, offset, iter)
        throw_qsub(split, offset, iter)

        # exit()

def validate(total=15, parallel=15, iter=1, incomp="False", reverse="False"):
    mach = "validate"#0.1
    dir = 1
    if incomp:
        pass
    if reverse:
        dir = -1
    split = parallel * iter
    split_num = int(total / split)
    for i in range(split_num)[::dir]:
        offset = i * split
        gen_fsource(offset, split, mach = mach, iter=iter)
        compile(split, offset, iter)
        throw_qsub(split, offset, iter)
        # exit()

def invMk2(total=16000, split=3200, debug=False, iter=2):
    if debug:
        total = 5
        split = 1
    split_num = int(total / split)
    for i in range(split_num):
        offset = i * split
        gen_fsource(offset, split, mach = "invMk2", iter=iter, app = "rhoCentralFoam")
        compile(split, offset, iter)
        throw_qsub(split, offset, iter)



if __name__ == '__main__':
    restart()
    # invMk2(debug = False)
    
    # validate()
    exit()
    restart()
    exit()
    incomp = False
    if not incomp:
        total = 28800   # CV
        total = 19200   # CV_append
    else:
        total = 9600

    parallel = 240
    iter = 20
    # total = 178 #
    #qsub = throw_qsub(split, offset)
    #compile(split, offset)
    #exit()
    #gen_fsource(offset=0, split=240)
    #exit()
    main(total=total, parallel=parallel, iter=iter, incomp=incomp, reverse = False)
