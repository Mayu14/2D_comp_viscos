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
    return shape + "_aoa" + str(angle) + "_ma" + str(ma) + "_re" + str(reynolds)
    

def gen_if_block(casename, i):
    block = "if "
    if not (i == 0):
        block = 'else ' + block
    block += '(my_rank == ' + str(i) + ') then\n'
    block += "\tdirectory = '" + casename + "'\n"
    return block

def gen_if_block_all(offset, split, old_style=False, incomp="", from_csv=False, csvName="no_finish.csv"):
    shapeHeader = 'NACA'
    ifblock = ""
    i = 0
    if not from_csv:
        if not incomp == "incomp":
            reynoldsList = [1000, 10000]
            machList = [0.5, 0.8, 1.1, 1.4]
        else:
            reynoldsList = [100, 1000, 10000, 100000]
            machList = ["incomp"]
            
        angleList = [-5.0, 0.0, 5.0, 10.0, 15.0, 20.0]

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
                                        ifblock += gen_if_block(case, i - offset)
                                    i += 1
    
        else:
            i3_list = [12, 22, 32, 42, 52, 62, 72, 82]
            naca45 = [True, False]
    
            for naca4 in naca45:
                if naca4:
                    i2_len = 1
                    i1_list = [0, 2, 5, 7, 9]
                    i2_list = [0, 1, 3, 4, 6, 7]
                else:
                    i2_len = 2
                    i1_list = [1, 2, 3, 4, 5]
                    i2_list = [10, 21, 30, 41, 50]
    
                for i1 in i1_list:
                    for i2 in i2_list:
                        if ((i1 == 0) and (i2 == 0)) or (i1 != 0):
                            for i34 in i3_list:
                                int_4 = str(i1) + str(i2).zfill(i2_len) + str(i34)
                                shape = shapeHeader + int_4
                                for angle in angleList:
                                    for reynolds in reynoldsList:
                                        for mach in machList:
                                            if ((i >= offset) and (i < split + offset)):
                                                case = gen_casename(shape, angle, mach, reynolds)
                                                ifblock += gen_if_block(case, i - offset)
                                            i += 1
    else:
        import csv
        with open(csvName) as f:
            reader = csv.reader(f)
            for row in reader:
                shape = shapeHeader + row[0]
                angle = float(row[1])
                mach = row[2]
                reynolds = int(row[3])
                case = gen_casename(shape, angle, mach, reynolds)
                print(case)
                exit()
                ifblock += gen_if_block(case, i)
                i += 1
                
    ifblock += 'end if\n'
    return ifblock

def gen_fsource(offset=960, split=960, use_offset=False, path='', mach=1.4, gen_if_block_all=gen_if_block_all):
    program_head = """
program runOF
    use mpi
    !$ use omp_lib
    implicit none

    integer :: ierr, my_rank, PETOT, iUnit, access
    integer :: Total_threads
    character(len=1000) :: chd, bmd, shm, exm, crp, rcf, rcf_timeout, command, directory
    character(len=1000) :: cmdmsg, cFileName
    """
    if use_offset:
        program_offset = "integer :: cmdstat, exitstat, offset = " + str(offset) + '\n'
    else:
        program_offset = "integer :: cmdstat, exitstat, offset = " + str(0) + '\n'
    program_mid = """
    double precision :: time_start, time_remaining, time_elapsed, maxCo

    time_remaining = 70000.0d0
    maxCo = 0.9d0

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, PETOT, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
    time_start = MPI_WTIME()
    my_rank = my_rank + offset
    directory = set_directory(my_rank)
    iUnit = 2 * (my_rank + 10)
    chd = "cd "//trim(adjustl(directory))//"/; "
    bmd = "blockMesh > log.blockMesh 2>&1; "
    shm = "snappyHexMesh -overwrite > log.snappyHexMesh 2>&1; "
    exm = "extrudeMesh > log.extrudeMesh 2>&1; "
    crp = "createPatch -overwrite > log.createPatch 2>&1; "
    """
    
    program_mid_12 = """
    command = trim(adjustl(chd))//trim(adjustl(bmd))//trim(adjustl(shm))//trim(adjustl(exm))//trim(adjustl(crp))
    if (access(trim(adjustl(directory))//"/log.extrudeMesh", " ") /= 0) then
        call execute_command_line(trim(adjustl(command)), wait=.true., exitstat=exitstat, cmdstat=cmdstat, cmdmsg=cmdmsg)
    end if
    time_elapsed = MPI_WTIME()
    time_remaining = time_remaining - (time_elapsed - time_start)

    rcf_timeout = set_timeout(rcf, time_remaining)
    command = trim(adjustl(chd))//trim(adjustl(rcf_timeout))"""
    rcf = '\trcf = "rhoSimpleFoam > log.rhoSimpleFoam 2>&1; "\n'

    dir = '50000'
    if mach == "incomp":
        rcf = '\trcf = "simpleFoam > log.simpleFoam 2>&1; "\n'

    restartcheck = '\nif (access(trim(adjustl(directory))//"/' + dir + '/U", " ") /= 0) then\n'
    program_mid2 = """
        call execute_command_line(trim(adjustl(command)), wait=.true., exitstat=exitstat, cmdstat=cmdstat, cmdmsg=cmdmsg)
        exitstat = re_run(exitstat, chd, rcf, directory, time_remaining, my_rank, maxCo)
    else
        exitstat = 1
    end if

    write(cFileName, '("", i1)') exitstat
    cFileName = trim(adjustl(cFileName))//"/"//trim(adjustl(directory))

    open(unit = iUnit, file = trim(adjustl(cFileName)), status = 'unknown')
        write(iUnit, *) ""
    close(iUnit)
    call MPI_FINALIZE(ierr)

    stop
contains
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
    rewrite = ""
    if not (mach == "incomp"):
        rewrite = """
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
        
    program_mid23 = """
        else
            exit_status = 0	!正常終了時
        end if

        return
    end function re_run

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

    function set_directory(my_rank) result(directory)
        integer, intent(in) :: my_rank
        character(len=1000) :: directory
    """
    program_tail ="""
        return
    end function set_directory

end program runOF
    """
    ifblock = gen_if_block_all(offset, split, incomp = mach)
    program = program_head + program_offset + program_mid + rcf + program_mid_12 + restartcheck + program_mid2 + rewrite + program_mid23 + ifblock + program_tail
    fname = path + get_program_name(split, offset) + '.f90'

    with open(fname, 'w') as f:
        f.write(program)

def throw_qsub(split, offset):
    qsub_head = """
#!/bin/bash
#$ -P GR06APR18
#$ -cwd
"""
    qsub_mid = """
#$ -l h_rt=20:00:00
#$ -sync n
"""
    if split > 1:
        jobclass = "#$ -jc dma.LS"
        para = '#$ -pe impi_pslots ' + str(split) + "\n"
    else:
        jobclass = "#$ -jc single"
        para = ""
    name = '#$ -N OF' + str(offset).zfill(5) + "\n"
    qsub_tail = """
. /etc/profile.d/modules.sh
module load intel/2018.2.046
module load openfoam/6_impi
"""
    fname = "auto_" + get_program_name(split, offset) + '.sh'
    run = 'mpirun ./' + get_program_name(split, offset)
    qsub = qsub_head + jobclass + qsub_mid + name + para + qsub_tail + run

    with open(fname, 'w') as f:
        f.write(qsub)
    cmd = 'qsub ' + fname
    run_unix(cmd)

def get_program_name(split, offset):
    name = 'runf90_' + str(split) + "_" + str(offset).zfill(5)
    if split > 1:
        return name
    else:
        return "sr_" + 'runf90_' + str(split) + "_" + str(offset).zfill(5)

def compile(split, offset, path='', run=False):
    p_name = get_program_name(split, offset)
    chd = ''
    if not path == '':
        chd = 'cd ' + path + ';'

    cmd = chd + 'mpiifort -xhost -qopenmp ' + p_name + '.f90 -qopenmp -o ' + p_name
    if run:
        cmd += '; ./' + p_name
    run_unix(cmd)

def main(total=19200, incomp=False):
    mach = 0.5
    if incomp:
        mach = "incomp"
    split = 960#1
    split_num = int(total / split)
    for i in range(split_num):
        offset = i * split
        gen_fsource(offset, split, mach = "incomp")
        compile(split, offset)
        throw_qsub(split, offset)

if __name__ == '__main__':
    incomp = True
    if not incomp:
        total = 19200
    else:
        total = 9600
    # total = 178 #
    #qsub = throw_qsub(split, offset)
    #compile(split, offset)
    #exit()
    #gen_fsource(offset=0, split=240)
    #exit()
    main(total, incomp)
