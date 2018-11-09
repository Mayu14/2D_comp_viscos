program CompareVTK
    use StructVar_Mod
    use LoopVar_Mod
    implicit none
    character(len=64) :: cFile1,cFile2
    character :: cAnnotate
    integer :: iUnit
    integer, allocatable :: Point(:,:)
    integer :: iMaxPoint
    double precision, allocatable :: Quantity1(:,:), Quantity2(:,:)


    iDim=2
    write(6,*) "Please input 1st file name what you want to compare."
    read(5,*) cFile1
    !cFile1 = "CX52VC23.vtk"

    write(6,*)  "Please input 2nd file name what you want to compare."
    read(5,*) cFile2
    !cFile2 = "CY52VC25.vtk"

    open(unit=1, file=trim(adjustl(cFile1)), status='unknown')
    open(unit=2, file=trim(adjustl(cFile2)), status='unknown')
        do iLoop=1,4
            do iUnit=1,2
                read(iUnit,*) cAnnotate
            end do
        end do

        allocate(Point(2,3)) !unit,XYZ

        do iUnit=1,2
            read(iUnit,*) cAnnotate,(Point(iUnit,iVariable),iVariable=1,3)
        end do

        iMaxPoint = maxval(maxval(Point,2),1)

        do iLoop=1, sum(Point,2)+3
            iUnit=1, 2
                read(iUnit,*) cAnnotate
            end do
        end do

        do

stop
end program

