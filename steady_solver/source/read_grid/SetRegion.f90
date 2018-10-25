!***********************************/
!	Name:MayuGrid2�`���̊i�q���ǂݍ��݃��\�b�h
!	Alias:GetCellAndEdgeCoordinate
!	Description:�o�[�W����0.9�ɑΉ�
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.10.??
!	Update:2018.10.25
!	Other:
!***********************************/
subroutine Set_Region(Conf, Geom, CC, CE)
    use StructVar_Mod_mk2
    implicit none

    type(Configulation), intent(in) :: Conf
    type(Geometry), intent(inout) :: Geom
    type(CellCenter), intent(inout) :: CC
    type(CellEdge), intent(inout) :: CE

    ! �i�q�_���ǂݍ���
    call read_mayugrid2(Conf%GridFileName, Geom)
    ! ���̒i�K�łǂ̕ϐ����ǂ̑傫���ŕK�v���������邽�߈�Ă�allocate���Ă���
    call MallocAll(Conf, Geom, CC, CE)
    ! �Z���E�E�ʂ̒��S���W���Z�o
    call GetCellAndEdgeCoordinate(Geom)
    ! ���R�r�s�񎮂��v�Z
    call GetGridJacobian(Geom, CC, CE)

    return
end subroutine Set_Region
