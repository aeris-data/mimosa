! Subroutine:
!    readecmr
!
! Purpose:
!    This function reads ECMWF isobaric files
!
! Input:
!    none
!
! Output:
!    none
!
! Restrictions:
!    none
!
! Comments:
!    none
!
! History
!    MIMOSA was written by A. Hauchecorne
!    201004   M.-A. DROUIN/LPC2E conversion into FORTRAN 95
!    201103   M.-A. DROUIN/LPC2E add comments
!
!======================================================================
subroutine readecmr(ian,mois,jour,iech,ul,vl,tl,nx,ny,np)

    implicit none

    !======================================================================
    ! Declarations & Initialisation
    !======================================================================

    ! Input variables
    !------------------------------------------------------------------
    integer, intent(in) :: ian
    integer, intent(in) :: iech
    integer, intent(in) :: jour
    integer, intent(in) :: mois
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    integer, intent(in) :: np

    ! Output variables
    !------------------------------------------------------------------
    real, dimension(nx,ny,np), intent(out) :: ul
    real, dimension(nx,ny,np), intent(out) :: vl
    real, dimension(nx,ny,np), intent(out) :: tl

    ! Local variables
    !------------------------------------------------------------------
    integer :: ip
    integer :: ix
    integer :: iy

    !======================================================================
    ! Program
    !======================================================================

    do ip = np,1,-1
        read(UNIT=21, FMT='(/////)')
        do iy = 1,ny
            read(UNIT=21, FMT='(10F7.2)') (tl(ix,iy,ip), ix = 1,nx)
        enddo
    enddo

    do ip = np,1,-1
        read(UNIT=21, FMT='(/////)')
        do iy = 1,ny
            read(UNIT=21, FMT='(10F7.2)') (ul(ix,iy,ip), ix = 1,nx)
        enddo
    enddo

    do ip = np,1,-1
        read(UNIT=21, FMT='(/////)')
        do iy = 1,ny
            read(UNIT=21, FMT='(10F7.2)') (vl(ix,iy,ip), ix = 1,nx)
        enddo
    enddo

    return

end
