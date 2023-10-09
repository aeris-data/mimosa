! Subroutine:
!    diffusion
!
! Purpose:
!    none
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
!    Following the ADI Cranck Nicholson Method
!	 Ref : Numerical Recipes F77, pag. 846
!
! History
!    MIMOSA was written by A. Hauchecorne
!    201004   M.-A. DROUIN/LPC2E conversion into FORTRAN 95
!    201103   M.-A. DROUIN/LPC2E add comments
!
!======================================================================
subroutine diffusion(gtr,gtrn,gind,gind00,gx0,gy0,ng,alpha)

    implicit none

!======================================================================
! Declarations & Initialisation
!======================================================================

    ! Input variables
    !------------------------------------------------------------------
    integer,                                   intent(in) :: ng
    integer(kind=1), dimension(-ng:ng,-ng:ng), intent(in) :: gind
    integer(kind=1), dimension(-ng:ng,-ng:ng), intent(in) :: gind00
    real, dimension(-ng:ng,-ng:ng),            intent(in) :: gx0
    real, dimension(-ng:ng,-ng:ng),            intent(in) :: gy0
    real,                                      intent(in) :: alpha

    !Input/Output variables
    !------------------------------------------------------------------
    real, dimension(-ng:ng,-ng:ng), intent(inout) :: gtr
    real, dimension(-ng:ng,-ng:ng), intent(inout) :: gtrn

    ! Local variables
    !------------------------------------------------------------------
    integer          :: ig
    integer          :: jg
    integer          :: x1
    integer          :: x2
    integer          :: y1
    integer          :: y2
    double precision :: dnx
    double precision :: dny

!======================================================================
! Program
!======================================================================

    do ig = -ng,ng
        do jg = -ng,ng
            gtrn(ig,jg)=gtr(ig,jg)
        enddo
    enddo

    do jg = -ng,ng
        do ig = -ng,ng
            x1 = nint(gx0(ig,jg)-.5)
            y1 = nint(gy0(ig,jg)-.5)

            if (x2 > ng) then
                x2 = ng
            endif
            if (x2 < -ng) then
                x2 = -ng
            endif
            if (y2 > ng) then
                y2 = ng
            endif
            if (y2 < -ng) then
                y2 = -ng
            endif

            if ( (gind(ig,jg)     == 1) .and. &
                 (gind00(x1,y1)   == 1) .and. &
                 (gind00(x1,y1-1) == 1) .and. &
                 (gind00(x1,y1+1) == 1) .and. &
                 (gind(ig+1,jg)   == 1) .and. &
                 (gind(ig-1,jg)   == 1)) then

                ! Definition of delta's
                !------------------------------------------------------
                dnx = gtr(ig+1,jg)-2*gtr(ig,jg)+gtr(ig-1,jg)
                dny = gtr(x1,y1+1)-2*gtr(x1,y1)+gtr(x1,y1-1)

                ! Calcul
                !------------------------------------------------------
                gtrn(ig,jg) = gtrn(ig,jg)+alpha*(dnx+dny)

            endif
        enddo
    enddo

    do ig=-ng,ng
        do jg=-ng,ng
            gtr(ig,jg)=gtrn(ig,jg)
        enddo
    enddo

    return

end subroutine diffusion
