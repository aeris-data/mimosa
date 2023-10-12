! Subroutine:
!    lisse2d
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
!    none
!
! History
!    MIMOSA was written by A. Hauchecorne
!    201004   M.-A. DROUIN/LPC2E conversion into FORTRAN 95
!    201103   M.-A. DROUIN/LPC2E add comments
!    201111   M.-A. DROUIN/LPC2E correct declaration of npoids and sum
!                                array in case nlis2D is greater than 15
!
!======================================================================
subroutine lisse2D(ng,dgpv,dgpvlis,nlis2D,gind0)

    implicit none
    
    !======================================================================
    ! Declarations & Initialisation
    !======================================================================

    ! Input variables
    !------------------------------------------------------------------
    integer,                                   intent(in) :: ng
    integer,                                   intent(in) :: nlis2D
    integer(kind=1), dimension(-ng:ng,-ng:ng), intent(in) :: gind0
    real, dimension(-ng:ng,-ng:ng),            intent(in) :: dgpv

    ! Output variables
    !------------------------------------------------------------------
    real, dimension(-ng:ng,-ng:ng), intent(out) :: dgpvlis

    ! Local variables
    !------------------------------------------------------------------
    integer                                :: ip
    integer                                :: jl
    integer                                :: il
    integer                                :: jg
    integer                                :: ig
    integer                                :: igl
    integer                                :: jgl
    integer                                :: ig2jg2
    integer                                :: nlis
    integer                                :: ndlis
    integer, dimension(-nlis2D:nlis2D,-nlis2D:nlis2D,0:nlis2D) :: npoids
    integer, dimension(0:nlis2D)           :: sum
    real                                   :: absig
    real                                   :: absjg
    real                                   :: rlismax
    real                                   :: r

    !======================================================================
    ! Program
    !======================================================================

    ! Filling of NPOIDS array
    !------------------------------------------------------------------
    do ip = 0,nlis2D

        sum(ip)=0
        do jl = -ip,ip
            do il = -ip,ip

                npoids(il,jl,ip) = 1+ip*ip-il*il-jl*jl
                if (npoids(il,jl,ip) < 0) then
                    npoids(il,jl,ip) = 0
                endif
                sum(ip) = sum(ip)+npoids(il,jl,ip)
            enddo
        enddo
    enddo

    do jg = -ng,ng
        do ig = -ng,ng

            ig2jg2 = ig*ig+jg*jg
            dgpvlis(ig,jg) = 0.
            if (ig2jg2 <= ng*ng) then

                absig = real(abs(ig))
                absjg = real(abs(jg))
                rlismax = (-absig-absjg+sqrt(2.*real(ng*ng)-(absig-absjg)**2))/2.
                nlis = min(int(rlismax-1.e-3),nlis2D)
                r = sqrt(real(ig2jg2))
                ndlis = int(real(nlis2D)*(sqrt(2.)-1.))+1
                nlis = min(int(real(ng-ndlis))-r*1.000001,real(nlis2D))
                nlis = max(nlis,0)

                do jgl = -nlis,nlis
                    do igl = -nlis,nlis

                        if (abs(ig+igl) > ng .or. abs(jg+jgl) > ng) then
                            print *,'abs', ig, jg, igl, jgl, nlis
                        endif

                        dgpvlis(ig,jg) = dgpvlis(ig,jg)+npoids(igl,jgl,nlis) &
                        * dgpv(ig+igl,jg+jgl)
                    enddo
                enddo

                dgpvlis(ig,jg) = dgpvlis(ig,jg)/sum(nlis)

            endif
        enddo
    enddo

    return

end subroutine lisse2d
