module interfaces_mod

    interface

        ! Functions
        !--------------------------------------------------------------------------------------------------------------------------
        real function ainterp2D(a,glong,glat,nx,ny,ig,jg,ng,paslong,latmin,paslat)

            ! Input variables
            integer,                        intent(in)  :: nx
            integer,                        intent(in)  :: ny
            integer,                        intent(in)  :: ig
            integer,                        intent(in)  :: jg
            integer,                        intent(in)  :: ng
            real, dimension(-ng:ng,-ng:ng), intent(in)  :: glong
            real, dimension(-ng:ng,-ng:ng), intent(in)  :: glat
            real, dimension(nx,ny),         intent(in)  :: a
            real,                           intent(in)  :: paslong
            real,                           intent(in)  :: paslat
            real,                           intent(in)  :: latmin

        end function ainterp2D

        integer(kind=1) function indgrille(i,j,gind0,ng)

            ! Input variables
            integer,                                   intent(in)  :: i
            integer,                                   intent(in)  :: j
            integer,                                   intent(in)  :: ng
            integer(kind=1), dimension(-ng:ng,-ng:ng), intent(in)  :: gind0

        end function indgrille

        integer function njour(j, m, ia)

            ! Input variables
            integer, intent(in) :: ia
            integer, intent(in) :: j
            integer, intent(in) :: m

        end function njour

        real function poids(x,i,ib)

            ! Input variables
            integer, intent(in) :: i
            integer, intent(in) :: ib
            real, intent(in) :: x

        end function poids

        ! Subroutines
        !--------------------------------------------------------------------------------------------------------------------------
        subroutine deplacement(gind0,gind,gu,gv,gx,gy,ng,nh)

            ! Input variables
            integer,                                   intent(in) :: ng
            integer,                                   intent(in) :: nh
            integer(kind=1), dimension(-ng:ng,-ng:ng), intent(in) :: gind0
            real,            dimension(-ng:ng,-ng:ng), intent(in) :: gu
            real,            dimension(-ng:ng,-ng:ng), intent(in) :: gv

            ! Input/Output variables
            real,            dimension(-ng:ng,-ng:ng), intent(inout) :: gx
            real,            dimension(-ng:ng,-ng:ng), intent(inout) :: gy
            integer(kind=1), dimension(-ng:ng,-ng:ng), intent(inout) :: gind

        end subroutine

        subroutine diffusion(gtr,gtrn,gind,gind00,gx0,gy0,ng,alpha)

            ! Input variables
            integer,                                   intent(in) :: ng
            integer(kind=1), dimension(-ng:ng,-ng:ng), intent(in) :: gind
            integer(kind=1), dimension(-ng:ng,-ng:ng), intent(in) :: gind00
            real, dimension(-ng:ng,-ng:ng),            intent(in) :: gx0
            real, dimension(-ng:ng,-ng:ng),            intent(in) :: gy0
            real,                                      intent(in) :: alpha

            !Input/Output variables
            real, dimension(-ng:ng,-ng:ng), intent(inout) :: gtr
            real, dimension(-ng:ng,-ng:ng), intent(inout) :: gtrn

        end subroutine diffusion

        subroutine initgrid(gind,gind0,gx,gy,glat,glong,ng,hem,ndeg)

            ! Input variables
            integer, intent(in) :: ng
            integer, intent(in) :: ndeg
            real,    intent(in) :: hem

            ! Output variables
            integer(kind=1), dimension(-ng:ng,-ng:ng), intent(out) :: gind
            integer(kind=1), dimension(-ng:ng,-ng:ng), intent(out) :: gind0
            real,            dimension(-ng:ng,-ng:ng), intent(out) :: gx
            real,            dimension(-ng:ng,-ng:ng), intent(out) :: gy
            real,            dimension(-ng:ng,-ng:ng), intent(out) :: glat
            real,            dimension(-ng:ng,-ng:ng), intent(out) :: glong


        end subroutine initgrid

        subroutine lisse2D(ng,dgpv,dgpvlis,nlis2D,gind0)

            ! Input variables
            integer,                                   intent(in) :: ng
            integer,                                   intent(in) :: nlis2D
            integer(kind=1), dimension(-ng:ng,-ng:ng), intent(in) :: gind0
            real, dimension(-ng:ng,-ng:ng),            intent(in) :: dgpv

            ! Output variables
            real, dimension(-ng:ng,-ng:ng), intent(out) :: dgpvlis

        end subroutine lisse2D

        subroutine njouri(j,m,ia,nj)

            ! Input variables
            integer, intent(in) :: nj

            ! Output variables
            integer, intent(out) :: ia
            integer, intent(out) :: j
            integer, intent(out) :: m

        end subroutine njouri

        subroutine pvp(theta0, nx, ny, np, np0, pl, ul, vl, tl, thetal, u, v, t, p, dtdp, pv, alatmin, alatmax)

            ! Input variables
            integer, intent(in) :: nx
            integer, intent(in) :: ny
            integer, intent(in) :: np
            integer, intent(in) :: np0
            real, intent(in) :: theta0
            real, intent(in) :: alatmin
            real, intent(in) :: alatmax
            real, dimension(nx,ny,np), intent(in) :: ul
            real, dimension(nx,ny,np), intent(in) :: vl
            real, dimension(nx,ny,np), intent(in) :: tl
            real, dimension(nx,ny,np), intent(in) :: pl

            ! Output variables
            real, dimension(nx,ny), intent(out) :: u
            real, dimension(nx,ny), intent(out) :: v
            real, dimension(nx,ny), intent(out) :: t
            real, dimension(nx,ny), intent(out) :: p
            real, dimension(nx,ny), intent(out) :: pv
            real, dimension(nx,ny), intent(out) :: dtdp
            real, dimension(nx,ny,np), intent(out) :: thetal

        end subroutine pvp

        subroutine raccord(gpv,gind,gind0,gxhem,gyhem,ng,lathem,ndeg,nlis2D)

            ! Input variables
            integer, intent(in) :: lathem
            integer, intent(in) :: ndeg
            integer, intent(in) :: ng
            integer, intent(in) :: nlis2D
            integer(kind=1), dimension(-ng:ng,-ng:ng), intent(in) :: gind0
            real, dimension(-ng:ng,-ng:ng), intent(in) :: gxhem
            real, dimension(-ng:ng,-ng:ng), intent(in) :: gyhem

            ! Input/Output variables
            real, dimension(-ng:ng,-ng:ng,2), intent(inout) :: gpv
            integer(kind=1), dimension(-ng:ng,-ng:ng,2), intent(inout) :: gind

        end subroutine raccord

        subroutine readecmr(ian,mois,jour,iech,ul,vl,tl,nx,ny,np)

            ! Input variables
            integer, intent(in) :: ian
            integer, intent(in) :: iech
            integer, intent(in) :: jour
            integer, intent(in) :: mois
            integer, intent(in) :: nx
            integer, intent(in) :: ny
            integer, intent(in) :: np

            ! Output variables
            real, dimension(nx,ny,np), intent(out) :: ul
            real, dimension(nx,ny,np), intent(out) :: vl
            real, dimension(nx,ny,np), intent(out) :: tl

        end subroutine readecmr

        subroutine readgrib(Sisobaric,nomSAap,pl,tl,ul,vl,nx,ny,np)

            ! Input variables
            integer, intent(in) :: nx
            integer, intent(in) :: ny
            integer, intent(in) :: np
            character(len=20), intent(in) :: NomSAap
            character(len=7), intent(in) :: Sisobaric

            ! Output variables
            real, dimension(nx,ny,np), intent(out) :: pl
            real, dimension(nx,ny,np), intent(out) :: ul
            real, dimension(nx,ny,np), intent(out) :: vl
            real, dimension(nx,ny,np), intent(out) :: tl

        end subroutine readgrib

        subroutine regrid2( gpv, gpv0, gpvn, gpoids, gx, gy, gind, gind0, ng, nhrelax, nhmod, indrelax, dgpv, dgpvlis, nlis2D)

            ! Input variables
            integer, intent(in) :: indrelax
            integer, intent(in) :: ng
            integer, intent(in) :: nhrelax
            integer, intent(in) :: nhmod
            integer, intent(in) :: nlis2D
            real, dimension(-ng:ng,-ng:ng), intent(in) :: gpv0
            integer(kind=1), dimension(-ng:ng,-ng:ng), intent(in) :: gind0

            ! Output variables
            integer(kind=2), dimension(-ng:ng,-ng:ng), intent(out) :: gpoids
            real, dimension(-ng:ng,-ng:ng), intent(out) :: dgpv
            real, dimension(-ng:ng,-ng:ng), intent(out) :: dgpvlis

            ! Input/Output variables
            real, dimension(-ng:ng,-ng:ng), intent(inout) :: gpv
            real, dimension(-ng:ng,-ng:ng), intent(inout) :: gpvn
            real, dimension(-ng:ng,-ng:ng), intent(inout) :: gx
            real, dimension(-ng:ng,-ng:ng), intent(inout) :: gy
            integer(kind=1), dimension(-ng:ng,-ng:ng) :: gind

        end subroutine regrid2

        subroutine xyinterp(inter,wxp,wyp,xxp,yyp,xxinterp,yyinterp)

            ! Input variables
            real, intent(in) :: xxp
            real, intent(in) :: yyp
            real, intent(in) :: wxp
            real, intent(in) :: wyp

            ! Output variables
            integer, intent(out) :: inter
            real, intent(out) :: xxinterp
            real, intent(out) :: yyinterp

        end subroutine xyinterp

    end interface

end module interfaces_mod
