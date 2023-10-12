! Subroutine:
!    readgrib
!
! Purpose:
!    This function reads ECMWF GRIB1/2 files
!
! Input:
!    none
!
! Output:
!    none
!
! Restrictions:
!    This function is compatible with GRIB2 only if
!   you use at least version 1.9.9 of grib_api library
!
! Comments:
!    The GRIB file should be created with ECMWF script provided
!   with MIMOSA
!
! History
!    201004   M.-A. DROUIN/LPC2E creation
!    201103   M.-A. DROUIN/LPC2E add comments
!    201109   M.-A. DROUIN/LPC2E corrected a memory leak due when using grib_api
!                                all the message loaded in memory were not released
!
!======================================================================
subroutine readgrib(Sisobaric,nomSAap,pl,tl,ul,vl,nx,ny,np)

     use grib_api

    implicit none

    !======================================================================
    ! Declarations & Initialisation
    !======================================================================

    ! Input variables
    !------------------------------------------------------------------
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    integer, intent(in) :: np
    character(len=20), intent(in) :: NomSAap
    character(len=7), intent(in) :: Sisobaric

    ! Output variables
    !------------------------------------------------------------------
    real, dimension(nx,ny,np), intent(out) :: pl
    real, dimension(nx,ny,np), intent(out) :: ul
    real, dimension(nx,ny,np), intent(out) :: vl
    real, dimension(nx,ny,np), intent(out) :: tl

    ! Local variables
    !------------------------------------------------------------------
	integer :: i
	integer :: j
	integer :: k
	integer :: ifile
	integer :: iret
	integer :: igrib
	integer :: is_missing
	integer :: tmp
    integer :: idx
	integer :: PVPresent
	integer :: nb_pv

    real, dimension(nx,ny)          :: lnspl
    real, dimension(np)             :: aaa
    real, dimension(np)             :: bbb
    real, dimension(np)             :: pres
	real, dimension(:), allocatable	:: pv
	real, dimension(:), allocatable	:: temp

	real :: tmp1

    !======================================================================
    ! Program
    !======================================================================

	call grib_open_file(ifile, Sisobaric//nomSAap, 'r', iret)
	if (iret /= GRIB_SUCCESS) then
		print *, "ERROR"
	endif

	call grib_new_from_file(ifile, igrib, iret)
	if (iret == GRIB_END_OF_FILE) then
		write(*,*) 'ERROR: end of file detected'
	endif
	if (iret /= GRIB_SUCCESS) then
		write(*,*) 'ERROR: file cannot be opened'
	endif

    ! Getting constants to convert model levels to pressure
	call grib_get(igrib, 'PVPresent', PVPresent)
	if (PVPresent == 1) then
		call grib_get_size(igrib,'pv', nb_pv)
		allocate(pv(nb_pv))
		call grib_get(igrib, 'pv', pv)
		if (nb_pv /= 2*(np+1)) then
			print *, 'nb_pv = ', nb_pv
			print *, 'np = ', np
			stop
		end if
		do i=1,np
			aaa(i) = (pv(i) + pv(i+1))*0.5
			bbb(i) = (pv(np+1+i) + pv(np+i+2))*0.5
		end do

		deallocate(pv)

	else
		print *, "there is no PV values in your GRIB message"
		stop
	endif

    ! Getting lnsp values
	call grib_get_int(igrib, 'level', tmp, iret)
	call grib_get_int(igrib, 'indicatorOfParameter', tmp, iret)
	call grib_get_size(igrib, 'values', tmp, iret)

	allocate(temp(tmp))
	call grib_get_real4_array(igrib, 'values', temp, iret)
	call grib_release(igrib)

	do i=1,ny
		lnspl(1:nx, i) = temp(((ny-i)*nx)+1:(ny-i)*nx+nx)
	end do

    ! Iteration to get values for each level of T, U and V
	do k =np,1,-1
		call grib_new_from_file(ifile,igrib, iret)
		call grib_get_real4_array(igrib, 'values', temp, iret)
		call grib_release(igrib)
		do i=1,ny
			tl(1:nx,i,k) = temp(((ny-i)*nx)+1:(ny-i)*nx+nx)
		end do
		call grib_new_from_file(ifile,igrib, iret)
		call grib_get_real4_array(igrib, 'values', temp, iret)
		call grib_release(igrib)
		do i=1,ny
			ul(1:nx,i,k) = temp(((ny-i)*nx)+1:(ny-i)*nx+nx)
		end do
		call grib_new_from_file(ifile,igrib, iret)
		call grib_get_real4_array(igrib, 'values', temp, iret)
		call grib_release(igrib)
		do i=1,ny
			vl(1:nx,i,k) = temp(((ny-i)*nx)+1:(ny-i)*nx+nx)
		end do
	end do

	deallocate(temp)
	!call grib_release(igrib)
	call grib_close_file(ifile)

    ! Creating the pressure array
	do k=1,np
		do j=1,ny
			do i=1,nx
				pl(i,j,k) = (aaa(np-k) + bbb(np-k)*exp(lnspl(i,j)))*0.01
			end do
		end do
	end do

end subroutine readgrib
