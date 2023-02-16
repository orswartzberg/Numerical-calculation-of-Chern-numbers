program chernsph
	use, intrinsic :: iso_fortran_env
	use triangulation
	use s2field
	implicit none
! 	integer(INT32), parameter :: nfc=8,nedg=12,maxgen=30,samples=1000
! 	real(REAL64), parameter :: tol=1.e-1_REAL64
! 	integer(INT32), parameter :: ranseed=22112010
	integer(INT32), parameter :: nfc=8,nedg=12
	integer(INT32) :: maxgen,samples,ranseed
	real(REAL64) :: tol
	namelist/parameters/maxgen,samples,ranseed,tol

!	type(vertex), pointer :: vl
	type(segment), dimension(nedg) :: edges
	type(triangle), dimension(nfc), target :: faces
	type(triangle), pointer :: fc
	integer(INT32) :: j,r
	real(REAL64), dimension(s) :: ch,ch1
!	integer(INT32) :: ntr(maxgen)
!	real(REAL64) :: err(maxgen-1)

	read(5,nml=parameters)
	print*,'Chern numbers on S^2, matrix-element cf = ','G',', corr length =',1._REAL64
	print*,'matrix size =',s,', stopping tolerance =',tol,', random seed =',ranseed
	call init_random_seed(ranseed)
	do r=1,samples
		call draw
!		nsa=0; nsf=0;
		call initocta(faces,edges)
! 		ch1=0
! 		do f=1,nfc
! 			ch1=ch1+faces(f)%flux()
! 		enddo
!			print*,ch1/(2*pi)
! 		err=0
! 		ntr=0
		do j=1,nfc
			fc=>faces(j)
			call adapt(fc,tol,maxgen)				
		enddo
		
		ch=0
		do j=1,nfc
			ch=ch+faces(j)%rflux()
		enddo
		ch=ch/(2*pi)
		print*,nint(ch)!,maxval(abs(ch-nint(ch))),sum(nint(ch))
! 		do j=1,maxgen-1
! 			if(ntr(j).eq.0) exit
! 			print*,j,ntr(j),ntr(j)/(4._REAL64)**j,2**j*err(j)
! 		enddo
		call free_v()
		do j=1,nedg
			call edges(j)%free()
		enddo
		call free_s()
		do j=1,nfc
			call faces(j)%free()
		enddo
!		print*,'nsa=',nsa,' nsf=',nsf
	enddo
contains
	recursive subroutine adapt(fc,tol,maxgen)
		use, intrinsic :: iso_fortran_env
		use triangulation
		implicit none
		integer(INT32), intent(in) :: maxgen
		real(REAL64), intent(in) :: tol
		type(triangle), pointer, intent(inout) :: fc
		integer(INT32) :: k
		type(triangle), pointer :: ch
		real(REAL64), dimension(s) :: sf2,sd

!		ntr(fc%gen)=ntr(fc%gen)+1
		if (fc%gen<maxgen) then
			call fc%refine()
			sf2=0
			do k=1,4
				sf2=sf2+fc%children(k)%flux()**2
			enddo
!			print*,'gen',fc%gen,sqrt(sf2/4-(fc%flux()/4)**2)
			sd=sqrt(sf2/4-(fc%flux()/4)**2)
!			if (maxval(sd)>err(fc%gen)) err(fc%gen)=maxval(sd)
			if (maxval(min(2._REAL64**fc%gen,1/abs(fc%flux()))*sd)>tol) then
				do k=1,4
					ch=>fc%children(k)
					call adapt(ch,tol,maxgen)
				enddo
			endif
		endif
	end subroutine adapt
end program
SUBROUTINE init_random_seed(s)
        implicit none
        integer, intent(in) :: s
        INTEGER :: i, n, clock
        INTEGER, DIMENSION(:), ALLOCATABLE :: seed

        CALL RANDOM_SEED(size = n)
        ALLOCATE(seed(n))

        seed=0
        seed(1)=s
        CALL RANDOM_SEED(PUT = seed)

        DEALLOCATE(seed)
END SUBROUTINE




