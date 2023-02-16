module s2field
	use, intrinsic :: iso_fortran_env
	implicit none
	real(REAL64), parameter :: pi=4*atan(1._REAL64)
	integer(INT32), parameter :: s=20
!	integer(INT32), parameter :: lmax=6			! Gaussian, scale=1
!	integer(INT32), parameter :: lmax=10			! Gaussian, scale=1/2
!	integer(INT32), parameter :: lmax=7			! Lorentzian, scale=1
!	integer(INT32), parameter :: lmax=1			! linear, scale=0.4
	integer(INT32), parameter :: lmax=3			! Gaussian, scale=3
	integer(INT32) :: par=0
	complex(REAL64), dimension(s,s,0:lmax,-lmax:lmax), target :: a1jklm,a2jklm
	complex(REAL64), dimension(:,:,:,:), pointer :: ajklm

contains
! 	subroutine draw
! 		use, intrinsic :: iso_fortran_env
! 		implicit none
! 		complex(REAL64), parameter :: i=(0._REAL64,1._REAL64)
! 	!	integer(INT32) :: j,k
! 		ajklm=0
! 		ajklm(1,1,0,0)=1._REAL64
! 		ajklm(1,2,0,0)=1._REAL64+1.0001_REAL64*i
! 		ajklm(2,1,0,0)=1._REAL64-1.0001_REAL64*i
! 		ajklm(2,2,0,0)=-1._REAL64
! 		ajklm(1,1,1,0)=1
! 		ajklm(1,2,1,1)=1
! 		ajklm(1,2,1,-1)=-i
! 		ajklm(2,1,1,1)=1
! 		ajklm(2,1,1,-1)=i
! 		ajklm(2,2,1,0)=-1
! 	end subroutine draw	

	subroutine draw
		use, intrinsic :: iso_fortran_env
		implicit none
		complex(REAL64), parameter :: i=(0._REAL64,1._REAL64)
!		real(REAL64), parameter :: cl(0:lmax)=								& 	! Gaussian, scale=1
!			(/ 0.4323323583816937,0.4060058497098381,0.1316325433592779,0.0259191791413436		&
!							,0.003665965774607803,0.0004045151299969687,0.00003658720669512191/)
!		real(REAL64), parameter :: cl(0:lmax)=								& 	! Gaussian, scale=1/2
!			(/ 0.12495806717150061, 0.28140724810680356, 0.27303127572403063, 0.178812179732226, 	&
!				0.0891288919057357, 0.035886115409867574, 0.012111857670659067, 			&
!				0.0035161456575747846, 0.0008949640630890462, 0.00020270519991702138, 0.00004134154306531507/)
!		real(REAL64), parameter :: cl(0:lmax)=								& 	! Lorentzian, scale=1
!		       (/0.128075, 0.0988727, 0.050585, 0.0232456, 0.0101642, 0.00431868, 0.00180107, 0.000741346/)
!	        real(REAL64), parameter :: cl(0:lmax)=								& 	! linear, scale=0.4
!		       (/0.36, 0.16/)
		real(REAL64), parameter :: cl(0:lmax)=								& 	! Gaussian, scale=3
			(/ 0.896682, 0.0995494, 0.00368571, 0.0000818887/)
		real(REAL64), dimension(s,s,0:lmax,-lmax:lmax) :: unif1, unif2
		integer(INT32) :: j,k,l
		
		if (par .eq. 0) then
			call random_number(unif1)
			call random_number(unif2)
			a1jklm=sqrt(-2*log(unif1))*cos(2*pi*unif2)
			a2jklm=sqrt(-2*log(unif1))*sin(2*pi*unif2)
			forall(l=0:lmax)
				a1jklm(:,:,l,:)=a1jklm(:,:,l,:)*sqrt(4*pi*cl(l)/(2*l+1))
				a2jklm(:,:,l,:)=a2jklm(:,:,l,:)*sqrt(4*pi*cl(l)/(2*l+1))
			end forall
			forall (j=1:s,k=1:s)
				a1jklm(j,k,:,:)=((1._REAL64+i)*a1jklm(j,k,:,:)+(1._REAL64-i)*a1jklm(k,j,:,:))/2._REAL64
				a2jklm(j,k,:,:)=((1._REAL64+i)*a1jklm(j,k,:,:)+(1._REAL64-i)*a2jklm(k,j,:,:))/2._REAL64
			end forall
			ajklm=>a1jklm
			par=1
		else
			ajklm=>a2jklm
			par=0
		endif
!		print*,'ajklm',ajklm(1,1,0,0),ajklm(1,1,1,-1),ajklm(1,1,1,0),ajklm(1,1,1,1)
	end subroutine draw	

! 		r2h=((1._REAL64+i)*r+(1._REAL64-i)*transpose(r))/2._REAL64
	
	subroutine sample(s,xyz,h)
		use, intrinsic :: iso_fortran_env
		implicit none
		integer(INT32), intent(in) :: s
		real(REAL64), dimension(3), intent(in) :: xyz
		complex(REAL64), dimension(s,s), intent(out) :: h
!		complex(REAL64), parameter :: i=(0._REAL64,1._REAL64)
		integer(INT32) :: l,m
		real(REAL64) :: theta,phi
		
		theta=atan2(sqrt(xyz(1)**2+xyz(2)**2),xyz(3))
		phi=atan2(xyz(2),xyz(1))
		h=0	
		do l=0,lmax
			do m=-l,l
!				print*,'lma',l,m,ajklm(:,:,l,m)
				h=h+ajklm(:,:,l,m)*real_sph_har_y(l,m,theta,phi)
			enddo
		enddo
! 		print*,'xyz',xyz
! 		print*,'h',h!h(1,1),h(1,2),h(2,1),h(2,2)
	end subroutine sample

	function real_sph_har_y(l,m,theta,phi)
		use, intrinsic :: iso_fortran_env
		IMPLICIT NONE
		INTEGER(INT32), INTENT(IN) :: l,m
		REAL(REAL64), INTENT(IN) :: theta,phi
		REAL(REAL64) :: real_sph_har_y
	
		integer(INT32) :: j
		real(REAL64) :: nf
	
	!  	stop 'stop1'
	!  	print*,'here'
	!   	write(11,*) legendre_p(l,m,cos(theta))
	!  	write(11,*) legendre_p(1,0,0.5_REAL64)
	!  	write(11,*) legendre_p(1,1,0.5_REAL64)
	!  	write(11,*) legendre_p(2,1,0.5_REAL64)
		nf=(2*l+1)/(4*pi)
	!  	write(11,*) 'here'
	 ! 	stop 'stop2'
		do j=l-abs(m)+1,l+abs(m)
			nf=nf/j
		enddo
		real_sph_har_y=sqrt(nf)*legendre_p(l,abs(m),cos(theta))
	!	print*,real_sph_har_y
		select case (m)
		case (:-1)
			real_sph_har_y=sqrt(2._REAL64)*real_sph_har_y*sin(abs(m)*phi)
		case(1:)
			real_sph_har_y=sqrt(2._REAL64)*real_sph_har_y*cos(m*phi)
		end select
	end function real_sph_har_y

	FUNCTION legendre_p(l,m,x)
		use, intrinsic :: iso_fortran_env
		IMPLICIT NONE
		INTEGER(INT32), INTENT(IN) :: l,m
		REAL(REAL64), INTENT(IN) :: x
		REAL(REAL64) :: legendre_p

		INTEGER(INT32) :: ll
		REAL(REAL64) :: pll,pmm,pmmp1,somx2

		if (.not.(m >= 0 .and. m <= l .and. abs(x) <= 1.0)) stop 'legendre_p args'
		pmm=1.0
		if (m > 0) then
			somx2=sqrt((1.0_REAL64-x)*(1.0_REAL64+x))
			pmm=product(arth(1.0_REAL64,2.0_REAL64,m))*somx2**m
			if (mod(m,2) == 1) pmm=-pmm
		end if
		if (l == m) then
			legendre_p=pmm
		else
			pmmp1=x*(2*m+1)*pmm
			if (l == m+1) then
				legendre_p=pmmp1
			else
				do ll=m+2,l
					pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
					pmm=pmmp1
					pmmp1=pll
				end do
				legendre_p=pll
			end if
		end if
	contains
	!reference implementation (for the scalar case) is definitional only, and neither parallelized nor optimized for roundoff error
		function arth(first,increment,n)
			use, intrinsic :: iso_fortran_env
			implicit none
			real(REAL64), INTENT(IN) :: first,increment 
			INTEGER(INT32), INTENT(IN) :: n
			real(REAL64), DIMENSION(n) :: arth
			INTEGER(INT32) :: k
			if (n > 0) arth(1)=first 
			do k=2,n
				arth(k)=arth(k-1)+increment 
			end do
		end function arth
	END FUNCTION legendre_p

end module
