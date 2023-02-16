module triangulation
	use, intrinsic :: iso_fortran_env
	use s2field
	implicit none

!	integer nva
!	integer nsa,nsf
	
	type vertex
		real(REAL64), dimension(3) :: xyz
		complex(REAL64), dimension(s,s) :: eigenvectors
		type(vertex), pointer :: previous		! Only for memory managment
	contains
		procedure :: make=>make_vertex
	end type
	type(vertex), pointer :: vl
	
	type segment
		type(vertex), pointer :: endpoint1, endpoint2
		real(REAL64) :: link(s)
		type(segment), pointer, dimension(:) :: children
		type(segment), pointer :: previous		! Only for memory managment
	contains
		procedure :: make=>make_segment
 		procedure :: bisect
 		procedure :: rlink
 		procedure :: free=>free_segment
	end type
	type(segment), pointer :: sl
	
	type edge
		type(segment), pointer :: sp
		integer(INT32) :: dir ! -1 <=> reversed
	contains
		procedure :: start
		procedure :: finish
	end type
	
	type triangle
		type(edge), dimension(3) :: edges	
		type(triangle), pointer, dimension(:) :: children
		integer(INT32) :: gen
!		type(triangle), pointer :: parent
	contains
		procedure :: make=>make_triangle
 		procedure :: refine
		procedure :: flux
 		procedure :: rflux
 		procedure :: free=>free_triangle
	end type
				
contains
!	vertex methods and memory deallocation
	subroutine make_vertex(v,xyz)
		use, intrinsic :: iso_fortran_env
		implicit none

		class(vertex) :: v
		real(REAL64), dimension(3), intent(in) :: xyz

		complex(REAL64), dimension(s,s) :: h
		integer(INT32) :: info,iwork(10*s+3)
		real(REAL64) :: eval(s),rwork(3*(s+3)*s)
		complex(REAL64) :: work(2*(s+2)*s)

!		nva=nva+1
		v%xyz=xyz
		call sample(s,xyz,h)
!		print*,'h=',h
		call zheevd('V','U',s,h,s,eval,work,2*(s+2)*s,rwork,3*(s+3)*s,iwork,10*s+3,info)
		if (info.ne.0) then 
			print*, 'zheevd error, info=',info
			stop
		endif
!		print*,'eval=',eval
		v%eigenvectors=h
!		print*,v%eigenvectors(:,1)
	end subroutine make_vertex
	
	subroutine free_v
		type (vertex), pointer :: vp
		do
			if (.not. associated(vl)) exit
			vp=>vl%previous
			deallocate(vl)
			vl=>vp
		enddo
	end subroutine free_v

!	segment methods and memory deallocation
	subroutine make_segment(e,ep1,ep2)
		use, intrinsic :: iso_fortran_env
		implicit none
		
		class(segment) :: e
		type(vertex), pointer, intent(in) :: ep1,ep2
		integer(INT32) :: j
		
		e%endpoint1=>ep1
		e%endpoint2=>ep2
		e%children=>null()
!  		print*,'endpoint1',e%endpoint1%xyz
!  		print*,'endpoint2',e%endpoint2%xyz
		do j=1,s
			e%link(j)=aimag(log(dot_product(e%endpoint1%eigenvectors(:,j),e%endpoint2%eigenvectors(:,j))))
!  			print*,'link(',j,') =',e%link(j)
		enddo
	end subroutine make_segment
	
	subroutine bisect(seg)	! do nothing if already bisected
		use, intrinsic :: iso_fortran_env
		implicit none
		class(segment) :: seg
		
		real(REAL64), dimension(3) :: mp_xyz
		type(vertex), pointer :: midpoint
		
		if (.not. associated(seg%children)) then
			allocate(midpoint)
			midpoint%previous=>vl
			vl=>midpoint
			mp_xyz=(seg%endpoint1%xyz+seg%endpoint2%xyz)/2
			mp_xyz=mp_xyz/norm2(mp_xyz)
			call midpoint%make(mp_xyz)
!			nsa=nsa+1
			allocate(seg%children(2))
			call seg%children(1)%make(seg%endpoint1,midpoint)
			call seg%children(2)%make(midpoint,seg%endpoint2)
		endif
	end subroutine bisect
		
	recursive function rlink(seg)
		use, intrinsic :: iso_fortran_env
		implicit none
		class(segment) :: seg
		real(REAL64), dimension(s) :: rlink
		
		if (associated(seg%children)) then
			rlink=seg%children(1)%rlink()+seg%children(2)%rlink()
		else
			rlink=seg%link
		endif
	end function rlink

	recursive subroutine free_segment(s)
		use, intrinsic :: iso_fortran_env
		implicit none
		
		class(segment) :: s
		integer(INT32) :: j
		
!		print*,'associated(s%children)=',associated(s%children)
		if (associated(s%children)) then
			do j=1,2
				call s%children(j)%free()
			enddo
!			nsf=nsf+1
			deallocate(s%children)
		endif
	end subroutine free_segment
	
	subroutine free_s
		type (segment), pointer :: sp
		do
			if (.not. associated(sl)) exit
!			nsf=nsf+1
			sp=>sl%previous
			call sl%free()
			deallocate(sl)
			sl=>sp
		enddo
	end subroutine free_s

!	edge methods
	function start(s)
		implicit none
		class(edge) :: s
		type(vertex) :: start
		
		select case (s%dir)
		case (1)
			start=s%sp%endpoint1
		case (-1)
			start=s%sp%endpoint2
		end select
	end function start			
	function finish(s)
		implicit none
		class(edge) :: s
		type(vertex) :: finish
		
		select case (s%dir)
		case (1)
			finish=s%sp%endpoint2
		case (-1)
			finish=s%sp%endpoint1
		end select
	end function finish

!	triangle methods
	subroutine make_triangle(t,gen,edge1,edge2,edge3,dir1,dir2,dir3)
		use, intrinsic :: iso_fortran_env
		implicit none
		
		class(triangle) :: t
		type(segment), target, intent(in) :: edge1,edge2,edge3
		integer(INT32), intent(in) :: gen,dir1,dir2,dir3
		integer(INT32) :: j
		type(vertex), dimension(3) :: vs,vf
		
		t%edges(1)%sp=>edge1
		t%edges(2)%sp=>edge2
		t%edges(3)%sp=>edge3
		t%edges(1)%dir=dir1
		t%edges(2)%dir=dir2
		t%edges(3)%dir=dir3
		do j=1,3
			vs(j)=t%edges(j)%start()
			vf(j)=t%edges(j)%finish()
! 			print*,'t%edges(',j,'))%start()',vs(j)%xyz
! 			print*,'t%edges(',j,'))%finish()',vf(j)%xyz
		enddo
		if (.not.(all(vf(1)%xyz.eq.vs(2)%xyz) .and. 			&
				  all(vf(2)%xyz.eq.vs(3)%xyz) .and. 			&
				  all(vf(3)%xyz.eq.vs(1)%xyz))) stop 'triangle edges do not match'
		t%children=>null()
!		t%parent=>null()
		t%gen=gen
! 		do j=1,s
!  			print*,'flux(',j,') =',t%flux(j)
! 		enddo
	end subroutine make_triangle

	subroutine refine(tri)
		use, intrinsic :: iso_fortran_env
		implicit none
		class(triangle), target :: tri
		
		integer(INT32) :: j,k,l
		type(segment), pointer :: internal
		
		if (.not. associated(tri%children)) then
			do j=1,3
! 				print*,j
! 				print*,tri%edges(j)%sp%endpoint1%xyz
! 				print*,tri%edges(j)%sp%endpoint2%xyz
				call tri%edges(j)%sp%bisect()
			enddo
			allocate(tri%children(4))
			do j=1,3
				k=1+modulo(j,3)
				l=1+modulo(k,3)
!				print*,j,k,l
				allocate(internal)
				call internal%make(tri%edges(k)%sp%children(2)%endpoint1,		&
							tri%edges(l)%sp%children(1)%endpoint2)
				internal%previous=>sl
				sl=>internal
!			nsa=nsa+1
				call tri%children(j)%make(tri%gen+1,tri%edges(k)%sp%children(2+(tri%edges(k)%dir-1)/2),	&
										tri%edges(l)%sp%children(1-(tri%edges(l)%dir-1)/2),		&
										internal,tri%edges(k)%dir,tri%edges(l)%dir,-1)
!				print*,'after'
!				tri%children(j)%parent=>tri
			enddo
! 			do j=1,3
! 				print*,j
! 				print*,tri%children(j)%edges(3)%sp%endpoint1%xyz
! 				print*,tri%children(j)%edges(3)%sp%endpoint1%xyz
! 				print*,tri%children(j)%edges(3)%dir
! 			enddo
! 			call tri%children(4)%make(tri%children(1)%edges(3)%sp,tri%children(2)%edges(3)%sp,		&
! 						tri%children(3)%edges(3)%sp,tri%edges(1)%dir,tri%edges(2)%dir,tri%edges(3)%dir)
			call tri%children(4)%make(tri%gen+1,tri%children(1)%edges(3)%sp,tri%children(2)%edges(3)%sp,		&
						tri%children(3)%edges(3)%sp,1,1,1)
!			tri%children(4)%parent=>tri
		else
			stop 'triangle already refined'
		endif
	end subroutine refine

	function flux(tri)
		use, intrinsic :: iso_fortran_env
		implicit none
		class(triangle) :: tri
		real(REAL64), dimension(s) :: flux
		integer(INT32) :: j

		flux=0
		do j=1,3
			flux=flux+tri%edges(j)%dir*tri%edges(j)%sp%rlink()
		enddo
		flux=modulo(flux+pi,2*pi)-pi
	end function flux

	recursive function rflux(tri)
		use, intrinsic :: iso_fortran_env
		implicit none
		class(triangle) :: tri
		real(REAL64), dimension(s) :: rflux
		
		integer(INT32) :: j

		if (associated(tri%children)) then
			rflux=0
			do j=1,4
				rflux=rflux+tri%children(j)%rflux()
			enddo
		else
			rflux=tri%flux()
		endif
	end function rflux
		
	recursive subroutine free_triangle(t)
		use, intrinsic :: iso_fortran_env
		implicit none
		
		class(triangle) :: t
		integer(INT32) :: j
		if (associated(t%children)) then
			do j=1,4
				call t%children(j)%free()
			enddo
			deallocate(t%children)
		endif
	end subroutine free_triangle

! Init
	subroutine initocta(faces,edges)
		use, intrinsic :: iso_fortran_env
		implicit none

		type(triangle), dimension(8), intent(out) :: faces
		type(segment), dimension(12), intent(out) :: edges
		type vp
			type(vertex), pointer :: v
		end type
		type(vp), dimension(6) :: vertices
			 
		real(REAL64), dimension(3) :: xyz=0._REAL64
		integer(INT32) :: j,k,l
		
		do j=1,6
			allocate(vertices(j)%v)
			if (j==1) then
				vertices(j)%v%previous=>null()
			else
				vertices(j)%v%previous=>vertices(j-1)%v
			endif
		enddo		
		vl=>vertices(6)%v
		do j=1,3
			xyz(j)=1._REAL64
			call vertices(2*j-1)%v%make(xyz)
			xyz(j)=-1._REAL64
			call vertices(2*j)%v%make(xyz)
			xyz(j)=0._REAL64
		enddo
		
		do j=1,4
			select case (j)
			case (1)
				k=1
				l=3
			case (2)
				k=3
				l=2
			case (3)
				k=2
				l=4
			case (4)
				k=4
				l=1
			end select
!			print*,'j=',j,' k=',k
			call edges(j)%make(vertices(5)%v,vertices(k)%v)
			call edges(4+j)%make(vertices(6)%v,vertices(k)%v)
			call edges(8+j)%make(vertices(k)%v,vertices(l)%v)
		enddo
		sl=>null()
		
		do j=1,4
			call faces(j)%make(1,edges(j),edges(8+j),edges(1+modulo(j,4)),1,1,-1)
			call faces(j+4)%make(1,edges(4+j),edges(5+modulo(j,4)),edges(8+j),-1,1,-1)
		enddo
	end subroutine initocta
end module
