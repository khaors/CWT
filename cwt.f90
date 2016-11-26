program cwt
implicit none;
!========================================================================================
!                 PROGRAM CWT (Continuous Wavelet Transform)       
!
! This program calculates the Continuous Wavelet Transform of a 1D signal (Time series)
! using the Morlet or Mexican Hat Mother Wavelet. This is a demonstrative program in the  
! sense that the implementation of the Continuous Waveler Transform is straightforward and 
! uses a simple integration routine (trapezoidal method).
!
! To compile using the console and gfortran:
!
! gfortran cwt.f90 -o cwt
!
! (c) 2006. Oscar Garcia-Cabrejo
!
!========================================================================================
character(len=*),parameter :: name='CWT'
real(kind=8),parameter :: version=0.1d0
character(len=30) :: signalfl,outputfl,parfl,arg
integer :: i,j,k,unit1,unit2,ndat,window_size,ierr,wavelet_type
real(kind=8),allocatable:: signal(:),time_(:) 
real(kind=8),dimension(50) :: var
real(kind=8) :: al,bl,ll,ci_pos
real(kind=8),allocatable:: scalogram(:,:),signal3(:),res2(:)
complex,allocatable:: res1(:),signal2(:)
complex :: i1
!
! Say hello to the user and get the name of the parameter file
!
write(*,*)
write(*,4) name,version
write(*,*)
call getarg(1,arg)
if(arg(1:1) .ne. ' ') then
  parfl=trim(adjustl(arg))
else
write(*,*) 'parameter file?'
read(*,'(a30)') parfl
if(parfl(1:1) .eq. ' ') parfl='cwt.par';
endif
write(*,*) 'parameter file= ',trim(parfl);
!
! Read parameter file
!
write(*,*)
write(*,*) 'Reading parameter file...'
write(*,*)
call read_parameters(parfl,signalfl,outputfl,window_size,wavelet_type);
write(*,*)
write(*,*) 'Reading parameter file...finished'
write(*,*)
!
! Define output and signal data files
!
outputfl=trim(adjustl(outputfl))//trim(adjustl(signalfl))//'.out';
signalfl=trim(adjustl(signalfl))//'.dat';
write(*,*) trim(adjustl(signalfl)),trim(adjustl(outputfl)),window_size
!
! Define input and output units
!
unit1=1;
unit2=2;
!
! Read signal
!
write(*,*)
write(*,*) 'Reading input data...'
write(*,*)
open(unit1,file=signalfl,status='unknown',action='read',access='sequential');
read(unit1,*) ndat
!
! Allocate memory
!
allocate(signal(ndat),stat=ierr);
allocate(time_(ndat),stat=ierr);
allocate(signal2(ndat),stat=ierr);
allocate(signal3(ndat),stat=ierr);
allocate(res1(ndat),stat=ierr);
allocate(res2(ndat),stat=ierr);
!
! Read input data and extract signal to be analyzed
!
do i=1,ndat
   read(unit1,*,err=90) (var(j),j=1,2);   
   signal(i)=var(2);
   time_(i)=var(1);
enddo
close(unit1);
write(*,*)
write(*,*) 'Reading input data...finished'
write(*,*)
!
! Define width parameter of the Morlet Wavelet
!
ll=5.0d0;
ll=2.5d0;;
allocate(scalogram(ndat,window_size),stat=ierr)
scalogram=0.0d0;
!
! Calculate tranform
!
do j=1,window_size
   if((int(j/10)*10)==j) write(*,*) 'currently on scale(m)= ',j
   do i=1,ndat
      al=dble(j);bl=dble(i);
      select case(wavelet_type)
      case(0)!Morlet Wavelet
      res1=morletwavelet(time_,al,bl,ll);
      signal2=signal*res1;
      scalogram(i,j)=abs((1.0/sqrt(al))*trapezoidal_complex(time_,signal2));      
      !
      case(1)!Mexican hat
      res2=mexicanhatwavelet(time_,al,bl);
      signal3=signal*res2;
      scalogram(i,j)=abs((1.0/sqrt(al))*trapezoidal_real(time_,signal3));
      end select            
   enddo
enddo
!
! Open output file
!
open(unit2,file=outputfl,status='unknown',action='write');
do i=1,ndat
   !write(unit1,5) signal(i),abs(res1(i)),signal2(i) 
   write(unit2,5) (scalogram(i,j),j=1,window_size);
enddo
close(unit2)
!
! Release memory
!
deallocate(time_,signal,res1,res2,signal2,scalogram,signal3); 
!
! finish program
!
stop
4 format(a,1x,f10.3)
5 format(100f15.6)
90 stop 'Error in data file'
!
 contains
!========================================================================================
  function morletwavelet(t,a,b,l) result(mw)
!========================================================================================
  real(kind=8),parameter :: pi =4.0d0*atan(1.0d0)
  complex,parameter :: uni=(0.0,1.0)
  real(kind=8),dimension(:),intent(in) :: t
  complex,dimension(size(t)) :: mw
  real(kind=8),intent(inout) :: a,b,l
!
  real(kind=8) :: p1,p2
  complex,dimension(size(t)) :: p3,p4 
!
  p1=pi**(-.25);
  p2=(a*l)**(-.5);
  p3=exp(-uni*2*pi*((t-b)/a));
  p4=exp(-0.5*(((t-b)/(a*l)))**2);
  mw=(1.0d0/dsqrt(1.4406d0))*(p1*p2*p3*p4);
!
  end function morletwavelet
!
!========================================================================================
  function mexicanhatwavelet(t,a,b) result(mw)
!========================================================================================
  real(kind=8),dimension(:),intent(in) :: t
  real(kind=8),intent(inout) :: a,b
  real(kind=8),dimension(size(t)) :: mw
!
  real(kind=8),dimension(size(t)) :: p1,p2,u2
!
  u2=((t-b)/a)**2
  p1=1.0-u2;
  p2=exp(-.5*u2);
  mw=p1*p2
!
  end function mexicanhatwavelet
!========================================================================================
  function trapezoidal_complex(t,f) result(i)
!========================================================================================
  complex :: i
  real(kind=8),dimension(:),intent(inout) :: t
  complex,dimension(:),intent(inout) :: f
!
  integer :: ndat,ncol,nrow
!
  ndat=size(t,1);nrow=size(f,1);
  if(ndat .ne. nrow) stop 'error in matrix'
  i=0.0;
  i=0.5*(f(1)+2.0*sum(f(2:ndat-1))+f(ndat));
!
  end function trapezoidal_complex
!========================================================================================
  function trapezoidal_real(t,f) result(i)
!========================================================================================
  complex(kind=8) :: i
  real(kind=8),dimension(:),intent(inout) :: t
  real(kind=8),dimension(:),intent(inout) :: f
!
  integer :: ndat,ncol,nrow
!
  ndat=size(t,1);nrow=size(f,1);
  if(ndat .ne. nrow) stop 'error in matrix'
  i=0.0;
  i=0.5d0*(f(1)+2.0d0*sum(f(2:ndat-1))+f(ndat));
!
  end function trapezoidal_real  
!========================================================================================
  subroutine read_parameters(parfl,inputfl,outputfl,window_size,wavelet_type)
!========================================================================================
  character(len=*),intent(inout) :: parfl
  character(len=*),intent(inout) :: inputfl,outputfl
  integer,intent(inout) :: window_size,wavelet_type
!
  integer :: unit_par
  logical :: check_exists
!
  unit_par=1;
!
  inquire(file=trim(parfl),exist=check_exists);
  if(.not. check_exists) then
    write(*,*) 'error: the parameter file ',trim(parfl),' does not exist';
    stop
  endif
!
  write(*,*) 'reading parameter file...'
  open(unit_par,file=trim(parfl),status='unknown')
  read(unit_par,'(a)') inputfl
  read(unit_par,'(a)') outputfl
  read(unit_par,*) window_size
  read(unit_par,*) wavelet_type !0=morlet,1=mexican_hat
  close(unit_par)
  write(*,*) 'reading parameter file...ok'
!
  end subroutine read_parameters
!  
end program cwt