program FFTcf
use omp_lib
use myomp
use p2grid
use head
implicit none

integer(4)::LL(3)
real(4),allocatable::pos(:,:)
real(4),allocatable::ran(:,:)
real(8),allocatable::gden(:,:,:)
real(8),allocatable::rden(:,:,:)
complex(4),allocatable::vc1(:,:,:)
complex(4),allocatable::vc2(:,:,:)

real(4)::boxinfo(2,3)

integer(4)::i,j,k,dir

integer(kind(ngal))::pid

real(8),allocatable::cf8(:,:)
real(4)::fs
real(4)::npeff
logical(4)::debug=.false.

call omp_set_num_threads(NUM_THREADS)
LL=[L,L,L]
boxinfo(1:2,1)=[0.,box]
boxinfo(1:2,2)=[0.,box]
boxinfo(1:2,3)=[0.,box]

call memo(0)
call readdata

fs=L/box
call ompshift(pos,3,ngal,[0.,0.,0.],boxinfo)
if (flag_r) call ompshift(ran,3,nran,[0.,0.,0.],boxinfo)

call ompscale(pos,3,ngal,[fs,fs,fs])
if (flag_r) call ompscale(ran,3,nran,[fs,fs,fs])

call massassign(ngal,pos,L,gden,painter)  ! get delta_g
if (flag_r) call massassign(nran,ran,L,rden,painter)  ! get delta_r

if (flag_r) gden=gden-rden

vc1=cmplx(gden,0.0)
call fftc2c_inplace(vc1,LL)  ! get delta(k)
 
vc1=vc1*conjg(vc1)/float(L)**3  ! get P(\vec{k})

call ifftc2c_inplace(vc1,LL)  ! get \xi(\vec{r})

if (debug) then
write(*,*) vc1(1,1,1)
do i=2,L/2
  write(*,*) vc1(i,1,1),vc1(L+2-i,1,1)
enddo
write(*,*) vc1(L/2+1,1,1)
endif

call calcxiell(vc1,LL,cf8,rbin,outname)

call memo(1)

contains

subroutine memo(command)
implicit none
integer(4)::command
if (command.eq.0) then
  write(*,*) 'memo for pos:',(float(ngal)*3)*4/1024.**3,'G'
  allocate(pos(3,ngal))
  write(*,*) 'memo for gden:',(float(L)**3)*8/1024.**3,'G'
  allocate(gden(L,L,L))
  if (flag_r) then
    write(*,*) 'memo for ran:',(float(nran)*3)*4/1024.**3,'G'
    allocate(ran(3,nran))
    write(*,*) 'memo for rden:',(float(L)**3)*8/1024.**3,'G'
    allocate(rden(L,L,L))
  endif
  write(*,*) 'memo for vc1:',(float(L)**3)*8/1024.**3,'G'
  allocate(vc1(L,L,L))
  write(*,*) 'memo for vc2:',(float(L)**3)*8/1024.**3,'G'
  allocate(vc2(L,L,L))

  allocate(cf8(9,rbin))
elseif(command.eq.1) then
  deallocate(pos)
  deallocate(gden)
  if (flag_r) then
    deallocate(ran)
    deallocate(rden)
  endif
  deallocate(vc1)
  deallocate(vc2)

  deallocate(cf8)
endif
endsubroutine memo

subroutine readdata
implicit none
integer(4)::colmax
real(4),allocatable::line(:)
! read in galaxy
colmax=maxval(gxyz)
allocate(line(colmax))
write(*,*) 'reading: ',trim(infg)
open(31,file=infg,status='old',form='formatted')
do pid=1,nghead
  read(31,*)
enddo
do pid=1,ngal
  read(31,*) line
  pos(1:3,pid)=line(gxyz(1:3))
enddo
deallocate(line)

! read in random
if (flag_r) then
colmax=maxval(rxyz)
allocate(line(colmax))
write(*,*) 'reading: ',trim(infr)
open(31,file=infr,status='old',form='formatted')
do pid=1,nrhead
  read(31,*)
enddo
do pid=1,nran
  read(31,*) line
  ran(1:3,pid)=line(rxyz(1:3))
enddo
deallocate(line)
endif
endsubroutine readdata

subroutine calcxiell(vc,LL,cf8,rbin,filename)
implicit none
integer(4)::LL(3)
complex(4)::vc(LL(1),LL(2),LL(3))
integer(4)::rbin
character(*)::filename
real(8)::cf8(9,rbin)
real(4)::ri,rj,rk,rmode,mu
real(4)::pl0,pl2,pl4
integer(4)::bin
real(4)::dr
real(4)::rbasic
real(4)::factor

dr=(rmax-rmin)/rbin
rbasic=box/L
cf8=0.
write(*,*) 'calculating correlation function'
write(*,*) 'rmin,rmax=',rmin,rmax
write(*,*) 'rbin=',rbin
write(*,*) 'dr=',dr
!$omp parallel do default(private) shared(vc,LL,rbasic,rmin,dr,rbin) &
!$omp reduction(+:cf8) schedule(guided)
do k=1,LL(3)
  if (k.le.LL(3)/2+1) then
    rk=k-1
  else
    rk=k-1-LL(3)
  endif
  rk=rbasic*rk
  do j=1,LL(2)
    if (j.le.LL(2)/2+1) then
      rj=j-1
    else
      rj=j-1-LL(2)
    endif
    rj=rbasic*rj
    do i=1,LL(1)
      if (i.le.LL(1)/2+1) then
        ri=i-1
      else
        ri=i-1-LL(1)
      endif
      ri=rbasic*ri
      if (i.eq.1.and.j.eq.1.and.k.eq.1) cycle
      if (i.eq.LL(1)/2+1.or.j.eq.LL(2)/2+1.or.k.eq.LL(3)/2+1) then
        factor=2
      else
        factor=1
      endif
      rmode=sqrt(ri*ri+rj*rj+rk*rk)
      mu=rk/rmode
      bin=ceiling((rmode-rmin)/dr)
      if (bin.lt.1.or.bin.gt.rbin) cycle
      pl0=1.*vc(i,j,k)
      pl2=0.5*(3*mu**2-1)*vc(i,j,k)
      pl4=0.125*(35*mu**4-30*mu**2+3)*vc(i,j,k)
      cf8(1,bin)=cf8(1,bin)+rmode*factor
      cf8(2,bin)=cf8(2,bin)+rmode**2*factor
      cf8(3,bin)=cf8(3,bin)+pl0*factor
      cf8(4,bin)=cf8(4,bin)+pl0**2*factor
      cf8(5,bin)=cf8(5,bin)+pl2*factor
      cf8(6,bin)=cf8(6,bin)+pl2**2*factor
      cf8(7,bin)=cf8(7,bin)+pl4*factor
      cf8(8,bin)=cf8(8,bin)+pl4**2*factor
      cf8(9,bin)=cf8(9,bin)+1.*factor
    enddo
  enddo
enddo
!$omp end parallel do

cf8(1,:)=cf8(1,:)/cf8(9,:)
cf8(2,:)=cf8(2,:)/cf8(9,:)-cf8(1,:)**2
cf8(3,:)=cf8(3,:)/cf8(9,:)
cf8(4,:)=cf8(4,:)/cf8(9,:)-cf8(3,:)**2
cf8(5,:)=cf8(5,:)/cf8(9,:)
cf8(6,:)=cf8(6,:)/cf8(9,:)-cf8(5,:)**2
cf8(7,:)=cf8(7,:)/cf8(9,:)
cf8(8,:)=cf8(8,:)/cf8(9,:)-cf8(7,:)**2

if (usecenter) then
  do i=1,rbin
    cf8(1,i)=rmin+(i-0.5)*dr
  enddo
  cf8(2,:)=0.
endif

write(*,*) 'writing:',trim(filename)
open(32,file=filename,form='formatted',status='replace')
do i=1,rbin
  write(32,'(9e14.6)') cf8(1:9,i)
  if (debug) write(*,*) real(cf8(1:9:2,i))
enddo
close(32)
endsubroutine calcxiell


endprogram FFTcf

include 'mkl_dfti.f90'
include 'intel_3d_fft.f90'
