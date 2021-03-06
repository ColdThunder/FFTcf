
module head
implicit none

real(4),parameter::box=1000.  ! box size in Mpc/h
integer(4),parameter::L=512   ! cells per dimension

character(512)::outname='./output/2pcf-rsd.txt'

character(512)::infg='./data/ELG-wpmax-v3-snap97-redshift0.99_dens2.dat'  ! galaxy catalogue name
integer(4)::gxyz(3)=[1,2,4]     ! columns for x y z
integer(4)::nghead=0            ! number of lines for header
integer(4)::ngal=2500000          ! number of lines for galaxy

logical(4)::flag_r=.false.      ! use random or not
character(512)::infr='./data/randomsx10_seed0.txt'    ! random catalgue name
integer(4)::rxyz(3)=[1,2,3]       ! columns for x y z
integer(4)::nrhead=0              ! number of lines for header
integer(4)::nran=32340460        ! number of lines for random

integer(4)::painter=2   ! 2 for cic

real(4)::rmin=0.
real(4)::rmax=200.
integer(4)::rbin=100
logical(4)::usecenter=.true.

integer(4)::NUM_THREADS=32

logical(4)::flag_xi2d=.true.
character(512)::outname2='./output/2pcf2d.txt'
real(4)::umin=0.
real(4)::umax=1.
integer(4)::ubin=12

endmodule head
