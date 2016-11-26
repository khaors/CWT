FF = gfortran
#
OFLAGS= -g -pg  -Wall -fall-intrinsics -fbounds-check -o 
#
cwt: cwt.f90
	$(FF) $^ $(OFLAGS) $@
