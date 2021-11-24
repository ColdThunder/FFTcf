########################################################################
#                                Makefile
#
EXEC = 3dcf.out

SRCS = 3dcf.f90

OBJS = $(SRCS:.f90=.o)

#MKL_PATH = -L/opt/intel/mkl/lib/

#MKL_INCL = -I/opt/intel/mkl/include/

MODSRCS = ./head.f90 ./myomp.f90 ./p2grid.f90

MODS = $(MODSRCS:.f90=.mod)

INCL =

LIBS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

OPTS = -fPIC -shared-intel -mcmodel=large -qopenmp

FC   = ifort

$(EXEC):$(SRCS) $(MODSRCS)
	$(FC) $(SRCS) $(MODSRCS) $(MKL_PATH) $(MKL_INCL) $(LIBS) $(OPTS) -o $(EXEC) 

.PHONY:clean
clean:
	-rm -rf $(EXEC) $(MODS)
########################################################################

