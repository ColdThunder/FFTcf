########################################################################
#                                Makefile
#
EXEC = 3dcf.out

SRCS = 3dcf.f90

OBJS = $(SRCS:.f90=.o)

MODSRCS = ./head.f90 ./myomp.f90 ./p2grid.f90

MODS = $(MODSRCS:.f90=.mod)

INCL = -I./

LIBS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

OPTS = -fPIC -shared-intel -mcmodel=large -qopenmp

FC   = ifort

$(EXEC):$(SRCS) $(MODS)
	$(FC) $(OPTS) $(SRCS) $(MODSRCS) $(INCL) $(LIBS) -o $(EXEC) 

$(MODS):$(MODSRCS)
	$(FC) -c $(OPTS) $(MODSRCS) $(INCL) $(LIBS)

.PHONY:clean
clean:
	-rm -rf $(EXEC) $(MODS)
########################################################################

