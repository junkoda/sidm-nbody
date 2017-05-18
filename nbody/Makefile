# % make
all: sidm-gadget

EXEC     = sidm-gadget
CC       = mpicc
CFLAGS  := -O2
LIBS    := -lm

# non-standard library path if necessary
FFTW2_DIR ?= #e.g. /Users/jkoda/Research/opt/gcc/fftw2
GSL_DIR ?=  #e.g. /Users/jkoda/Research/opt/gcc/gsl

DIR_PATH   = $(FFTW2_DIR) $(GSL_DIR)


#
# Options
#
#OPT  +=  -DPERIODIC        # enables periodic boundaries conditions
#OPT  +=  -DDIAG            # outputs performance diagnostic for tree walks
#OPT  +=  -DRECTREECONS     # switches on recursive computation of multipole moments
OPT   +=  -DINLINE          # compiler accepts the `inline' statement
OPT   +=  -DBMAX            # enables very conservative node-opening 
OPT   +=  -DSIDM            # enables Self-Interacting Dark Matter
OPT   +=  -DREFLECTIONBOUNDARY
                            # spherical reflection boundary for a isolated halo
#OPT  +=  -DSCATTERLOG      # record scattering event
#OPT  +=  -DFINDNBRLOG     
#OPT  +=  -DNOSCATTER       # no velocity update due to scatter (for code test)
OPT   += -DRANDOM_GSL       # use GSL for random number generator

OPT   += -DCROSS_SECTION_TYPE=0
# CROSS_SECTION_TYPE
# 0: Hard sphere (velocity independent, isotropic direction)
# 1: Maxwellian  (sigma ~ v^-1, isotropic direction)
# 2: Yukawa-like (isotropic direction)
# 3: Power law   (sigma ~ v^-n)
# 4: Yukawa interaction

CFLAGS  += $(OPT)
CFLAGS  += $(foreach dir, $(DIR_PATH), -I$(dir)/include)
LIBS    += $(foreach dir, $(DIR_PATH), -L$(dir)/lib)

OBJS    := main.o run.o predict.o begrun.o endrun.o global.o
OBJS    += timestep.o init.o restart.o  io.o sfr.o
OBJS    += accel.o read_ic.o cooling.o dissolvegas.o
OBJS    += system.o allocate.o density.o
OBJS    += gravtree.o timeline.o hydra.o veldisp.o
OBJS    += domain.o allvars.o potential.o
OBJS    += forcetree.o ewald.o read_ic_cluster.o
OBJS    += sidm.o reflection.o sidm_rand.o

ifneq (,$(findstring DRANDOM_GSL, $(CFLAGS)))
  LIBS +=  -lgsl -lgslcblas
endif


$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): Makefile

# cc -MM *.c
accel.o: accel.c allvars.h tags.h proto.h forcetree.h cooling.h
allocate.o: allocate.c allvars.h tags.h proto.h forcetree.h cooling.h
allvars.o: allvars.c allvars.h tags.h
begrun.o: begrun.c allvars.h tags.h proto.h forcetree.h cooling.h \
  sidm_rand.h
cooling.o: cooling.c allvars.h tags.h proto.h forcetree.h cooling.h
density.o: density.c allvars.h tags.h proto.h forcetree.h cooling.h
dissolvegas.o: dissolvegas.c allvars.h tags.h proto.h forcetree.h \
  cooling.h
domain.o: domain.c allvars.h tags.h proto.h forcetree.h cooling.h \
  domain.h
endrun.o: endrun.c allvars.h tags.h proto.h forcetree.h cooling.h
ewald.o: ewald.c
forcetree.o: forcetree.c allvars.h tags.h proto.h forcetree.h cooling.h
global.o: global.c allvars.h tags.h proto.h forcetree.h cooling.h
gravtree.o: gravtree.c allvars.h tags.h proto.h forcetree.h cooling.h
hydra.o: hydra.c allvars.h tags.h proto.h forcetree.h cooling.h
init.o: init.c allvars.h tags.h proto.h forcetree.h cooling.h
io.o: io.c allvars.h tags.h proto.h forcetree.h cooling.h
main.o: main.c allvars.h tags.h proto.h forcetree.h cooling.h
potential.o: potential.c allvars.h tags.h proto.h forcetree.h cooling.h
predict.o: predict.c allvars.h tags.h proto.h forcetree.h cooling.h
read_ic.o: read_ic.c allvars.h tags.h proto.h forcetree.h cooling.h
read_ic_cluster.o: read_ic_cluster.c allvars.h tags.h proto.h forcetree.h \
  cooling.h
reflection.o: reflection.c allvars.h tags.h
restart.o: restart.c allvars.h tags.h proto.h forcetree.h cooling.h
run.o: run.c allvars.h tags.h proto.h forcetree.h cooling.h
sfr.o: sfr.c allvars.h tags.h proto.h forcetree.h cooling.h
sidm.o: sidm.c allvars.h tags.h proto.h forcetree.h cooling.h sidm.h \
  sidm_rand.h
sidm_rand.o: sidm_rand.c
system.o: system.c allvars.h tags.h proto.h forcetree.h cooling.h
timeline.o: timeline.c allvars.h tags.h proto.h forcetree.h cooling.h
timestep.o: timestep.c allvars.h tags.h proto.h forcetree.h cooling.h
veldisp.o: veldisp.c allvars.h tags.h proto.h forcetree.h cooling.h


.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)
