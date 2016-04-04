#  There are a number of code options that need to be switched on/off
#  on compile-time.
#
#  From the following list, please select the desired options:

#OPT   +=  -DPERIODIC       # enables periodic boundaries conditions
#OPT   +=  -DCOOLING        # radiative heating and cooling
#OPT   +=  -DDIAG           # outputs performance diagnostic for tree walks
OPT   +=  -DINLINE         # compiler accepts the `inline' statement
#OPT   +=  -DT3E            # needed on the T3E
#OPT   +=  -DRECTREECONS    # switches on recursive computation of multipole moments
OPT   +=  -DBMAX           # enables very conservative node-opening 
#OPT   +=  -DVELDISP        # selects computation of velocity dispersion/density for dark matter

#OPT   +=  -DMOREPARAMS    # enables additional parameters in parameterfile for star formation sector
#OPT   +=  -DSFR           # switches on star formation
#OPT   +=  -DCLOUDS        # switches on multi-phase model for star formation
#OPT   +=  -DSTELLARAGE    # Switches on computation of mean time of stellar formation redshift/age
#OPT13   =  -DSTOCHASTICSFR # switches on stochastic star formation events as opposed to continuous treatment in standard scheme
OPT   +=  -DSIDM
OPT   +=  -DREFLECTIONBOUNDARY
OPT   +=  -DSCATTERLOG
OPT   +=  -DFINDNBRLOG
OPT   +=  -DNOSCATTER
OPT   +=  -DRANDOM_GSL     
OPT   +=  -DCROSS_SECTION_TYPE=4

# 0: Hard sphere (velocity independent)
# 1: Maxwellian (sigma ~ v^-1)
# 2: Yukawa-like
# 3: Power law (sigma ~ v^-n)
# 4: Yukawa interaction

CC       =   mpicc   # sets the C-compiler
OPTIMIZE =   -O2 # -xT # -Wall  

EXEC   = P-Gadget

OBJS   = main.o  run.o  predict.o begrun.o endrun.o global.o  \
	 timestep.o  init.o restart.o  io.o sfr.o \
	 accel.o   read_ic.o  cooling.o dissolvegas.o \
	 system.o  allocate.o  density.o \
	 gravtree.o timeline.o hydra.o veldisp.o \
	 domain.o  allvars.o potential.o \
         forcetree.o  ewald.o read_ic_cluster.o sidm.o reflection.o\
         sidm_rand.o

INCL   = allvars.h  proto.h  forcetree.h tags.h cooling.h domain.h

ifeq ($(SYSTYPE),TACC)
CC      =  mpicc
OPTIMIZE += -xT
GSL_INC =  -I$(TACC_GSL_INC)
GSL_LIB =  -L$(TACC_GSL_LIB)
endif

ifeq ($(SYSTYPE),GALILEO)
CC       =  mpicc   
GSL_INC  = -I$(GSL_INC_GCC)
GSL_LIB  = -L$(GSL_LIB_GCC)
endif

ifeq ($(SYSTYPE),)
CC       = mpicc
endif

OPTIONS =  $(OPTIMIZE) $(GADGET_FLAGS) $(OPT) 
CFLAGS =   $(OPTIONS) $(GSL_INC)
LIBS   =   -lm -lgsl -lgslcblas $(GSL_LIB)

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 


.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)







