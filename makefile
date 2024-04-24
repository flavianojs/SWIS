SHELL = /bin/sh
####################################################
####################################################
#             SHE makefile
####################################################
#
# to compile use the following options:
#
# make          		       -  default compilation using mpif90
#
# make PLATFORM=(location) -  uses specific compiler options
# Configured locations:
# uff, iff, osx, juropa, juropatest, jureca, juqueen
#
# make DEBUG=debug         -  uses debug flags
#
# make PARALLEL=omp        -  uses openmp
#
# make FILE=(filename.exe) -  creates executable (filename.exe)
#
# make clean    			     -  removes *.o *.mod *__genmod* *.exe files
#
# make cleanall  			     -  also removes dependency files .dep and executable
#
# make recompile  			   -  removes all objects and modules and recompile
#
####################################################
# Filename and folders configuration               #
####################################################
SRCDIR = ./source
BINDIR = ./bin
OBJDIR = ./build$(addprefix /,$(PLATFORM))
ifdef FILE
	FILENAME = $(BINDIR)/$(FILE)
else
	FILENAME = $(BINDIR)/main.exe
endif

####################################################
# Suffixes used in this makefile                   #
####################################################
.SUFFIXES:
.SUFFIXES: .f90 .o .mod .dep

####################################################
# Source files                                     #
####################################################
SRC =$(wildcard $(SRCDIR)/*.f90)

####################################################
# Objects to compile                               #
####################################################
OBJ = $(addprefix $(OBJDIR)/,$(notdir $(SRC:.f90=.o)))

####################################################
# Dependency files                                #
####################################################
DEP = $(OBJ:.o=.dep)

####################################################
# Performance tool                                #
####################################################
ifeq ($(PERFORM),scalasca)
	PREP = scorep
endif

#=======================================================================
#============================ DEFAULT VALUES ===========================
#=======================================================================

####################################################
#  Compiler                                        #
####################################################
FC = ifort
#FC = mpif90

####################################################
#  Libraries                                       #
####################################################
#LLIBS =-mkl
LLIBS = -qmkl

####################################################
#  Flags                                           #
####################################################
#FFLAGS =-O3 -xHost -heap-arrays -fPIC
# FFLAGS =
# FFLAGS =-mcmodel=medium
FFLAGS = -fopenmp -O3 -xHost

####################################################
#  Preprocessor                                    #
####################################################
# CPP = -fpp
CPP = 

####################################################
#  Debugger Flags                                  #
####################################################
ifeq ($(DEBUG),debug)
FFLAGS =-CB -check all -check uninit -ftrapuv -debug all -traceback -g -warn all -O0
endif

####################################################
#  Module and Include folder                       #
####################################################
#FFLAGS += -module $(OBJDIR)/ -I$(OBJDIR)/

####################################################
#  Parallelization                                 #
####################################################
# ifeq ($(PARALLEL),omp)
LLIBS +=-qopenmp
#FFLAGS +=-qopenmp
FFLAGS +=-fopenmp
# endif

#=======================================================================
#========================== LOCATION SPECIFIC ==========================
#=======================================================================

####################################################
#                    PSI mpc2976                   #
####################################################
ifeq ($(PLATFORM),psi)
# Compiler
FC = ifort
# Preprocessor
CPP =
# Libraries
LLIBS = -qmkl
#Flags
FFLAGS = -fopenmp -O3 -xHost
endif
####################################################
#                       OSX                        #
####################################################
ifeq ($(PLATFORM),osx)
# Compiler
FC = mpif90
# Preprocessor
CPP = -fpp
# Libraries
LLIBS =-mkl -L$(HOME)/lib -lkibe
#Flags
FFLAGS =-O3 -xHost -qoverride-limits #-guide -parallel # -vec-threshold0
#Debugger
ifeq ($(DEBUG),debug)
FFLAGS =-CB -check all -check uninit -ftrapuv -debug all -traceback -g -warn all -O0
endif
FFLAGS += -module $(OBJDIR)/ -I$(OBJDIR)/
# Parallelization
ifeq ($(PARALLEL),omp)
LLIBS +=-qopenmp
FFLAGS +=-qopenmp
endif
endif
####################################################
#                       Ubuntu                     #
####################################################
ifeq ($(PLATFORM),ubu)
# Compiler
#FC = gfortran
FC = ifort
# Preprocessor
CPP = 
# Libraries
#LLIBS = -fopenmp -ffree-line-length-1000 -L/usr/lib -llapack -L/usr/lib -lblas
LLIBS = -fopenmp -L/usr/lib -llapack -L/usr/lib -lblas -fPIE
#Flags
#FFLAGS =-fopenmp -ffree-line-length-1000 -L/usr/lib -llapack -L/usr/lib -lblas
FFLAGS =-fopenmp -L/usr/lib -llapack -L/usr/lib -lblas -fPIE
endif
####################################################
#                       UFF                        #
####################################################
ifeq ($(PLATFORM),uff)
# Compiler
FC = mpif90
# Preprocessor
CPP = -fpp -D _UFF
# Libraries
LLIBS =-mkl -static-intel -L$(HOME)/lib -lnag
#Flags
FFLAGS =-O3 -xSSE4.2
#Debugger
ifeq ($(DEBUG),debug)
FFLAGS =-CB -check all -check uninit -ftrapuv -debug all -traceback -g -warn all -O0
endif
FFLAGS += -module $(OBJDIR)/ -I$(OBJDIR)/
# Parallelization
ifeq ($(PARALLEL),omp)
LLIBS +=-openmp
FFLAGS +=-openmp
endif
endif
####################################################
#                       IFF                        #
####################################################
ifeq ($(PLATFORM),iff)
# Compiler
FC = mpiifort
# Preprocessor
CPP = -fpp
# Libraries
LLIBS =-mkl -L$(HOME)/lib -lnag
# Flags
FFLAGS =-O3 -xHost
# Debugger
ifeq ($(DEBUG),debug)
FFLAGS =-CB -check all -check uninit -ftrapuv -debug all -traceback -g -warn all -O0
endif
FFLAGS += -module $(OBJDIR)/ -I$(OBJDIR)/
# Parallelization
ifeq ($(PARALLEL),omp)
LLIBS +=-openmp
FFLAGS +=-openmp
endif
endif

FFLAGS += -module $(OBJDIR)/ -I$(OBJDIR)/

####################################################
# Linking                                          #
####################################################
all: $(FILENAME)

$(FILENAME): $(OBJ)
	@echo Creating executable $(addprefix $(BINDIR)/,$(notdir $@)) $(and $(strip $(DEBUG) $(PARALLEL) $(PERFORM)), with $(strip $(DEBUG) $(PARALLEL) $(PERFORM)))
	@$(FC) $^ -o $@ $(LLIBS)

$(OBJDIR)/%.o $(OBJDIR)/%.mod: $(SRCDIR)/%.f90
	@echo Compiling file $(addprefix $(SRCDIR)/,$(notdir $<)) to $(addprefix $(OBJDIR)/,$(notdir $(patsubst %.f90,%.o,$<))) $(and $(strip $(DEBUG) $(PARALLEL) $(PERFORM)), with $(strip $(DEBUG) $(PARALLEL) $(PERFORM)))
	@$(FC) $(FFLAGS) $(CPP) -c $< -o $(addprefix $(OBJDIR)/,$(notdir $(patsubst %.f90,%.o,$<)))

####################################################
# Dependencies                                     #
####################################################
# $(OBJDIR)/%.dep: $(SRCDIR)/%.f90 f90_mod_deps.py
# 	@echo Building dependency of file $(addprefix $(SRCDIR)/,$(notdir $<))
# 	@./f90_mod_deps.py -o $@ -d "(mod_.*)" -D "$(OBJDIR)/\1.mod" -m "(.*)" -M "$(OBJDIR)/\1.mod" -O "$(OBJDIR)/\1.o" $<

ifeq ($(filter $(MAKECMDGOALS),clean cleanall),)
-include $(DEP)
endif

####################################################
# Clean up                                         #
####################################################
cleandep:
	@echo Removing dependency files...
	@rm -f $(DEP)

cleandebug:
	@echo Removing debug files...
	@rm -f $(OBJDIR)/*genmod*
	@rm -f *prep.opari*
	@rm -f *.i90
	@rm -f $(OBJDIR)/*.optrpt

cleanobj:
	@echo Removing object files...
	@rm -f $(OBJDIR)/*.o

cleanmod:
	@echo Removing module files...
	@rm -f $(OBJDIR)/*.mod

cleanexe:
	@echo Removing executable $(FILENAME)...
	@rm -f $(FILENAME)

clean: cleandebug cleanobj cleanmod

cleanall: cleandep cleandebug cleanobj cleanmod cleanexe

recompile: clean all
