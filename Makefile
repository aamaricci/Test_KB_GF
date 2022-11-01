##$ COMPILER: suppprted compilers are ifort, gnu >v4.7 or use mpif90
FC=mpif90

##$ PLATFORM: supported platform are intel, gnu 
PLAT=gnu

##$ SET THE TARGET DIRECTORY WHERE TO PUT THE EXECUTABLE (default if $HOME/.bin in the PATH)
DIREXE=test

##$ CHOOSE THE DRIVER CODE:
EXE=testOP_DERIVED3

##$ SET INCLUDE AND LINK OPTIONS USING pkg-config
INCARGS=$(shell pkg-config --cflags dmft_tools scifor)
LIBARGS=$(shell pkg-config --libs   dmft_tools scifor)

SFINC=$(shell pkg-config --cflags scifor)
SFLIB=$(shell pkg-config --libs   scifor)


ifeq ($(PLAT),intel)
FFLAG=-O2 -ftz
OFLAG=-O3 -ftz
DFLAG=-p -O0 -g -fpe0 -warn -warn errors -debuEg extended -traceback -check all,noarg_temp_created
FPPFLAG =-fpp
endif

ifeq ($(PLAT),gnu)
FFLAG = -O2 -ffree-line-length-none
DFLAG = -O0 -p -g -fimplicit-none -Wsurprising  -Waliasing -fwhole-file -fcheck=all -pedantic -fbacktrace -ffree-line-length-none
OFLAG = -O3 -ffast-math -march=native -funroll-all-loops -fno-protect-parens -flto -ffree-line-length-none
FPPFLAG =-cpp
endif



##$ Extends the implicit support of the Makefile to .f90 files

.SUFFIXES: .f90

OBJS =  KB_VARS_GLOBAL.o \
	KB_CONTOUR.o \
	KB_AUX.o \
	KB_GF_COMMON.o \
	KB_GF_SUM.o \
	KB_GF_CONVOLUTE.o \
	KB_GF_VIE.o \
	KB_GF_DYSON.o \
	KB_GF_FREE.o \
	KB_GF_MAIN.o \
	KB_LIB.o


all: all libkbgf.a compile

lib: lib libkbgf.a

debug: debug libkbgf.a compile


lib: FPPFLAG+=-D_

all: FPPFLAG+=-D_

debug: FFLAG=$(DFLAG)
debug: FPPFLAG+=-D_


libkbgf.a:$(OBJS)
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@
	@echo 'LIBRARY IS DONE.........'
	@echo ''

compile: $(OBJS)
	@echo " ..................... compile ........................... "
	$(FC) $(FFLAG) $(INCARGS) $(OBJS) $(EXE).f90 -o $(DIREXE)/$(EXE)$(BRANCH) $(LIBARGS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)$(BRANCH)


.f90.o:	
	$(FC) $(FFLAG) $(SFINC) -c $<


clean: 
	@echo "Cleaning:"
	@rm -fv *.mod *.o *~ *.a





