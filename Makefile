#*****************************************************
#  General Makefile
#
#*****************************************************

#**************************************************************************
# Macros defining the Fortran, C/C++ compiler and linker.

CC = mpifort
LINK = mpifort
#FFLAGS = -g -fbacktrace -ffpe-trap=zero,overflow,underflow
FFLAGS = -O3 
BLAS_UBUNTU = -L/usr/lib/x86_64-linux-gnu -lblas
#**************************************************************************
# List of .o files that EXENAME depends on.  Edit as appropriate for MP.

OBJS = \
	DataStruct/NumConst.o DataStruct/PhysConst.o DataStruct/Pgrid.o \
	DataStruct/Data.o \
	Func/Timer.o Func/Transpose.o Func/Fldmger.o Func/Ptclmger.o  \
	Func/FFT.o Func/Bessel.o Func/Filter.o Func/Utility.o \
        Appl/BPM_NLLThick.o Appl/CCL.o Appl/CCDTL.o Appl/DTL.o Appl/SC.o \
	Appl/DriftTube.o Appl/Quadrupole.o Appl/ConstFoc.o Appl/SolRF.o \
	Appl/Sol.o Appl/Dipole.o Appl/Multipole_NLL.o Appl/EMfld.o Appl/TWS.o \
	Appl/NonlinearLens_NLLThick_Reversible.o Appl/BeamLineElem_NLLThick.o \
	Appl/CompDom.o Appl/BeamBunch_NLLThick_SymplecticSC.o Appl/sym2dsolver.o \
        Appl/Field_2D.o Appl/Distribution.o \
	Contrl/Input.o Contrl/Output.o Contrl/AccSimulator_NLLThick_SympSC.o Contrl/main.o

OBJS1 = \
	base.o beam_mod.o book.o cex.o \
        comm1.o elem.o genm.o hamdrift.o hmltn12345.o liea.o \
	liea_mod.o math.o parallel_mod.o pcmap.o pure.o sbend.o set_pscale.o sss.o trac.o

OBJS2 = \
	NumConst.o PhysConst.o Pgrid.o Data.o \
        Timer.o Transpose.o Fldmger.o Ptclmger.o FFT.o Bessel.o Filter.o Utility.o \
	BPM_NLLThick.o CCL.o CCDTL.o DTL.o SC.o DriftTube.o Quadrupole.o ConstFoc.o \
	SolRF.o Sol.o Dipole.o Multipole_NLL.o EMfld.o TWS.o \
	NonlinearLens_NLLThick_Reversible.o BeamLineElem_NLLThick.o CompDom.o \
	BeamBunch_NLLThick_SymplecticSC.o sym2dsolver.o Field_2D.o Distribution.o \
	Input.o Output.o AccSimulator_NLLThick_SympSC.o main.o	
#**************************************************************************
# Change this line if you don't like 'a.out'.

EXENAME = xmain

#************************************************************************
# disable predefined suffixes and define your own set of allowed suffixes
 .SUFFIXES:
 .SUFFIXES: .o .f .F .c .f90 .F90

#*************************************************************************
# inference rules (how to compile object files that have no explicit rules)
#  $* = name part of target
#  $@ = full target name
#  $< = dependent name

.f90.o:
	$(CC) -c  $(FFLAGS) $<

#**************************************************************************
# Rules for building EXENAME from OBJS and OBJS from your source.

$(EXENAME):  $(OBJS) 
	$(LINK) -o $(EXENAME) $(OBJS1) $(OBJS2) $(BLAS_UBUNTU)
	
all: $(EXENAME)
#************************************************************************
# if you wish to compile a certain object with different flags
# or in some special way, then specify the target & dependency explicitly
# the following line say Timer.o is depended on Timer.f90
#Timer.o: Timer.f90
#	$(CC) -c -O3 Timer.f90

	cp  AccSimulator_NLLThick_SympSC.o main.o Input.o Output.o Filter.o Utility.o Contrl
	cp  BPM_NLLThick.o CCL.o CCDTL.o DTL.o SC.o DriftTube.o Quadrupole.o \
	    ConstFoc.o BeamLineElem_NLLThick.o BeamBunch_NLLThick_SymplecticSC.o \
            Field_2D.o CompDom.o \
	    Multipole_NLL.o sym2dsolver.o Distribution.o SolRF.o Sol.o Dipole.o TWS.o \
	    NonlinearLens_NLLThick_Reversible.o EMfld.o Appl
	cp  Timer.o Transpose.o Fldmger.o Ptclmger.o FFT.o Bessel.o Func	
	cp  NumConst.o PhysConst.o Data.o Pgrid.o DataStruct 
#***********************************************************************
clean:
	-rm *.o $(EXENAME) *.mod Appl/*.o Func/*.o Contrl/*.o DataStruct/*.o
