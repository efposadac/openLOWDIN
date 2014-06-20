TOPDIR=../../
include $(TOPDIR)/CONFIG

TARGET=lowdin-MolecularMechanics.x
LDR=$(FC)
INC=-I$(TOPDIR)src/core/
INCI=-I$(TOPDIR)src/ints/
LDPATH=$(TOPDIR)lib/lowdincore.a
FLAGS=$(FCFLAGS)
LIBS=$(ALL_LIBS)

OBJECTS=AromaticityFinder.o AtomTypeUFF.o MMCommons.o MMFunctions.o MolecularMechanics.o RingFinder.o  

$(TARGET) : $(OBJECTS)
	$(LDR) $(FLAGS) -o $(TARGET) $(OBJECTS) $(LDPATH) $(LIBS)
	mv $(TARGET) $(TOPDIR)bin/

AromaticityFinder.o:	AromaticityFinder.f90 MMCommons.o RingFinder.o 
	$(FC) -c $(INC) $(INCI) $(FLAGS) AromaticityFinder.f90

AtomTypeUFF.o:	AtomTypeUFF.f90 MMCommons.o AromaticityFinder.o 
	$(FC) -c $(INC) $(INCI) $(FLAGS) AtomTypeUFF.f90

MMCommons.o:	MMCommons.f90 
	$(FC) -c $(INC) $(INCI) $(FLAGS) MMCommons.f90

MMFunctions.o:	MMFunctions.f90 AtomTypeUFF.o 
	$(FC) -c $(INC) $(INCI) $(FLAGS) MMFunctions.f90

MolecularMechanics.o:	MolecularMechanics.f90 MMFunctions.o 
	$(FC) -c $(INC) $(INCI) $(FLAGS) MolecularMechanics.f90

RingFinder.o:	RingFinder.f90 MMCommons.o 
	$(FC) -c $(INC) $(INCI) $(FLAGS) RingFinder.f90

clean:: 
	rm -f *.o
	rm -f *.mod
	rm makefile
