TOPDIR=..

include $(TOPDIR)/CONFIG

SUBDIRS = ints HF scf integralsTransformation MP2 PT cosmo CalcProp output CI MM DFTB
ALLDIRS = $(SUBDIRS) core

LDR=$(FC)
INC=-I$(TOPDIR)/src/core/
LDPATH=$(TOPDIR)/lib/lowdincore.a
FLAGS=$(FCFLAGS)
LIBS=$(ALL_LIBS)

TARGET=lowdin.x
OBJECTS=lowdin.x.o

default::
	(cd core && $(MAKE) $(JOBS)) || exit 1;\

	for dir in $(SUBDIRS); \
	do \
	  (cd $${dir} && $(MAKE) $(JOBS)) || exit 1; \
	done

	$(FC) -c $(INC) $(FLAGS) lowdin.x.f90

	$(LDR) $(FLAGS) -o $(TARGET) $(OBJECTS) $(LDPATH) $(LIBS)
	mv $(TARGET) $(TOPDIR)/bin/

clean::
	for dir in $(ALLDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(DODEPENDOPT) clean) || exit 1; \
	  done


