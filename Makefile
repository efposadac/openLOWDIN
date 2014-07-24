TOPDIR=.
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(TOPDIR)/CONFIG

SUBDIRS = src
ALLSUBDIRS = $(SUBDIRS) doc bin

default::
	for dir in $(SUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(JOBS)) || exit 1; \
	  done

doc::
	cd $(TOPDIR)/doc && doxygen Doxyfile


install:: bin/lowdin bin/lowdin.x
	mkdir -p ${HOME}/.lowdin2
	cp -rf $(TOPDIR)/bin/lowdinvars.sh ${HOME}/.lowdin2/
	cp -rf lib ${HOME}/.lowdin2/
	mkdir -p ${HOME}/.lowdin2/bin
	cp -rf $(TOPDIR)/bin/*.x ${HOME}/.lowdin2/bin
	cp -rf $(TOPDIR)/bin/lowdin $(PREFIX)/lowdin2

uninstall:: bin/lowdin bin/lowdin.x
	rm -rf ${HOME}/.lowdin2
	rm -rf $(PREFIX)/lowdin2

clean::
	for dir in $(SUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(DODEPENDOPT) clean) || exit 1; \
	  done

distclean::
	find . -name "*.o" -exec rm -f {} \;
	find . -name "*.mod" -exec rm -f {} \;
	find . -name "*~" -exec rm -f {} \;
	find . -name "*.x" -exec rm -f {} \;
	find . -name "*.a" -exec rm -f {} \;
	rm -rf $(TOPDIR)/doc/html
	rm -rf $(TOPDIR)/doc/latex
