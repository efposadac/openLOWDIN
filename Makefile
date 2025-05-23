TOPDIR=.
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(TOPDIR)/CONFIG

ifeq ($(strip $(PREFIX)),)
	PREFIX="/usr/local"
endif
ifeq ($(strip $(EXENAME)),)
	EXENAME="openlowdin"
endif

SUBDIRS = utilities src
ALLSUBDIRS = $(SUBDIRS) doc bin

default::
	for dir in $(SUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(JOBS)) || exit 1; \
	  done

doc::
	cd $(TOPDIR)/doc && doxygen Doxyfile

install:: bin/lowdin bin/lowdin.x
	if [ -e $(PREFIX)/.$(EXENAME) ]; then \
		rm -rf $(PREFIX)/.$(EXENAME)/bin \
		rm -rf $(PREFIX)/.$(EXENAME)/lib ; fi
	mkdir -p $(PREFIX)/.$(EXENAME)
	cp -rf $(TOPDIR)/bin/lowdinvars.sh $(TOPDIR)
	$(SED) -i  's|PREFIX|$(PREFIX)|g' $(TOPDIR)/lowdinvars.sh
	$(SED) -i  's|EXENAME|$(EXENAME)|g' $(TOPDIR)/lowdinvars.sh
	$(SED) -i  's|SCRATCH_DIR|$(SCRATCH)|g' $(TOPDIR)/lowdinvars.sh
	$(SED) -i "s|COMMIT_ID|$(shell git --no-pager log -1 --pretty=format:"%H")|g" $(TOPDIR)/lowdinvars.sh
	$(SED) -i 's|COMPILATION_DATE|$(shell date)|g' $(TOPDIR)/lowdinvars.sh
	cp -rf $(TOPDIR)/lowdinvars.sh $(PREFIX)/.$(EXENAME)/
	cp -rf lib/ $(PREFIX)/.$(EXENAME)/lib/	
	if [ -e utilities/erkale/build/erkale/basis ]; then \
		cp -rf utilities/erkale/build/erkale/basis $(PREFIX)/.$(EXENAME)/lib/erkaleBasis ; fi
	mkdir -p $(PREFIX)/.$(EXENAME)/bin
	cp -rf $(TOPDIR)/bin/*.x $(PREFIX)/.$(EXENAME)/bin
	if [ -e utilities/erkale/erkale/bin/erkale_loc ]; then \
		cp -rf utilities/erkale/erkale/bin/erkale_fchkpt utilities/erkale/erkale/bin/erkale_loc $(PREFIX)/.$(EXENAME)/bin ; fi
	cp -rf $(TOPDIR)/bin/lowdin $(TOPDIR)
	$(SED) -i  's|PREFIX|$(PREFIX)|g' $(TOPDIR)/lowdin
	$(SED) -i  's|EXENAME|$(EXENAME)|g' $(TOPDIR)/lowdin
	cp -rf $(TOPDIR)/lowdin $(PREFIX)/$(EXENAME)
	cp -rf $(TOPDIR)/lowdin $(PREFIX)/.$(EXENAME)/$(EXENAME)
	rm -rf $(TOPDIR)/lowdin
	rm -rf $(TOPDIR)/lowdinvars.sh

uninstall:: bin/lowdin bin/lowdin.x
	rm -rf $(PREFIX)/.$(EXENAME)
	rm -rf $(PREFIX)/$(EXENAME)

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
	find . -name "*.i90" -exec rm -f {} \;
	find . -name "*.i" -exec rm -f {} \;
	rm -rf $(TOPDIR)/doc/html
	rm -rf $(TOPDIR)/doc/latex

test::
	cd $(TOPDIR)/test/; \
	sh clean.sh $(EXENAME); \
	sh runtest.sh $(EXENAME) && sh clean.sh $(EXENAME) 


