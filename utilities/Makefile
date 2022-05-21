TOPDIR=..

include $(TOPDIR)/CONFIG

SUBDIRS = aduw gepol molden2aim jadamilu
ALLDIRS = $(SUBDIRS) 

default::

	for dir in $(SUBDIRS); \
	do \
	  (cd $${dir} && $(MAKE) $(JOBS)) || exit 1; \
	done

clean::
	for dir in $(ALLDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(DODEPENDOPT) clean) || exit 1; \
	  done


