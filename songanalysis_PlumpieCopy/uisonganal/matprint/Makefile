
prefix = /usr/local
libdir = ${prefix}/lib
matlabdir = ${libdir}/matlab/
matlabsitedir = ${matlabdir}/site-packages/
matlabpkgdir = ${matlabsitedir}/uisonganal/

MEX = /usr/bin/mex

SUBDIRS =

matlab_files = \
	pgsetup.m \
	print2.m \
	printopt2.m \
	uiprint.m

mex_targets =

CLEAN_FILES =

all_targets = \
	${mex_targets}

all: all-subdirs $(all_targets)

all-subdirs: dummy
	@if test x"$(SUBDIRS)" != x ; then \
	  for i in $(SUBDIRS); do \
	    make -C $$i ; \
	  done ; \
	fi

install: dummy
	install -d ${matlabpkgdir}/matprint
	for i in ${matlab_files}; do \
	  install -m644 $$i ${matlabpkgdir}/matprint; \
	done
	for i in ${mex_targets}; do \
	  install -m755 $$i ${matlabpkgdir}/matprint; \
	done

clean-subdirs: dummy
	@for i in $(SUBDIRS); do \
	  make -C $$i clean; \
	done

clean: clean-subdirs
	-rm -f *.o
	rm -f $(CLEAN_FILES)

dummy:

