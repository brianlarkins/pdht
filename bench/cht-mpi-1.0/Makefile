SUBDIRS= source \
	 test_source 
all:
	list='$(SUBDIRS)'; \
	for subdir in $$list; do \
	(cd $$subdir && $(MAKE)) \
	|| exit 1; \
	done; \
	echo; \
	echo 'Updated all targets in directories: $(SUBDIRS)'; \
	echo;
check:
	list='$(SUBDIRS)'; \
	for subdir in $$list; do \
	(cd $$subdir && $(MAKE) check) \
	|| exit 1; \
	done; \
	echo; \
	echo 'Successfully completed all tests in $(SUBDIRS)'; \
	echo;
clean:
	list='$(SUBDIRS)'; \
	for subdir in $$list; do \
	(cd $$subdir && $(MAKE) clean) \
	done;

isis:
	$(MAKE) BLASFLAGS='/bubo/home/h20/eliasr/acml/installation/gfortran64/lib/libacml.a  /bubo/home/h20/eliasr/acml/installation/gfortran64/lib/libacml_mv.a -lgfortran'
