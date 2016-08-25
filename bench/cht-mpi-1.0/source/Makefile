SUBDIRS= services \
	mpi_wrapper \
	utilities

TMPOBJDIR = tmp_objdir_for_make

all:
	list='$(SUBDIRS)'; \
	for subdir in $$list; do \
	(cd $$subdir && $(MAKE)) \
	|| exit 1; \
	done; \
	echo; \
	echo 'Updated all targets in directories: $(SUBDIRS)'; \
	echo; \
	$(MAKE) libs

libs: libcht.a 
	

libcht.a: services/libservice.a mpi_wrapper/libmpiwrapper.a utilities/libutils.a
	(rm -rf $(TMPOBJDIR) && mkdir $(TMPOBJDIR) && \
	 cd $(TMPOBJDIR) && \
	 ar x ../services/libservice.a && \
	 ar x ../utilities/libutils.a && \
	 ar x ../mpi_wrapper/libmpiwrapper.a) \
	|| exit 1; \
	ar rcs libcht.a $(TMPOBJDIR)/*.o;

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
	rm -rf *.a $(TMPOBJDIR)
