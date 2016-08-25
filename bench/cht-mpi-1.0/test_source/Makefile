
SUBDIRS=  \
	 test_comp_comm_overlap \
	 test_matrix

all_tests:
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
	done; 
	@echo; 
	@echo 'Successfully completed all tests in $(SUBDIRS)'; 
	@echo;
	@stars="***************************************************" ; \
	printf "\n$$stars\n All tests finished OK. This is great! \n$$stars\n\n"

clean:
	list='$(SUBDIRS)'; \
	for subdir in $$list; do \
        (cd $$subdir && $(MAKE) clean) \
	done;
