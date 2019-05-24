
SUBDIRS = stencil1d stencil3d smith-waterman cholesky
#SUBDIRS = stencil1d stencil3d smith-waterman
     
subdirs:
	for dir in $(SUBDIRS); do \
    	$(MAKE) -C $$dir; \
    done

clean:
	for dir in $(SUBDIRS); do \
    	$(MAKE) -C $$dir clean; \
    done
