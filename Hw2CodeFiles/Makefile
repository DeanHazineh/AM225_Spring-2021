# Include the configuration file from the parent directory
include ../config.mk

# Lists of files to be built
objs = sol_IRK5.o sol_rk4d.o sol_rk4classic.o 
src=$(patsubst %.o,%.cc,$(objs) $(mg_objs))
execs= driver_1p driver_6p driver_7p driver_KuramotoODE

all: $(execs)

# Include the file dependencies
-include Makefile.dep

# A Makefile target to refresh the dependency file
depend:
	$(cxx) $(iflags) -MM $(src) >Makefile.dep

# A Makefile target to remove all the built files
clean:
	rm -f $(objs) $(execs)

%.o: %.cc
	$(cxx) $(cflags) $(iflags) -c $<

driver_1p: driver_1p.cc $(objs)
	$(cxx) $(cflags) -o $@ $^

driver_6p: driver_6p.cc $(objs)
	$(cxx) $(cflags) -o $@ $^

driver_7p: driver_7p.cc $(objs)
	$(cxx) $(cflags) -o $@ $^

driver_KuramotoODE: driver_KuramotoODE.cc $(objs)
	$(cxx) $(cfalgs) -o $@ $^

.PHONY: clean depend
