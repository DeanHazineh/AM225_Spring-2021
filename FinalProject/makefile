# Include the configuration file from the parent directory
include ../config.mk

# Lists of files to be built
objs = FDFD.o file_output.o
src=$(patsubst %.o,%.cc,$(objs) $(mg_objs))
execs = FDFD_test1

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

FDFD_test1: FDFD_test1.cc $(objs)
	$(cxx) $(cflags)  -o $@ $^ $(lp_lflags)


.PHONY: clean depend
