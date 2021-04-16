# Include the configuration file 	
include ../config.mk

# Lists of files to be built
objs= ho_transport.o LaxWendroffSolver.o
src=$(patsubst %.o,%.cc,$(objs))
execs= ho_solve LaxWendroffDriver

all: $(execs)

# Include the file dependencies
-include Makefile.dep

# A Makefile target to refresh the dependency file
depend:
	$(cxx) -MM $(src) >Makefile.dep

# A Makefile target to remove all the built files
clean:
	rm -f $(objs) $(execs)

%.o: %.cc
	$(cxx) $(cflags) -c $< $(ode_iflags)

LaxWendroffDriver: LaxWendroffDriver.cc $(objs)
	$(cxx) $(cflags) -o $@ $^ $(lp_lflags)

ho_solve: ho_solve.cc $(objs)
	$(cxx) $(cflags) -o $@ $^ $(lp_lflags)


.PHONY: clean depend
